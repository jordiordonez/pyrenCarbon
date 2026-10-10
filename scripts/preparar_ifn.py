"""Construeix el conjunt de calibratge de creixement i mortalitat a partir de l'IFN2 i l'IFN3 de Lleida.

Ús:
    python scripts/preparar_ifn.py DIR_IFN

DIR_IFN ha de contenir els fitxers públics del MITECO descomprimits:
    ifn2_25/PIESMA25.DBF, ifn2_25/DATEST25.DBF  (IFN2, Lleida)
    ifn3p25/Ifn3p25.mdb                          (IFN3, Lleida, camp)
    sig25/Sig_25.mdb                             (IFN3, Lleida, SIG)
Cal tenir `mdb-export` (mdbtools) al PATH.

Escriu data/ifn_lleida_arbres.csv: un registre per pi (roig o negre) de l'IFN2 en parcel·les
remesurades a l'IFN3 sense tallades, amb el seu estat el 1990, la competència de la parcel·la i
el resultat el 2001 (viu amb diàmetre nou, o mort).
"""

import io
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from dbfread import DBF

DESTI = Path(__file__).resolve().parent.parent / "data" / "ifn_lleida_arbres.csv"
ESPECIES = {"21": "Pinus sylvestris", "22": "Pinus uncinata"}


def factor_expansio(dbh_cm: pd.Series) -> pd.Series:
    """Peus/ha que representa cada arbre segons les parcel·les concèntriques de l'IFN."""
    radi = np.select([dbh_cm < 12.5, dbh_cm < 22.5, dbh_cm < 42.5], [5.0, 10.0, 15.0], 25.0)
    return pd.Series(1e4 / (np.pi * radi**2), index=dbh_cm.index)


def llegir_mdb(ruta: Path, taula: str) -> pd.DataFrame:
    sortida = subprocess.run(["mdb-export", str(ruta), taula], check=True, capture_output=True, text=True).stdout
    return pd.read_csv(io.StringIO(sortida), dtype=str)


def llegir_dbf(ruta: Path) -> pd.DataFrame:
    return pd.DataFrame(iter(DBF(str(ruta), encoding="latin-1")))


def arbres_ifn2(dir_ifn: Path) -> pd.DataFrame:
    a = llegir_dbf(dir_ifn / "ifn2_25" / "PIESMA25.DBF")
    a["especie"] = a["ESPECIE"].str.strip()
    d1 = pd.to_numeric(a["DIAMETRO1"].str.strip(), errors="coerce")
    d2 = pd.to_numeric(a["DIAMETRO2"].str.strip(), errors="coerce")
    a["dbh"] = (d1 + d2) / 20.0
    a["ht"] = pd.to_numeric(a["ALTURA"].str.strip().str.replace(",", "."), errors="coerce")
    a = a[a["dbh"] >= 7.5].rename(columns={"ESTADILLO": "parcela", "NUMORDEN": "ordre"})
    a["f"] = factor_expansio(a["dbh"])
    a["g"] = np.pi * (a["dbh"] / 200.0) ** 2 * a["f"]
    return a[["parcela", "ordre", "especie", "dbh", "ht", "f", "g"]]


def competencia(arbres: pd.DataFrame) -> pd.DataFrame:
    """AB i N de la parcel·la, i BAL (AB dels arbres més grossos) de cada arbre."""
    arbres = arbres.sort_values(["parcela", "dbh"], ascending=[True, False]).copy()
    acumulada = arbres.groupby("parcela")["g"].cumsum()
    arbres["bal"] = acumulada - arbres["g"]
    parcela = arbres.groupby("parcela").agg(ab=("g", "sum"), n=("f", "sum"))
    return arbres.join(parcela, on="parcela")


def main(dir_ifn: str) -> None:
    dir_ifn = Path(dir_ifn)
    ifn2 = competencia(arbres_ifn2(dir_ifn))

    parceles2 = llegir_dbf(dir_ifn / "ifn2_25" / "DATEST25.DBF")
    parceles2 = parceles2.rename(columns={"ESTADILLO": "parcela"})
    parceles2["any_ifn2"] = 1900 + pd.to_numeric(parceles2["ANO"], errors="coerce")
    parceles2["altitud"] = pd.to_numeric(parceles2["ALTITUD1"], errors="coerce") * 100.0
    parceles2["utm_x"] = pd.to_numeric(parceles2["COORDEX"], errors="coerce") * 1000.0
    parceles2["utm_y"] = pd.to_numeric(parceles2["COORDEY"], errors="coerce") * 1000.0

    mdb = dir_ifn / "ifn3p25" / "Ifn3p25.mdb"
    parceles3 = llegir_mdb(mdb, "PCParcelas")
    parceles3 = parceles3[(parceles3["Cla"] == "A") & (parceles3["Subclase"] == "1")]
    parceles3 = pd.DataFrame(
        {"parcela": parceles3["Estadillo"].astype(int), "any_ifn3": parceles3["Ano"].astype(int)}
    )

    ifn3 = llegir_mdb(mdb, "PCMayores")
    ifn3 = ifn3[(ifn3["Cla"] == "A") & (ifn3["Subclase"] == "1")].copy()
    ifn3["parcela"] = ifn3["Estadillo"].astype(int)
    # OrdenIf3 = 000: arbre de l'IFN2 desaparegut (tallat); 888: mort; 999: no localitzat.
    parceles_tallades = set(ifn3.loc[ifn3["OrdenIf3"] == "000", "parcela"])
    remesurats = ifn3[ifn3["OrdenIf2"].str.fullmatch(r"\d{3}") & ~ifn3["OrdenIf2"].isin(["000", "999"])].copy()
    remesurats = remesurats[remesurats["OrdenIf3"] != "999"]
    remesurats["ordre"] = remesurats["OrdenIf2"].astype(int)
    remesurats["mort"] = remesurats["OrdenIf3"] == "888"
    remesurats["dbh_ifn3"] = (remesurats["Dn1"].astype(float) + remesurats["Dn2"].astype(float)) / 20.0
    remesurats["especie_ifn3"] = remesurats["Especie"].str.lstrip("0")
    remesurats = remesurats[["parcela", "ordre", "mort", "dbh_ifn3", "especie_ifn3"]]

    dades = ifn2[ifn2["especie"].isin(ESPECIES)].merge(remesurats, on=["parcela", "ordre"], how="inner")
    dades = dades[~dades["parcela"].isin(parceles_tallades)]
    dades = dades[dades["mort"] | (dades["especie_ifn3"] == dades["especie"])]
    dades = dades.merge(parceles2[["parcela", "any_ifn2", "altitud", "utm_x", "utm_y"]], on="parcela")
    dades = dades.merge(parceles3, on="parcela")
    dades["anys"] = dades["any_ifn3"] - dades["any_ifn2"]
    dades = dades[dades["anys"].between(8, 14)]
    dades["id_dbh"] = np.where(dades["mort"], np.nan, (dades["dbh_ifn3"] - dades["dbh"]) / dades["anys"])
    # Errors de mesura evidents: decreixements o creixements de més de 2 cm/any.
    dades = dades[dades["mort"] | dades["id_dbh"].between(-0.2, 2.0)]
    dades["especie"] = dades["especie"].map(ESPECIES)
    columnes = ["parcela", "ordre", "especie", "dbh", "ht", "f", "bal", "ab", "n", "altitud", "utm_x", "utm_y",
                "any_ifn2", "any_ifn3", "anys", "mort", "dbh_ifn3", "id_dbh"]
    dades[columnes].round(4).to_csv(DESTI, index=False)
    resum = dades.groupby("especie").agg(arbres=("ordre", "size"), parceles=("parcela", "nunique"), morts=("mort", "sum"))
    print(resum.to_string())


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])
