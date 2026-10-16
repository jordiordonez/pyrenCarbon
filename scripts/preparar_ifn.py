"""Construeix el conjunt de calibratge de creixement i mortalitat amb les remesures de l'IFN a Lleida.

Ús:
    python scripts/preparar_ifn.py DIR_IFN

DIR_IFN ha de contenir els fitxers públics del MITECO descomprimits:
    ifn2_25/PIESMA25.DBF, ifn2_25/DATEST25.DBF  (IFN2, Lleida)
    ifn3p25/Ifn3p25.accdb (o .mdb)               (IFN3, Lleida, camp)
    ifn4_25/Ifn4_Cataluna.accdb                  (IFN4, Catalunya, camp)
Cal tenir `mdb-export` (mdbtools) al PATH, i data/ifn_parceles_comarca.csv (scripts/assignar_comarca.py).

Dues remesures, cadascuna amb la mateixa lògica:
    IFN2-IFN3  estat inicial IFN2 (1989-90), final IFN3 (2000-01); interval de 8 a 14 anys
    IFN3-IFN4  estat inicial IFN3 (2000-01), final IFN4 (2014-16); interval de 8 a 16 anys
Només es fan servir les parcel·les de les comarques de COMARQUES (comarca del node de malla de l'IFN3).

Escriu data/ifn_lleida_arbres.csv: un registre per pi (roig o negre) viu a la mesura inicial en una
parcel·la remesurada sense tallades, amb el seu estat inicial, la competència de la parcel·la i el
resultat a la mesura final (viu amb diàmetre nou, o mort). Per als morts, l'agent causant que anota l'IFN
(camp Agente de l'IFN3 i l'IFN4) i el grup de causa (vegeu CAUSES). Columnes de període:
    periode   "IFN2-IFN3" o "IFN3-IFN4"
    comarca   comarca de la parcel·la
    agent     codi d'agent causant de la mort (camp Agente; buit si l'arbre és viu o no en té)
    causa     grup de causa de la mort (CAUSES); buit si l'arbre és viu
    any_ifn2  any de la mesura inicial; any_ifn3 i dbh_ifn3: any i diàmetre de la mesura final. Els noms
              es mantenen per compatibilitat: a les files IFN3-IFN4 corresponen a l'IFN3 i l'IFN4.
Una parcel·la remesurada en tots dos períodes manté el mateix número (`parcela`) a totes dues.
"""

from __future__ import annotations

import io
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from dbfread import DBF

DATA = Path(__file__).resolve().parent.parent / "data"
DESTI = DATA / "ifn_lleida_arbres.csv"
COMARCA = DATA / "ifn_parceles_comarca.csv"
ESPECIES = {"21": "Pinus sylvestris", "22": "Pinus uncinata"}
# Comarques pirinenques de Lleida més semblants a Andorra (acord de la reunió del 15/10/2026, opció c).
COMARQUES = ("Pallars Sobirà", "Alt Urgell", "Alta Ribagorça", "Cerdanya")
INTERVALS = {"IFN2-IFN3": (8, 14), "IFN3-IFN4": (8, 16)}
COLUMNES = ["periode", "comarca", "parcela", "ordre", "especie", "dbh", "ht", "f", "bal", "ab", "n", "altitud",
            "utm_x", "utm_y", "any_ifn2", "any_ifn3", "anys", "mort", "dbh_ifn3", "id_dbh", "agent", "causa"]
# Grups de causa de mort a partir de l'agent causant de l'IFN (camp Agente; documentador de l'IFN4, annex 8).
# «fons» és la mortalitat que calibra el model (competència, causa desconeguda o sense agent anotat); la resta
# són pertorbacions, que l'app tracta com a escenaris a part.
CAUSES = {
    "fons": (100, 200, 316),  # sense dany, causa desconeguda, dominància (competència); també sense agent
    "foc": (421, 413),  # incendi, llamp
    "sequera": (412, 414),  # sequera, gelada
    "neu_vent": (410, 411, 415),  # neu, vent, calamarsa
    "biotic": (300, 310, 311, 312, 313, 314, 315),  # biòtic genèric, fongs, insectes, vesc, epífits, fauna, bestiar
    # maquinària, aprofitament de fusta, persones, abiòtic genèric, despreniments, erosió
    "altres": (320, 321, 322, 400, 422, 423),
}
CAUSA_PER_AGENT = {agent: causa for causa, agents in CAUSES.items() for agent in agents}


def factor_expansio(dbh_cm: pd.Series) -> pd.Series:
    """Peus/ha que representa cada arbre segons les parcel·les concèntriques de l'IFN."""
    radi = np.select([dbh_cm < 12.5, dbh_cm < 22.5, dbh_cm < 42.5], [5.0, 10.0, 15.0], 25.0)
    return pd.Series(1e4 / (np.pi * radi**2), index=dbh_cm.index)


def llegir_mdb(ruta: Path, taula: str) -> pd.DataFrame:
    sortida = subprocess.run(["mdb-export", str(ruta), taula], check=True, capture_output=True, text=True).stdout
    return pd.read_csv(io.StringIO(sortida), dtype=str)


def llegir_dbf(ruta: Path) -> pd.DataFrame:
    return pd.DataFrame(iter(DBF(str(ruta), encoding="latin-1")))


def base_ifn3(dir_ifn: Path) -> Path:
    """Base de camp de l'IFN3 de Lleida (el MITECO la distribueix com a .accdb; abans, .mdb)."""
    for nom in ("Ifn3p25.accdb", "Ifn3p25.mdb"):
        if (dir_ifn / "ifn3p25" / nom).exists():
            return dir_ifn / "ifn3p25" / nom
    raise FileNotFoundError(dir_ifn / "ifn3p25" / "Ifn3p25.accdb")


def causa_mort(agent: pd.Series, mort: pd.Series) -> pd.Series:
    """Grup de causa (CAUSES) de cada arbre mort; un mort sense agent o amb un codi no previst compta com a «fons»."""
    causa = pd.to_numeric(agent, errors="coerce").map(CAUSA_PER_AGENT).fillna("fons")
    return causa.where(mort.astype(bool), None)


def ordre_valid(s: pd.Series) -> pd.Series:
    """Ordre de peu real (001-887), no un codi especial."""
    return s.fillna("").str.fullmatch(r"\d{3}") & ~s.isin(["000", "888", "999"])


def competencia(arbres: pd.DataFrame) -> pd.DataFrame:
    """AB i N de la parcel·la, i BAL (AB dels arbres més grossos) de cada arbre."""
    arbres = arbres.sort_values(["parcela", "dbh"], ascending=[True, False]).copy()
    acumulada = arbres.groupby("parcela")["g"].cumsum()
    arbres["bal"] = acumulada - arbres["g"]
    parcela = arbres.groupby("parcela").agg(ab=("g", "sum"), n=("f", "sum"))
    return arbres.join(parcela, on="parcela")


def amb_factor(a: pd.DataFrame) -> pd.DataFrame:
    a = a[a["dbh"] >= 7.5].copy()
    a["f"] = factor_expansio(a["dbh"])
    a["g"] = np.pi * (a["dbh"] / 200.0) ** 2 * a["f"]
    return a[["parcela", "ordre", "especie", "dbh", "ht", "f", "g"]]


def altitud_ifn2(dir_ifn: Path) -> pd.Series:
    """Altitud de la parcel·la a l'IFN2 (ALTITUD1, en centenars de metres), per número d'estadillo."""
    d2 = llegir_dbf(dir_ifn / "ifn2_25" / "DATEST25.DBF")
    return pd.Series(pd.to_numeric(d2["ALTITUD1"], errors="coerce").to_numpy() * 100.0,
                     index=pd.to_numeric(d2["ESTADILLO"], errors="coerce")).groupby(level=0).first()


# ---------------------------------------------------------------- IFN2 -> IFN3

def arbres_ifn2(dir_ifn: Path) -> pd.DataFrame:
    a = llegir_dbf(dir_ifn / "ifn2_25" / "PIESMA25.DBF")
    a["especie"] = a["ESPECIE"].str.strip()
    d1 = pd.to_numeric(a["DIAMETRO1"].str.strip(), errors="coerce")
    d2 = pd.to_numeric(a["DIAMETRO2"].str.strip(), errors="coerce")
    a["dbh"] = (d1 + d2) / 20.0
    a["ht"] = pd.to_numeric(a["ALTURA"].str.strip().str.replace(",", "."), errors="coerce")
    return amb_factor(a.rename(columns={"ESTADILLO": "parcela", "NUMORDEN": "ordre"}))


def periode_ifn23(dir_ifn: Path) -> tuple[pd.DataFrame, dict]:
    ifn2 = competencia(arbres_ifn2(dir_ifn))

    parceles2 = llegir_dbf(dir_ifn / "ifn2_25" / "DATEST25.DBF")
    parceles2 = parceles2.rename(columns={"ESTADILLO": "parcela"})
    parceles2["any_ifn2"] = 1900 + pd.to_numeric(parceles2["ANO"], errors="coerce")
    parceles2["altitud"] = pd.to_numeric(parceles2["ALTITUD1"], errors="coerce") * 100.0
    parceles2["utm_x"] = pd.to_numeric(parceles2["COORDEX"], errors="coerce") * 1000.0
    parceles2["utm_y"] = pd.to_numeric(parceles2["COORDEY"], errors="coerce") * 1000.0

    mdb = base_ifn3(dir_ifn)
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
    remesurats["agent"] = pd.to_numeric(remesurats["Agente"], errors="coerce").where(remesurats["mort"])
    remesurats = remesurats[["parcela", "ordre", "mort", "dbh_ifn3", "especie_ifn3", "agent"]]

    dades = ifn2[ifn2["especie"].isin(ESPECIES)].merge(remesurats, on=["parcela", "ordre"], how="inner")
    dades = dades[~dades["parcela"].isin(parceles_tallades)]
    dades = dades[dades["mort"] | (dades["especie_ifn3"] == dades["especie"])]
    dades = dades.merge(parceles2[["parcela", "any_ifn2", "altitud", "utm_x", "utm_y"]], on="parcela")
    dades = dades.merge(parceles3, on="parcela")
    recompte = {
        "parceles_remesurades": set(ifn3["parcela"]),
        "parceles_tallades": parceles_tallades,
        "arbres_tallats": ifn3.loc[ifn3["OrdenIf3"] == "000", "parcela"],
        "arbres_morts": ifn3.loc[ifn3["OrdenIf3"] == "888", "parcela"],
        "arbres_no_trobats": ifn3.loc[ifn3["OrdenIf3"] == "999", "parcela"],
    }
    return dades.assign(periode="IFN2-IFN3"), recompte


# ---------------------------------------------------------------- IFN3 -> IFN4

def arbres_ifn3(dir_ifn: Path) -> pd.DataFrame:
    """Peus vius de l'IFN3 (totes les espècies, per a la competència) a les parcel·les de classe A."""
    a = llegir_mdb(base_ifn3(dir_ifn), "PCMayores")
    a = a[(a["Cla"] == "A") & ordre_valid(a["OrdenIf3"])].copy()
    a = a[a["Calidad"].fillna("") != "6"]  # Calidad 6 = peu mort dempeus
    a["parcela"] = a["Estadillo"].astype(int)
    a["ordre"] = a["OrdenIf3"].astype(int)
    a["especie"] = a["Especie"].str.lstrip("0")
    a["dbh"] = (a["Dn1"].astype(float) + a["Dn2"].astype(float)) / 20.0
    a["ht"] = pd.to_numeric(a["Ht"], errors="coerce")
    return amb_factor(a)


def periode_ifn34(dir_ifn: Path) -> tuple[pd.DataFrame, dict]:
    """Codis de l'IFN4 (PCMayores): OrdenIf3 = ordre del peu a l'IFN3 (000 = peu nou; 888/999 = peu ja mort o
    no trobat a l'IFN3); OrdenIf4 = ordre a l'IFN4 (888 = mort; 000 = peu desaparegut, sense mesures, que es
    tracta com a tallat; 999 = no localitzat)."""
    ifn3 = competencia(arbres_ifn3(dir_ifn))
    mdb = base_ifn3(dir_ifn)
    p3 = llegir_mdb(mdb, "PCParcelas")
    p3 = p3[p3["Cla"] == "A"]
    p3 = pd.DataFrame({"parcela": p3["Estadillo"].astype(int), "any_ifn2": p3["Ano"].astype(int)})
    p3 = p3.drop_duplicates("parcela")
    m3 = llegir_mdb(mdb, "PCDatosMap")
    m3 = pd.DataFrame({"parcela": m3["Estadillo"].astype(int),
                       "utm_x": pd.to_numeric(m3["CoorX"], errors="coerce"),
                       "utm_y": pd.to_numeric(m3["CoorY"], errors="coerce")}).drop_duplicates("parcela")
    # Altitud: la mateixa font que a l'IFN2-IFN3 (la base SIG de l'IFN3 només dona classes de 200 m). Les
    # parcel·les noves de l'IFN3 sense altitud a l'IFN2 es descarten.
    m3["altitud"] = m3["parcela"].map(altitud_ifn2(dir_ifn))
    p3 = p3.merge(m3, on="parcela", how="left")

    ifn4 = dir_ifn / "ifn4_25" / "Ifn4_Cataluna.accdb"
    p4 = llegir_mdb(ifn4, "PCParcelas")
    p4 = p4[(p4["Provincia"] == "25") & (p4["Cla"] == "A") & (p4["Subclase"] == "1")]
    p4 = pd.DataFrame({"parcela": p4["Estadillo"].astype(int), "any_ifn3": p4["Ano"].astype(int)})

    a = llegir_mdb(ifn4, "PCMayores")
    a = a[(a["Provincia"] == "25") & (a["Cla"] == "A") & (a["Subclase"] == "1")].copy()
    a["parcela"] = a["Estadillo"].astype(int)
    tallades = set(a.loc[ordre_valid(a["OrdenIf3"]) & (a["OrdenIf4"] == "000"), "parcela"])
    r = a[ordre_valid(a["OrdenIf3"]) & (a["OrdenIf4"] != "999") & (a["OrdenIf4"] != "000")].copy()
    r["ordre"] = r["OrdenIf3"].astype(int)
    r["mort"] = (r["OrdenIf4"] == "888") | (r["Calidad"].fillna("") == "6")
    r["dbh_ifn3"] = (pd.to_numeric(r["Dn1"], errors="coerce") + pd.to_numeric(r["Dn2"], errors="coerce")) / 20.0
    r.loc[r["mort"], "dbh_ifn3"] = np.nan  # com a l'IFN2-IFN3: els morts no tenen diàmetre final
    r["especie_ifn3"] = r["Especie"].str.lstrip("0")
    r["agent"] = pd.to_numeric(r["Agente"], errors="coerce").where(r["mort"])
    r = r[["parcela", "ordre", "mort", "dbh_ifn3", "especie_ifn3", "agent"]]

    dades = ifn3[ifn3["especie"].isin(ESPECIES)].merge(r, on=["parcela", "ordre"], how="inner")
    dades = dades[~dades["parcela"].isin(tallades)]
    dades = dades[dades["mort"] | (dades["especie_ifn3"] == dades["especie"])]
    dades = dades.merge(p3[["parcela", "any_ifn2", "altitud", "utm_x", "utm_y"]], on="parcela")
    dades = dades.merge(p4, on="parcela")
    dades = dades[dades["altitud"].notna()]
    codis = a[ordre_valid(a["OrdenIf3"])]
    recompte = {
        "parceles_remesurades": set(p4["parcela"]) & set(p3["parcela"]),
        "parceles_tallades": tallades,
        "arbres_tallats": codis.loc[codis["OrdenIf4"] == "000", "parcela"],
        "arbres_morts": codis.loc[(codis["OrdenIf4"] == "888") | (codis["Calidad"].fillna("") == "6"), "parcela"],
        "arbres_no_trobats": codis.loc[codis["OrdenIf4"] == "999", "parcela"],
    }
    return dades.assign(periode="IFN3-IFN4"), recompte


# ---------------------------------------------------------------- Conjunt de calibratge

def comarques() -> pd.Series:
    return pd.read_csv(COMARCA).set_index("parcela")["comarca"]


def construir(dir_ifn: str | Path) -> tuple[pd.DataFrame, dict]:
    """Dades de calibratge de les dues remesures i, per període, el recompte del filtratge dins les comarques."""
    dir_ifn = Path(dir_ifn)
    comarca = comarques()
    dins = set(comarca[comarca.isin(COMARQUES)].index)
    trossos, recomptes = [], {}
    for periode, funcio in [("IFN2-IFN3", periode_ifn23), ("IFN3-IFN4", periode_ifn34)]:
        dades, recompte = funcio(dir_ifn)
        dades["comarca"] = dades["parcela"].map(comarca)
        dades = dades[dades["parcela"].isin(dins)]
        anys_min, anys_max = INTERVALS[periode]
        dades["anys"] = dades["any_ifn3"] - dades["any_ifn2"]
        dades = dades[dades["anys"].between(anys_min, anys_max)]
        dades["id_dbh"] = np.where(dades["mort"], np.nan, (dades["dbh_ifn3"] - dades["dbh"]) / dades["anys"])
        # Errors de mesura evidents: decreixements o creixements de més de 2 cm/any.
        dades = dades[dades["mort"] | dades["id_dbh"].between(-0.2, 2.0)]
        dades["especie"] = dades["especie"].map(ESPECIES)
        dades["causa"] = causa_mort(dades["agent"], dades["mort"])
        trossos.append(dades)
        remesurades = recompte["parceles_remesurades"] & dins
        recomptes[periode] = {
            "parceles_remesurades": len(remesurades),
            "parceles_tallades": len(recompte["parceles_tallades"] & remesurades),
            "arbres_tallats": int(recompte["arbres_tallats"].isin(dins).sum()),
            "arbres_morts": int(recompte["arbres_morts"].isin(dins).sum()),
            "arbres_no_trobats": int(recompte["arbres_no_trobats"].isin(dins).sum()),
            "interval_admes": [anys_min, anys_max],
        }
    dades = pd.concat(trossos, ignore_index=True)
    dades = dades.sort_values(["periode", "parcela", "ordre"])[COLUMNES].round(4)
    dades["agent"] = dades["agent"].astype("Int64")
    return dades, recomptes


def main(dir_ifn: str) -> None:
    dades, recomptes = construir(dir_ifn)
    dades.to_csv(DESTI, index=False)
    resum = dades.groupby(["periode", "especie"]).agg(
        arbres=("ordre", "size"), parceles=("parcela", "nunique"), morts=("mort", "sum"),
        anys_min=("anys", "min"), anys_max=("anys", "max"))
    print(resum.to_string())
    print(dades.groupby(["comarca", "periode"])["parcela"].nunique().unstack().to_string())
    morts = dades[dades["mort"]]
    print("Morts per causa:")
    print(morts.groupby(["periode", "especie", "causa"]).size().unstack(fill_value=0).to_string())
    print("Morts sense agent:", morts.groupby(["periode", "especie"])["agent"].apply(lambda a: int(a.isna().sum()))
          .to_dict())
    print("Parcel·les a tots dos períodes:", int(dades.groupby("parcela")["periode"].nunique().eq(2).sum()))
    for periode, r in recomptes.items():
        print(periode, r)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])
