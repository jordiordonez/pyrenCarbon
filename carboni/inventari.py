"""Llistes d'arbres de partida: parcel·les de l'INF-AND, parcel·les anàlogues o inventaris propis.

Una llista d'arbres és un DataFrame amb una fila per arbre:
    sp      codi d'espècie de l'INF-AND (Ps, Pu, Aa...)
    dbh     diàmetre normal (cm)
    f       peus per hectàrea que representa l'arbre
    ht_obs  alçada mesurada (m), opcional
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from carboni.equacions import ESPECIES, alcada, habitat_de_codi

ASPECTES_VIUS = (1, 2)
FRACCIO_PINAR = 0.7


@dataclass(frozen=True)
class Rodal:
    arbres: pd.DataFrame
    habitat: str
    altitud: float
    descripcio: str


def carregar_inf1and(directori: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    parceles = pd.read_csv(directori / "parceles.csv", encoding="utf-8-sig")
    arbres = pd.read_csv(directori / "arbres.csv", encoding="utf-8-sig")
    arbres = arbres[arbres["asp"].isin(ASPECTES_VIUS)].copy()
    radi = arbres["parcela"].map(parceles.set_index("parcela")["radi"])
    arbres["f"] = 1e4 / (np.pi * radi**2)
    return parceles, arbres


def resum_parceles(parceles: pd.DataFrame, arbres: pd.DataFrame) -> pd.DataFrame:
    a = arbres.assign(g=np.pi * (arbres["dbh"] / 200.0) ** 2 * arbres["f"])
    totals = a.groupby("parcela").agg(n=("f", "sum"), ab=("g", "sum"))
    per_especie = a.groupby(["parcela", "sp"])["g"].sum().unstack(fill_value=0.0)
    fraccions = per_especie.div(per_especie.sum(axis=1), axis=0)
    totals["especie_dominant"] = fraccions.idxmax(axis=1)
    totals["fraccio_dominant"] = fraccions.max(axis=1)
    totals["dq"] = np.sqrt(4.0 * totals["ab"] / (np.pi * totals["n"])) * 100.0
    info = parceles.set_index("parcela")[["habitat", "estrat", "lat_4258", "lon_4258", "elev_el", "radi"]]
    resum = info.join(totals, how="inner").reset_index()
    resum["habitat_hd"] = resum["habitat"].map(habitat_de_codi)
    resum = resum.rename(columns={"lat_4258": "lat", "lon_4258": "lon", "elev_el": "altitud"})
    resum["altitud_estimada"] = resum["altitud"].isna()
    mediana_habitat = resum.groupby("habitat")["altitud"].transform("median")
    resum["altitud"] = resum["altitud"].fillna(mediana_habitat).fillna(resum["altitud"].median())
    return resum


def pinars(resum: pd.DataFrame, sp: str | None = None) -> pd.DataFrame:
    seleccio = resum[resum["especie_dominant"].isin(ESPECIES) & (resum["fraccio_dominant"] >= FRACCIO_PINAR)]
    if sp is not None:
        seleccio = seleccio[seleccio["especie_dominant"] == sp]
    return seleccio


def rodal_de_parcela(resum: pd.DataFrame, arbres: pd.DataFrame, codi: str) -> Rodal:
    info = resum.set_index("parcela").loc[codi]
    llista = arbres.loc[arbres["parcela"] == codi, ["sp", "dbh", "f", "hm"]].rename(columns={"hm": "ht_obs"})
    altitud = f"{info['altitud']:.0f} m" + (" (estimada)" if info["altitud_estimada"] else "")
    return Rodal(
        arbres=llista.reset_index(drop=True),
        habitat=info["habitat_hd"],
        altitud=float(info["altitud"]),
        descripcio=f"Parcel·la {codi} de l'INF-AND ({info['habitat']}, {altitud})",
    )


def rodal_analog(resum: pd.DataFrame, arbres: pd.DataFrame, sp: str, habitat: str, altitud: float,
                 n: float, ab: float) -> Rodal:
    """Pren la distribució diamètrica de la parcel·la andorrana més semblant i l'escala a (n, ab)."""
    candidates = pinars(resum, sp)
    mateix_habitat = candidates[candidates["habitat_hd"] == habitat]
    if len(mateix_habitat) >= 5:
        candidates = mateix_habitat
    distancia = np.hypot(np.log(candidates["n"] / n), np.log(candidates["ab"] / ab))
    codi = candidates.loc[distancia.idxmin(), "parcela"]
    llista = arbres.loc[(arbres["parcela"] == codi) & (arbres["sp"] == sp), ["sp", "dbh", "f"]].copy()
    n_ref = llista["f"].sum()
    ab_ref = (np.pi * (llista["dbh"] / 200.0) ** 2 * llista["f"]).sum()
    llista["f"] *= n / n_ref
    llista["dbh"] *= np.sqrt(ab / (ab_ref * n / n_ref))
    llista["ht_obs"] = np.nan
    return Rodal(
        arbres=llista.reset_index(drop=True),
        habitat=habitat,
        altitud=altitud,
        descripcio=f"Distribució diamètrica de la parcel·la {codi} de l'INF-AND, escalada a les dades introduïdes",
    )


CODIS_PER_NOM = {nom.lower(): codi for codi, nom in ESPECIES.items()} | {
    "pi roig": "Ps", "pi negre": "Pu", "ps": "Ps", "pu": "Pu",
}


def rodal_de_csv(taula: pd.DataFrame, habitat: str, altitud: float, superficie_m2: float | None) -> Rodal:
    """Inventari propi. Columnes: especie, dbh (cm), i peus_ha o bé una superfície de parcel·la; ht opcional."""
    columnes = {c.lower().strip(): c for c in taula.columns}
    if "especie" not in columnes or "dbh" not in columnes:
        raise ValueError("El fitxer ha de tenir com a mínim les columnes 'especie' i 'dbh'.")
    especie = taula[columnes["especie"]].astype(str).str.strip()
    sp = especie.str.lower().map(CODIS_PER_NOM).fillna(especie)
    llista = pd.DataFrame({"sp": sp, "dbh": pd.to_numeric(taula[columnes["dbh"]], errors="coerce")})
    if "peus_ha" in columnes:
        llista["f"] = pd.to_numeric(taula[columnes["peus_ha"]], errors="coerce")
    elif superficie_m2:
        llista["f"] = 1e4 / superficie_m2
    else:
        raise ValueError("Cal una columna 'peus_ha' o indicar la superfície de la parcel·la.")
    llista["ht_obs"] = pd.to_numeric(taula[columnes["ht"]], errors="coerce") if "ht" in columnes else np.nan
    llista = llista.dropna(subset=["dbh", "f"])
    llista = llista[llista["dbh"] >= 7.5]
    if not llista["sp"].isin(ESPECIES).any():
        raise ValueError("L'inventari no conté cap pi roig (Ps) ni pi negre (Pu).")
    return Rodal(arbres=llista.reset_index(drop=True), habitat=habitat, altitud=altitud,
                 descripcio="Inventari propi")


def factor_alcada_local(rodal: Rodal) -> float:
    """Ajust de la corba alçada–diàmetre a les alçades mesurades del rodal, si n'hi ha."""
    a = rodal.arbres
    mesurats = a[a["sp"].isin(ESPECIES) & a["ht_obs"].notna() & (a["ht_obs"] > 1.3)]
    if len(mesurats) < 3:
        return 1.0
    predita = alcada(mesurats["sp"], mesurats["dbh"], rodal.habitat)
    return float(np.clip(np.median(mesurats["ht_obs"] / predita), 0.7, 1.4))
