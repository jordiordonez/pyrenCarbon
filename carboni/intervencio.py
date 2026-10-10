"""Aclarides aplicades a una llista d'arbres."""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pandas as pd

from carboni.equacions import ESPECIES
from carboni.inventari import Rodal

TIPUS = {
    "baix": "Per baix (es tallen primer els peus més prims)",
    "dalt": "Per dalt (es tallen primer els peus més gruixuts)",
    "sistematica": "Sistemàtica (la mateixa proporció a totes les mides)",
}


def area_basimetrica(arbres: pd.DataFrame) -> float:
    return float((np.pi * (arbres["dbh"] / 200.0) ** 2 * arbres["f"]).sum())


def aclarir(rodal: Rodal, tipus: str, fraccio_ab: float) -> tuple[Rodal, pd.DataFrame]:
    """Extreu `fraccio_ab` de l'àrea basimètrica de pins. Retorna el rodal resultant i els arbres extrets."""
    if tipus not in TIPUS:
        raise ValueError(f"Tipus d'aclarida desconegut: {tipus}")
    if not 0.0 <= fraccio_ab < 1.0:
        raise ValueError("La fracció d'àrea basimètrica a extreure ha d'estar entre 0 i 1.")
    arbres = rodal.arbres.copy().reset_index(drop=True)
    pins = arbres["sp"].isin(ESPECIES).to_numpy()
    g = np.pi * (arbres["dbh"].to_numpy() / 200.0) ** 2 * arbres["f"].to_numpy()
    objectiu = fraccio_ab * g[pins].sum()
    f = arbres["f"].to_numpy()
    extret = np.zeros(len(arbres))

    if tipus == "sistematica":
        extret[pins] = f[pins] * fraccio_ab
    else:
        candidats = np.where(pins)[0]
        ordre = np.argsort(arbres["dbh"].to_numpy()[candidats], kind="stable")
        if tipus == "dalt":
            ordre = ordre[::-1]
        pendent = objectiu
        for i in candidats[ordre]:
            if pendent <= 0:
                break
            g_peu = g[i] / f[i]
            peus = min(f[i], pendent / g_peu)
            extret[i] = peus
            pendent -= peus * g_peu

    restants = arbres.assign(f=f - extret)
    restants = restants[restants["f"] > 1e-9].reset_index(drop=True)
    extrets = arbres.assign(f=extret)
    extrets = extrets[extrets["f"] > 1e-9].reset_index(drop=True)
    return replace(rodal, arbres=restants), extrets
