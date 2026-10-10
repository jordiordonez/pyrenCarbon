"""Equacions d'arbre individual per a pi roig (Ps) i pi negre (Pu).

Totes provenen de l'informe de carboni del 1r Inventari Nacional Forestal d'Andorra
(Beltran et al. 2023, https://doi.org/10.5281/zenodo.12519202):
- alçada–diàmetre pròpies de l'INF-AND, per hàbitat xeròfil (x) o mesòfil (m) (taula 5)
- biomassa aèria total, IEFC Catalunya (taula 6)
- biomassa subterrània total, INIA Espanya (taula 7)
- fracció de carboni de 0,5 t C per tona de biomassa
"""

from __future__ import annotations

import numpy as np
import pandas as pd

ESPECIES = {"Ps": "Pinus sylvestris", "Pu": "Pinus uncinata"}
NOMS_COMUNS = {"Ps": "Pi roig", "Pu": "Pi negre"}

ALCADA_DIAMETRE = {
    ("Ps", "x"): (2.2485, 0.5067),
    ("Ps", "m"): (2.9159, 0.4898),
    ("Pu", "x"): (4.3740, 0.2806),
    ("Pu", "m"): (3.2135, 0.3974),
}
BIOMASSA_AERIA = {
    "Ps": (0.05539335, 2.04694079, 0.55768901),
    "Pu": (0.05370209, 2.03016707, 0.59328472),
}
BIOMASSA_SUBTERRANIA = {
    "Ps": (0.01088662, 2.62841),
    "Pu": (0.01147806, 2.7558),
}
FRACCIO_CARBONI = 0.5
CO2_PER_C = 44.0 / 12.0
HABITATS = {"x": "xeròfil", "m": "mesòfil"}


def habitat_de_codi(codi: str | None) -> str:
    """'x' o 'm' a partir del codi d'hàbitat de l'INF-AND (PRX, PNM...); mesòfil per defecte."""
    if isinstance(codi, str) and codi.upper().endswith("X"):
        return "x"
    return "m"


def alcada(sp: pd.Series, dbh: pd.Series, habitat: str) -> np.ndarray:
    a = sp.map(lambda s: ALCADA_DIAMETRE[(s, habitat)][0]).to_numpy(dtype=float)
    b = sp.map(lambda s: ALCADA_DIAMETRE[(s, habitat)][1]).to_numpy(dtype=float)
    return a * np.power(dbh.to_numpy(dtype=float), b)


def biomassa_aeria_kg(sp: pd.Series, dbh: pd.Series, ht: pd.Series) -> np.ndarray:
    p = np.array([BIOMASSA_AERIA[s] for s in sp]).reshape(-1, 3)
    return p[:, 0] * np.power(dbh.to_numpy(dtype=float), p[:, 1]) * np.power(ht.to_numpy(dtype=float), p[:, 2])


def biomassa_subterrania_kg(sp: pd.Series, dbh: pd.Series) -> np.ndarray:
    p = np.array([BIOMASSA_SUBTERRANIA[s] for s in sp]).reshape(-1, 2)
    return p[:, 0] * np.power(dbh.to_numpy(dtype=float), p[:, 1])
