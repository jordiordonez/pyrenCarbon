"""Descarrega la base de dades del 1r Inventari Nacional Forestal d'Andorra (BD-INF1-AND).

Font: Andorra Recerca + Innovació i CTFC, https://doi.org/10.5281/zenodo.11388445 (CC BY-SA 4.0).
"""

import urllib.request
from pathlib import Path

REGISTRE = "https://zenodo.org/records/11388445/files"
FITXERS = ("parceles.csv", "arbres.csv")
DESTI = Path(__file__).resolve().parent.parent / "data" / "inf1and"


def main() -> None:
    DESTI.mkdir(parents=True, exist_ok=True)
    for nom in FITXERS:
        urllib.request.urlretrieve(f"{REGISTRE}/{nom}?download=1", DESTI / nom)
        print(f"{nom}: {(DESTI / nom).stat().st_size} bytes")


if __name__ == "__main__":
    main()
