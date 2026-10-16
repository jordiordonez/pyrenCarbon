"""Assigna la comarca a cada parcel·la de l'IFN de Lleida a partir del node de malla de l'IFN3.

Ús:
    python scripts/assignar_comarca.py DIR_IFN

DIR_IFN ha de contenir:
    ifn3p25/Ifn3p25.accdb (o .mdb)          IFN3, Lleida, camp (taules PCDatosMap i PCParcelas)
    icgc_divisions_administratives/          Divisions administratives v2r2 de l'ICGC (comarques i
                                             municipis 1:5.000, shapefile)
Cal tenir `mdb-export` (mdbtools) al PATH i shapely, pyshp i pyproj (requirements-dev.txt).

Mètode: les coordenades del node de malla de l'IFN3 (PCDatosMap, CoorX/CoorY, ED50 / UTM 31N,
EPSG:23031) es transformen a ETRS89 / UTM 31N (EPSG:25831) i s'intersequen amb els polígons de les
comarques de l'ICGC. Si un punt queda fora de Catalunya, es fa servir la comarca actual del municipi
del codi INE de l'IFN3. Les parcel·les de l'IFN4 mantenen el número d'estadillo de l'IFN3, de manera
que la mateixa taula serveix per a les dues remesures.

Escriu data/ifn_parceles_comarca.csv: una fila per parcel·la de classe A de l'IFN3 de Lleida.
"""

import sys
from pathlib import Path

import pandas as pd
import shapefile
from pyproj import Transformer
from shapely.geometry import Point, shape
from shapely.strtree import STRtree

sys.path.insert(0, str(Path(__file__).resolve().parent))
from preparar_ifn import base_ifn3, llegir_mdb  # noqa: E402

DESTI = Path(__file__).resolve().parent.parent / "data" / "ifn_parceles_comarca.csv"
PER_COORDENADES = "coordenades IFN3 x límits comarcals ICGC"


def capa(dir_icgc: Path, nivell: str) -> shapefile.Reader:
    ruta = sorted(dir_icgc.glob(f"divisions-administratives-v2r2-{nivell}-5000-*.shp"))[-1]
    return shapefile.Reader(str(ruta), encoding="utf-8")


def main(dir_ifn: str) -> None:
    dir_ifn = Path(dir_ifn)
    mdb = base_ifn3(dir_ifn)
    mapa = llegir_mdb(mdb, "PCDatosMap")
    mapa = mapa[mapa["Clase"] == "A"].copy()
    mapa["parcela"] = mapa["Estadillo"].astype(int)
    mapa = mapa.drop_duplicates("parcela").sort_values("parcela")
    pc = llegir_mdb(mdb, "PCParcelas")
    ine = pc[pc["Cla"] == "A"].assign(parcela=lambda d: d["Estadillo"].astype(int)).drop_duplicates("parcela")
    ine = ine.set_index("parcela")["INE"]

    dir_icgc = dir_ifn / "icgc_divisions_administratives"
    rc = capa(dir_icgc, "comarques")
    comarques = [(r["NOMCOMAR"], shape(s.__geo_interface__)) for r, s in zip(rc.records(), rc.shapes())]
    arbre = STRtree([g for _, g in comarques])
    comarca_municipi = {r["CODIMUNI"][:5]: r["NOMCOMAR"] for r in capa(dir_icgc, "municipis").records()}

    transf = Transformer.from_crs("EPSG:23031", "EPSG:25831", always_xy=True)
    files = []
    for fila in mapa.itertuples():
        x, y = float(fila.CoorX), float(fila.CoorY)
        punt = Point(*transf.transform(x, y))
        dins = [i for i in arbre.query(punt) if comarques[i][1].covers(punt)]
        if dins:
            nom, poligon = comarques[dins[0]]
            files.append({"parcela": fila.parcela, "comarca": nom, "x_ed50": x, "y_ed50": y,
                          "dist_limit_m": round(poligon.boundary.distance(punt)),
                          "metode": PER_COORDENADES})
        else:
            codi = ine.get(fila.parcela, fila.INE)
            files.append({"parcela": fila.parcela, "comarca": comarca_municipi.get(f"25{codi}"), "x_ed50": x,
                          "y_ed50": y, "dist_limit_m": None,
                          "metode": "codi INE del municipi (coordenada fora de Catalunya)"})
    taula = pd.DataFrame(files)
    taula.to_csv(DESTI, index=False)
    print(f"{len(taula)} parcel·les; sense comarca: {int(taula['comarca'].isna().sum())}; "
          f"per codi INE: {int((taula['metode'] != PER_COORDENADES).sum())}")
    print(taula["comarca"].value_counts().to_string())


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])
