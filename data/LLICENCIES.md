# Llicències de les dades

La llicència MIT del repositori cobreix només el codi. Les dades d'aquesta carpeta mantenen la llicència de la seva font.

| Fitxers | Origen | Llicència |
|---|---|---|
| `inf1and/parceles.csv`, `inf1and/arbres.csv` | BD-INF1-AND, Andorra Recerca + Innovació i CTFC, [doi:10.5281/zenodo.11388445](https://doi.org/10.5281/zenodo.11388445) | CC BY-SA 4.0 |
| `calibracio.json` | Coeficients ajustats amb l'IFN de Lleida; el límit d'autoaclarida usa també la BD-INF1-AND | CC BY-SA 4.0 |
| `ifn_lleida_arbres.csv`, `validacio_arbres.csv`, `validacio_parceles.csv` | Derivats del Segon, Tercer i Quart Inventari Forestal Nacional, província de Lleida (MITECO). IFN2 i IFN3: Banco de Datos de la Naturaleza; IFN4: base de dades de camp de Catalunya (`ifn4_cataluna`), [Cuarto Inventario Forestal Nacional](https://www.miteco.gob.es/es/biodiversidad/temas/inventarios-nacionales/inventario-forestal-nacional/cuarto_inventario.html) | Reutilització amb citació de la font (MITECO) |
| `ifn_parceles_comarca.csv` | Comarca de cada parcel·la: coordenades del node de malla de l'IFN3 (MITECO) creuades amb els límits comarcals de les Divisions administratives v2r2 de l'Institut Cartogràfic i Geològic de Catalunya (ICGC, [icgc.cat](https://www.icgc.cat)) | Coordenades: citació del MITECO; límits comarcals: CC BY 4.0 (ICGC) |

Els coeficients de les equacions d'alçada, biomassa i carboni de `carboni/equacions.py` provenen de l'informe de carboni de l'INF1-AND ([doi:10.5281/zenodo.12519202](https://doi.org/10.5281/zenodo.12519202), CC BY-SA 4.0), que s'hi cita.
