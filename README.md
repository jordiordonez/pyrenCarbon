# pyrenCarbon

**Carboni a les pinedes d'Andorra.** Aplicació Streamlit per estimar l'estoc i la captació de CO₂ de les pinedes de **pi roig** (*Pinus sylvestris*) i **pi negre** (*Pinus uncinata*) d'Andorra, i l'efecte d'una aclarida, en un horitzó de fins a 80 anys.

- **Punt de partida**: qualsevol pineda del 1r Inventari Nacional Forestal d'Andorra, unes dades de rodal (densitat i àrea basimètrica) o un inventari propi en CSV.
- **Creixement arbre a arbre**: creixement diametral, mortalitat i autoaclarida calibrats amb les remesures IFN2→IFN3 i IFN3→IFN4 de l'Inventari Forestal Nacional a les comarques pirinenques de Lleida (Pallars Sobirà, Alt Urgell, Alta Ribagorça i Cerdanya): 19.603 registres de pins (11.546 de pi roig i 8.057 de pi negre) en 784 parcel·les, 243 de les quals es fan servir en tots dos períodes.
- **Mortalitat de fons i pertorbacions**: el model de mortalitat s'ajusta només amb els pins morts per competència o causa desconeguda (agent causant de l'IFN); el foc, la sequera, la neu i el vent, les plagues i altres pertorbacions es tracten a part, amb tres escenaris: sense pertorbacions, risc mitjà (per defecte) o un episodi puntual.
- **Biomassa i carboni** amb les mateixes equacions que l'informe de carboni de l'INF1-AND (alçades pròpies d'Andorra, IEFC, INIA).
- **Validació**: validació creuada per parcel·les i contrast amb l'informe de l'INF1-AND. Els detalls són a [docs/metode.md](docs/metode.md) i a les pestanyes *Model i validació* i *Com funciona* de l'app.

## Execució local

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
streamlit run app.py
```

## Reproduir el calibratge

```bash
pip install -r requirements-dev.txt
python scripts/descarregar_inf1and.py          # BD-INF1-AND des de Zenodo
python scripts/assignar_comarca.py DIR_IFN      # (opcional) refà data/ifn_parceles_comarca.csv
python scripts/preparar_ifn.py DIR_IFN          # IFN2, IFN3 i IFN4 de Lleida (MITECO), cal mdbtools
python scripts/calibrar.py                      # ajust, validació creuada i data/calibracio.json
pytest
```

`DIR_IFN` ha de contenir, descomprimits, els fitxers de parcel·les de l'IFN2 de Lleida (`ifn2_25/PIESMA25.DBF`, `ifn2_25/DATEST25.DBF`), la base de dades de camp de l'IFN3 de Lleida (`ifn3p25/Ifn3p25.accdb`) i la de l'IFN4 de Catalunya (`ifn4_25/Ifn4_Cataluna.accdb`), descarregables del web del MITECO. `scripts/assignar_comarca.py` necessita a més la capa de comarques i municipis de les divisions administratives de l'ICGC (`icgc_divisions_administratives/`) i les dependències geogràfiques de `requirements-dev.txt` (shapely, pyshp, pyproj).

El calibratge fa servir les dues remesures juntes amb un terme de període (`p34`); el model de creixement té efectes aleatoris niats de parcel·la i de parcel·la×període. `data/calibracio.json` desa els coeficients del període més recent (IFN3→IFN4), que són els que fa servir l'app, i a part el coeficient del període i el seu error estàndard, els components de variància i la prova de raó de versemblances de l'efecte parcel·la×període. `scripts/preparar_ifn.py` hi afegeix l'agent causant de cada pi mort i el grup de causa (constant `CAUSES`); `scripts/calibrar.py` ajusta la mortalitat només amb les morts de fons i desa les pertorbacions al bloc `pertorbacions`.

## Estructura

```
app.py                       interfície Streamlit
carboni/equacions.py         alçada, biomassa i carboni per arbre
carboni/inventari.py         rodals de partida (INF-AND, rodal anàleg, CSV)
carboni/creixement.py        simulador: creixement, mortalitat de fons, pertorbacions i autoaclarida
carboni/intervencio.py       aclarides
data/inf1and/                BD-INF1-AND (parcel·les i arbres)
data/ifn_lleida_arbres.csv   conjunt de calibratge derivat de l'IFN2, l'IFN3 i l'IFN4 de Lleida (dos períodes)
data/ifn_parceles_comarca.csv  comarca de cada parcel·la de l'IFN de Lleida (node de malla de l'IFN3)
data/calibracio.json         coeficients, pertorbacions i mètriques de validació
scripts/                     descàrrega, preparació i calibratge
docs/metode.md               metodologia
```

## Fonts i llicències

- **BD-INF1-AND**: Beltran Barba, M., Busquets Olivé, E., Piqué Nicolau, M., Domènech Ferrés, M. i Travesset-Baro, O. (2023). *1r Inventari Nacional Forestal (INF1) del Principat d'Andorra. Base de dades.* Andorra Recerca + Innovació i CTFC. [doi:10.5281/zenodo.11388445](https://doi.org/10.5281/zenodo.11388445). CC BY-SA 4.0.
- **Equacions d'alçada, biomassa i carboni**: Beltran Barba, M. et al. (2023). *1r INF del Principat d'Andorra. Estimació de l'estoc i l'embornal de carboni total als boscos d'Andorra.* [doi:10.5281/zenodo.12519202](https://doi.org/10.5281/zenodo.12519202). CC BY-SA 4.0.
- **Segon, Tercer i Quart Inventari Forestal Nacional**, província de Lleida (IFN4: fitxer de Catalunya). Ministerio para la Transición Ecológica y el Reto Demográfico, Banco de Datos de la Naturaleza i [Cuarto Inventario Forestal Nacional](https://www.miteco.gob.es/es/biodiversidad/temas/inventarios-nacionales/inventario-forestal-nacional/cuarto_inventario.html).
- **Límits comarcals**: Divisions administratives v2r2, Institut Cartogràfic i Geològic de Catalunya (ICGC). CC BY 4.0.

## Llicència

El codi es distribueix amb llicència [MIT](LICENSE). Les dades de `data/` mantenen la llicència de la seva font: les derivades de la BD-INF1-AND, CC BY-SA 4.0; les derivades de l'IFN, reutilització amb citació del MITECO; la taula de comarques, CC BY 4.0 (ICGC). El detall és a [data/LLICENCIES.md](data/LLICENCIES.md).
