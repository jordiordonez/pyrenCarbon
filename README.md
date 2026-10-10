# pyrenCarbon

**Carboni a les pinedes d'Andorra.** Aplicació Streamlit per estimar l'estoc i la captació de CO₂ de les pinedes de **pi roig** (*Pinus sylvestris*) i **pi negre** (*Pinus uncinata*) d'Andorra, i l'efecte d'una aclarida, en un horitzó de fins a 80 anys.

- **Punt de partida**: qualsevol pineda del 1r Inventari Nacional Forestal d'Andorra, unes dades de rodal (densitat i àrea basimètrica) o un inventari propi en CSV.
- **Creixement arbre a arbre**: creixement diametral, mortalitat i autoaclarida calibrats amb 15.536 pins remesurats de l'Inventari Forestal Nacional a Lleida.
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
python scripts/preparar_ifn.py DIR_IFN          # IFN2 i IFN3 de Lleida (MITECO), cal mdbtools
python scripts/calibrar.py                      # ajust, validació creuada i data/calibracio.json
pytest
```

`DIR_IFN` ha de contenir, descomprimits, el fitxer de parcel·les de l'IFN2 de Lleida (`ifn2_25/PIESMA25.DBF`, `ifn2_25/DATEST25.DBF`) i les bases de dades de camp i SIG de l'IFN3 de Lleida (`ifn3p25/Ifn3p25.mdb`, `sig25/Sig_25.mdb`), descarregables del Banc de Dades de la Naturalesa del MITECO.

## Estructura

```
app.py                       interfície Streamlit
carboni/equacions.py         alçada, biomassa i carboni per arbre
carboni/inventari.py         rodals de partida (INF-AND, rodal anàleg, CSV)
carboni/creixement.py        simulador: creixement, mortalitat i autoaclarida
carboni/intervencio.py       aclarides
data/inf1and/                BD-INF1-AND (parcel·les i arbres)
data/ifn_lleida_arbres.csv   conjunt de calibratge derivat de l'IFN2 i l'IFN3 de Lleida
data/calibracio.json         coeficients i mètriques de validació
scripts/                     descàrrega, preparació i calibratge
docs/metode.md               metodologia
```

## Fonts i llicències

- **BD-INF1-AND**: Beltran Barba, M., Busquets Olivé, E., Piqué Nicolau, M., Domènech Ferrés, M. i Travesset-Baro, O. (2023). *1r Inventari Nacional Forestal (INF1) del Principat d'Andorra. Base de dades.* Andorra Recerca + Innovació i CTFC. [doi:10.5281/zenodo.11388445](https://doi.org/10.5281/zenodo.11388445). CC BY-SA 4.0.
- **Equacions d'alçada, biomassa i carboni**: Beltran Barba, M. et al. (2023). *1r INF del Principat d'Andorra. Estimació de l'estoc i l'embornal de carboni total als boscos d'Andorra.* [doi:10.5281/zenodo.12519202](https://doi.org/10.5281/zenodo.12519202). CC BY-SA 4.0.
- **Segon i Tercer Inventari Forestal Nacional**, província de Lleida. Ministerio para la Transición Ecológica y el Reto Demográfico, Banco de Datos de la Naturaleza.

## Llicència

El codi es distribueix amb llicència [MIT](LICENSE). Les dades de `data/` mantenen la llicència de la seva font: les derivades de la BD-INF1-AND, CC BY-SA 4.0; les derivades de l'IFN, reutilització amb citació del MITECO. El detall és a [data/LLICENCIES.md](data/LLICENCIES.md).
