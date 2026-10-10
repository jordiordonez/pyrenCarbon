from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from carboni.creixement import carregar_models, simular, sdi
from carboni.equacions import CO2_PER_C, FRACCIO_CARBONI, alcada, biomassa_aeria_kg
from carboni.intervencio import aclarir, area_basimetrica
from carboni.inventari import (
    carregar_inf1and,
    pinars,
    resum_parceles,
    rodal_analog,
    rodal_de_csv,
    rodal_de_parcela,
)

DATA = Path(__file__).resolve().parent.parent / "data"


@pytest.fixture(scope="module")
def inf():
    parceles, arbres = carregar_inf1and(DATA / "inf1and")
    return resum_parceles(parceles, arbres), arbres


@pytest.fixture(scope="module")
def models():
    return carregar_models(DATA / "calibracio.json")


@pytest.fixture(scope="module")
def rodal_pu(inf):
    resum, arbres = inf
    codi = pinars(resum, "Pu").sort_values("parcela")["parcela"].iloc[0]
    return rodal_de_parcela(resum, arbres, codi)


def test_inventari_nomes_arbres_vius(inf):
    _, arbres = inf
    assert arbres["asp"].isin([1, 2]).all()


def test_hi_ha_pinars_de_les_dues_especies(inf):
    resum, _ = inf
    assert len(pinars(resum, "Pu")) >= 80
    assert len(pinars(resum, "Ps")) >= 60
    assert pinars(resum)["altitud"].notna().all()


def test_estoc_per_estrat_coherent_amb_informe_ari(inf, models):
    """Carboni aeri mitjà per estrat dins del ±10 % de la taula 14 de l'informe de l'INF1-AND."""
    resum, arbres = inf
    informe = {"PNM": 68.3, "PNX": 70.2, "PRM": 86.8, "PRX": 66.6}
    seleccio = pinars(resum)
    for habitat, valor in informe.items():
        codis = seleccio.loc[seleccio["habitat"] == habitat, "parcela"]
        cat = [simular(rodal_de_parcela(resum, arbres, c), models, 0).serie["cat"].iloc[0] for c in codis]
        assert np.mean(cat) == pytest.approx(valor, rel=0.10), habitat


def test_simulacio_creix_i_perd_peus(rodal_pu, models):
    serie = simular(rodal_pu, models, 30).serie
    assert len(serie) == 31
    assert serie["dq"].is_monotonic_increasing
    assert serie["densitat"].iloc[-1] < serie["densitat"].iloc[0]
    assert serie["carboni"].iloc[-1] > serie["carboni"].iloc[0]
    assert (serie["increment_brut"].iloc[1:] > 0).all()
    assert serie["co2"].to_numpy() == pytest.approx(serie["carboni"].to_numpy() * CO2_PER_C)


def test_autoaclarida_no_supera_el_maxim(rodal_pu, models):
    dens = replace(rodal_pu, arbres=rodal_pu.arbres.assign(f=rodal_pu.arbres["f"] * 3))
    serie = simular(dens, models, 40).serie
    assert serie["sdi"].iloc[1:].max() <= serie["sdi"].iloc[0] * 1.001


def test_aclarida_per_baix_treu_els_prims(rodal_pu):
    abans = area_basimetrica(rodal_pu.arbres)
    restant, extrets = aclarir(rodal_pu, "baix", 0.3)
    assert area_basimetrica(restant.arbres) == pytest.approx(abans * 0.7, rel=1e-6)
    assert extrets["dbh"].max() <= restant.arbres.loc[restant.arbres["sp"] == "Pu", "dbh"].min() + 1e-9 or len(
        extrets) == 0


def test_aclarida_sistematica_mante_el_dq(rodal_pu):
    restant, _ = aclarir(rodal_pu, "sistematica", 0.4)
    dq = lambda a: np.sqrt((a["dbh"] ** 2 * a["f"]).sum() / a["f"].sum())  # noqa: E731
    assert dq(restant.arbres) == pytest.approx(dq(rodal_pu.arbres))


def test_rodal_analog_reprodueix_n_i_ab(inf):
    resum, arbres = inf
    rodal = rodal_analog(resum, arbres, "Ps", "m", 1700.0, 900.0, 35.0)
    assert rodal.arbres["f"].sum() == pytest.approx(900.0)
    assert area_basimetrica(rodal.arbres) == pytest.approx(35.0)


def test_rodal_de_csv():
    taula = pd.DataFrame({"especie": ["Pinus uncinata", "pi negre", "Pu", "Aa"], "DBH": [12, 25.5, 40, 30]})
    rodal = rodal_de_csv(taula, "m", 2000.0, superficie_m2=400.0)
    assert list(rodal.arbres["sp"]) == ["Pu", "Pu", "Pu", "Aa"]
    assert rodal.arbres["f"].iloc[0] == pytest.approx(25.0)
    with pytest.raises(ValueError):
        rodal_de_csv(pd.DataFrame({"especie": ["Pu"], "dbh": [20]}), "m", 2000.0, superficie_m2=None)


def test_biomassa_arbre_individual():
    sp, dbh = pd.Series(["Ps"]), pd.Series([30.0])
    ht = pd.Series(alcada(sp, dbh, "m"))
    esperat = 0.05539335 * 30**2.04694079 * ht.iloc[0] ** 0.55768901
    assert biomassa_aeria_kg(sp, dbh, ht)[0] == pytest.approx(esperat)
    assert FRACCIO_CARBONI == 0.5


def test_sdi_reineke():
    assert sdi(1000, 25.0) == pytest.approx(1000)


def test_creixement_validat_sense_biaix_greu():
    import json

    calibracio = json.loads((DATA / "calibracio.json").read_text(encoding="utf-8"))
    for especie, p in calibracio["especies"].items():
        v = p["validacio"]
        assert abs(v["parceles"]["biaix"]) < 0.5, especie
        assert v["parceles"]["r2"] > 0.9, especie
        assert v["mortalitat"]["morts_predits"] == pytest.approx(v["mortalitat"]["morts_observats"], rel=0.05)
