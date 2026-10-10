from pathlib import Path

from streamlit.testing.v1 import AppTest

APP = str(Path(__file__).resolve().parent.parent / "app.py")


def iniciar():
    return AppTest.from_file(APP, default_timeout=120).run()


def test_parcela_inf_amb_aclarida():
    at = iniciar()
    assert not at.exception
    assert len(at.metric) == 4
    assert at.metric[2].label == "Captació neta amb aclarida"
    assert at.success or at.error


def test_canvi_de_parcela_i_tipus_d_aclarida():
    at = iniciar()
    at.selectbox[0].set_value(at.selectbox[0].options[10]).run()
    assert not at.exception
    for tipus in ("dalt", "sistematica"):
        at.selectbox[1].set_value(tipus).run()
        assert not at.exception


def test_dades_de_rodal_sense_aclarida():
    at = iniciar()
    at.radio[0].set_value("Dades del rodal").run()
    assert not at.exception
    at.toggle[0].set_value(False).run()
    assert not at.exception
    assert at.metric[2].label == "Carboni aeri avui"


def test_rodal_molt_dens_avisa_d_extrapolacio():
    at = iniciar()
    at.radio[0].set_value("Dades del rodal").run()
    at.number_input[2].set_value(85.0).run()
    assert not at.exception
    assert any("extrapola" in w.value for w in at.warning)


def test_inventari_csv_sense_fitxer_demana_fitxer():
    at = iniciar()
    at.radio[0].set_value("Inventari propi (CSV)").run()
    assert not at.exception
    assert at.info
