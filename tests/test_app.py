from pathlib import Path

from streamlit.testing.v1 import AppTest

APP = str(Path(__file__).resolve().parent.parent / "app.py")


def iniciar():
    return AppTest.from_file(APP, default_timeout=120).run()


def text_markdown(at) -> str:
    return "\n".join(m.value for m in at.markdown)


def test_parcela_inf_amb_aclarida():
    at = iniciar()
    assert not at.exception
    text = text_markdown(at)
    assert "CO₂ emmagatzemat avui" in text
    assert "amb aclarida" in text
    assert "CO₂ retirat en l'aclarida" in text
    assert "respecte de l'escenari sense aclarida" in text


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
    text = text_markdown(at)
    assert "Àrea basimètrica avui" in text
    assert "CO₂ retirat" not in text


def test_rodal_molt_dens_avisa_d_extrapolacio():
    at = iniciar()
    at.radio[0].set_value("Dades del rodal").run()
    at.number_input[2].set_value(85.0).run()
    assert not at.exception
    assert "Projecció parcialment extrapolada" in text_markdown(at)
    assert any("s'extrapola" in m.value for m in at.markdown)


def test_inventari_csv_sense_fitxer_demana_fitxer():
    at = iniciar()
    at.radio[0].set_value("Inventari propi (CSV)").run()
    assert not at.exception
    assert at.info


def test_sliders_amb_pas_d_un():
    at = iniciar()
    pes, anys = at.slider[0], at.slider[1]
    assert (pes.value, pes.step, anys.value, anys.step) == (30, 1, 20, 1)
    pes.set_value(17).run()
    assert not at.exception
    at.slider[1].set_value(80).run()
    assert not at.exception
    assert "horitzó <strong>80 anys</strong>" in text_markdown(at)
    assert "del 17 % de l'AB" in text_markdown(at)


def test_pertorbacions_tres_modes():
    at = iniciar()
    radio = next(r for r in at.radio if r.label == "Escenari de pertorbacions")
    assert radio.value == "risc_mitja"
    assert "pertorbacions: risc mitjà" in text_markdown(at)
    assert "Interval:" in text_markdown(at)
    assert any("cas pessimista" in c.value for c in at.caption)
    radio.set_value("cap").run()
    assert not at.exception
    assert "pertorbacions: sense pertorbacions" in text_markdown(at)
    assert "Interval:" not in text_markdown(at)
    next(r for r in at.radio if r.label == "Escenari de pertorbacions").set_value("episodi").run()
    assert not at.exception
    tipus = next(s for s in at.selectbox if s.label == "Tipus d'episodi")
    tipus.set_value("sequera").run()
    assert not at.exception
    assert "episodi de sequera l'any 10" in text_markdown(at)
    any_ep = next(s for s in at.slider if s.label == "Any de l'episodi")
    any_ep.set_value(3).run()
    assert not at.exception
    assert "l'any 3" in text_markdown(at)
