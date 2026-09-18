"""Write Structure uses the standard structure selection and reads old flowcharts."""

import read_structure_step


def _legacy(structures, configurations="all"):
    return {
        "file": {"value": "out.sdf", "units": None},
        "structures": {"value": structures, "units": None},
        "configurations": {"value": configurations, "units": None},
    }


def test_defaults_are_the_current_configuration():
    P = read_structure_step.WriteStructureParameters()
    assert P["source systems"].value == "current"
    assert P["source configurations"].value == "current"
    assert "structures" not in P and "configurations" not in P


def test_legacy_translation():
    P = read_structure_step.WriteStructureParameters()
    P.from_dict(_legacy("current configuration"))
    assert (P["source systems"].value, P["source configurations"].value) == (
        "current",
        "current",
    )
    P.from_dict(_legacy("current system", "all"))
    assert (P["source systems"].value, P["source configurations"].value) == (
        "current",
        "all",
    )
    P.from_dict(_legacy("all systems", "optimized"))
    assert P["source systems"].value == "all"
    assert P["source configurations"].value == "name is"
    assert P["source configuration name"].value == "optimized"
    P.from_dict(_legacy("$dimers"))
    assert P["source systems"].value == "$dimers"
    assert P["file"].value == "out.sdf"


def test_description(monkeypatch):
    node = read_structure_step.WriteStructure()
    node._id = (1,)
    P = node.parameters.values_to_dict()
    P["file"] = "all.sdf"
    P["source systems"] = "all"
    P["source configurations"] = "all"
    text = " ".join(node.description_text(P).split())
    assert "All configurations of every system will be written to all.sdf." in text
    P["ignore missing"] = "no"
    assert "error if the selection matches no structures" in " ".join(
        node.description_text(P).split()
    )
