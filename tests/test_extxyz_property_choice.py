"""Choosing which property becomes REF_energy / REF_stress in an extxyz file.

A quantity averaged over an MD trajectory is stored alongside its statistics,
named with a comma. A trailing '*' in the search pattern matches those too, so
taking the last match selected a statistic instead of the quantity.

The names below are the ones a LAMMPS OPLS-AA run actually wrote (job 4785).
"""

from read_structure_step.formats.extxyz.extxyz import _latest_property


class FakeProperties:
    """Just the `list` method `_latest_property` uses."""

    def __init__(self, keys):
        self._keys = keys

    def list(self, pattern):
        stem = pattern[:-1] if pattern.endswith("*") else pattern
        return [k for k in self._keys if k.startswith(stem)]


LAMMPS_MD = FakeProperties(
    [
        "potential energy#LAMMPS#oplsaa+",
        "potential energy, stderr#LAMMPS#oplsaa+",
        "potential energy, tau#LAMMPS#oplsaa+",
        "potential energy, inefficiency#LAMMPS#oplsaa+",
    ]
)


def test_the_quantity_is_chosen_not_its_statistics():
    """The regression: 'inefficiency' was picked -- dimensionless, so writing
    it raised ValueError on Q_('')."""
    assert _latest_property(LAMMPS_MD, "potential energy*") == (
        "potential energy#LAMMPS#oplsaa+"
    )


def test_stderr_is_not_chosen():
    """The quieter danger: stderr *is* an energy in the right units, so it
    would have been written as REF_energy with nothing to show for it."""
    keys = [
        "potential energy, stderr#LAMMPS#oplsaa+",
        "potential energy#LAMMPS#oplsaa+",
    ]
    assert "stderr" not in _latest_property(FakeProperties(keys), "potential energy*")


def test_stress_too():
    props = FakeProperties(
        [
            "stress#LAMMPS#oplsaa+",
            "stress, stderr#LAMMPS#oplsaa+",
            "stress, tau#LAMMPS#oplsaa+",
            "stress, inefficiency#LAMMPS#oplsaa+",
        ]
    )
    assert _latest_property(props, "stress*") == "stress#LAMMPS#oplsaa+"


def test_several_models_keeps_the_most_recent():
    """The original intent of taking the last match: several models, not
    several statistics."""
    props = FakeProperties(
        [
            "total energy#LAMMPS#oplsaa+",
            "total energy, stderr#LAMMPS#oplsaa+",
            "total energy#MACE#water",
        ]
    )
    assert _latest_property(props, "total energy*") == "total energy#MACE#water"


def test_no_match_returns_none():
    assert _latest_property(LAMMPS_MD, "enthalpy of formation*") is None


def test_a_bare_name_without_a_model_still_works():
    props = FakeProperties(["energy", "energy, stderr"])
    assert _latest_property(props, "energy*") == "energy"


def test_the_writer_is_still_the_registered_writer():
    """The helper was first added between @register_writer and the function it
    decorates, so the decorator bound the helper instead: writing any extxyz
    file then died with

        TypeError: _latest_property() got an unexpected keyword argument 'extension'

    The unit tests above all passed, because they import the helper directly
    and never go through the registry.
    """
    import read_structure_step.formats.extxyz.extxyz  # noqa: F401 - registers it
    from read_structure_step.formats.registries import REGISTERED_WRITERS

    registered = REGISTERED_WRITERS[".extxyz"]["function"]
    assert registered.__name__ == "write_extxyz", registered.__name__
