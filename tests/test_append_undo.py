# -*- coding: utf-8 -*-

"""Write Structure's appends are safe to repeat (phase 5 resume)."""

from pathlib import Path
from types import SimpleNamespace

from read_structure_step.write_structure import WriteStructure

undo = WriteStructure._undo_earlier_appends


def test_first_time_records_sizes(tmp_path):
    step = SimpleNamespace(directory=str(tmp_path / "1"))
    existing = tmp_path / "all.sdf"
    existing.write_text("frame 0\n")
    new = tmp_path / "new.sdf"
    undo(step, [existing, new])
    existing.write_text("frame 0\nframe 1\n")  # the step appends
    new.write_text("frame 1\n")
    # The step runs again in the same directory: its appends are undone.
    undo(step, [existing, new])
    assert existing.read_text() == "frame 0\n"
    assert not new.exists()


def test_other_directory_is_another_step(tmp_path):
    """Each loop iteration has its own directory, so appends accumulate."""
    target = tmp_path / "all.sdf"
    for i in range(3):
        step = SimpleNamespace(directory=str(tmp_path / f"iter_{i}" / "1"))
        undo(step, [target])
        with open(target, "a") as fd:
            fd.write(f"frame {i}\n")
    assert target.read_text() == "frame 0\nframe 1\nframe 2\n"


def test_target_files_split():
    names = WriteStructure._target_files(Path("/j/out.sdf.gz"), 2, 5)
    assert [n.name for n in names] == ["out_1.sdf.gz", "out_3.sdf.gz", "out_5.sdf.gz"]
    assert WriteStructure._target_files(Path("/j/out.sdf"), "all", 5) == [
        Path("/j/out.sdf")
    ]
