# -*- coding: utf-8 -*-

"""Write Structure's appends are safe to repeat (phase 5 resume)."""

from pathlib import Path
from types import SimpleNamespace

import pytest
import seamm

from read_structure_step.write_structure import WriteStructure

undo = WriteStructure._undo_earlier_appends


class FakeCheckpointer:
    def __init__(self, run_id):
        self.run_id = run_id


@pytest.fixture(autouse=True)
def a_run():
    """Write Structure inside a checkpointed run, as in a job."""
    seamm.checkpoint.set_checkpointer(FakeCheckpointer("run-1"))
    yield
    seamm.checkpoint.set_checkpointer(None)


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


def test_another_run_starts_afresh(tmp_path):
    """A rerun from the top (a new run) must not cut back to sizes the previous
    run recorded: its iterations' outputs may differ in length."""
    step = SimpleNamespace(directory=str(tmp_path / "iter_1" / "1"))
    target = tmp_path / "all.sdf"
    undo(step, [target])
    target.write_text("frame from run 1\n")
    # Run 2, from the top in the same directory; the file has grown since
    seamm.checkpoint.set_checkpointer(FakeCheckpointer("run-2"))
    target.write_text("frame from run 1\nand more from run 1's later steps\n")
    undo(step, [target])
    assert target.read_text().endswith("later steps\n")  # left as it is


def test_without_a_checkpoint_nothing_is_cut(tmp_path):
    seamm.checkpoint.set_checkpointer(None)
    step = SimpleNamespace(directory=str(tmp_path / "1"))
    target = tmp_path / "all.sdf"
    undo(step, [target])
    target.write_text("frame\n")
    undo(step, [target])
    assert target.read_text() == "frame\n"
