import joblib
import pytest
from pydefect.analyzer.grids import Grids
from pymatgen.core import Lattice, Structure
import numpy as np
from pymatgen.io.vasp import Chgcar
from vise.tests.helpers.assertion import assert_dataclass_almost_equal


@pytest.fixture
def grids():
    return Grids(lattice=Lattice.cubic(10),
                 dim=(1, 1, 5),
                 distance_data=np.array([[[0.0, 2.0, 4.0, 4.0, 2.0]]]))


@pytest.fixture
def chgcar():
    struc = Structure(lattice=Lattice.cubic(10), species=["H"], coords=[[0]*3])
    data = {"total": np.array([[[0.0, 2.0, 4.0, 6.0, 8.0]]]),
            "diff": np.array([[[0.0, 1.0, 2.0, 3.0, 4.0]]])}
    return Chgcar(struc, data=data)


def test_grids_joblib_roundtrip(tmpdir, grids):
    tmpdir.chdir()
    print(tmpdir)
    with open("tmp.joblib", mode="wb") as f:
        joblib.dump(grids, f, compress=3)
    with open("tmp.joblib", mode="rb") as f:
        actual = joblib.load(f)
    assert_dataclass_almost_equal(actual, grids)


def test_grids_np_save_load_roundtrip(tmpdir, grids):
    tmpdir.chdir()
    print(tmpdir)
    grids.dump()
    actual = grids.from_file()
    assert_dataclass_almost_equal(actual, grids)


def test_grids_from_chgcar(grids, chgcar):
    actual = Grids.from_chgcar(chgcar)
    assert_dataclass_almost_equal(actual, grids)


def test_shift_distance_data(grids):
    actual = grids.shifted_distance_data(center=[0, 0, 1])
    expected = np.array([[[2.0, 0.0, 2.0, 4.0, 4.0]]])
    np.testing.assert_array_almost_equal(actual, expected)


def test_shift_distance_data2():
    grids = Grids(lattice=Lattice.cubic(10),
                  dim=(2, 2, 2),
                  distance_data=np.array([[[0.0, 5.0], [5.0, 7.07]],
                                          [[5.0, 7.07], [7.07, 8.66]]]))
    actual = grids.shifted_distance_data(center=[1, 1, 1])
    expected = np.array([[[8.66, 7.07], [7.07, 5.0]],
                         [[7.07, 5.0], [5.0, 0.0]]])
    np.testing.assert_array_almost_equal(actual, expected)


def test_spherical_dist(grids):
    # distance_data=np.array([[[0.0, 2.0, 4.0, 4.0, 2.0]]]))
    actual = grids.spherical_dist(data=np.array([[[0.0, 0.0, 0.0, 0.0, 1.0]]]),
                                  center=[0, 0, 1],
                                  distance_bins=np.array([0.0, 2.5, 5.0]))
    #  Divide by 2 since there are 2 points at 4.0 distance. volume=1000.
    expected = [0.0, 1.0 / 2 / 1000]
    assert actual == expected


"""
TODO
- Check how to revolve numpy array.
- Add distances_data
- Add _calc_histogram(chgcar, distances_data, center) method

DONE
- Add defect_center_idxs
- Add defect_center_coords property
"""