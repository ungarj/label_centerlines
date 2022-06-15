import fiona
import os
import pytest
from shapely.geometry import shape

TESTDATA_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), "testdata")


@pytest.fixture
def alps_shape():
    with fiona.open(os.path.join(TESTDATA_DIR, "alps.geojson"), "r") as src:
        return shape(next(iter(src))["geometry"])
