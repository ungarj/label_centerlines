from click.testing import CliRunner

from label_centerlines import __version__, get_centerline
from label_centerlines.cli import main


def test_cli():
    runner = CliRunner()
    result = runner.invoke(main, ["--version"])
    assert result.exit_code == 0
    assert __version__ in result.output


def test_centerline(alps_shape):
    cl = get_centerline(alps_shape)
    assert cl.is_valid
    assert cl.geom_type == "LineString"
