import click
import fiona
import logging
from shapely.geometry import shape, mapping
import tqdm

from label_centerlines import __version__, get_centerline

formatter = logging.Formatter('%(asctime)s %(levelname)s %(name)s %(message)s')
stream_handler = logging.StreamHandler()
stream_handler.setFormatter(formatter)
stream_handler.setLevel(logging.INFO)
logging.getLogger().addHandler(stream_handler)

logger = logging.getLogger(__name__)


@click.command()
@click.version_option(version=__version__, message='%(version)s')
@click.argument('input_path')
@click.argument('output_path')
@click.option(
    '--segmentize_maxlen',
    type=float,
    help="Maximum segment length for polygon borders. (default: 0.5)",
    default=0.5
)
@click.option(
    "--max_points",
    type=int,
    help="Number of points per geometry allowed before simplifying. "
         "(default: 3000)",
    default=3000
)
@click.option(
    "--simplification",
    type=float,
    help="Simplification threshold. "
         "(default: 0.05)",
    default=0.05
)
@click.option(
    "--smooth",
    type=int,
    help="Smoothness of the output centerlines. "
         "(default: 5)",
    default=5
)
@click.option(
    "--output_driver",
    type=click.Choice(['GeoJSON', 'GPKG']),
    help="Output format. "
         "(default: 'GeoJSON')",
    default="GeoJSON"
)
@click.option(
    "--debug",
    is_flag=True,
    help="show debug log output"
)
def main(
    input_path, output_path, segmentize_maxlen, max_points, simplification,
    smooth, output_driver, debug
):
    if debug:
        logging.getLogger("label_centerlines").setLevel(logging.DEBUG)
        stream_handler.setLevel(logging.DEBUG)
    with fiona.open(input_path, "r") as src:
        with fiona.open(
            output_path,
            "w",
            schema=dict(src.schema.copy(), geometry="LineString"),
            crs=src.crs,
            driver=output_driver
        ) as dst:
            for feature in tqdm.tqdm(src):
                centerline = get_centerline(shape(feature["geometry"]))
                if centerline is None:
                    logger.error(
                        "centerline could not be extracted from feature %s",
                        feature["properties"]
                    )
                    continue
                dst.write({
                    'properties': feature['properties'],
                    'geometry': mapping(centerline)
                })
