import click
import concurrent.futures
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
        log_level = logging.DEBUG
    else:
        log_level = logging.INFO
    logging.getLogger("label_centerlines").setLevel(log_level)
    stream_handler.setLevel(log_level)

    with fiona.open(input_path, "r") as src:
        with fiona.open(
            output_path,
            "w",
            schema=dict(src.schema.copy(), geometry="LineString"),
            crs=src.crs,
            driver=output_driver
        ) as dst:
            with concurrent.futures.ProcessPoolExecutor() as executor:
                shapes = [shape(f["geometry"]) for f in src]
                for centerline, feature in tqdm.tqdm(
                    zip(executor.map(
                        get_centerline,
                        shapes,
                        [segmentize_maxlen for _ in range(len(src))],
                        [max_points for _ in range(len(src))],
                        [simplification for _ in range(len(src))],
                        [smooth for _ in range(len(src))]
                    ), src),
                    disable=debug,
                    total=len(src)
                ):
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
