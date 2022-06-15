import click
import concurrent.futures
from contextlib import ExitStack
import fiona
import logging
from shapely.geometry import shape, mapping
import time
import tqdm

from label_centerlines import __version__, get_centerline
from label_centerlines.exceptions import CenterlineError


class TqdmHandler(logging.StreamHandler):
    """Custom handler to avoid log outputs to interfere with progress bar."""

    def __init__(self):
        logging.StreamHandler.__init__(self)

    def emit(self, record):
        msg = self.format(record)
        tqdm.tqdm.write(msg)


formatter = logging.Formatter("%(asctime)s %(levelname)s %(name)s %(message)s")
stream_handler = TqdmHandler()
stream_handler.setFormatter(formatter)
stream_handler.setLevel(logging.INFO)
logging.getLogger().addHandler(stream_handler)

logger = logging.getLogger(__name__)


@click.command()
@click.version_option(version=__version__, message="%(version)s")
@click.argument("input_path")
@click.argument("output_path")
@click.option(
    "--segmentize_maxlen",
    type=float,
    help="Maximum segment length for polygon borders. (default: 0.5)",
    default=0.5,
)
@click.option(
    "--max_points",
    type=int,
    help="Number of points per geometry allowed before simplifying. " "(default: 3000)",
    default=3000,
)
@click.option(
    "--simplification",
    type=float,
    help="Simplification threshold. " "(default: 0.05)",
    default=0.05,
)
@click.option(
    "--smooth",
    type=int,
    help="Smoothness of the output centerlines. " "(default: 5)",
    default=5,
)
@click.option(
    "--max_paths",
    type=int,
    help="Number of longest paths used to create the centerlines. " "(default: 5)",
    default=5,
)
@click.option(
    "--output_driver",
    type=click.Choice(["GeoJSON", "GPKG"]),
    help="Output format. " "(default: 'GeoJSON')",
    default="GeoJSON",
)
@click.option("--verbose", is_flag=True, help="show information on processed features")
@click.option("--debug", is_flag=True, help="show debug log messages")
def main(
    input_path,
    output_path,
    segmentize_maxlen,
    max_points,
    simplification,
    smooth,
    max_paths,
    output_driver,
    verbose,
    debug,
):
    """
    Read features, convert to centerlines and write to output.

    Multipart features (MultiPolygons) from input will be converted to
    singlepart features, i.e. all output features written will be LineString
    geometries, not MultiLineString geometries.
    """
    # set up logger
    log_level = logging.DEBUG if debug else logging.INFO
    logging.getLogger("label_centerlines").setLevel(log_level)
    stream_handler.setLevel(log_level)

    with ExitStack() as es:
        # set up context managers for fiona & process pool
        src = es.enter_context(fiona.open(input_path, "r"))
        dst = es.enter_context(
            fiona.open(
                output_path,
                "w",
                schema=dict(src.schema.copy(), geometry="LineString"),
                crs=src.crs,
                driver=output_driver,
            )
        )
        executor = es.enter_context(concurrent.futures.ProcessPoolExecutor())

        tasks = (
            executor.submit(
                _feature_worker,
                feature,
                segmentize_maxlen,
                max_points,
                simplification,
                smooth,
                max_paths,
            )
            for feature in src
        )
        for task in tqdm.tqdm(
            concurrent.futures.as_completed(tasks), disable=debug, total=len(src)
        ):
            # output is split up into parts of single part geometries to meet
            # GeoPackage requirements
            for part in task.result():
                feature, elapsed = part
                if "geometry" in feature:
                    dst.write(feature)
                else:
                    logger.error(
                        "centerline could not be extracted from feature %s",
                        feature["properties"],
                    )
                if verbose:
                    tqdm.tqdm.write("%ss: %s" % (elapsed, feature["properties"]))


def _feature_worker(
    feature, segmentize_maxlen, max_points, simplification, smooth, max_paths
):
    try:
        start = time.time()
        centerline = get_centerline(
            shape(feature["geometry"]),
            segmentize_maxlen,
            max_points,
            simplification,
            smooth,
            max_paths,
        )
    except CenterlineError:
        return [(dict(properties=feature["properties"]), round(time.time() - start, 3))]
    finally:
        elapsed = round(time.time() - start, 3)

    if centerline.geom_type == "LineString":
        return [(dict(feature, geometry=mapping(centerline)), elapsed)]
    elif centerline.geom_type == "MultiLineString":
        return [
            (dict(feature, geometry=mapping(subgeom)), elapsed)
            for subgeom in centerline.geoms
        ]
