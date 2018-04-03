from itertools import combinations
import logging
import networkx as nx
from networkx.exception import NetworkXNoPath
import numpy as np
from scipy.spatial import Voronoi
from scipy.ndimage import filters
from shapely.geometry import LineString, MultiLineString, Point, MultiPoint

logger = logging.getLogger(__name__)


def get_centerline(
    geom,
    segmentize_maxlen=0.5,
    max_points=3000,
    simplification=0.05,
    smooth_sigma=5
):
    """
    Return centerline from geometry.

    Parameters:
    -----------
    geom : shapely Polygon or MultiPolygon
    segmentize_maxlen : Maximum segment length for polygon borders.
        (default: 0.5)
    max_points : Number of points per geometry allowed before simplifying.
        (default: 3000)
    simplification : Simplification threshold.
        (default: 0.05)
    smooth_sigma : Smoothness of the output centerlines.
        (default: 5)
    """
    logger.debug("geometry type %s", geom.geom_type)
    if geom.geom_type == "MultiPolygon":
        logger.debug("MultiPolygon found with %s sub-geometries", len(geom))
        centerlines = []
        for subgeom in geom:
            centerline = get_centerline(
                subgeom, segmentize_maxlen, max_points, simplification,
                smooth_sigma
            )
            if centerline is not None:
                logger.debug("subgeometries:")
                logger.debug(centerline)
                centerlines.append(centerline)

        logger.debug("create MultiLineString")
        return MultiLineString(centerlines)

    elif geom.geom_type == "Polygon":
        # segmentized Polygon outline
        outline = _segmentize(
            geom.exterior,  # if len(geom.interiors) else geom.outline,
            segmentize_maxlen
        )
        logger.debug("outline: %s", outline)

        # simplify segmentized geometry if necessary and get points
        outline_points = outline.coords
        simplification_updated = simplification
        while len(outline_points) > max_points:
            # If geometry is too large, apply simplification until geometry
            # is simplified enough (indicated by the "max_points" value)
            simplification_updated += simplification
            outline_points = outline.simplify(simplification_updated).coords
        logger.debug("simplification used: %s", simplification_updated)
        logger.debug("simplified points: %s", MultiPoint(outline_points))

        # calculate Voronoi diagram and convert to graph but only use points
        # from within the original polygon
        vor = Voronoi(outline_points)
        graph = _graph_from_voronoi(vor, geom)
        logger.debug(
            "voronoi diagram: %s", _multilinestring_from_voronoi(vor, geom)
        )

        # determine longest path between all end nodes from graph
        end_nodes = _get_end_nodes(graph)
        if len(end_nodes) < 2:
            logger.debug("Polygon has too few points")
            return None
        logger.debug("get longest path from %s end nodes", len(end_nodes))
        longest_paths = _get_longest_paths(end_nodes, graph)
        if logger.getEffectiveLevel() <= 10:
            logger.debug("longest paths:")
            for path in longest_paths:
                logger.debug(LineString(vor.vertices[path]))

        # get least curved path from the five longest paths and convert to
        # centerline
        best_path = _get_least_curved_path(longest_paths[:5], vor.vertices)

        logger.debug("smooth linestring")
        centerline = _smooth_linestring(
            LineString(vor.vertices[best_path]), smooth_sigma)
        logger.debug("centerline: %s", centerline)
        return centerline

    else:
        raise TypeError(
            "Geometry type must be Polygon or MultiPolygon, not %s" %
            geom.geom_type
        )


def _segmentize(geom, max_len):
    """Interpolate points on segments if they exceed maximum length."""
    points = []
    p_xy = None
    for xy in geom.coords:
        if p_xy is not None:
            line_segment = LineString([p_xy, xy])
            points.extend([
                line_segment.interpolate(max_len * i).coords[0]
                for i in range(int(line_segment.length / max_len))
            ])
        p_xy = xy
        points.append(xy)
    return LineString(points)


def _smooth_linestring(linestring, smooth_sigma):
    """Use a gauss filter to smooth out the LineString coordinates."""
    return LineString(
        zip(
            np.array(filters.gaussian_filter1d(linestring.xy[0], smooth_sigma)),
            np.array(filters.gaussian_filter1d(linestring.xy[1], smooth_sigma))
        )
    )


def _get_longest_paths(nodes, graph, maxnum=5):
    """Return longest paths of all possible paths between a list of nodes."""
    def _gen_paths_distances():
        for node1, node2 in combinations(nodes, r=2):
            try:
                # path = nx.shortest_path(graph, node1, node2, "weight")
                yield nx.single_source_dijkstra(
                    graph, node1, node2, 1000000, "weight"
                )
            except NetworkXNoPath:
                continue

    return [
        x for (y, x) in sorted(_gen_paths_distances(), reverse=True)
    ][:maxnum]


def _get_least_curved_path(paths, vertices):
    """Return path with smallest angles."""
    angle_sums = []
    for path in paths:
        path_angles = _get_path_angles(path, vertices)
        angle_sum = abs(sum(path_angles))
        angle_sums.append(angle_sum)
    return [x for (y, x) in sorted(zip(angle_sums, paths))][0]


def _get_path_angles(path, vertices):
    """Return all angles between edges from path."""
    angles = []
    for index, point in enumerate(path):
        if index > 0 and index < len(path) - 1:
            previous_point = vertices[path[index - 1]]
            current_point = vertices[point]
            next_point = vertices[path[index + 1]]
            angles.append(
                _get_angle(
                    (previous_point, current_point),
                    (current_point, next_point)
                )
            )
    return angles


def _get_angle(edge1, edge2):
    """Return angle between edges."""
    v1 = edge1[0] - edge1[1]
    v2 = edge2[0] - edge2[1]
    return np.degrees(np.math.atan2(np.linalg.det([v1, v2]), np.dot(v1, v2)))


def _get_end_nodes(graph):
    """Return list of nodes with just one neighbor node."""
    return [i for i in graph.nodes() if len(list(graph.neighbors(i))) == 1]


def _graph_from_voronoi(vor, geometry):
    """Return networkx.Graph from Voronoi diagram within geometry."""
    graph = nx.Graph()
    for x, y, dist in _yield_ridge_vertices(vor, geometry, dist=True):
        graph.add_nodes_from([x, y])
        graph.add_edge(x, y, weight=dist)
    return graph


def _multilinestring_from_voronoi(vor, geometry):
    """Return MultiLineString geometry from Voronoi diagram."""
    return MultiLineString([
        LineString([
            Point(vor.vertices[[x, y]][0]),
            Point(vor.vertices[[x, y]][1])
        ])
        for x, y in _yield_ridge_vertices(vor, geometry)
    ])


def _yield_ridge_vertices(vor, geometry, dist=False):
    """Yield Voronoi ridge vertices within geometry."""
    for x, y in vor.ridge_vertices:
        if x >= 0 and y >= 0:
            point1 = Point(vor.vertices[[x, y]][0])
            point2 = Point(vor.vertices[[x, y]][1])
            # Eliminate all points outside our geometry.
            if point1.within(geometry) and point2.within(geometry):
                if dist:
                    yield x, y, point1.distance(point2)
                else:
                    yield x, y
