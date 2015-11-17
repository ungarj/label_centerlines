#!/usr/bin/env python

import os
import sys
import argparse
import fiona
from shapely.geometry import shape, LineString, MultiLineString, Point, mapping
from shapely.wkt import loads
import ogr
from scipy.spatial import Voronoi
import networkx as nx
from itertools import combinations
import numpy as np
from scipy.ndimage import filters

MAXPOINTS = 1000

def main(args):

    parser = argparse.ArgumentParser()
    parser.add_argument("input_shp", type=str)
    parser.add_argument("output_shp", type=str)
    parser.add_argument("segmentize_maxlen", type=float)
    parser.add_argument("smooth", type=float)
    parsed = parser.parse_args(args)
    input_shp = parsed.input_shp
    output_shp = parsed.output_shp
    segmentize_maxlen = parsed.segmentize_maxlen
    smooth_sigma = parsed.smooth


    with fiona.open(input_shp, "r") as regions:
        out_schema = regions.schema.copy()
        out_schema['geometry'] = "LineString"
        with fiona.open(
            output_shp,
            "w",
            schema=out_schema,
            crs=regions.crs,
            driver="GeoJSON"
            ) as out_centerlines:
            for region in regions:
                geom = shape(region['geometry'])
                # print "calculating", region['properties']['name']
                try:
                    centerlines_geom = get_centerlines_from_geom(
                        geom,
                        segmentize_maxlen,
                        smooth_sigma
                        )
                except:
                    raise
                if centerlines_geom:
                    centerline = {
                        'properties': region['properties'],
                        'geometry': mapping(centerlines_geom)
                        }
                    out_centerlines.write(centerline)
                else:
                    print region['properties']['name'], "FAILED"

def get_centerlines_from_geom(
    geometry,
    segmentize_maxlen=0.1,
    smooth_sigma=5
    ):
    """
    Returns centerline (for Polygon) or centerlines (for MultiPolygons) as
    LineString or MultiLineString geometries.
    """
    if geometry.geom_type == "MultiPolygon":
        out_centerlines = MultiLineString([
            get_centerlines_from_geom(subgeom, segmentize_maxlen)
            for subgeom in geometry
            if get_centerlines_from_geom(subgeom, segmentize_maxlen) != None
            ])
        return out_centerlines
    else:

        # Convert Polygon to Linestring.
        if len(geometry.interiors) > 0:
            boundary = geometry.exterior
        else:
            boundary = geometry.boundary

        # Convert to OGR object and segmentize.
        ogr_boundary = ogr.CreateGeometryFromWkb(boundary.wkb)
        ogr_boundary.Segmentize(segmentize_maxlen)
        segmentized = loads(ogr_boundary.ExportToWkt())

        # Get points.
        points = segmentized.coords

        tolerance = 0.1
        while len(points) > MAXPOINTS:
            # print "WARNING: geometry too large, simplifying."
            tolerance += 0.1
            # Simplify geometry.
            simplified = boundary.simplify(tolerance)
            # Get points.
            points = simplified.coords

        # Calculate Voronoi diagram.
        vor = Voronoi(points)

        # Convert to networkx graph.
        graph = graph_from_voronoi(vor, geometry)

        # Get end nodes from graph.
        end_nodes = get_end_nodes(graph)

        if len(end_nodes) == 0:
            return None

        # Get longest path.
        longest_path = get_longest_path(
            end_nodes,
            graph,
            vor.vertices
            )
        centerline = LineString(vor.vertices[longest_path])

        # Smooth out geometry.
        centerline_smoothed = smooth_linestring(centerline, smooth_sigma)

        out_centerline = centerline_smoothed

        return out_centerline


def smooth_linestring(linestring, smooth_sigma):
    """
    Uses a gauss filter to smooth out the LineString coordinates.
    """
    smooth_x = np.array(filters.gaussian_filter1d(
        linestring.xy[0],
        smooth_sigma)
        )
    smooth_y = np.array(filters.gaussian_filter1d(
        linestring.xy[1],
        smooth_sigma)
        )
    smoothed_coords = np.hstack((smooth_x, smooth_y))
    smoothed_coords = zip(smooth_x, smooth_y)
    linestring_smoothed = LineString(smoothed_coords)
    return linestring_smoothed


def get_longest_path(nodes, graph, vertices):
    """
    Returns longest path of all possible paths between a list of nodes.
    """
    paths = []
    distances = []
    possible_paths = list(combinations(nodes, r=2))
    for node1, node2 in possible_paths:
        try:
            path = nx.shortest_path(graph, node1, node2, "weight")
        except Exception,e:
            path = []
        if len(path)>1:
            distance = get_path_distance(path, graph)
            paths.append(path)
            distances.append(distance)
    paths_sorted = [x for (y,x) in sorted(zip(distances,paths),reverse=True)]
    longest_path = paths_sorted[0]
    return longest_path


def get_path_distance(path, graph):
    """
    Returns weighted path distance.
    """
    distance = 0
    for i,w in enumerate(path):
        j=i+1
        if j<len(path):
            distance += round(graph.edge[path[i]][path[j]]["weight"], 6)
    return distance


def get_end_nodes(graph):
    """
    Returns list of nodes with just one neighbor node.
    """
    nodelist = [
        i
        for i in graph.nodes_iter()
        if len(graph.neighbors(i))==1
        ]
    return nodelist


def graph_from_voronoi(vor, geometry):
    """
    Creates a networkx graph out of all Voronoi ridge vertices which are inside
    the original geometry.
    """
    graph = nx.Graph()
    for i in vor.ridge_vertices:
        if i[0]>-1 and i[1]>-1:
            point1 = Point(vor.vertices[i][0])
            point2 = Point(vor.vertices[i][1])
            if point1.within(geometry) and point2.within(geometry):
                dist = point1.distance(point2)
                graph.add_nodes_from([i[0], i[1]])
                graph.add_edge(i[0], i[1], weight=dist)
    return graph


if __name__ == "__main__":
        main(sys.argv[1:])
