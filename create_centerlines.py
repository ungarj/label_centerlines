#!/usr/bin/env python

# Author:  Joachim Ungar <joachim.ungar@eox.at>
#
#-------------------------------------------------------------------------------
# Copyright (C) 2015 EOX IT Services GmbH
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies of this Software or works derived from this Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#-------------------------------------------------------------------------------

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
from progressbar import ProgressBar

def main(args):

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "input_shp",
        type=str,
        help="input polygons"
        )
    parser.add_argument(
        "output_geojson",
        type=str,
        help="output centerlines"
        )
    parser.add_argument(
        "--segmentize_maxlen",
        type=float,
        help="maximum length used when segmentizing polygon borders",
        default=0.5
        )
    parser.add_argument(
        "--max_points",
        type=int,
        help="number of points per geometry allowed before simplifying",
        default=1000
        )
    parser.add_argument(
        "--simplification",
        type=float,
        help="value which increases simplification when necessary",
        default=0.1
        )
    parser.add_argument(
        "--smooth",
        type=int,
        help="smoothness of the output centerlines",
        default=5
        )
    parsed = parser.parse_args(args)
    input_shp = parsed.input_shp
    output_geojson = parsed.output_geojson
    segmentize_maxlen = parsed.segmentize_maxlen
    max_points = parsed.max_points
    simplification = parsed.simplification
    smooth_sigma = parsed.smooth


    with fiona.open(input_shp, "r") as regions:
        out_schema = regions.schema.copy()
        out_schema['geometry'] = "LineString"
        with fiona.open(
            output_geojson,
            "w",
            schema=out_schema,
            crs=regions.crs,
            driver="GeoJSON"
            ) as out_centerlines:
            pbar = ProgressBar(maxval=len(regions)).start()
            progress_counter = 0
            failed_features = []
            print len(regions), "features"
            for region in regions:
                progress_counter += 1
                pbar.update(progress_counter)
                geom = shape(region['geometry'])
                try:
                    centerlines_geom = get_centerlines_from_geom(
                        geom,
                        segmentize_maxlen=segmentize_maxlen,
                        max_points=max_points,
                        simplification=simplification,
                        smooth_sigma=smooth_sigma
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
                    failed_features.append(region['properties']['name'])
            pbar.finish()
            if len(failed_features) > 0:
                print "failed:"
                for feature in failed_features:
                    print feature


def get_centerlines_from_geom(
    geometry,
    segmentize_maxlen=0.5,
    max_points=1000,
    simplification=0.1,
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

        # Simplify segmentized geometry if necessary. This step is required
        # as huge geometries slow down the centerline extraction significantly.
        tolerance = simplification
        while len(points) > max_points:
            # If geometry is too large, apply simplification until geometry
            # is simplified enough (indicated by the "max_points" value)
            tolerance += simplification
            simplified = boundary.simplify(tolerance)
            points = simplified.coords

        # Calculate Voronoi diagram.
        vor = Voronoi(points)

        # The next three steps are the most processing intensive and probably
        # not the most efficient method to get the skeleton centerline. If you
        # have any recommendations, I would be very happy to know.

        # Convert to networkx graph.
        graph = graph_from_voronoi(vor, geometry)

        # Get end nodes from graph.
        end_nodes = get_end_nodes(graph)

        if len(end_nodes) < 2:
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
    paths_sorted = [x for (y,x) in sorted(zip(distances, paths), reverse=True)]
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
            # Eliminate all points outside our geometry.
            if point1.within(geometry) and point2.within(geometry):
                dist = point1.distance(point2)
                graph.add_nodes_from([i[0], i[1]])
                graph.add_edge(i[0], i[1], weight=dist)
    return graph


if __name__ == "__main__":
        main(sys.argv[1:])
