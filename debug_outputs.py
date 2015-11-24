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
from shapely.geometry import shape, mapping
from progressbar import ProgressBar

from src_create_centerlines import get_centerlines_from_geom

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
        default=3000
        )
    parser.add_argument(
        "--simplification",
        type=float,
        help="value which increases simplification when necessary",
        default=0.005
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
            debug_outputs = {
                "original_points": [],
                "segmentized_points": [],
                "voronoi": [],
                "centerline": []
                }
            print len(regions), "features"
            for region in regions:
                progress_counter += 1
                pbar.update(progress_counter)
                geom = shape(region['geometry'])
                from src_create_centerlines import debug_output
                try:
                    centerlines_geom = get_centerlines_from_geom(
                        geom,
                        segmentize_maxlen=segmentize_maxlen,
                        max_points=max_points,
                        simplification=simplification,
                        smooth_sigma=smooth_sigma,
                        debug=True
                        )
                except:
                    raise
                if centerlines_geom:
                    try:
                        for step in debug_outputs:
                            debug_outputs[step].append(debug_output[step])
                    except:
                        pass
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

    out_schema = {
        "geometry": "GeometryCollection",
        "properties": {
            "type": "str"
            }
        }
    with fiona.open(
        "temp/debug.geojson",
        "w",
        schema=out_schema,
        crs=regions.crs,
        driver="GeoJSON"
        ) as debug_file:

        for step, output in debug_outputs.iteritems():
            for item in output:
                debug_feature = {
                    "geometry": mapping(item),
                    "properties": {
                        "type": step
                        }
                    }
                debug_file.write(debug_feature)




if __name__ == "__main__":
    main(sys.argv[1:])
