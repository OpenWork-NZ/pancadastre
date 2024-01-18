#! /usr/local/bin/python3

from xml.etree import ElementTree as ET
import rdflib
from rdflib.namespace import RDFS, RDF, GEO
import json

from datetime import datetime
from pyproj import Transformer, Geod

from math import pi, isclose, ceil, sqrt, isnan
tau = 2*pi

from data import *
from landxml import *
from canon import *
from jsonfg import *
from rdf import *
from wa import *
from reports import *

#class ReducedVertPos(Observation) for 3D?

## Processors
# NOTE: Only interpolates in 2D, not 3D.
def interpCurves(data, delta):
  def inner(subdata):
    for geom in subdata.geom:
      segments = []
      for segment in geom.segments:
        if isinstance(segment, Cubic):
          points = segment.asObs().interpolatedPath()
          subsegs = IDList()
          for i in range(1, len(points)):
            subsegs.append(Line(segment.id + str(i), points[i-1], points[i]))
          subsegs.id = segment.id
          segments.append(subsegs)
          continue
        elif not isinstance(segment, Curve):
          segments.append(segment)
          continue
        interp = interpCurve(segment.start, segment.mid, segment.end)


        #east0, north0 = segment.start
        #east1, north1 = segment.end
        d = segment.start - segment.end
        dist = sqrt(d.coord1*d.coord1 + d.coord2*d.coord2 + (d.coord3*d.coord3 if d.coord3 is not None else 0))
        subdivisions = ceil(dist/delta)

        subsegs = IDList()
        prev = interp(0)
        for t in range(subdivisions):
          end = interp((t+1)/subdivisions)
          subsegs.append(Line(t, prev, end))
          prev = end
        subsegs.id = segment.id # Attach extra property to preserve original ID.
        segments.append(subsegs)
      geom.segments = segments

  for parcel in data.parcels: inner(parcel)
  for feature in data.features: inner(feature)
  for observations in data.survey.observationGroups.values():
    for observation in observations:
      if isinstance(observation, ReducedArcObservation) or isinstance(observation, ArcByChord):
        segment = observation.geom
        interp = interpCurve(segment.start, segment.mid, segment.end)
        #east0, north0 = segment.start
        #east1, north1 = segment.end
        d = segment.start - segment.end
        dist = sqrt(d.coord1*d.coord1 + d.coord2*d.coord2 + (d.coord3*d.coord3 if d.coord3 is not None else 0))
        if isnan(dist):
          print("WARNING: Failed to compute distance between the 2 points!",
              (east0, north0), (east1, north1))
          dist = 1
        subdivisions = ceil(dist/delta)

        subsegs = IDList()
        prev = interp(0)
        for t in range(subdivisions):
          end = interp((t+1)/subdivisions)
          subsegs.append(Line(t, prev, end))
          prev = end
        subsegs.id = segment.id
        observation.geom = subsegs

def interpCurve(a, m, b):
  # See https://observablehq.com/@jrus/circle-arc-interpolation
  b_m = b - m
  m_a = m - a
  ab_m = a.coord1*b_m.coord1 - a.coord2*b_m.coord2, a.coord1*b_m.coord2 + a.coord2*b_m.coord1 # 3rd dimension?
  ab_m = Point(None, None, None, ab_m[0], ab_m[1], a.coord3, a.objID)
  bm_a = b.coord1*m_a.coord1 - b.coord2*m_a.coord2, b.coord1*m_a.coord2 + b.coord2*m_a.coord1 # 3rd dimension?
  bm_a = Point(None, None, None, bm_a[0], bm_a[1], b.coord3, b.objID)

  def inner(t):
    num = ab_m.mul(1-t) + bm_a.mul(t)
    den = b_m.mul(1-t) + m_a.mul(t)
    dist2 = den.coord1*den.coord1 + den.coord2*den.coord2
    if den.coord3 is not None: dist2 += den.coord3*den.coord3
    x = (num.coord1*den.coord1 + num.coord2*den.coord2)/dist2
    y = (num.coord2*den.coord1 - num.coord1*den.coord2)/dist2
    z = a.coord3 # 3rd dimension?
    return Point(None, None, None, x, y, z, str(a.objID) + ':' + str(t))
  return inner

## Commandline interface

if __name__ == "__main__":
  import argparse

  argparser = argparse.ArgumentParser()
  argparser.add_argument('-X', '--LANDXML', help="LandXML input file", nargs='+')
  argparser.add_argument('-V', '--vocab', help="Jurisdictional vocabulary to align to")
  argparser.add_argument('-C', '--CSDM', help="JSON-Topology input file", nargs='+')
  argparser.add_argument('-W', '--WA', help="Western Australia CSD input format (dramatically incomplete!)")
  argparser.add_argument('--interpolate', help="Interpolate curves, possibly specifying length in meters of each segment (default 1m)",
      const=1, type=int, nargs='?')
  argparser.add_argument('--epsg', help="Overwrite the projection being used.")
  argparser.add_argument('-j', '--jsonfg', help="JSONfg output file", const="?.json", nargs='?')
  argparser.add_argument('-k', '--jsonfgparcel', help="JSONfg output file for parcels", const="?.parcels.json", nargs='?')
  argparser.add_argument('-r', '--rdf', help="RDF syntax to output", const="ttl", nargs='?')
  argparser.add_argument('-o', '--output', help="RDF file to output to", const="?.rdf", nargs='?')
  argparser.add_argument('-c', '--csdm', help="CSDM JSON output file", const="?.json", nargs='?')
  argparser.add_argument('-s', '--summary', help="Textual summary of parsed data", const="?.txt", nargs='?')
  argparser.add_argument('-g', '--geojson', help="GeoJSON output file", const="?.json", nargs='?')
  argparser.add_argument('-p', '--geojsonparcel', help="GeoJSON output file for parcels", const="?.parcels.json", nargs='?')
  argparser.add_argument('-e', '--errorlog', help="Log of any encountered errors", const="?.log", nargs='?')
  args = argparser.parse_args()

  wrote_output = False
  def export(infile, data):
    global wrote_output
    from os import path
    filename = path.splitext(path.basename(infile))[0]
    
    if data.projection is None or data.projection.horizontal is None:
      print("WARNING! No projection system specified! Defaulting to wgs84")
      data.projection = Projection('wgs84')

    if args.epsg: data.projection = Projection(args.epsg)
    if args.interpolate:
      interpCurves(data, args.interpolate)

    if args.jsonfg:
      with open(args.jsonfg.replace("?", filename), "w") as f: exportJSONfg(data, f)
      wrote_output = True
    if args.jsonfgparcel:
      with open(args.jsonfgparcel.replace("?", filename), "w") as f: exportJSONfgParcel(data, f)
      wrote_output = True
    # Extra condition suppresses GeoJSON files for empty geometry counts.
    group = data.survey.observationGroups
    if args.geojson and (len(data.monuments) +
        sum(1 for parcel in data.parcels
            for geom in parcel.geom for seg in geom.segments) +
        sum(1 for _ in data.survey.observationGroups.items()
            for i, observation in enumerate(group))
        > 0):
      with open(args.geojson.replace("?", filename), "w") as f:
        exportGeoJSON(data, f)
      wrote_output = True
    if args.geojsonparcel and (len(data.monuments) +
        sum(1 for parcel in data.parcels
            for geom in parcel.geom for seg in geom.segments) +
        sum(1 for _ in data.survey.observationGroups.items()
            for i, observation in enumerate(group))
        > 0):
      with open(args.geojsonparcel.replace("?", filename), "w") as f:
        exportGeoJSONParcel(data, f)
      wrote_output = True
    if args.output or args.rdf:
      with open((args.output or "?.rdf").replace("?", filename), "w") as f:
        exportRDF(data, f, format=args.rdf or "turtle")
      wrote_output = True
    if args.csdm:
      with open(args.csdm.replace("?", filename), "w") as f: exportCSDM(data, f)
      wrote_output = True
    if args.summary:
      with open(args.summary.replace("?", filename), "w") as f: f.write(exportSummary(data))

  def redirect_errors(infile):
    import sys
    from os import path
    filename = path.splitext(path.basename(infile))[0]

    if args.errorlog:
      sys.stdout = open(args.errorlog.replace("?", filename), "w")
      print(infile)

  read_input = False
  for source in args.LANDXML or []:
    redirect_errors(source)
    try:
      with open(source) as f: export(source, importLandXML(f))
    except KeyError as err:
      print("Unresolvable link!", err)
    except Exception as err:
      print(err)
    read_input = True
  for source in args.CSDM or []:
    redirect_errors(source)
    try:
      with open(source) as f: export(source, importCSDM(f))
    except KeyError as err:
      print("Unresolvable link or missing property!", err)
    except Exception as err:
      print(err)
    read_input = True
  for source in args.WA or []:
    redirect_errors(source)
    try:
      with open(source) as f: export(source, importWA_CSD(f))
    except KeyError as err:
      print("Unresolvable link!", err)
    except Exception as err:
      print(err)
    read_input = True

  if not read_input:
    print("Please specify at least one input file!")
    argparser.print_help()
