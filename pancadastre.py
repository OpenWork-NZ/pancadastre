#! /usr/local/bin/python3

from xml.etree import ElementTree as ET
import rdflib
from rdflib.namespace import RDFS, RDF, GEO
import json

from datetime import datetime
from pyproj import Transformer, Geod

from math import pi, isclose, ceil, sqrt
tau = 2*pi

from data import *
from landxml import *
from canon import *
from jsonfg import *
from rdf import *

#class ReducedVertPos(Observation) for 3D?

## Processors
class IDList(list):
  def __init__(self):
    self.id = ""
def interpCurves(data, delta):
  geod = Geod(ellps='WGS84')
  trans = Transformer.from_crs(data.projection.horizontal, 'wgs84')
  def inner(subdata):
    for geom in subdata.geom:
      segments = []
      for segment in geom.segments:
        if not isinstance(segment, Curve):
          segments.append(segment)
          continue
        interp = interpCurve(segment.start, segment.mid, segment.end)

        long0, lat0 = trans.transform(*segment.start)
        long1, lat1 = trans.transform(*segment.end)
        _a, _b, dist = geod.inv(long0, lat0, long1, lat1)
        subdivisions = ceil(dist/delta)

        subsegs = IDList()
        prev = interp(0)
        for t in range(subdivisions):
          end = interp((t+1)/subdivisions)
          subsegs.append(Line(t, Point(None, None, None, prev[1], prev[0]), Point(None, None, None, end[1], end[0])))
          prev = end
        subsegs.id = segment.id # Attach extra property to preserve original ID.
        segments.append(subsegs)
      geom.segments = segments

  for parcel in data.parcels: inner(parcel)
  for feature in data.features: inner(feature)

def interpCurve(a, m, b):
  b_m = b.easting - m.easting, b.northing - m.northing
  m_a = m.easting - a.easting, m.northing - a.northing
  ab_m = a.easting*b_m[0] - a.northing*b_m[1], a.easting*b_m[1] + a.northing*b_m[0]
  bm_a = b.easting*m_a[0] - b.northing*m_a[1], b.easting*m_a[1] + b.northing*m_a[0]

  def inner(t):
    num = ab_m[0]*(1-t) + bm_a[0]*t, ab_m[1]*(1-t) + bm_a[1]*t
    den = b_m[0]*(1-t) + m_a[0]*t, b_m[1]*(1-t) + m_a[1]*t
    dist2 = den[0]*den[0] + den[1]*den[1]
    return (num[0]*den[0] + num[1]*den[1])/dist2, (num[1]*den[0] - num[0]*den[1])/dist2
  return inner

## Commandline interface

if __name__ == "__main__":
  import argparse

  argparser = argparse.ArgumentParser()
  argparser.add_argument('-X', '--LANDXML', help="LandXML input file", nargs='+')
  argparser.add_argument('-V', '--vocab', help="Jurisdictional vocabulary to align to")
  argparser.add_argument('-C', '--CSDM', help="JSON-Topology input file", nargs='+')
  argparser.add_argument('--interpolate', help="Interpolate curves, possibly specifying length in meters of each segment (default 1m)",
      const=1, type=int, nargs='?')
  argparser.add_argument('--epsg', help="Overwrite the projection being used.")
  argparser.add_argument('-j', '--jsonfg', help="JSONfg output file", const="?.json", nargs='?')
  argparser.add_argument('-r', '--rdf', help="RDF syntax to output", const="ttl", nargs='?')
  argparser.add_argument('-o', '--output', help="RDF file to output to", const="?.rdf", nargs='?')
  argparser.add_argument('-c', '--csdm', help="CSDM JSON output file", const="?.json", nargs='?')
  argparser.add_argument('-s', '--summary', help="Textual summary of parsed data", const="?.txt", nargs='?')
  argparser.add_argument('-e', '--errorlog', help="Log of any encountered errors", const="?.log", nargs='?')
  args = argparser.parse_args()

  wrote_output = False
  def export(infile, data):
    global wrote_output
    from os import path
    filename = path.splitext(path.basename(infile))[0]

    if args.epsg: data.projection = Projection(args.epsg)
    if args.interpolate:
      interpCurves(data, args.interpolate)

    if args.jsonfg:
      with open(args.jsonfg.replace("?", filename), "w") as f: exportJSONfg(data, f)
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

  if not read_input:
    print("Please specify at least one input file!")
    argparser.print_help()
