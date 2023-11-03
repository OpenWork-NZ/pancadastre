import rdflib
from rdflib.namespace import RDFS, RDF, GEO
from data import *

def exportWKTGeom(data):
  def pt(point):
    # Latitude,Longitude order
    return str(pt.northing) + " " + str(pt.easting)
  rets = []
  for seg in segments:
    if isinstance(seg, Line):
      rets.append("LINESTRING( " + pt(seg.start) + ", " + pt(seg.end) + " )")
    elif isinstance(seg, Curve):
      # I'm confused by the spec: There's discussion about curves, but I don't see where the syntax is specified.
      rets.append("CURVE( " + ptr(seg.start) + ", " + pt(seg.mid) + ", " + pt(seg.end) + " )")
    else:
      raise Exception("Unknown geom type: " ++ str(type(seg)))
  return "GEOMETRYCOLLECTION( " + ", ".join(rets) + " )"

def exportRDF(data, file, format="turtle"):
  g = rdflib.Graph()
  surv = rdflib.Namespace("https://data.surroundaustralia.com/def/csdm/surveyfeats/")
  geom = rdflib.Namespace("https://data.surroundaustralia.com/def/csdm/geometryprim/")
  commonpatterns = rdflib.Namespace("https://data.surroundaustralia.com/def/csdm/commonpatterns/")
  geo = rdflib.Namespace("http://www.opengis.net/ont/geosparql#")
  self = rdflib.Namespace(".")

  for i, point in enumerate(data.points):
    oPt = self["pt" + str(i)]
    oGeom = rdflib.BNode()
    g.add((oPt, RDF.type, surv.SurveyPoint))
    g.add((oPt, GEO.hasGeometry, oGeom))
    g.add((oGeom, RDF.type, geom.Point))
    # FIXME: Translate this literal to a URL... Next 3 properties have been remodelled.
    g.add((oGeom, geom.srsName, rdflib.Literal(point.projection or data.projection.horizontal)))
    g.add((oGeom, geom.srsDimension, rdflib.Literal(2)))
    g.add((oGeom, geom.pos, rdflib.Literal(str(point.northing) + " " + str(point.easting))))
    g.add((oGeom, geo.asWKT, rdflib.Literal("POINT( " + str(point.northing) + " " + str(point.easting) + " )")))
    
    g.add((oPt, surv.hasPointQuality, rdflib.Literal(point.observation.distanceQuality)))
    g.add((oPt, surv.hasProvenance, rdflib.Literal(point.state)))
    g.add((oPt, surv.hasPurpose, rdflib.Literal(point.purpose)))
    g.add((oPt, RDFS.label, rdflib.Literal(point.name)))

  for i, marker in enumerate(data.monuments):
    oMrk = self["m" + str(i)]
    oGeom = rdflib.BNode()
    pt = marker.point
    g.add((oMrk, RDF.type, surv.CadastralMark))
    g.add((oMrk, GEO.hasGeometry, oGeom))
    g.add((oGeom, RDF.type, geom.Point))
    g.add((oGeom, geom.srsName, rdflib.Literal(pt.projection or data.projection.horizontal)))
    g.add((oGeom, geom.srsDimension, rdflib.Literal(2)))
    g.add((oGeom, geom.pos, rdflib.Literal(str(point.northing) + " " + str(point.easting))))
    g.add(())

    g.add((oGeom, surv.hasPointQuality, rdflib.Literal(point.observation.distanceQuality)))
    g.add((oGeom, surv.hasProvenance, rdflib.Literal(point.state)))
    g.add((oGeom, surv.hasPurpose, rdflib.Literal(point.purpose)))
    
    g.add((oMrk, surv.monumentForm, rdfLib.Literal(marker.type)))
    g.add((oMrk, surv.monumentCondition, rdfLib.Literal(marker.condition)))
    g.add((oMrk, surv.markLabel, rdfLib.Literal(marker.name)))
    # surv:hasComment?
    
    g.add((oMrk, RDFS.label, marker.name))
    # TODO: What boundary marks? How'm I storing that?

  # FIXME: This respect seems to have been redesigned...
  for group, observations in data.survey.observationGroups.items():
    # TODO: Additional observation collection properties?
    oGrp = self[group]
    g.add((oGrp, RDF.type, surv.vectorObservationCollection))
    for i, observation in enumerate(observations):
      oObs = self[group + "-" + str(i)]
      g.add((oGrp, surv.inputRawObservation, oObs))
      g.add((oObs, RDF.type, surv.vectorDetermination))

      oDist = rdflib.BNode()
      g.add((oObs, surv.hasDistance, oDist))
      g.add((oDist, RDF.type, surv.Distance))
      # FIXME: Quality measure?
      g.add((oDist, commonpatterns.value, rdflib.Literal(observation.horizDistance)))
      g.add(oObs, surv.distType, rdflib.Literal(observation.distType)) # FIXME: Use URIs instead of literals.

      oHNgl = rdflib.BNode()
      g.add((oObs, surv.hasHorizontalAngle, oHNgl))
      g.add((oHNgl, RDF.type, surv.Angle))
      g.add((oHNgl, commonpatterns.value, rdflib.Literal(obs.azimuth))) # Is this right?

      # Vertical angle?
      g.add((oObs, surv.resultTime, rdflib.Literal(obs.date)))
      g.add((oObs, surv.instrumentHeight, rdflib.Literal(obs.setup.height)))
      g.add((oObs, surv.targetHeight, rdflib.Literal(obs.targetSetup.height)))

  # Surveyed vectors?

  g.serialize(file, format=format)

def exportSummary(data):
  ret = ""

  if data.projection.horizontal == data.projection.vertical:
    ret += "Projection:\t" + str(data.projection.horizontal) + "\n"
  else:
    ret += "Projection:\t" + str(data.projection.horizontal) + ", " + str(data.projection.vertical) + "\n"
  if data.survey.metadata.name is not None:   ret += "Name:\t" + str(data.survey.metadata.name) + "\n"
  if data.survey.metadata.firm is not None:   ret += "Firm:\t" + str(data.survey.metadata.firm) + "\n"
  if data.survey.metadata.surveyor is not None:ret+= "Surveyor:\t" + str(data.survey.metadata.surveyor) + "\n"
  if data.survey.metadata.type is not None:   ret += "Type:\t" + str(data.survey.metadata.type) + "\n"
  if data.survey.metadata.jurisdiction is not None:
    ret += "Jurisdiction:\t" + str(data.survey.metadata.jurisdiction) + "\n"
  if data.survey.metadata.format is not None: ret += "Format:\t" + str(data.survey.metadata.format) + "\n"
  if data.survey.metadata.headOfPower is not None:
    ret += "Head of Power:\t" + str(data.survey.metadata.headOfPower) + "\n"
  ret += "\n"

  ret += "Monuments:\t" + str(len(data.monuments)) + "\n"
  for monument in data.monuments or []:
    ret += "\t" + str(monument.name) + "\t" + str(monument.point) + "\t" + str(monument.state) + "\t" + str(monument.condition) + "\n"
  ret += "Points:\t" + str(len(data.points)) + "\n"
  for point in data.points or []:
    ret += "\t" + str(point.name) + "\t" + str(point.state) + "\t" + str(point) + "\n"
  ret += "Parcels:\t" + str(len(data.parcels)) + "\n"
  for parcel in data.parcels or []:
    ret += "\t" + str(parcel.name) + "\t" + str(parcel.type) + "\t" + str(parcel.area) + "\t" + str(parcel.center) + "\t" + \
        str(parcel.state) + "\t" + str(parcel.klass) + "\t" + str(parcel.desc) + "\n"
  ret += "Instruments:\t" + str(len(data.survey.instruments)) + "\n"
  for instrument in data.survey.instruments or []:
    ret += "\t" + str(instrument.stationName or instrument.id) + "\t" + str(instrument.height) + "\t" + str(instrument.point) + "\n"
  ret += "Observations:\t" + str(len([y for x in data.survey.observationGroups.values() for y in x])) + "\n"
  for label, observations in data.survey.observationGroups.items():
    ret += "\t\t" + label + "\n"
    for observation in observations:
      ret += "\t" + str(observation.name) + "\t" + str(observation.date) + "\t" + str(observation.purpose) + "\t"
      if isinstance(observation, ReducedObservation):
        ret += str(observation.setup.stationName or observation.setup.id) + "\t"
        ret += str(observation.targetSetup.stationName or observation.targetSetup.id) + "\t"
        ret += str(observation.azimuth) + "\t" + str(observation.horizDist) + "\n"
      elif isinstance(observation, ReducedArcObservation):
        ret += str(observation.setup.stationName or observation.setup.id) + "\t"
        ret += str(observation.targetSetup.stationName or observation.targetSetup.id) + "\t"
        ret += str(observation.chordAzimuth) + "\t" + str(observation.length) + "\t" + str(observation.radius) + "\n"
      elif isinstance(observation, RedHorizPos):
        ret += str(observation.northing) + "," + str(observation.easting) + "\n"
      else:
        ret += "Unsupported observation type: " + repr(type(observation)) + "\n"

  return ret
