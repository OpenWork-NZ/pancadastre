import rdflib
from rdflib.namespace import RDFS, RDF, GEO
from data import *

def exportWKTGeom(data):
  def pt(point):
    # Latitude,Longitude order
    return str(pt.coord1) + " " + str(pt.coord2)
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
    g.add((oGeom, geom.pos, rdflib.Literal(str(point.coord1) + " " + str(point.coord2))))
    g.add((oGeom, geo.asWKT, rdflib.Literal("POINT( " + str(point.coord1) + " " + str(point.coord2) + " )")))
    
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
    g.add((oGeom, geom.pos, rdflib.Literal(str(point.coord1) + " " + str(point.coord2))))
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
