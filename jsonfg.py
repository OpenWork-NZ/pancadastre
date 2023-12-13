from data import *
from pyproj import Transformer
import json

# TODO: Tweak according to surveyfeatures.json, circular_arc.json
def exportJSONfg(data, file = None, isGeoJSON = False):
  def simplifyPoly(trans, lines):
    lines = [line for line in lines if line is not None]
    if len(lines) == 2:
      return [transform(lines[0].start, trans),
          transform(lines[0].end, trans), transform(lines[1].end, trans)]
    elif len(lines) == 1:
      return [transform(lines[0].start, trans), transform(lines[0].end, trans)]
    elif len(lines) == 0:
      return []

    lines = list(lines)
    for i, line in enumerate(lines):
      if i == len(lines) - 1:
        if line.end in list(lines[i-1]):
          line.start, line.end = line.end, line.start
        elif line.start not in list(lines[i-1]):
          print("WARNING: Non-contiguous line segments (last)!", line, lines)
        # Otherwise, it's fine!
      else:
        if line.start in list(lines[i+1]):
          line.start, line.end = line.end, line.start
        elif line.end not in list(lines[i+1]):
          print("WARNING: Non-contiguous line segments!", line, lines)
        # Otherwise, it's fine!

    points = []
    for line in lines:
      points.append(transform(line.start, trans))
      points.append(transform(line.end, trans))
    i = 1
    while i < len(points):
      if points[i] == points[i-1]: points.pop(i)
      else: i += 1
    return list(points)

  def exportGeom(parcel, proj):
    trans = Transformer.from_crs(fileproj, proj)
    return {
      'type': 'Polygon',
      'coordinates': [simplifyPoly(trans, Geom.flatten(geom.segments)) for geom in parcel.geom]
    }
  def exportObs(obs, proj):
    trans = Transformer.from_crs(fileproj, proj)
    obs.populateProperties() # Make them more legible!
    if isinstance(obs, ReducedObservation) or (
        isinstance(obs, ReducedArcObservation) and obs.geom is None):
      return {
          'type': "LineString",
          'featureType': 'observation',
          'coordinates': [transform(obs.setupPoint, trans),
              transform(obs.targetPoint, trans)],
        }
    elif isinstance(obs, ReducedArcObservation) and obs.geom is not None:
      return {
        'type': "LineString",
        'featureType': 'observation',
        'coordinates': simplifyPoly(trans, Geom.flatten([obs.geom])),
      }
    elif isinstance(obs, RedHorizPos):
      return {
        'type': "Point",
        'featureType': 'observation',
        'coordinates': transform(obs, trans),
      }
    elif isinstance(obs, SubtendedAngle):
      return {
        #'type': "Point",
        #'featureType': 'subtendedAngle',
        #'coordinates': transform(obs.setupPt, trans) # Redundant info!
      }
    elif isinstance(obs, CircleByCenter):
      obs.properties["radius"] = obs.radius
      return {
        #'type': "Point",
        #'featureType': 'circle',
        #'coordinates': transform(obs.centerSetup.point, trans) # Redundant info!
      }
    elif isinstance(obs, CubicSplineObservation):
      return {
        'type': "LineString",
        'featureType': 'cubicSpline',
        'coordinates': [transform(pt, trans) for pt in obs.interpolatedPath()]
      }
    else:
      print("Unexpected observation type!", type(obs))
  fileproj = data.projection.horizontal
  if fileproj[:5] == "epsg:": fileproj = fileproj[5:]
  trans = Transformer.from_crs(fileproj, 'wgs84')
  def transform(pt, transformer = trans):
    crs = transformer.target_crs
    if crs.name.lower().startswith("wgs 84"):
      return list(reversed(transformer.transform(pt.coord1, pt.coord2)))
    else:
      return list(transformer.transform(pt.coord1, pt.coord2))
 
  ret = {
    '$schema': 'schema.json',
    'type': 'FeatureCollection',
    'horizontalCRS': data.projection.horizontal, # What to do with vertical?
    'id': getattr(data.survey.metadata, 'objID', None),
    'name': data.survey.metadata.name,
    'purpose': data.survey.metadata.purpose,
    'referencedCSDs': [{
        'id': ref.id,
        'name': ref.name,
        'adminUnit': {
          'href': ref.href,
          'rel': ref.rel,
          'role': ref.role
        },
        'bearingRotation': ref.bearing,
        'time': str(ref.time)
      } for ref in data.referencedCSDs],
    # TODO: Capture data to transfer
    'links': [],
    'provenance': {},
    'features': [{
        'id': getattr(monument.point, 'objID', None) or i,
        'type': "Feature",
        'featureType': "SurveyPoint", # Either CadastralMark, Boundary, or Geodetic
        'geometry': {
          'type': "Point",
          'coordinates': transform(monument.point) if monument.point is not None else []
        },
        'place': {
          'type': "Point",
          'coordinates': [monument.point.coord1, monument.point.coord2] if monument.point is not None else []
        },
        'properties': monument.properties
      } for i, monument in enumerate(data.monuments)] + [{
        'id': observation.name or 'obs' + str(i),
        'type': "Feature",
        'featureType': "sosa:ObservationCollection",
        'geometry': exportObs(observation, 'wgs84'),
        'place': exportObs(observation, fileproj),
        'properties': getattr(observation, 'properties', {})
      } for groupID, group in data.survey.observationGroups.items() for i, observation in enumerate(group)]
  }
  if file is not None: json.dump(ret, file, indent=4)
  else: return ret

def exportJSONfgParcel(data, file = None, isGeoJSON = False):
  def simplifyPoly(trans, lines):
    lines = [line for line in lines if line is not None]
    if len(lines) == 2:
      return [transform(lines[0].start, trans),
          transform(lines[0].end, trans), transform(lines[1].end, trans)]
    elif len(lines) == 1:
      return [transform(lines[0].start, trans), transform(lines[0].end, trans)]
    elif len(lines) == 0:
      return []

    lines = list(lines)
    for i, line in enumerate(lines):
      if i == len(lines) - 1:
        if line.end in list(lines[i-1]):
          line.start, line.end = line.end, line.start
        elif line.start not in list(lines[i-1]):
          print("WARNING: Non-contiguous line segments (last)!", line, lines)
        # Otherwise, it's fine!
      else:
        if line.start in list(lines[i+1]):
          line.start, line.end = line.end, line.start
        elif line.end not in list(lines[i+1]):
          print("WARNING: Non-contiguous line segments!", line, lines)
        # Otherwise, it's fine!

    points = []
    for line in lines:
      points.append(transform(line.start, trans))
      points.append(transform(line.end, trans))
    i = 1
    while i < len(points):
      if points[i] == points[i-1]: points.pop(i)
      else: i += 1
    return list(points)

  def exportGeom(parcel, proj):
    trans = Transformer.from_crs(fileproj, proj)
    return {
      'type': 'Polygon',
      'coordinates': [simplifyPoly(trans, Geom.flatten(geom.segments)) for geom in parcel.geom]
    }
  def exportObs(obs, proj):
    trans = Transformer.from_crs(fileproj, proj)
    obs.populateProperties() # Make them more legible!
    if isinstance(obs, ReducedObservation) or (
        isinstance(obs, ReducedArcObservation) and obs.geom is None):
      return {
          'type': "LineString",
          'featureType': 'observation',
          'coordinates': [transform(obs.setupPoint, trans),
              transform(obs.targetPoint, trans)],
        }
    elif isinstance(obs, ReducedArcObservation) and obs.geom is not None:
      return {
        'type': "LineString",
        'featureType': 'observation',
        'coordinates': simplifyPoly(trans, Geom.flatten([obs.geom])),
      }
    elif isinstance(obs, RedHorizPos):
      return {
        'type': "Point",
        'featureType': 'observation',
        'coordinates': transform(obs, trans),
      }
    elif isinstance(obs, SubtendedAngle):
      return {
        'type': "Point",
        'featureType': 'subtendedAngle',
        'coordinates': transform(obs.setupPt, trans)
      }
    elif isinstance(obs, CircleByCenter):
      obs.properties["radius"] = obs.radius
      return {
        'type': "Point",
        'featureType': 'circle',
        'coordinates': transform(obs.centerSetup.point, trans)
      }
    elif isinstance(obs, CubicSplineObservation):
      return {
        'type': "LineString",
        'featureType': 'cubicSpline',
        'coordinates': [transform(pt, trans) for pt in obs.interpolatedPath()]
      }
    else:
      print("Unexpected observation type!", type(obs))
  fileproj = data.projection.horizontal
  if fileproj[:5] == "epsg:": fileproj = fileproj[5:]
  trans = Transformer.from_crs(fileproj, 'wgs84')
  def transform(pt, transformer = trans):
    crs = transformer.target_crs
    if crs.name.lower().startswith("wgs 84"):
      return list(reversed(transformer.transform(pt.coord1, pt.coord2)))
    else:
      return list(transformer.transform(pt.coord1, pt.coord2))
 
  ret = {
    '$schema': 'schema.json',
    'type': 'FeatureCollection',
    'horizontalCRS': data.projection.horizontal, # What to do with vertical?
    'id': getattr(data.survey.metadata, 'objID', None),
    'name': data.survey.metadata.name,
    'purpose': data.survey.metadata.purpose,
    'referencedCSDs': [{
        'id': ref.id,
        'name': ref.name,
        'adminUnit': {
          'href': ref.href,
          'rel': ref.rel,
          'role': ref.role
        },
        'bearingRotation': ref.bearing,
        'time': str(ref.time)
      } for ref in data.referencedCSDs],
    # TODO: Capture data to transfer
    'links': [],
    'provenance': {},
    'features': [{
        'id': getattr(monument.point, 'objID', None) or i,
        'type': "Feature",
        'featureType': "SurveyPoint", # Either CadastralMark, Boundary, or Geodetic
        'geometry': {
          'type': "Point",
          'coordinates': transform(monument.point) if monument.point is not None else []
        },
        'place': {
          'type': "Point",
          'coordinates': [monument.point.coord1, monument.point.coord2] if monument.point is not None else []
        },
        'properties': monument.properties
      } for i, monument in enumerate(data.monuments)] + [{
        'type': 'Feature',
        'id': parcel.oid or i,
        'featureType': 'parcel',
        'links': [],
        'time': '', # Appears to be missing from LandXML
        'geometry': exportGeom(parcel, 'wgs84'),
        'place': exportGeom(parcel, fileproj),
        'properties': parcel.properties
      } for i, parcel in enumerate(data.parcels)]
  }
  if file is not None: json.dump(ret, file, indent=4)
  else: return ret

def exportGeoJSON(data, file):
  data = exportJSONfg(data, isGeoJSON = True)
  del data["horizontalCRS"]
  for datum in data["features"]:
    del datum["place"]
  json.dump(data, file, indent=4)

def exportGeoJSONParcel(data, file):
  data = exportJSONfgParcel(data, isGeoJSON = True)
  del data["horizontalCRS"]
  for datum in data["features"]:
    del datum["place"]
  json.dump(data, file, indent=4)
