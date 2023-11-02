from data import *
from pyproj import Transformer

# TODO: Tweak according to surveyfeatures.json, circular_arc.json
def exportJSONfg(data, file = None):
  import json
  def simplifyPoly(trans, lines):
    if len(lines) == 2:
      return [list(trans.transform(*lines[0].start)), list(trans.transform(*lines[0].end)), list(trans.transform(*lines[1].end))]
    elif len(lines) == 1:
      return [list(trans.transform(*lines[0].start)), list(trans.transform(*lines[0].end))]
    elif len(lines) == 0:
      return []

    lines = list(lines)
    for i, line in enumerate(lines):
      if i == len(lines) - 1:
        if line.end in list(lines[i-1]):
          line.start, line.end = line.end, line.start
        elif line.start not in list(lines[i-1]):
          print("Non-contiguous line segments!", line, lines)
        # Otherwise, it's fine!
      else:
        if line.start in list(lines[i+1]):
          line.start, line.end = line.end, line.start
        elif line.end not in list(lines[i+1]):
          print("Non-contiguous line segments!", line, lines)
        # Otherwise, it's fine!

    points = []
    for line in lines:
      points.append(trans.transform(*line.start))
      points.append(trans.transform(*line.end))
    i = 1
    while i < len(points):
      if points[i] == points[i-1]: points.pop(i)
      else: i += 1
    return list(map(list, points))

  def exportGeom(parcel, proj):
    trans = Transformer.from_crs(fileproj, proj, always_xy = proj.lower() == 'wgs84')
    return {
      'type': 'Polygon',
      'coordinates': [simplifyPoly(trans, Geom.flatten(geom.segments)) for geom in parcel.geom]
    }
  def exportObs(obs, proj):
    trans = Transformer.from_crs(fileproj, proj, always_xy = proj.lower() == 'wgs84')
    if isinstance(obs, ReducedObservation) or (
        isinstance(obs, ReducedArcObservation) and obs.geom is None):
      return {
          'type': "LineString",
          'featureType': 'observation',
          'coordinates': [trans.transform(*obs.setupPoint), trans.transform(*obs.targetPoint)],
          'properties': {
            # What should go in here?
          }
        }
    elif isinstance(obs, ReducedArcObservation) and obs.geom is not None:
      return {
        'type': "LineString",
        'featureType': 'observation',
        'coordinates': simplifyPoly(trans, Geom.flatten([obs.geom])),
        'properties': {
          # What should go in here?
        }
      }
    elif isinstance(obs, RedHorizPos):
      return {
        'type': "Point",
        'featureType': 'observation',
        'coordinates': list(trans.transform(obs.northing, obs.easting)),
        'properties': {
          # What should go in here?
        }
      }
    else:
      print("Unexpected observation type!", type(obs))
  fileproj = data.projection.horizontal
  if fileproj[:5] == "epsg:": fileproj = fileproj[5:]
  trans = Transformer.from_crs(fileproj, 'wgs84', always_xy = True)
  
  lines = set()
  observedVecs = []
  for i, seg in enumerate(seg for parcel in data.parcels for geom in parcel.geom for seg in geom.segments):
    if isinstance(seg, list):
      x = {
        'id': getattr(seg, "id", "curve" + str(i)),
        'type': "Feature",
        # Is this the VectorPurpose value?
        'featureType': 'boundary', # FIXME: This will need to support other feature types. What's the logic here?
        'geometry': {
          'type': "LineString",
          'coordinates': simplifyPoly(trans, seg),
          'properties': {'comment': ""}
        },
        'place': {
          'type': "LineString",
          'coordinates': simplifyPoly(Transformer.from_crs(fileproj, fileproj), seg),
          'properties': {'comment': ""}
        }
      }
    else:
      if seg in lines: continue
      lines.add(seg)
      x = {
        'id': seg.id or i,
        'type': "Feature",
        # Is this the VectorPurpose value?
        'featureType': 'boundary', # FIXME: This will need to support other feature types. What's the logic here?
        'geometry': {
          'type': "LineString",
          'coordinates': [trans.transform(*seg.start), trans.transform(*seg.end)],
          'properties': {
            'comment': ""
           }
        },
        'place': {
          'type': "LineString",
          'coordinates': [list(seg.start), list(seg.end)],
          'properties': {
            'comment': ""
          }
        }
      }
    observedVecs.append(x)
 
  ret = {
    '$schema': 'schema.json',
    'type': 'FeatureCollection',
    'horizontalCRS': data.projection.horizontal, # What to do with vertical?
    'id': getattr(data.survey.metadata, 'objID', None),
    'name': data.survey.metadata.name,
    'purpose': data.survey.metadata.purpose,
    # TODO: Capture data to transfer
    'links': [],
    'provenance': {},
    'features': [{
        'id': getattr(monument.point, 'objID', None) or i,
        'type': "Feature",
        'featureType': "SurveyPoint", # Either CadastralMark, Boundary, or Geodetic
        'geometry': {
          'type': "Point",
          'coordinates': list(trans.transform(*monument.point)) if monument.point is not None else []
        },
        'place': {
          'type': "Point",
          'coordinates': [monument.point.northing, monument.point.easting] if monument.point is not None else []
        },
        'properties': monument.properties
      } for i, monument in enumerate(data.monuments)] + observedVecs + [{
        'type': 'Feature',
        'id': parcel.name or i,
        'featureType': 'boundary',
        'links': [],
        'time': '', # Appears to be missing from LandXML
        'place': exportGeom(parcel, 'wgs84'),
        'geometry': exportGeom(parcel, fileproj),
        'properties': parcel.properties
      } for i, parcel in enumerate(data.parcels)] + [{
        'type': "Feature",
        'featureType': "sosa:ObservationCollection",
        'geometry': exportObs(observation, 'wgs84'),
        'place': exportObs(observation, fileproj)
      } for groupID, group in data.survey.observationGroups.items() for i, observation in enumerate(group)]
  }
  if file is not None: json.dump(ret, file, indent=4)
  else: return ret
