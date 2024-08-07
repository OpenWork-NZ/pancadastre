from data import *
from pyproj import Transformer
from pyproj.crs import CompoundCRS
import json

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
    if all(isinstance(face, Face) for geom in parcel.geom for face in geom.segments):
      return {
        'type': 'Polyhedron',
        'coordinates': [
          [transform(line.start, trans) for line in face]
          for geom in parcel.geom
          for segment in geom.segments
          for face in segment.faces]
      }
    else:
      return {
        'type': 'Polygon',
        'coordinates': [simplifyPoly(trans, Geom.flatten(geom.segments)) for geom in parcel.geom]
      }
  def exportObs(obs, proj):
    trans = Transformer.from_crs(fileproj, proj)
    obs.populateProperties() # Make them more legible!
    if isinstance(obs, ReducedObservation) or (
        isinstance(obs, ReducedArcObservation) and obs.geom is None) or (
        isinstance(obs, ArcByChord) and obs.geom is None):
      return {
          'type': "LineString",
          'featureType': 'observation',
          'coordinates': [transform(obs.setupPoint, trans),
              transform(obs.targetPoint, trans)],
        }
    elif (isinstance(obs, ReducedArcObservation) or isinstance(obs, ArcByChord)) and obs.geom is not None:
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
    elif isinstance(obs, ArcByChord):
      return
    else:
      print("Unexpected observation type!", type(obs))
  fileproj = data.projection.horizontal
  if fileproj[:5] == "epsg:": fileproj = fileproj[5:]
  vertproj = data.projection.vertical
  if vertproj and vertproj != data.projection.horizontal:
    if vertproj[:5] == "epsg:": vertproj = vertproj[5:]
    fileproj = CompoundCRS(data.projection.horizontal + ';' + vertproj, components=[fileproj, vertproj])

    dest = CompoundCRS('wgs84;' + vertproj, components=['wgs84', vertproj])
    trans = Transformer.from_crs(fileproj, dest)
  else:
    trans = Transformer.from_crs(fileproj, 'wgs84')
  def transform(pt, transformer = trans):
    crs = transformer.target_crs
    if crs.name.lower().startswith("wgs 84"):
      ret = list(transformer.transform(pt.coord1, pt.coord2, pt.coord3))
      return [ret[1], ret[0], ret[2]] if len(ret) == 3 else [ret[1], ret[0]]
    else:
      return list(transformer.transform(pt.coord1, pt.coord2, pt.coord3))
 
  ret = {
    '$schema': 'schema.json',
    'type': 'FeatureCollection',
    'coordRefSys': data.projection.horizontal,
    'verticalCRS': data.projection.vertical,
    'id': getattr(data.survey.metadata, 'objID', None),
    'name': data.survey.metadata.name,
    'purpose': data.survey.metadata.purpose,
    'featureType': "CSD",
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
          'coordinates': list(monument.point) if monument.point is not None else []
        },
        'properties': monument.properties
      } for i, monument in enumerate(data.monuments)] + [{
        'id': observation.name or ('obs' + str(i)),
        'type': "Feature",
        'featureType': "sosa:ObservationCollection",
        'geometry': exportObs(observation, 'wgs84'),
        'place': exportObs(observation, fileproj),
        'properties': getattr(observation, 'properties', {})
      } for groupID, group in data.survey.observationGroups.items() for i, observation in enumerate(group)] + [{
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
        #elif line.start not in list(lines[i-1]):
        #  print("WARNING: Non-contiguous line segments (last)!", line, lines)
        # Otherwise, it's fine!
      else:
        if line.start in list(lines[i+1]):
          line.start, line.end = line.end, line.start
        #elif line.end not in list(lines[i+1]):
        #  print("WARNING: Non-contiguous line segments!", line, lines)
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
    if all(isinstance(face, Face) for geom in parcel.geom for face in geom.segments):
      return {
        'type': 'Polyhedron',
        'coordinates': [
          [transform(line.start, trans) for line in face]
          for geom in parcel.geom
          for segment in geom.segments
          for face in segment.faces]
      }
    else:
      return {
        'type': 'Polygon',
        'coordinates': [simplifyPoly(trans, Geom.flatten(geom.segments)) for geom in parcel.geom]
      }
  fileproj = data.projection.horizontal
  if fileproj[:5] == "epsg:": fileproj = fileproj[5:]
  vertproj = data.projection.vertical
  if vertproj and vertproj != data.projection.horizontal:
    if vertproj[:5] == "epsg:": vertproj = vertproj[5:]
    fileproj = CompoundCRS(data.projection.horizontal + ';' + vertproj, components=[fileproj, vertproj])

    dest = CompoundCRS('wgs84;' + vertproj, components=['wgs84', vertproj])
    trans = Transformer.from_crs(fileproj, dest)
  else:
    trans = Transformer.from_crs(fileproj, 'wgs84')
  def transform(pt, transformer = trans):
    crs = transformer.target_crs
    if crs.name.lower().startswith("wgs 84"):
      ret = list(transformer.transform(pt.coord1, pt.coord2, pt.coord3))
      return [ret[1], ret[0], ret[2]] if len(ret) == 3 else [ret[1], ret[0]]
    else:
      return list(transformer.transform(pt.coord1, pt.coord2, pt.coord3))

  ret = {
    '$schema': 'schema.json',
    'type': 'FeatureCollection',
    'horizontalCRS': data.projection.horizontal, # What to do with vertical?
    'id': getattr(data.survey.metadata, 'objID', None),
    'name': data.survey.metadata.name,
    'purpose': data.survey.metadata.purpose,
    'surveyType': data.survey.metadata.type,
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
          'coordinates': list(monument.point) if monument.point is not None else []
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
      } for i, parcel in enumerate(data.parcels)] + [{
        'id': tins.id or ('tin'+str(i)),
        'type': "FeatureCollection",
        'featureType': tins.type,
        'properties': tins.properties,
        'features': [{
          'id': tin.id or ('TIN'+str(j)),
          'type': "feature",
          'geometry': {
            'type': "TIN",
            'coordinates': [
                [transform(seg.start) for segs in face.faces for seg in segs]
                for faces in tin.faces for face in faces]
          },
          'place': {
            'type': "TIN",
            'coordinates': [
                [list(seg.start) for segs in face.faces for seg in segs]
                for faces in tin.faces for face in faces]
          }
        } for j, tin in enumerate(tins.features)]
      } for i, tins in enumerate(data.surfaceTINs)] + [{
        'id': curves.id or ("tic" + str(i)),
        'type': "FeatureCollection",
        'featureType': curves.type,
        'properties': curves.properties,
        'features': [{
          'id': curve.id or ('TIC' + str(j)),
          'type': "feature",
          'geometry': {
            'type': "MultiLineString",
            'coordinates': simplifyPoly(trans, curve.topology)
          },
          'place': {
            'type': "MultiLineString",
            'coordinates': simplifyPoly(Transformer.from_crs(fileproj, fileproj), curve.topology)
          }
        } for j, curve in enumerate(curves.features)]
      } for i, curves in enumerate(data.terrainIntersectionCurves)]
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
