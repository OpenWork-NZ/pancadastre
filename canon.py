from data import *
import json
from jsonfg import * # I question this...

def importCSDM(file):
  from copy import deepcopy
  data = json.load(file)
  assert data["type"] == "FeatureCollection"
  projection = Projection(data.get("horizontalCRS"))
  metadata = SurveyMetadata(data["name"], None, None, None, None, None, data["purpose"], None, None, None, None)
  
  instruments = []
  def instrument(*args):
    ret = InstrumentSetup(*args)
    instruments.append(ret)
    return ret
  
  monuments = []
  points = []
  observations = []
  measures = {}

  for group in data["vectorObservations"]:
    for observation in group["features"]:
      measure = Measure.fromProperties(group["properties"], observation["properties"])
      measures[observation["properties"]["hasFeatureOfInterest"]] = measure

  pointsIndex = {}
  for collection in data["points"]:
    for monument in collection["features"]:
      monumentedBy = monument["properties"]["monumentedBy"]
      if "state" not in monumentedBy or "form" not in monumentedBy or "condition" not in monumentedBy:
        if "monumentState" in monumentedBy:
          print("Outdated schema! Key 'monumentState' should be 'state'!")
        else:
          print("Invalid schema! Not attaching monument info!")
      point = Point(None, None, monumentedBy.get("state", monumentedBy.get("monumentState")), monument["place"]["coordinates"][0], monument["place"]["coordinates"][1], monument["id"])
      points.append(point)
      pointsIndex[monument["id"]] = deepcopy(monument)
      pointsIndex[monument["id"]][""] = point
      # The observations will be initialized later
      monuments.append(Monument((monument["properties"].get("name") or {}).get("label", monument["id"]), point, monumentedBy.get("state", monumentedBy.get("monumentState")), monumentedBy.get("form", monumentedBy.get("monumentForm")), monumentedBy.get("condition", monumentedBy.get("monumentCondition")), properties = monument["properties"]))

  def n(comment): return None # Indicates a TODO...
  def d(x):
    if isinstance(x, dict): return x["$ref"]
    else: return x

  segments = {}
  observationGroups = {} # TODO: Build!
  for group in data["observedVectors"]:
    observations = []
    for observation in group["features"]:
      # FIXME: Is it a bug if a point doesn't have a measure? Should this be reported?
      measure = measures.get(observation["id"]) or Measure.fromProperties({}, {})
      if observation["topology"]["type"].lower() == "linestring":
        for i in range(1, len(observation["topology"]["references"])):
          start = pointsIndex[d(observation["topology"]["references"][i-1])]
          end = pointsIndex[d(observation["topology"]["references"][i])]
          # Distance comes from vector observations, we might want to split that class out!
          obs = ReducedObservation(None, instrument(start["id"], (start["properties"].get("name") or {}).get("label", start["id"]), None, start[""]), instrument(end["id"], (end["properties"].get("name") or {}).get("label", end["id"]), None, end[""]), measure.azimuth, measure.dist, measure.equipment, measure.distType, measure.azimuthType, observation["id"], start["time"], start["properties"], measure = measure)
          observations.append(obs)

          geom = Line(observation["id"], start[""], end[""])
          segments[observation["id"]] = geom
      elif observation["topology"]["type"].lower() == "arc":
        start = pointsIndex[d(observation["topology"]["references"][0])]
        mid = pointsIndex[d(observation["topology"]["references"][1])]
        end = pointsIndex[d(observation["topology"]["references"][2])]
        # Distance, radius, etc comes from vector observations.
        obs = ReducedArcObservation(None, instrument(start["id"], (start["properties"].get("name") or {}).get("label", start["id"]), None, start[""]), instrument(end["id"], (end["properties"].get("name") or {}).get("label", end["id"]), None, end[""]), measure.azimuth, measure.radius, measure.dist, measure.is_clockwise, measure.equipment, measure.angleAccuracy, measure.arcType, observation["id"], start["time"], start["properties"], measure.distanceQuality, measure.distanceAccuracy, measure.angleQuality, measure = measure)
        observations.append(obs)

        geom = Curve(observation["id"], measure.is_clockwise, measure.radius, start[""], mid[""], end[""])
        segments[observation["id"]] = geom
      else:
        print("Unexpected observedVector topology-type: ", observation["topology"]["type"])
    observationGroups[group["id"]] = observations

  parcels = []
  for parcel in data.get("parcels", []):
    geoms = []
    properties = None # I'm not particularly keen on how this bit is structured.
    for geom in parcel["features"]:
      if properties is None: properties = geom["properties"]
      geoms.append(Geom(geom["id"],
        [segments.get(ref if isinstance(ref, str) else ref.get('$ref'))
          for ref in geom["topology"]["references"]]))
    parcels.append(Parcel.fromProperties(None, None, geoms, properties))
    
  return Cadastre(projection, {}, None, monuments, points, parcels, Survey(metadata, instruments, observationGroups))

def exportCSDM(data, file):
  import json
  def simplifyPoly(trans, lines):
    if len(lines) == 2:
      return [list(lines[0].start), list(lines[0].end), list(lines[1].end)]
    elif len(lines) == 1:
      return [list(lines[0].start), list(lines[0].end)]
    elif len(lines) == 0:
      return []

    lines = list(lines)
    for i, line in enumerate(lines):
      if i == len(lines) - 1:
        if line.end in list(lines[i-1]):
          line.start, line.end = line.end, line.start
        elif line.start not in list(lines[i-1]):
          print("Non-contiguous line segments!", lines)
        # Otherwise, it's fine!
      else:
        if line.start in list(lines[i+1]):
          line.start, line.end = line.end, line.start
        elif line.end not in list(lines[i+1]):
          print("Non-contiguous line segments!", lines)
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
    trans = Transformer.from_crs(fileproj, proj)
    return {
      'type': 'Linestring',
      'coordRefSys': proj,
      'coordinates': [simplifyPoly(trans, Geom.flatten(geom.segments)) for geom in parcel.geom]
    }
  fileproj = data.projection.horizontal
  if fileproj[:5] == "epsg:": fileproj = fileproj[5:]
  trans = Transformer.from_crs(fileproj, 'wgs84')
  def ref(id): return id

  lines = set()
  observedVecs = []
  for i, seg in enumerate(seg for parcel in data.parcels for geom in parcel.geom for seg in geom.segments):
    if isinstance(seg, list):
      x = {
        'id': getattr(seg, "id", "curve" + str(i)),
        'type': "Feature",
        # Is this the VectorPurpose value?
        'featureType': 'boundary', # FIXME: This will need to support other feature types. What's the logic here?
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
        'id': "observedVectors",
        'type': "FeatureCollection",
        'featureType': 'surv:ObservedVector',
        'features': [{
          'id': seg.id or i,
          'type': "Feature",
          # Is this the VectorPurpose value?
          'featureType': 'boundary', # FIXME: This will need to support other feature types. What's the logic here?
          'place': None,
          'topology': {
            'type': "LineString",
            'references': [ref(seg.start.objID), ref(seg.end.objID)], # TODO: Not all of these IDs would be dereferencable, find them!
            'properties': {
              'comment': ""
            }
          }
        }],
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
    'referencedCSDs': [], # TODO: Capture this data.
    'points': [{
      'id': 'surveymarks',
      'type': 'FeatureCollection',
      'featureType': 'csdm:SurveyMark',
      'features': [{
        'id': getattr(monument.point, 'objID', None) or i,
        'type': "Feature",
        'featureType': "SurveyPoint", # Could be "cadastralMark", etc?
        'time': monument.point.date,
        'place': {
          'type': "Point",
          'coordinates': [monument.point.northing, monument.point.easting] if monument.point is not None else []
        },
        'properties': monument.properties
      } for i, monument in enumerate(data.monuments)]
    }],
    'observedVectors': observedVecs,
    'parcels': [{
      'id': parcel.name or i,
      'type': 'FeatureCollection',
      'featureType': 'surv:PrimaryParcel',
      'properties': None,
      'features': [{
        'id': geom.name or j,
        'type': 'Feature',
        'place': None,
        'topology': {
          'type': 'Polygon',
          'references': [ref(seg.id) for seg in geom.segments]
        },
        'properties': parcel.properties
      } for j, geom in enumerate(parcel.geom)]
    } for i, parcel in enumerate(data.parcels)],
    'vectorObservations': [{
      'id': groupName,
      'type': "FeatureCollection",
      'featureType': 'sosa:ObservationCollection',
      'usedProcedure': 'surveyproc:calculation', # TODO: Fill in from obs
      'properties': {
        # FIXME: Need to remodel, where does this data come from?
#        'resultTime': '???',
#        'observedProperty': '???',
#        'baseSensor': '???',
#        'roverSensor': '???'
      },
      'features': [{
        'id': obs.name or (groupName + str(i)),
        'type': "Feature",
        'place': None,
        'properties': obs.properties
      } for i, obs in enumerate(group)]
    } for groupName, group in data.survey.observationGroups.items()],
    'features': exportJSONfg(data)["features"] # I (Adrian) question of the value of this...
  }
  json.dump(ret, file, indent=4)
