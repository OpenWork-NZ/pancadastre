from data import *
import json
from jsonfg import * # I question this...

def importCSDM(file):
  from copy import deepcopy
  ids = set() # For validation
  def isUID(id):
    if id in ids: print("Warning: Duplicate ID!", id)
    ids.add(id)
  
  data = json.load(file)
  assert data["type"] == "FeatureCollection"
  projection = Projection(data.get("compoundCRS", data.get("horizontalCRS")), data.get("compoundCRS", data.get("verticalDatum")))
  annotations = [Annotation(ann.get("role"), ann.get("id"), ann["description"]) for ann in data.get("annotations", [])]
  metadata = SurveyMetadata(data["name"], None, None, None, None, None, data["purpose"], None, None, None, annotations)
  
  instruments = []
  def instrument(*args):
    ret = InstrumentSetup(*args)
    instruments.append(ret)
    return ret

  referencedCSDs = [ReferencedCSD(ref["id"], ref["name"], ref.get("adminUnit", {}).get("href"), ref.get("adminUnit", {}).get("rel"), ref.get("adminUnit", {}).get("role"), ref["bearingRotation"], ref["time"]) for ref in data["referencedCSDs"]]

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
      isUID(monument["id"])
      monumentedBy = monument["properties"].get("monumentedBy")
      if monumentedBy is None:
        print("'monumentedBy' property is mandatory!", monument["id"], monument["properties"])
        monumentedBy = {}
      elif "state" not in monumentedBy or "form" not in monumentedBy or "condition" not in monumentedBy:
        if "monumentState" in monumentedBy:
          print("Outdated schema! Key 'monumentState' should be 'state'!")
        else:
          print("Invalid schema! Not attaching monument info!")
      place = monument.get("place", monument.get("topology", monument.get("geometry"))) # "geometry" doesn't seem right here...
      point = Point(None, None, monumentedBy.get("state", monumentedBy.get("monumentState")), place["coordinates"][0], place["coordinates"][1], place["coordinates"][2] if len(place["coordinates"]) >= 3 else None, monument["id"])
      points.append(point)
      pointsIndex[monument["id"]] = deepcopy(monument)
      pointsIndex[monument["id"]][""] = point
      # The observations will be initialized later
      monuments.append(Monument((monument["properties"].get("name") or {}).get("label", monument["id"]), point, monumentedBy.get("state", monumentedBy.get("monumentState")), monumentedBy.get("form", monumentedBy.get("monumentForm")), monumentedBy.get("condition", monumentedBy.get("monumentCondition")), properties = monument["properties"], klass = monument.get("featureType")))

  def n(comment): return None # Indicates a TODO...
  def d(x):
    if isinstance(x, dict): return x["$ref"]
    else: return x

  segments = {}
  observationGroups = {} # TODO: Build!
  observationsIndex = {}
  for group in data["observedVectors"]:
    observations = []
    for observation in group["features"]:
      if "id" in observation: isUID(observation["id"])
      # FIXME: Is it a bug if a point doesn't have a measure? Should this be reported?
      measure = measures.get(observation.get("id")) or Measure.fromProperties({}, {})
      if observation["topology"]["type"].lower() == "linestring":
        for i in range(1, len(observation["topology"]["references"])):
          start = pointsIndex[d(observation["topology"]["references"][i-1])]
          end = pointsIndex[d(observation["topology"]["references"][i])]
          # Distance comes from vector observations, we might want to split that class out!
          obs = ReducedObservation(None, instrument(start["id"], (start["properties"].get("name") or {}).get("label", start["id"]), None, start[""]), instrument(end["id"], (end["properties"].get("name") or {}).get("label", end["id"]), None, end[""]), measure.azimuth, measure.dist, measure.equipment, measure.distType, measure.azimuthType, observation["id"], start.get("time", data.get("time")), start["properties"], measure = measure)
          observations.append(obs)
          observationsIndex[observation["id"]] = obs

          geom = Line(observation["id"], start[""], end[""])
          segments[observation["id"]] = [geom]
      elif observation["topology"]["type"].lower() == "arc":
        start = pointsIndex[d(observation["topology"]["references"][0])]
        mid = pointsIndex[d(observation["topology"]["references"][1])]
        end = pointsIndex[d(observation["topology"]["references"][2])]
        # Distance, radius, etc comes from vector observations.
        obs = ReducedArcObservation(None, instrument(start["id"], (start["properties"].get("name") or {}).get("label", start["id"]), None, start[""]), instrument(end["id"], (end["properties"].get("name") or {}).get("label", end["id"]), None, end[""]), measure.azimuth, measure.radius, measure.dist, measure.is_clockwise, measure.equipment, measure.angleAccuracy, measure.arcType, observation["id"], start.get("time", data.get("time")), start["properties"], measure.distanceQuality, measure.distanceAccuracy, measure.angleQuality, measure = measure)
        observations.append(obs)
        observationsIndex[observation["id"]] = obs

        geom = Curve(observation["id"], measure.is_clockwise, measure.radius, start[""], mid[""], end[""])
        segments[observation["id"]] = [geom]
        obs.geom = geom
        obs.mid = mid
      elif observation["topology"]["type"].lower() == "arcwithcenter":
        # FIXME: center isn't mid?
        start = pointsIndex[d(observation["topology"]["references"][0])]
        end = pointsIndex[d(observation["topology"]["references"][1])]
        mid = pointsIndex[d(observation["topology"]["references"][2])]
        obs = ReducedArcObservation(None, instrument(start["id"], (start["properties"].get("name") or {}).get("label", start["id"]), None, start[""]), instrument(end["id"], (end["properties"].get("name") or {}).get("label", end["id"]), None, end[""]), measure.azimuth, measure.radius, measure.dist, measure.is_clockwise, measure.equipment, measure.angleAccuracy, measure.arcType, observation["id"], start.get("time", data.get("time")), start["properties"], measure.distanceQuality, measure.distanceAccuracy, measure.angleQuality, measure = measure)
        observations.append(obs)
        observationsIndex[observation["id"]] = obs

        geom = Curve.from_center(observation["id"], measure.is_clockwise, measure.radius, start[""], mid[""], end[""])
        segments[observation["id"]] = [geom]
        obs.geom = geom
        obs.center = mid
      elif observation["topology"]["type"].lower() == "subtendedangle":
        # Handle in a separate pass...
        pass
      elif observation["topology"]["type"].lower() == "circlebycenter":
        center = pointsIndex[d(observation["topology"]["references"][0])]
        obs = CircleByCenter(observation["id"], center.get("time", data.get("time")), None, instrument(center["id"], (center["properties"].get("name") or {}).get("label", center["id"]), None, center[""]), observation["topology"]["radius"], observation["properties"])
        observations.append(obs)
        observationsIndex[observation["id"]] = obs

        radius = observation["topology"]["radius"]
        geoms = []
        geoms.append(Curve.from_center(observation["id"] + "-cw", True, radius, center[""].offset1(-radius), center[""], center[""].offset1(radius)))
        geoms.append(Curve.from_center(observation["id"] + "-ccw", False, radius, center[""].offset1(-radius), center[""], center[""].offset1(radius)))
        segments[observation["id"]] = geoms
      elif observation["topology"]["type"].lower() == "cubicspline":
        def refs(feat):
          return [pointsIndex[d(pt)][""] for pt in feat["references"]]
        obsId = observation.get("id")
        topology = observation["topology"]
        obs = CubicSplineObservation(obsId, refs(topology["startTangentVector"]), refs(topology["endTangentVector"]), refs(topology), observation["properties"])
        observations.append(obs)

        if obsId is not None:
          observationsIndex[obsId] = obs
          geom = Cubic(obsId, refs(topology["startTangentVector"]), refs(topology["endTangentVector"]), refs(topology))
          segments[obsId] = [geoms]
      else:
        print("Unexpected observedVector topology-type: ", observation["topology"]["type"])
    observationGroups[group["id"]] = observations
  for group in data["observedVectors"]:
    for observation in group["features"]:
      if observation["topology"]["type"].lower() == "subtendedangle":
        refs = observation["topology"]["references"]
        setup = observationsIndex.get(d(refs[0]), pointsIndex.get(d(refs[0])))
        if setup is None: print("Warning: Failed to dereference subtended angle setup as either point or observation!")
        elif isinstance(setup, dict): setup = setup[""]
        backsight = observationsIndex[d(refs[1])]
        target = observationsIndex.get(d(refs[2]), pointsIndex.get(d(refs[2])))
        if target is None: print("Warning: Failed to dereference subtended angle target as either point or observation!")
        elif isinstance(target, dict): target = target[""]
        obs = SubtendedAngle(observation["id"], start.get("time", data.get("time")), None, setup, backsight, target, observation["properties"])
        observationGroups[group["id"]].append(obs)

  parcels = []
  for parcel in data.get("parcels", []):
    for geom in parcel["features"]:
      isUID(geom["id"])
      geoms = []
      refs = geom["topology"]["references"]
      if len(refs) == 1 and isinstance(refs[0], list): refs = refs[0] # Handle double-nesting.
      geoms.append(Geom(geom["id"],
        [segment for ref in refs for segment in segments.get(ref if isinstance(ref, str) else ref.get('$ref'), [])]))
      parcels.append(Parcel.fromProperties(None, None, geoms, geom["properties"], geom["id"], geom.get("featureType"))) # TODO: Differentiate primary vs secondary
    
  return Cadastre(projection, {}, None, monuments, points, parcels, Survey(metadata, instruments, observationGroups), referencedCSDs = referencedCSDs)

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
    'annotations': [{
      'role': ann.type,
      'description': ann.desc
    } for ann in data.survey.metadata.annotations],
    'referencedCSDs': [{
        'id': ref.id,
        'name': ref.name,
        'adminUnit': {
          'href': ref.href,
          'rel': ref.rel,
          'role': ref.role
        },
        'bearingRotation': ref.bearing,
        'time': ref.time
      } for ref in data.referencedCSDs],
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
          'coordinates': [monument.point.coord1, monument.point.coord2] if monument.point is not None else []
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
