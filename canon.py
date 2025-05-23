from data import *
import json
from jsonfg import * # I question this...
from math import sqrt

def importCSDM(file):
  from copy import deepcopy
  ids = set() # For validation
  def isUID(id):
    if id in ids: print("Warning: Duplicate ID!", id)
    ids.add(id)
  
  data = json.load(file)
  assert data["type"] == "FeatureCollection"
  projection = Projection(data.get("compoundCRS", data.get("horizontalCRS")), data.get("compoundCRS", data.get("verticalDatum")))
  annotations = [Annotation(ann.get("role"), ann.get("id"), ann.get("description")) for ann in data.get("annotations", [])]
  metadata = SurveyMetadata(data["name"], None, None, None, None, None, data["purpose"], None, None, None, annotations)
  
  instruments = []
  def instrument(*args):
    ret = InstrumentSetup(*args)
    instruments.append(ret)
    return ret

  referencedCSDs = [ReferencedCSD(ref.get("id"), ref["name"], ref.get("adminUnit", {}).get("href"), ref.get("adminUnit", {}).get("rel"), ref.get("adminUnit", {}).get("role"), ref["bearingRotation"], ref["time"]) for ref in data.get("referencedCSDs", [])]
  referencedDocs = [SupportingDocument(doc["title"], doc["href"], doc.get("role"), doc.get("rel")) for doc in data.get("supportingDocuments", [])]

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
      name = monument["properties"].get("name") or {}
      if isinstance(name, dict): name = name.get("label", monument["id"])
      monuments.append(Monument(name, point, monumentedBy.get("state", monumentedBy.get("monumentState")), monumentedBy.get("form", monumentedBy.get("monumentForm")), monumentedBy.get("condition", monumentedBy.get("monumentCondition")), properties = monument["properties"], klass = monument.get("featureType")))

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
      measure = measures.get(observation.get("id"))
      if measure is None:
        measure = Measure.fromProperties({}, {})
        measure.is_clockwise = None # Code could be cleaner, constructors set it to false.
        print("Expected to find measurements corresponding to observation", observation.get("id"), " but found none!")
      geoms = []
      if observation["topology"]["type"].lower() == "linestring":
        for i in range(1, len(observation["topology"]["references"])):
          start = pointsIndex[d(observation["topology"]["references"][i-1])]
          end = pointsIndex[d(observation["topology"]["references"][i])]
          # Distance comes from vector observations, we might want to split that class out!
          startname = start["properties"].get("name") or {}
          if isinstance(startname, dict): startname = startname.get("label", start["id"])
          endname = end["properties"].get("name") or {}
          if isinstance(endname, dict): endname = endname.get("label", end["id"])
          obs = ReducedObservation(None, instrument(start["id"], startname, None, start[""]), instrument(end["id"], endname, None, end[""]), measure.azimuth, measure.dist, measure.equipment, measure.distType, measure.azimuthType, observation["id"], start.get("time", data.get("time")), start["properties"], measure = measure)
          observations.append(obs)
          observationsIndex[observation["id"]] = obs

          geom = Line(observation["id"], start[""], end[""])
          geoms.append(geom)
        segments[observation["id"]] = geoms
      elif observation["topology"]["type"].lower() == "arc":
        start = pointsIndex[d(observation["topology"]["references"][0])]
        mid = pointsIndex[d(observation["topology"]["references"][1])]
        end = pointsIndex[d(observation["topology"]["references"][2])]
        # Distance, radius, etc comes from vector observations.
        startname = start["properties"].get("name") or {}
        if isinstance(startname, dict): startname = startname.get("label", start["id"])
        endname = end["properties"].get("name") or {}
        if isinstance(endname, dict): endname = endname.get("label", end["id"])
        isClockwise = measure.is_clockwise
        if isClockwise is None and "orientation" in observation["topology"]:
          isClockwise = observation["topology"]["orientation"] == "cw"
        obs = ReducedArcObservation(None, instrument(start["id"], startname, None, start[""]), instrument(end["id"], endname, None, end[""]), measure.azimuth, measure.radius, measure.dist, isClockwise, measure.equipment, measure.angleAccuracy, measure.arcType, observation["id"], start.get("time", data.get("time")), start["properties"], measure.distanceQuality, measure.distanceAccuracy, measure.angleQuality, measure = measure)
        observations.append(obs)
        observationsIndex[observation["id"]] = obs

        geom = Curve(observation["id"], isClockwise, measure.radius, start[""], mid[""], end[""])
        segments[observation["id"]] = [geom]
        obs.geom = geom
        obs.mid = mid
      elif observation["topology"]["type"].lower() == "arcwithcenter":
        # FIXME: center isn't mid?
        start = pointsIndex[d(observation["topology"]["references"][0])]
        end = pointsIndex[d(observation["topology"]["references"][1])]
        mid = pointsIndex[d(observation["topology"]["references"][2])]
        startname = start["properties"].get("name") or {}
        if isinstance(startname, dict): startname = startname.get("label", start["id"])
        endname = end["properties"].get("name") or {}
        isClockwise = measure.is_clockwise if measure.is_clockwise is not None else observation["topology"]["orientation"] == "cw"
        obs = ReducedArcObservation(None, instrument(start["id"], startname, None, start[""]), instrument(end["id"], endname, None, end[""]), measure.azimuth, measure.radius, measure.dist, isClockwise, measure.equipment, measure.angleAccuracy, measure.arcType, observation["id"], start.get("time", data.get("time")), start["properties"], measure.distanceQuality, measure.distanceAccuracy, measure.angleQuality, measure = measure)
        observations.append(obs)
        observationsIndex[observation["id"]] = obs

        geom = Curve.from_center(observation["id"], isClockwise, measure.radius, start[""], mid[""], end[""])
        segments[observation["id"]] = [geom]
        obs.geom = geom
        obs.center = mid
      elif observation["topology"]["type"].lower() == "subtendedangle":
        # Handle in a separate pass...
        pass
      elif observation["topology"]["type"].lower() == "circlebycenter":
        center = pointsIndex[d(observation["topology"]["references"][0])]
        centername = center["properties"].get("name") or {}
        if isinstance(centername, dict): centername = centername.get("label", center["id"])
        obs = CircleByCenter(observation["id"], center.get("time", data.get("time")), None, instrument(center["id"], centername, None, center[""]), observation["topology"]["radius"], observation["properties"])
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
        obs = CubicSplineObservation(obsId, refs(topology.get("startTangentVector")), refs(topology.get("endTangentVector")), refs(topology), observation["properties"])
        observations.append(obs)

        if obsId is not None:
          observationsIndex[obsId] = obs
          geom = Cubic(obsId, refs(topology.get("startTangentVector")), refs(topology.get("endTangentVector")), refs(topology))
          segments[obsId] = [geom]
      elif observation["topology"]["type"].lower() == "arcbychord":
        obsId = observation.get("id")
        pointA = pointsIndex[d(observation["topology"]["references"][0])][""]
        pointB = pointsIndex[d(observation["topology"]["references"][1])][""]
        radius = observation["topology"]["radius"]
        isClockwise = observation["topology"].get("orientation", "cw") == "cw"
        if measure.is_clockwise is not None: isClockwise = measure.is_clockwise

        obs = ArcByChord(obsId, pointA, pointB, radius, isClockwise, observation.get("properties"))
        observations.append(obs)

        # https://stackoverflow.com/questions/24928317/draw-a-curved-line-with-given-radius-and-two-locations
        dist = (pointA - pointB).dist2D()/2
        saggita = radius - sqrt(radius*radius - dist*dist)
        mid = (pointA + pointB).div()
        vec = pointB - pointA
        if not isClockwise: # FIXME: Do have these branches the right way around?
          vec.coord1, vec.coord2 = -vec.coord2, vec.coord1
        else:
          vec.coord1, vec.coord2 = vec.coord2, -vec.coord1
        pointC = vec.div(vec.dist2D()).mul(saggita) + mid
        geom = Curve(obsId, isClockwise, radius, pointA, pointC, pointB)
        segments[obsId] = [geom]
        obs.geom = geom
      else:
        print("Unexpected observedVector topology-type: ", observation["topology"]["type"])
    observationGroups[group["id"]] = observations

  for group in data["faces"]:
    for face in group["features"]:
      isUID(face["id"])
      assert face["topology"]["type"] == "Polygon"
      segments[face["id"]] = [Face(
          face["id"],
          face.get('type', "Feature"),
          [[segment for ref in refs for segment in segments[d(ref)]] for refs in face["topology"]["references"]]
        )]

  for group in data["observedVectors"]:
    for observation in group["features"]:
      if observation["topology"]["type"].lower() == "subtendedangle":
        refs = observation["topology"]["references"]
        if d(refs[0]) in observationsIndex and d(refs[0]) not in pointsIndex:
          print("WARNING! Subtended angle setup referencing an observation instead of point!")
        setup = pointsIndex[d(refs[0])][""]
        if d(refs[1]) in pointsIndex and d(refs[1]) not in observationsIndex:
          print("WARNING! Subtended angle backsight line references a point instead of an observation!")
        backsight = observationsIndex[d(refs[1])]
        if d(refs[2]) in pointsIndex and d(refs[2]) not in observationsIndex:
          print("WARNING! Subtended angle target line references a point instead of an observation!")
        target = observationsIndex[d(refs[2])]
        obs = SubtendedAngle(observation["id"], start.get("time", data.get("time")), None, setup, backsight, target, observation["properties"])
        observationGroups[group["id"]].append(obs)

  surfaceTINs = []
  for tins in data.get("surfaceTin", []):
    isUID(tins["id"])
    assert tins["type"].lower() == "featurecollection"
    features = []
    for tin in tins["features"]:
      isUID(tin["id"])
      assert tin["type"].lower() == "feature"
      assert tin["topology"]["type"].upper() == "TIN"
      faces = [[segments[x][0] for x in face] for face in tin["topology"]["references"]]
      assert all(isinstance(face, Face) for x in faces for face in x)
      features.append(SurfaceTIN(tin["id"], faces))
    surfaceTINs.append(SurfaceTINs(tins["id"], tins["featureType"], features, tins.get("properties", {})))

  terrainIntersectionCurves = []
  for curves in data.get("terrainIntersectionCurve", []):
    isUID(curves["id"])
    assert curves["type"].lower() == "featurecollection"
    features = []
    for curve in curves["features"]:
      isUID(curve["id"])
      assert curve["type"].lower() == "feature"
      assert curve["topology"]["type"].lower() == "multilinestring"
      topology = [y for x in curve["topology"]["references"] for y in segments[x]]
      features.append(TerrainIntersectionCurve(curve["id"], topology))
    terrainIntersectionCurves.append(TerrainIntersectionCurves(curves["id"], curves["featureType"], features, curves.get("properties", {})))

  parcels = []
  for parcel in data.get("parcels", []):
    for geom in parcel["features"]:
      isUID(geom["id"])
      geoms = []
      refs = geom["topology"]["references"]
      if len(refs) == 1 and isinstance(refs[0], list): refs = refs[0] # Handle double-nesting.
      geoms.append(Geom(geom["id"],
        [segment for ref in refs for segment in segments.get(ref if isinstance(ref, str) else ref.get('$ref'), [])],
        geom.get("type", "Feature")))
      parcels.append(Parcel.fromProperties(None, None, geoms, geom["properties"], geom["id"], geom.get("featureType"))) # TODO: Differentiate primary vs secondary
    
  return Cadastre(projection, {}, None, monuments, points, parcels, Survey(metadata, instruments, observationGroups), referencedCSDs = referencedCSDs, supportingDocs = referencedDocs, surfaceTINs = surfaceTINs, terrainIntersectionCurves = terrainIntersectionCurves)

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
        'time': str(ref.time)
      } for ref in data.referencedCSDs],
    'points': [{
      'id': 'surveymarks',
      'type': 'FeatureCollection',
      'featureType': 'csdm:SurveyMark',
      'features': [{
        'id': getattr(monument.point, 'objID', None) or i,
        'type': "Feature",
        'featureType': "SurveyPoint", # Could be "cadastralMark", etc?
        'time': str(monument.point.date), # Right fix?
        'place': {
          'type': "Point",
          'coordinates': [monument.point.coord1, monument.point.coord2] + ([monument.point.coord3] if monument.point.coord3 is not None else []) if monument.point is not None else []
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
        'type': geom.type,
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
    'surfaceTin': [{
      'id': tins.id or ("tin"+str(i)),
      'type': "FeatureCollection",
      'featureType': tins.type,
      'properties': tins.properties,
      'features': [{
        'id': tin.id or ("TIN"+str(j)),
        'type': "feature",
        'geometry': None,
        'topology': {
          'type': "TIN",
          'references': [[y.id for y in x] for x in tin.faces]
        }
      } for j, tin in enumerate(tins.features)]
    } for i, tins in enumerate(data.surfaceTINs)],
    'terrainIntersectionCurve': [{
      'id': curves.id or ("tic"+str(i)),
      'type': "FeatureCollection",
      'featureType': curves.type,
      'properties': curves.properties,
      'features': [{
        'id': feature.id or ("TIC"+str(j)),
        'type': "feature",
        'geometry': None,
        'topology': {
          'type': "MultiLineString",
          'references': [x.id for x in feature.topology]
        }
      } for j, feature in enumerate(curves.features)]
    } for i, curves in enumerate(data.terrainIntersectionCurves)],
    'features': exportJSONfg(data)["features"] # I (Adrian) question of the value of this...
  }
  json.dump(ret, file, indent=4, default=str)
