#! /usr/local/bin/python3

from xml.etree import ElementTree as ET
import rdflib
from rdflib.namespace import RDFS, RDF, GEO
import json

from datetime import datetime
from pyproj import Transformer, Geod

from math import pi, isclose, ceil, sqrt
tau = 2*pi

## Intermediate Data Model

class Cadastre:
  def __init__(self, projection, features, units, monuments, points, parcels, survey, boundaries = []):
    self.projection = projection or Projection("epsg:4326")
    self.features = features
    self.units = units # Refactor to assign units directly against points, observations, etc?
    self.monuments = monuments
    self.points = points
    self.parcels = parcels
    self.survey = survey # Haven't actually used the multiplicity here... Merge it's fields in?
    self.boundaries = boundaries # Any survey boundary vectors

class App:
  def __init__(self, app, appversion, appAuthors):
    self.app = app
    self.appversion = appversion
    self.appAuthors = appAuthors

class Projection:
  def __init__(self, horizontal, vertical = None):
    self.horizontal = horizontal
    self.vertical = vertical if vertical is not None else horizontal
    # Support PyProj input format for custom projections -Proj4?

class Unit:
  def __init__(self, areaUnit, linearUnit, volumeUnit, temperatureUnit, pressureUnit, angularUnit, directionUnit):
    # Ask Andrew
    self.areaUnit = areaUnit
    self.linearUnit = linearUnit
    self.volumeUnit = volumeUnit
    self.temperatureUnit = temperatureUnit
    self.pressureUnit = pressureUnit
    self.angularUnit = angularUnit
    self.directionUnit = directionUnit

class Monument:
  def __init__(self, name, point, state, type, condition):
    self.name = name
    self.point = point
    self.state = state
    self.type = type
    self.condition = condition
    self.point.monuments.append(self)
    # Modelling remonumenting?

class Point:
  def __init__(self, survey, name, state, northing, easting, objID = None):
    self.survey = survey
    self.name = name
    self.state = state # State belongs on monument, not point?
    self.northing = northing
    self.easting = easting
    self.objID = objID
    self.instruments = []
    self.monuments = []
    self.segments = []

  @property
  def observations(self): return {x for setup in self.instruments for x in setup.observations}

  @property
  def observation(self):
    for instrument in self.instruments:
      for observation in instrument.observations:
        return observation
    return None

  @property
  def date(self): return self.observation.date if self.observation is not None else None
  @property
  def purpose(self): return self.observation.purpose if self.observation is not None else None

  def __iter__(self):
    yield self.easting
    yield self.northing

  def __eq__(self, other):
    return isclose(self.northing, other.northing) and isclose(self.easting, other.easting)

  def __repr__(self):
    return repr(self.easting) + "," + repr(self.northing)

  def __add__(self, other):
    return Point(self.survey, self.name, self.state, self.northing + other.northing, self.easting + other.easting, str(self.objID) + "+" + str(other.objID))
  def __sub__(self, other):
    return Point(self.survey, self.name, self.state, self.northing - other.northing, self.easting - other.easting, str(self.objID) + "-" + str(other.objID))
  def div(self, other = 2):
    return Point(self.survey, self.name, self.state, self.northing/other, self.easting/other, str(self.objID) + "/" + str(other))

class Parcel:
  def __init__(self, name, area, type, state, klass, format, center, geom, titleDoc, address=None, desc = None):
    self.name = name
    self.area = area
    self.type = type
    self.state = state
    self.klass = klass
    self.format = format
    self.center = center # Do we want?
    self.geom = geom
    self.titleDoc = titleDoc # Contains LandXML Title element or CSDM interestRef property, currently unused.
    self.address = address # Do we want?

class Address: # Ask Nick if there's a better address model
  def __init__(self, type, num, roadNameType, roadName, roadType, adminArea):
    self.type = type
    self.num = num
    self.roadNameType = roadNameType
    self.roadName = roadName
    self.roadType = roadType
    self.adminArea = adminArea

class Geom:
  def __init__(self, name, segments):
    self.name = name # Identifier
    self.segments = segments

  @staticmethod
  def flatten(segments):
    ret = []
    for segment in segments:
      if isinstance(segment, list): ret.extend(segment)
      else: ret.append(segment)
    return ret

class Segment:
  pass

class Line(Segment):
  def __init__(self, id, start, end, state = None):
    self.id = id
    self.start = start
    self.end = end
    self.state = state # NOTE: Are states on lines needed?
    self.start.segments.append(self)
    self.end.segments.append(self)

  @property
  def observation(self):
    for observation in self.start.observations:
      if self.start == observation.setupPoint and self.end == observation.targetPoint: return observation
      if self.end == observation.setupPoint and self.start == observation.targetPoint: return observation

  @property
  def observationProperties(self): return self.observation.properties
    
  def __repr__(self):
    return "Line " + str(self.id) + "(" + repr(self.start) + " -> " + repr(self.end) + ")"

  def __iter__(self):
    yield self.start
    yield self.end

  def __hash__(self):
    return hash(self.id)
  def __eq__(self, other):
    return self.id == other.id

class Curve(Segment):
  def __init__(self, id, is_clockwise, radius, start, mid, end, state = None):
    self.id = id
    self.is_clockwise = is_clockwise
    self.radius = radius
    self.start = start
    self.end = end
    self.mid = mid
    self.state = state # NOTE: Are states on curves needed? As per lines?
    self.start.segments.append(self)
    self.end.segments.append(self)
    self.mid.segments.append(self)

  @staticmethod
  def from_center(id, is_clockwise, r, start, center, end, state = None):
    # Do we need to consider is_clockwise in this formula?
    if isclose(r, 0):
      dist = center - start
      r = sqrt(dist.northing*dist.northing + dist.easting*dist.easting)
    a, b = (start, end) if is_clockwise else (end, start)
    ab = b - a
    lab = sqrt(ab.easting*ab.easting + ab.northing*ab.northing)
    uab = ab.div(lab)
    mab = (a + b).div()
    f = r - sqrt(r*r - lab*lab/4)
    mid = Point(center.survey, center.name, center.state, mab.northing + uab.easting*f, mab.easting - uab.northing*f, center.objID)
    return Curve(id, is_clockwise, r, start, mid, end, state)

  @property
  def observation(self):
    for observation in self.start.observations:
      if self.start == observation.setupPoint and self.end == observation.targetPoint: return observation
      if self.end == observation.setupPoint and self.start == observation.targetPoint: return observation

  @property
  def observationProperties(self): return self.observation.properties

  def __iter__(self):
    yield self.start
    yield self.end

  def __hash__(self): return hash(self.id)
  def __eq__(self, other): return self.id == other.id

class Feature:
  def __init__(self, desc, name, geom, properties = {}):
    self.desc = desc
    self.name = name
    self.geom = geom
    self.properties = properties

class Survey:
  def __init__(self, metadata, instruments, observationGroups):
    self.metadata = metadata
    self.instruments = instruments
    self.observationGroups = observationGroups

class SurveyMetadata:
  def __init__(self, name, firm, surveyor, type, jurisdiction, format, purpose, headOfPower, adminAreas, personnel, annotations):
    # Ask Andrew what's needed here
    self.name = name
    self.firm = firm
    self.surveyor = surveyor
    self.type = type
    self.jurisdiction = jurisdiction
    self.format = format
    self.purpose = purpose
    self.headOfPower = headOfPower
    self.adminAreas = adminAreas
    self.personnel = personnel
    self.annotations = annotations # Do these belong elsewhere?

class AdminArea:
  def __init__(self, type, name, code):
    self.type = type
    self.name = name
    self.code = code

class Personnel:
  def __init__(self, name, role, registerType, registerID):
    self.name = name
    self.role = role
    self.registerType = registerType
    self.registerID = registerID

class Annotation:
  def __init__(self, type, name, desc, parcel = None):
    self.type = type
    self.name = name
    self.desc = desc
    self.parcel = parcel

class InstrumentSetup:
  def __init__(self, id, stationName, height, point):
    self.id = id
    self.stationName = stationName
    self.height = height
    self.point = point
    self.point.instruments.append(self)
    self.observations = []

class Observation:
  pass

class ReducedObservation(Observation):
  def __init__(self, purpose, setup, targetSetup, azimuth, horizDist, equipment, distType, azimuthType, name, date, properties = {}):
    self.purpose = purpose
    self.setup = setup
    self.targetSetup = targetSetup
    self.azimuth = azimuth
    self.horizDist = horizDist
    self.equipment = equipment
    self.distType = distType
    self.azimuthType = azimuthType
    self.name = name
    self.date = date
    self.properties = properties
    self.setup.observations.append(self)
    if self.targetSetup is not None:
      self.targetSetup.observations.append(self)

  @property
  def setupPoint(self): return self.setup.point

  @property
  def targetPoint(self): return self.targetSetup.point

  @property
  def line(self):
    for segment in self.setupPoint.segments:
      if self.setupPoint == segment.start and self.targetPoint == segment.end: return segment
      if self.targetPoint == segment.start and self.setupPoint == segment.end: return segment

class ReducedArcObservation(Observation):
  def __init__(self, purpose, setup, targetSetup, chordAzimuth, radius, length, is_clockwise, equipmentUsed, arcLengthAccuracy, arcType, name, date, properties = {}):
    self.purpose = purpose
    self.setup = setup
    self.targetSetup = targetSetup
    self.chordAzimuth = chordAzimuth
    self.radius = radius
    self.length = length
    self.is_clockwise = is_clockwise
    self.equipmentUsed = equipmentUsed
    self.arcLengthAccuracy = arcLengthAccuracy
    self.arcType = arcType
    self.name = name
    self.date = date
    self.properties = properties
    self.setup.observations.append(self)
    if self.targetSetup is not None:
      self.targetSetup.observations.append(self)

  @property
  def line(self):
    for segment in self.setupPoint.segments:
      if self.setupPoint == segment.start and self.targetPoint == segment.end: return segment
      if self.targetPoint == segment.start and self.setupPoint == segment.end: return segment


class RedHorizPos(Observation): # Spell out "Reduced"
  def __init__(self, desc, name, objID, setup, date, horizDatum, northing, easting, horizFix, order, properties = {}):
    self.desc = desc
    self.name = name
    self.objID = objID
    self.setup = setup
    self.date = date
    self.horizDatum = horizDatum
    self.northing = northing
    self.easting = easting
    self.horizFix = horizFix
    self.order = order
    self.properties = properties
    self.setup.observations.append(self)
    self.line = None

#class ReducedVertPos(Observation) for 3D?

## Importers

def importLandXML(file):
  # Look at both CRSs, if survey header has different projection apply to places.
  file = ET.parse(file)

  app = None
  appEl = file.find("{http://www.landxml.org/schema/LandXML-1.2}Application")
  if appEl is not None:
    appAuthors = [author.get('createdBy') for author in appEl if author.get('createdBy')]
    app = App(appEl.get('name'), appEl.get('version'), appAuthors)

  # Capture both projections, does it need a datamodel change?
  # Or do we change projection to what the surveyer worked in (inner proj) as we parse this?
  # Discuss with Andrew, some jurisdictions may want it in a particular CRS.
  projEl = file.find('{http://www.landxml.org/schema/LandXML-1.2}CoordinateSystem')
  proj = None
  if projEl is not None and projEl.get('datum'):
    proj = Projection(projEl.get('horizontalDatum', projEl.get('datum')), projEl.get('datum'))
  elif projEl is not None and projEl.get('name'):
    proj = Projection(projEl.get('name'))

  #features = {feature.get('name'): feature.get('version')
  #            for feature in file.findall('FeatureDictionary')} # Do we actually need this? (Vic/ePlan profile). Rename!
  units = [Unit(l.get('areaUnit'), l.get('linearUnit'), l.get('volumeUnit'),
                l.get('temperatureUnit'), l.get('pressureUnit'), l.get('angularUnit'),
                l.get('directionUnit')) for l in file.find('{http://www.landxml.org/schema/LandXML-1.2}Units')]

  pointsIndex = dict()
  points = []
  for pointEl in file.find('{http://www.landxml.org/schema/LandXML-1.2}CgPoints'):
    coords = list(map(float, pointEl.text.split()))
    if len(coords) != 2: continue

    point = Point(pointEl.get('pntSurv'), pointEl.get('name'), pointEl.get('state'), # difference in ePlan.
                  coords[0], coords[1], pointEl.get('oID'))
    points.append(point)
    if pointEl.get('name'): pointsIndex[pointEl.get('name')] = point

  monuments = [Monument(l.get('name'), pointsIndex.get(l.get('pntRef')),
                        l.get('state'), l.get('type'), l.get('condition'))
              for l in file.find('{http://www.landxml.org/schema/LandXML-1.2}Monuments')] # Capture desc!

  def parseGeoms(el):
    geoms = [] # Parcel boundary lines, edges of a polygon.
    for geom in el.findall('{http://www.landxml.org/schema/LandXML-1.2}CoordGeom'):
      path = []
      name = geom.get('name', 'g')
      for i, segment in enumerate(list(geom)):
        if segment.tag == "{http://www.landxml.org/schema/LandXML-1.2}Line":
          start = segment.find('{http://www.landxml.org/schema/LandXML-1.2}Start')
          if start is not None and start.get('pntRef') and pointsIndex.get(start.get('pntRef')):
            start = pointsIndex.get(start.get('pntRef'))
          elif start.text.strip() != "":
            start = Point(None, None, None, *[float(x) for x in start.text.strip().split()])
          else:
            print("Invalid start reference")
            continue
          end = segment.find('{http://www.landxml.org/schema/LandXML-1.2}End')
          if end is not None and end.get('pntRef') and pointsIndex.get(end.get('pntRef')):
            end = pointsIndex.get(end.get('pntRef'))
          elif end.text.strip() != "":
            start = Point(None, None, None, *[float(x) for x in end.text.strip().split()])
          else:
            print("Invalid end reference")
            continue
          path.append(Line(segment.get('oID', name + str(i)), start, end, segment.get('state')))
        elif segment.tag == "{http://www.landxml.org/schema/LandXML-1.2}Curve":
          start = segment.find('{http://www.landxml.org/schema/LandXML-1.2}Start')
          if start is not None and start.get('pntRef') and pointsIndex.get(start.get('pntRef')):
            start = pointsIndex.get(start.get('pntRef'))
          elif start.text.strip() != "":
            start = Point(None, None, None, *[float(x) for x in start.text.strip().split()])
          else:
            print("invalid start reference")
            continue
          mid = segment.find('{http://www.landxml.org/schema/LandXML-1.2}Center')
          if mid is not None and mid.get('pntRef') and pointsIndex.get(mid.get('pntRef')):
            mid = pointsIndex.get(mid.get('pntRef'))
          elif mid.text.strip() != "":
            mid = Point(None, None, None, *[float(x) for x in mid.text.strip().split()])
          else:
            print("invalid mid reference")
            continue
          end = segment.find('{http://www.landxml.org/schema/LandXML-1.2}End')
          if end is not None and end.get('pntRef') and pointsIndex.get(end.get('pntRef')):
            end = pointsIndex.get(end.get('pntRef'))
          elif end.text.strip() != "":
            end = Point(None, None, None, *[float(x) for x in end.text.strip().split()])
          else:
            print("invalid end reference")
            continue
          path.append(Curve.from_center(segment.get('oID', name + str(i)), segment.get('rot') == 'cw',
                                      float(segment.get('radius', '0')), start, mid, end, segment.get('state')))
        else:
          print("Unsupported geometry type!", segment.tag)
      geoms.append(Geom(name, path))
    return geoms

  def i(x):
    if x is None: return {}
    else: return x

  parcels = []
  parcelsIndex = {}
  for parcel in i(file.find('{http://www.landxml.org/schema/LandXML-1.2}Parcels')):
    center = parcel.find('{http://www.landxml.org/schema/LandXML-1.2}Center') # CSDM model does not capture this, should we?
    if not center.get('pntRef'): center = {'pntRef': None}

    geoms = parseGeoms(parcel)
    # titleType is Vic only, want to preserve name.
    titles = {l.get('name'): l.get('titleType') for l in parcel.findall('{http://www.landxml.org/schema/LandXML-1.2}Title')}

    addressEl = parcel.find('{http://www.landxml.org/schema/LandXML-1.2}LocationAddress') # Vic ePlan only.
    address = None
    if addressEl:
      road = addressEl.find('{http://www.landxml.org/schema/LandXML-1.2}RoadName')
      adminAreaEl = addressEl.find('{http://www.landxml.org/schema/LandXML-1.2}AdministrativeArea')

      adminArea = None
      if adminAreaEl:
        adminArea = AdminArea(adminAreaEl.get('adminAreaType'), adminAreaEl.get('adminAreaName'),
                              adminAreaEl.get('adminAreaCode'))
      address = Address(addressEl.get('addressType'), addressEl.get('numberFirst'),
                        road.get('roadNameType'), road.get('roadName'), road.get('roadType'),
                        adminArea)

    parcel = Parcel(parcel.get('name'), parcel.get('area'), parcel.get('parcelType'),
                    parcel.get('state'), parcel.get('class'), parcel.get('parcelFormat'),
                    pointsIndex.get(center.get('pntRef')), geoms, titles, address,
                    parcel.get('desc'))
    parcels.append(parcel)
    parcelsIndex[parcel.name] = parcel

  features = dict() # Where did we see these? Failing to find the element again in sample files.
  for features in file.find('{http://www.landxml.org/schema/LandXML-1.2}PlanFeatures') or []:
    features[features.get('name')] = [Feature(l.get('desc'), l.get('name'), parseGeoms(l))
                                      for l in features.iter()]

  surveyEl = file.find('{http://www.landxml.org/schema/LandXML-1.2}Survey')
  survey = None
  if surveyEl:
    surveyHeaderEl = surveyEl.find('{http://www.landxml.org/schema/LandXML-1.2}SurveyHeader')
    # This is based on Vic ePlan, need to modify for LandOnline.
    surveyMeta = SurveyMetadata(surveyHeaderEl.get('name'), surveyHeaderEl.get('surveyFirm'),
                                surveyHeaderEl.get('surveyorReference'),
                                surveyHeaderEl.get('type'), surveyHeaderEl.get('jurisdiction'),
                                surveyHeaderEl.get('surveyFormat'),
                                # TODO: Capture county, desc, endTime, surveyor, surveyorFirm, surveyorReference, class
                                # LINZ has it as SurveyPurpose
                                i(surveyHeaderEl.find('{http://www.landxml.org/schema/LandXML-1.2}PurposeOfSurvey')).get('name'),
                                i(surveyHeaderEl.find('{http://www.landxml.org/schema/LandXML-1.2}HeadOfPower')).get('name'),
                                [AdminArea(l.get('adminAreaType'), l.get('adminAreaName'),
                                          l.get('adminAreaCode'))
                                  for l in surveyHeaderEl.findall('{http://www.landxml.org/schema/LandXML-1.2}AdministrativeArea')],
                                [Personnel(l.get('name'), l.get('role'), l.get('regType'),
                                          l.get('regNumber'))
                                  for l in surveyHeaderEl.findall('{http://www.landxml.org/schema/LandXML-1.2}Personnel')],
                                [Annotation(l.get('type'), l.get('name'), l.get('desc'),
                                            parcelsIndex.get(l.get('pclRef')))
                                  for l in surveyHeaderEl.findall('{http://www.landxml.org/schema/LandXML-1.2}Annotation')])

    instruments = []
    instrumentsIndex = dict()
    for el in surveyEl.findall('{http://www.landxml.org/schema/LandXML-1.2}InstrumentSetup'):
      instrument = InstrumentSetup(el.get('id'), el.get('stationName'), el.get('instrumentHeight'),
                        pointsIndex.get(el.find('{http://www.landxml.org/schema/LandXML-1.2}InstrumentPoint').get('pntRef')))
      instruments.append(instrument)
      instrumentsIndex[instrument.id] = instrument

    observationGroups = dict()
    for el in surveyEl.findall('{http://www.landxml.org/schema/LandXML-1.2}ObservationGroup'):
      observations = []
      for observation in list(el):
        # FIXME: skip duplicate IDs.
        # FIXME: target point may be a sub-element for LandOnline
        # Be resilient to missing purpose.
        if observation.tag == "{http://www.landxml.org/schema/LandXML-1.2}ReducedObservation": # All in LandOnline are this type.
          observations.append(ReducedObservation(observation.get('purpose'),
              instrumentsIndex.get(observation.get('setupID')),
              instrumentsIndex.get(observation.get('targetSetupID')),
              float(observation.get('azimuth')), float(observation.get('horizDistance', '0')),
              observation.get('equipmentUsed'), observation.get('distanceType'),
              observation.get('azimuthType'), observation.get('name'), observation.get('date'))) # Capture azimuthAccuracy, distanceAccuracy
              # If azimuthType or distanceType is "adopted" need to capture adoptedAzimuthSurvey, azimuthDistanceSurvey, & azimuthAdoptionFactor.
        elif observation.tag == "{http://www.landxml.org/schema/LandXML-1.2}ReducedArcObservation":
          observations.append(ReducedArcObservation(observation.get('purpose'),
              instrumentsIndex.get(observation.get('setupID')),
              instrumentsIndex.get(observation.get('targetSetupID')),
              float(observation.get('chordAzimuth')), float(observation.get('radius')),
              float(observation.get('length')), observation.get('rot') == 'cw',
              observation.get('equipmentUsed'), float(observation.get('arcLengthAccuracy')),
              observation.get('arcType'), observation.get('name'), observation.get('date')))
        elif observation.tag == "{http://www.landxml.org/schema/LandXML-1.2}RedHorizontalPosition":
          observation.append(RedHorizPos(observation.get('desc'),
              observation.get('name'), observation.get('oID'),
              instrumentsIndex[observation.get('setupID')],
              datetime.parse(observation.get('date'), '%Y-%m-%d'),
              observation.get('horizontalDatum'),
              float(observation.get('latitude')), float(observation.get('longitude')),
              observation.get('horizontalFix'), int(observation.get('order'))))
        else:
          print("Unsupported tag!", observation.tag)
      observationGroups[el.get('id')] = observations
    survey = Survey(surveyMeta, instruments, observationGroups)

  return Cadastre(proj, features, units, monuments, points, parcels, survey)

def importTopology(file):
  data = derefJSONLinks(json.load(file))
  assert data["type"] == "FeatureCollection"
  # TODO: provenance schema?
  surveyMeta = SurveyMetadata(data["name"], "TODO: firm", "TODO: surveyer", "TODO: type", None, "TODO: format", data["purpose"], None, None, None, None)
  proj = Projection(str(data["crs"])) if "crs" in data else None
  points = []

  monuments = []
  for feature in data["features"]:
    assert feature is not None
    assert feature["type"] == "Feature"
    assert feature["geometry"] is None or feature["geometry"]["type"] == "Point"

    pt = None
    if feature["geometry"] is not None:
      coords = feature["geometry"]["coordinates"]
      pt = Point(feature.get("fromSurvey"), feature["id"], None, coords[1], coords[0]) # state?
      points.append(pt)
    props = feature.get("properties", {"name": {}})
    monuments.append(Monument(props["name"].get("label"), pt,
        feature.get("ptQuality"), feature.get("featureType"), None)) # Condition? State?

  observations = {}
  for i, group in enumerate(data["vectorObservations"]):
    assert group["type"] == "FeatureCollection"
    assert group["featureType"] == "sosa:ObservationCollection"
    subObservations = []

    for observation in group["features"]:
      assert observation["type"] == "Feature"
      assert observation["geometry"] is None or observation["geometry"]["type"] == "LineString"

      props = observation.get("properties", {"hasResult": {}})
      kind = props.get("distanceType")
      item = None
      if kind == "icsmdistancetype:ellipsoidal":
        item = ReducedArcObservation("TODO: purpose", "TODO: setup", "TODO: targetSetup", props["hasResult"].get("angle"), "TODO: radius", props["hasResult"].get("distance"), "TODO: isClockwise", "TODO: equipmentUsed", props.get("angleAccuracy"), props.get("angleType"), props.get("hasFeatureOfInterest"))
      elif "pose" in props["hasResult"]:
        item = ReducedArcObservation("TODO: purpose", "TODO: setup", "TODO: targetSetup", props["hasResult"]["pose"].get("angles", {}).get("yaw"), "TODO: radius", props["hasResult"].get("distance"), "TODO: isClockwise", "TODO: equipmentUsed", "TODO angleArrucay", "TODO: angleType", props.get("hasFeatureOfInterest"))
      else:
        raise Exception("Unsupported observation type: " + str(kind))
      subObservations.append(observation)
    observations[group.get("properties", {}).get("sensor", i)] = subObservations

  def parseGeom(name, el):
    if el["type"].lower() == "linestring":
      geom = []
      lastPoint = None
      for i, edge in enumerate(el["references"]):
        if edge is None: continue
        assert edge["type"] == "Feature"
        geometry = edge.get("topology", edge.get("geometry", {"type": ""}))
        assert geometry["type"].lower() == "point"
        coord = geometry["coordinates"]
        point = Point(edge.get("fromSurvey"), edge["id"], None, coord[1], coord[0]) # state?
        points.append(point)
        if lastPoint: geom.append(Line(name + str(i), lastPoint, point))
        lastPoint = point
      return Geom(name, geom)
    elif el["type"].lower() == "polygon":
      geom = []
      for i, edge in enumerate(el["references"]):
        if edge is None: continue
        assert edge["type"] == "Feature"
        geom.extend(parseGeom(edge.get("id", i), edge).segments)
      return Geom(el.get("id", name), geom)
    elif el["type"].lower() == "feature":
      return parseGeom(el.get("id", name), el.get("topology", el.get("geometry", {})))
    else:
      raise Exception("Unsupported geometry type: " + el["type"])

  parcels = []
  for parcel in data["parcels"]:
    assert group["type"] == "FeatureCollection"
    for feature in parcel["features"]:
      assert feature["type"] == "Feature"
      geometry = feature.get("topology", feature.get("geometry", {"type": ""}))

      # NOTE: Many of these fields have been moved into the contained features.
      parcels.append(Parcel(feature["properties"]["name"]["label"], feature["properties"]["area"], parcel["featureType"], feature["properties"]["state"], feature["properties"]["class"], "TODO: format", "TODO: center", [parseGeom(feature["properties"]["name"]["label"], geometry)], feature["properties"]["interestRef"]))

  return Cadastre(proj, {}, None, monuments, points, parcels, Survey(None, [], observations))

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

## Exporters

# TODO: Tweak according to surveyfeatures.json, circular_arc.json
def exportJSONfg(data, file = None):
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
    trans = Transformer.from_crs(fileproj, proj)
    return {
      'type': 'Polygon',
      'coordinates': [simplifyPoly(trans, Geom.flatten(geom.segments)) for geom in parcel.geom]
    }
  fileproj = data.projection.horizontal
  if fileproj[:5] == "epsg:": fileproj = fileproj[5:]
  trans = Transformer.from_crs(fileproj, 'wgs84')
  
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
  # Alternate algorithm:
  # For reduced observation, select all lines in all parcels where the targetpoint == endpoint.
  # If not more than one select the objID for that line & that becomes hasFeatureOfInterest.
  # If more than one find the instrumentID & select the InstrumentSetup, get the pntRef from the instrumentPoint
  # Select from set of lines the line with the start == that.
  
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
        'properties': {
          'name': {
            'label': monument.name,
            'hasPart': [
              {'type': "Provenance", 'label': getattr(monument.point, 'state', None)},
              {'type': "Stamp", 'label': monument.type} # Is this right?
            ],
            # pointQuality?
            # purpose?
            'monumentForm': monument.type,
            'monumentCondition': monument.condition,
            'markState': monument.state
            # comment?
          }
        }
      } for i, monument in enumerate(data.monuments)] + observedVecs + [{
      'type': 'Feature',
      'id': parcel.name or i,
      'featureType': 'boundary',
      'links': [],
      'time': '', # Appears to be missing from LandXML
      'place': exportGeom(parcel, 'wgs84'),
      'geometry': exportGeom(parcel, fileproj),
      'properties': {
        'comment': ""
      } # Place miscallanea here!
    } for i, parcel in enumerate(data.parcels)]# + [{
#      'type': "Feature",
#      'featureType': "sosa:ObservationCollection",
#      'properties': {
#        'usedProcedure': "icsm:XXX"
#      } # FIXME: This seems incomplete...
#    } for groupID, group in data.survey.observationGroups.items() for i, observation in enumerate(group)]
  }
  if file is not None: json.dump(ret, file, indent=4)
  else: return ret

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
  
  def observationProperties(obs):
    if isinstance(obs, ReducedObservation):
      # Fields not echoed: setup & targetSetup, purpose, & azimuthType
      return {
          'hasFeatureOfInterest': obs.setupPoint.objID,
          'resultTime': obs.date,
          'hasResult': {
            'distance': obs.horizDist, # or obs.distance depending on type.
            'azimuth': obs.azimuth
          },
          'distanceType': obs.distType, # Reformat?
          'distanceProvenance': obs.setup.point.state, # Reformat?
          'equipmentUsed': obs.equipment, # Reformat?
          #'distanceQuality': 'A', # need to model quality!
          #'angleQuality': 'A', # need to model quality!
          #'distanceAccuracy': '???', # need to model accuracy!
          #'angleAccuracy': '???' # need to model accuracy!
         }
    elif isinstance(obs, ReducedArcObservation):
      # Fields not echoed: setup & targetSetup, & purpose
      return {
          'hasFeatureOfInterest': obs.setup.point.objID,
          'resultTime': obs.date,
          'hasResult': {
            'distance': obs.length,
            'angle': obs.chordAzimuth,
            'radius': obs.radius
          },
          'rot': "cw" if obs.is_clockwise else "ccw",
          'distanceType': 'icsmdistancetype:ellipsoidal', # What should this be?
          'distanceProvenance': obs.setup.point.state, # Reformat?
          'angleType': obs.arcType, # Reformat?
          'equipmentUsed': obs.equipmentUsed,
          #'distanceQuality': 'A', # need to model quality!
          #'angleQuality': 'A', # need to model quality!
          #'distanceAccuracy': '???', # need to model accuracy!
          'angleAccuracy': obs.arcLengthAccuracy
         }
    elif isinstance(obs, RedHorizPos):
      # Fields not echoed: objID, setup
      return {
          'hasFeatureOfInterest': obs.setup.point.objID, # FIXME: Parse the association!
          'resultTime': obs.date,
          'comment': obs.desc,
          'projection': obs.horizDatum,
          'horizontalFix': obs.horizFix,
          'hasResult': {
            'coord': [obs.northing, obs.easting]
          },
          'distanceType': 'icsmdistancetype:ellipsoidal', # from type
          'distanceProvenance': obs.setup.point.state, # Reformat?
          'angleType': obs.arcType, # depending on type, reformat?
          'order': obs.order,
          #'distanceQuality': 'A', # need to model quality!
          #'angleQuality': 'A', # need to model quality!
          #'distanceAccuracy': '???', # need to model accuracy!
          #'angleAccuracy': '???' # need to model accuracy!
         }
    else:
      print("Unexpected observation type!", type(obs))
      return {}
  
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
        'id': "observedVectors",
        'type': "FeatureCollection",
        'featureType': 'surv:ObservedVector',
        'features': [{
          'id': seg.id or i,
          'type': "Feature",
          # Is this the VectorPurpose value?
          'featureType': 'boundary', # FIXME: This will need to support other feature types. What's the logic here?
          'geometry': None,
          'topology': {
            'type': "LineString",
            'references': [ref(seg.start.objID), ref(seg.end.objID)], # TODO: Not all of these IDs would be dereferencable, find them!
            'properties': {
              'comment': ""
            }
          }
        }],
      }
    observedVecs.append(x)  # Alternate algorithm:
  # For reduced observation, select all lines in all parcels where the targetpoint == endpoint.
  # If not more than one select the objID for that line & that becomes hasFeatureOfInterest.
  # If more than one find the instrumentID & select the InstrumentSetup, get the pntRef from the instrumentPoint
  # Select from set of lines the line with the start == that.
  
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
        'geometry': {
          'type': "Point",
          'coordinates': list(trans.transform(*monument.point)) if monument.point is not None else []
        },
        'place': {
          'type': "Point",
          'coordinates': [monument.point.northing, monument.point.easting] if monument.point is not None else []
        },
        'properties': {
          'name': {
            'label': monument.name,
            'hasPart': [
              # TODO: What else is in "hasPart"?
              {'type': "Provenance", 'label': getattr(monument.point, 'state', None)},
              {'type': "Stamp", 'label': monument.type} # Is this right?
            ],
            # ptQuality?
            'purpose': monument.point.purpose,
            'comment': None, # TODO: Capture this data!
            'monumentedBy': {
              'monumentForm': monument.type, # TODO: Namespace
              'monumentCondition': monument.condition, # TODO: Namespace
              'monumentState': monument.state, # TODO: Namespace
            }
          }
        }
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
        'geometry': None,
        'topology': {
          'type': 'Polygon',
          'references': [ref(seg.id) for seg in geom.segments]
        },
        'properties': {
          'name': {
            'label': parcel.name,
            'hasPart': [] # TODO: What are the parts?
          },
          'area': parcel.area,
          'parcelType': parcel.type,
          'state': parcel.state,
          'class': parcel.klass,
          'interestRef': parcel.titleDoc
        }
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
        'geometry': None,
        'properties': observationProperties(obs)
      } for i, obs in enumerate(group)]
    } for groupName, group in data.survey.observationGroups.items()],
    'features': exportJSONfg(data)["features"] # I (Adrian) question of the value of this...
#    'features': [{
#      'type': 'Feature',
#      'id': parcel.name or i,
#      'featureType': 'boundary',
#      'links': [],
#      'time': '', # Appears to be missing from LandXML
#      'place': exportGeom(parcel, 'wgs84'),
#      'geometry': exportGeom(parcel, fileproj),
#      'properties': {
#        'comment': ""
#      } # Place miscallanea here!
#    } for i, parcel in enumerate(data.parcels)]# + [{
#      'type': "Feature",
#      'featureType': "sosa:ObservationCollection",
#      'properties': {
#        'usedProcedure': "icsm:XXX"
#      } # FIXME: This seems incomplete...
#    } for groupID, group in data.survey.observationGroups.items() for i, observation in enumerate(group)]
  }
  json.dump(ret, file, indent=4)

def exportRDF(data, file, format="turtle"):
  g = rdflib.Graph()
  surv = rdflib.Namespace("https://data.surroundaustralia.com/def/csdm/surveyfeatures/")
  geom = rdflib.Namespace("https://data.surroundaustralia.com/def/csdm/geom/")
  commonpatterns = rdflib.Namespace("https://data.surroundaustralia.com/def/csdm/commonpatterns/")
  self = rdflib.Namespace(".")

  for i, point in enumerate(data.points):
    oPt = self["pt" + str(i)]
    oGeom = rdflib.BNode()
    g.add((oPt, RDF.type, surv.surveyPoint))
    g.add((oPt, GEO.hasGeometry, oGeom))
    g.add((oGeom, RDF.type, geom.Point))
    # FIXME: Translate this literal to a URL...
    g.add((oGeom, geom.srsName, rdflib.Literal(point.projection or data.projection.horizontal)))
    g.add((oGeom, geom.srsDimension, rdflib.Literal(2)))
    g.add((oGeom, geom.pos, rdflib.Literal(str(point.northing) + " " + str(point.easting))))
    
    # TODO: Where does surv:hasPointQuality, surv:hasProvenance, & surv:hasPurpose come from?
    g.add((oPt, RDFS.label, rdflib.Literal(point.name)))

  for i, marker in enumerate(data.monuments):
    oMrk = self["m" + str(i)]
    oGeom = rdflib.BNode()
    pt = marker.point
    g.add((oMrk, RDF.type, surv.cadastralMark))
    g.add((oMrk, GEO.hasGeometry, oGeom))
    g.add((oGeom, RDF.type, geom.Point))
    g.add((oGeom, geom.srsName, rdflib.Literal(pt.projection or data.projection.horizontal)))
    g.add((oGeom, geom.srsDimension, rdflib.Literal(2)))
    g.add((oGeom, geom.pos, rdflib.Literal(str(point.northing) + " " + str(point.easting))))

    # TODO: Where does surv:hasPointQuality, surv:hasProvenance, & surv:hasPurpose come from?
    # TODO: How about surv:monumentForm, surv:monumentCondition, surv:markLabel, surv:markState, & surv:hasComment?
    g.add((oMrk, RDFS.label, marker.name))
    # TODO: What boundary marks? How'm I storing that?

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

      # Vertical hangle?
      g.add((oObs, surv.resultTime, rdflib.Literal(obs.date)))
      g.add((oObs, surv.instrumentHeight, rdflib.Literal(obs.setup.height)))
      g.add((oObs, surv.targetHeight, rdflib.Literal(obs.targetSetup.height)))

  # Surveyed vectors?

  g.serialize(file, format=format)

## Utils
def derefJSONLinks(root):
    index = {}
    def indexJSON(data):
        if isinstance(data, list):
            for datum in data: indexJSON(datum)
        elif isinstance(data, dict):
            if "id" in data: index["#" + str(data["id"])] = data
            for datum in data.values(): indexJSON(datum)
    indexJSON(root)

    def expand(data):
        if isinstance(data, list):
            return list(map(expand, data))
        elif isinstance(data, dict):
            if "$ref" in data:
                if data["$ref"] in index:
                    return index[data["$ref"]]
                elif data["$ref"][:2] == "#/":
                    ret = root
                    for prop in data["$ref"][2:].split("/"):
                        if prop in ret: ret = ret[prop]
                        elif isinstance(ret, list):
                            found = False
                            for item in ret:
                                if item.get("id") == prop:
                                    ret = item
                                    found = True
                                    break
                            if not found:
                                if prop.isdigit() and int(prop) < len(ret): ret = ret[int(prop)]
                                else:
                                    print("Failed to dereference path:",
                                            data["$ref"])
                                    return None
                        else:
                            print("Failed to dereference path:", data["$ref"])
                            return None
                    return ret
                else:
                    print("Failed to dereference:", data["$ref"])
                    return None

            for key, datum in data.items():
                data[key] = expand(datum)
            return data
        else:
            return data
    return expand(root)

## Commandline interface

if __name__ == "__main__":
  import argparse

  argparser = argparse.ArgumentParser()
  argparser.add_argument('-X', '--LANDXML', help="LandXML input file", nargs='+')
  argparser.add_argument('-V', '--vocab', help="Jurisdictional vocabulary to align to")
  argparser.add_argument('-T', '--TOPOLOGY', help="JSON-Topology input file", nargs='+')
  argparser.add_argument('--interpolate', help="Interpolate curves, possibly specifying length in meters of each segment (default 1m)",
      const=1, type=int, nargs='?')
  argparser.add_argument('--epsg', help="Overwrite the projection being used.")
  argparser.add_argument('-j', '--jsonfg', help="JSONfg output file", const="?.json", nargs='?')
  argparser.add_argument('-r', '--rdf', help="RDF syntax to output", const="ttl", nargs='?')
  argparser.add_argument('-o', '--output', help="RDF file to output to", const="?.rdf", nargs='?')
  argparser.add_argument('-c', '--csdm', help="CSDM JSON output file", const="?.json", nargs='?')
  #argparser.add_argument('-t', '--topology', help="JSON-Topology output file")
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

  read_input = False
  for source in args.LANDXML or []:
    with open(source) as f: export(source, importLandXML(f))
    read_input = True
  for source in args.TOPOLOGY or []:
    with open(source) as f: export(source, importTopology(f))
    read_input = True

  if not read_input:
    print("Please specify at least one input file!")
    argparser.print_help()
  elif not wrote_output:
    print("Please specify at least one output file!")
    argparser.print_help()
