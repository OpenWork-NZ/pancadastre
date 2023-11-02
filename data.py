from datetime import datetime
from math import isclose, sqrt

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
    self.areaUnit = areaUnit
    self.linearUnit = linearUnit
    self.volumeUnit = volumeUnit
    self.temperatureUnit = temperatureUnit
    self.pressureUnit = pressureUnit
    self.angularUnit = angularUnit
    self.directionUnit = directionUnit

class Monument:
  def __init__(self, name, point, state, type, condition, properties = None, oid = None):
    self.name = name
    self.point = point
    self.state = state
    self.type = type
    self.condition = condition
    self.oid = oid
    self.properties_ = properties
    self.point.monuments.append(self)
    # Modelling remonumenting?

  @property
  def properties(self):
    if self.properties_ is not None: return self.properties_
    return {
      'name': {
        'label': self.name,
        'hasPart': [
          # TODO: What else is in "hasPart"?
          {'type': "Provenance", 'label': getattr(self.point, 'state', None)},
          {'type': "Stamp", 'label': self.type} # Is this right?
        ]
      },
      'ptQuality': self.point.observation.distanceQuality
          if self.point.observation is not None else None,
      'purpose': self.point.purpose,
      'comment': None, # TODO: Capture this data!
      'monumentedBy': {
        'form': self.type, # TODO: Namespace
        'condition': self.condition, # TODO: Namespace
        'state': self.state, # TODO: Namespace
      }
    }

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
    yield self.northing
    yield self.easting

  def __eq__(self, other):
    return isclose(self.northing, other.northing) and isclose(self.easting, other.easting)

  def __str__(self):
    return str(self.northing) + "," + str(self.easting)

  def __repr__(self):
    return repr(self.northing) + "," + repr(self.easting)

  def __add__(self, other):
    return Point(self.survey, self.name, self.state, self.northing + other.northing, self.easting + other.easting, str(self.objID) + "+" + str(other.objID))
  def __sub__(self, other):
    return Point(self.survey, self.name, self.state, self.northing - other.northing, self.easting - other.easting, str(self.objID) + "-" + str(other.objID))
  def div(self, other = 2):
    return Point(self.survey, self.name, self.state, self.northing/other, self.easting/other, str(self.objID) + "/" + str(other))

class Parcel:
  def __init__(self, name, area, type, state, klass, format, center, geom, titleDoc, address=None, desc = None, properties = None):
    props = properties or {}
    self.name = name or props.get('name', {}).get('label')
    self.area = area or props.get('area')
    self.type = type or props.get('parcelType')
    self.state = state or props.get('state')
    self.klass = klass or props.get('class')
    self.desc = desc or props.get('comment')
    self.format = format
    self.center = center # Do we want?
    self.geom = geom
    self.titleDoc = titleDoc or props.get('interestRef') # Contains LandXML Title element or CSDM interestRef property.
    self.address = address # Do we want?
    self.properties = properties
    if properties is None: self.populateProperties()

  @staticmethod
  def fromProperties(format, center, geom, properties, address = None):
    return Parcel(None, None, None, None, None, format, center, geom, None, address = address, properties = properties)

  def populateProperties(self):
    self.properties['name'] = {'label': self.name, 'hasPart': []}
    self.properties['area'] = self.area
    self.properties['parcelType'] = self.type
    self.properties['state'] = self.state
    self.properties['class'] = self.klass
    self.properties['interestRef'] = self.titleDoc
    pass

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
    # Handle missing inputs
    if is_clockwise is None: is_clockwise = True # Better default?
    if r is None or isclose(r, 0):
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

class Measure:
  def __init__(self, equipment, featureOfInterest, dist, distType, azimuth, azimuthType, date, distanceQuality, distanceAccuracy, angleQuality, angleAccuracy, radius = None, is_clockwise = None, arcType = None, properties = None):
    props = properties or {}
    result = props.get("hasResult", {})
    self.equipment = equipment
    self.featureOfInterest = featureOfInterest or props.get("hasFeatureOfInterest")
    self.dist = dist or result.get("distance") or result.get("arcLength")
    self.distType = distType or props.get("distanceType")
    self.azimuth = azimuth or result.get("angle") or result.get("chordAngle")
    self.azimuthType = azimuthType # What should this be?
    self.date = date or props.get("resultTime")
    self.distanceQuality = distanceQuality or props.get("distanceQuality")
    self.distanceAccuracy = distanceAccuracy or props.get("distanceAccuracy")
    self.angleQuality = angleQuality or props.get("angleQuality")
    self.angleAccuracy = angleAccuracy or props.get("angleAccuracy")
    self.radius = radius or result.get("radius")
    self.is_clockwise = is_clockwise or result.get("orientation") == "cw"
    self.arcType = arcType or props.get("angleType")
    self.properties = props
    if properties is None: self.populateProperties()

  @staticmethod
  def fromProperties(equipment, properties):
    return Measure(equipment, None, None, None, None, None, None, None, None, None, None, properties = properties)

  def populateProperties(self):
    self.properties["hasFeatureOfInterest"] = self.featureOfInterest
    if not isinstance(self.properties["hasResult"], dict): self.properties["hasResult"] = {}
    self.properties["hasResult"]["distance"] = self.properties["hasResult"]["arcLength"] = self.dist
    self.properties["distanceType"] = self.distType
    self.properties["hasResult"]["angle"] = self.properties["hasResult"]["chordAngle"] = self.azimuth
    self.properties["resultTime"] = self.date
    self.properties["distanceQuality"] = self.distanceQuality
    self.properties["distanceAccuracy"] = self.distanceAccuracy
    self.properties["angleQuality"] = self.angleQuality
    self.properties["angleAccuracy"] = self.angleAccuracy
    self.properties["radius"] = self.radius
    self.properties["orientation"] = "cw" if self.is_clockwise else "ccw"
    self.properties["angleType"] = self.arcType

class Observation:
  pass

class ReducedObservation(Observation):
  def __init__(self, purpose, setup, targetSetup, azimuth, horizDist, equipment, distType, azimuthType, name, date, properties = None, distanceQuality = None, distanceAccuracy = None, angleQuality = None, angleAccuracy = None, measure = None):
    # TODO: Fully split out measurement properties from Reduced Observations...
    props = properties or {}
    self.purpose = purpose or props.get("purpose")
    self.setup = setup
    self.targetSetup = targetSetup
    self.azimuth = azimuth
    self.horizDist = horizDist
    self.equipment = equipment
    self.distType = distType
    self.azimuthType = azimuthType
    self.name = name or props.get("name", {}).get("label")
    self.date = date
    self.properties = properties or {}
    self.distanceQuality = distanceQuality
    self.distanceAccuracy = distanceAccuracy
    self.angleQuality = angleQuality
    self.angleAccuracy = angleAccuracy
    self.setup.observations.append(self)
    if self.targetSetup is not None:
      self.targetSetup.observations.append(self)
    if properties is None: self.populateProperties()
    self.measure = measure

  @property
  def setupPoint(self): return self.setup.point

  @property
  def targetPoint(self): return self.targetSetup.point

  @property
  def line(self):
    for segment in self.setupPoint.segments:
      if self.setupPoint == segment.start and self.targetPoint == segment.end: return segment
      if self.targetPoint == segment.start and self.setupPoint == segment.end: return segment

  def populateProperties(self):
    # Start leaving hasPart blank...
    self.properties['name'] = {'label': self.name, 'hasPart': []}
    self.properties['ptQuality'] = self.point.observation.distanceQuality
    self.properties['purpose'] = self.purpose
    self.properties['comment'] = None
    for monument in self.setupPoint.monuments:
      self.properties['monumentedBy'] = {
        'form': monument.type,
        'condition': monument.condition,
        'state': monument.state
      }
      break

class ReducedArcObservation(Observation):
  def __init__(self, purpose, setup, targetSetup, chordAzimuth, radius, length, is_clockwise, equipmentUsed, arcLengthAccuracy, arcType, name, date, properties = None, distanceQuality = None, distanceAccuracy = None, angleQuality = None, measure = None):
    props = properties or {}
    self.purpose = purpose or props.get("purpose")
    self.setup = setup
    self.targetSetup = targetSetup
    self.chordAzimuth = chordAzimuth or props.get("hasResult", {}).get("angle")
    self.radius = radius or props.get("hasResult", {}).get("radius")
    self.length = length or props.get("hasResult", {}).get("distance")
    self.is_clockwise = is_clockwise or props.get("rot") == 'cw'
    self.equipmentUsed = equipmentUsed or props.get("equipmentUsed")
    self.arcLengthAccuracy = arcLengthAccuracy or props.get("angleAccuracy")
    self.arcType = arcType or props.get("angleType")
    self.name = name
    self.date = date or (props.get("resulttime") and datetime.fromisoformat(props["resultTime"]))
    self.properties = properties or {}
    self.distanceQuality = distanceQuality or props.get("distanceQuality")
    self.distanceAccuracy = distanceAccuracy or props.get("distanceAccuracy")
    self.angleQuality = angleQuality or props.get("angleQuality")
    self.setup.observations.append(self)
    if self.targetSetup is not None:
      self.targetSetup.observations.append(self)
    if self.properties == {}: self.populateProperties()
    self.measure = measure
    self.geom = None # For the CSDM schema
    self.center = None
    self.mid = None

  @property
  def setupPoint(self): return self.setup.point

  @property
  def targetPoint(self): return self.targetSetup.point

  @property
  def line(self):
    for segment in self.setupPoint.segments:
      if self.setupPoint == segment.start and self.targetPoint == segment.end: return segment
      if self.targetPoint == segment.start and self.setupPoint == segment.end: return segment

  def populateProperties(self):
    self.properties['purpose'] = self.purpose
    self.properties['hasFeatureOfInterest'] = self.setup.point.objID
    self.properties['resultTime'] = self.date
    self.properties['hasResult'] = {
        'distance': self.length, 'angle': self.chordAzimuth, 'radius': self.radius
    }
    self.properties['rot'] = "cw" if obs.is_clockwise else "ccw"
    self.properties['distanceType'] = 'icsmdistancetype:ellipsoidal' # What should this be?
    self.properties['distanceProvenance'] = self.setup.point.state # Reformat?
    self.properties['angleType'] = self.arcType # Reformat?
    self.properties['equipmentUsed'] = self.equipmentUsed
    self.properties['distanceQuality'] = self.distanceQuality
    self.properties['angleQuality'] = self.angleQuality
    self.properties['distanceAccuracy'] = self.distanceAccuracy
    self.properties['angleAccuracy'] = self.arcLengthAccuracy

class RedHorizPos(Observation): # Spell out "Reduced"
  def __init__(self, desc, name, objID, setup, date, horizDatum, northing, easting, horizFix, order, properties = None, distanceQuality = None, distanceAccuracy = None, angleQuality = None, angleAccuracy = None):
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
    self.properties = properties or {}
    self.distanceQuality = distanceQuality
    self.distanceAccuracy = distanceAccuracy
    self.angleQuality = angleQuality
    self.angleAccuracy = angleAccuracy
    self.setup.observations.append(self)
    self.line = None
    if self.properties == {}: self.populateProperties()
    self.purpose = "" # Easier if all 3 have this property...

  def populateProperties(self):
    self.properties['hasFeatureOfInterest'] = self.setup.point.objID # FIXME: Parse association!
    self.properties['resultTime'] = self.date
    self.properties['comment'] = self.desc
    self.properties['projection'] = self.horizDatum
    self.properties['horizontalFix'] = self.horizFix
    self.properties['hasResult'] = {'coord': [self.northing, self.easting]}
    self.properties['distanceType'] = 'icsmdistancetype:ellipsoidal'
    self.properties['distanceProvenance'] = self.setup.point.state # Reformat?
    self.properties['angleType'] = self.arcType # Reformat?
    self.properties['order'] = self.order
    self.properties['distanceQuality'] = self.distanceQuality # TODO: Thesaurus
    self.properties['angleQuality'] = self.angleQuality
    self.properties['distanceAccuracy'] = self.distanceAccuracy
    self.properties['angleAccuracy'] = self.angleAccuracy
