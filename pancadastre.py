#! /usr/local/bin/python3

from xml.etree import ElementTree as ET
import rdflib
from rdflib.namespace import RDFS, RDF, GEO
import json

from datetime import datetime
from pyproj import Transformer, Geod

from math import pi, isclose, ceil, sqrt
tau = 2*pi

from data import *

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
              datetime.fromisoformat(observation.get('date')),
              observation.get('horizontalDatum'),
              float(observation.get('latitude')), float(observation.get('longitude')),
              observation.get('horizontalFix'), int(observation.get('order'))))
        else:
          print("Unsupported tag!", observation.tag)
      observationGroups[el.get('id')] = observations
    survey = Survey(surveyMeta, instruments, observationGroups)

  return Cadastre(proj, features, units, monuments, points, parcels, survey)

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
      point = Point(None, None, monument["properties"]["monumentedBy"]["state"], monument["geometry"]["coordinates"][0], monument["geometry"]["coordinates"][1], monument["id"])
      points.append(point)
      pointsIndex[monument["id"]] = deepcopy(monument)
      pointsIndex[monument["id"]][""] = point
      # The observations will be initialized later
      monuments.append(Monument(monument["properties"]["name"]["label"], point, monument["properties"]["monumentedBy"]["state"], monument["properties"]["monumentedBy"]["form"], monument["properties"]["monumentedBy"]["condition"], properties = monument["properties"]))

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
          obs = ReducedObservation(None, instrument(start["id"], start["properties"]["name"]["label"], None, start[""]), instrument(end["id"], end["properties"]["name"]["label"], None, end[""]), measure.azimuth, measure.dist, measure.equipment, measure.distType, measure.azimuthType, observation["id"], start["time"], start["properties"], measure = measure)
          observations.append(obs)

          geom = Line(observation["id"], start[""], end[""])
          segments[observation["id"]] = geom
      elif observation["topology"]["type"].lower() == "arc":
        start = pointsIndex[d(observation["topology"]["references"][0])]
        mid = pointsIndex[d(observation["topology"]["references"][1])]
        end = pointsIndex[d(observation["topology"]["references"][2])]
        # Distance, radius, etc comes from vector observations.
        obs = ReducedArcObservation(None, instrument(start["id"], start["properties"]["name"]["label"], None, start[""]), instrument(end["id"], end["properties"]["name"]["label"], None, end[""]), measure.azimuth, measure.radius, measure.dist, measure.is_clockwise, measure.equipment, measure.angleAccuracy, measure.arcType, observation["id"], start["time"], start["properties"], measure.distanceQuality, measure.distanceAccuracy, measure.angleQuality, measure = measure)
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

def exportLandXML(data, file):
  root = ET.Element("{http://www.landxml.org/schema/LandXML-1.2}LandXML") # FIXME: I forgot to capture the attributes on the root element
  if data.units:
    unitsL = ET.SubElement(root, "{http://www.landxml.org/schema/LandXML-1.2}Units")
    ET.SubElement(unitsL, "{http://www.landxml.org/schema/LandXML-1.2}Metric", {
        'areaUnit': data.units.areaUnit, 'linearUnit': data.units.linearUnit, 'volumeUnit': data.units.volumeUnit,
        'temperatureUnit': data.units.temparatureUnit, 'pressureUnit': data.units.pressureUnit,
        'angularUnit': data.units.angularUnit, 'directionUnit': data.units.directionUnit
    })
  ET.SubElement(root, "{http://www.landxml.org/schema/LandXML-1.2}CoordinateSystem", {
      # NOTE: Inconsistency in how this is denoted, outputting a mixed schema.
      'horizontalDatum': data.proj.horizontal, 'datum': data.proj.vertical, 'name': data.proj.horizontal
  })

  surveyL = ET.SubElement(root, "Survey")
  surveyHeaderL = ET.SubElement(surveyL, "SurveyHeader", {
    # TODO: endTime attribute? surveyorReference? class?
    'name': data.survey.metadata.name, 'surveyPurpose': data.survey.metadata.purpose,
    'surveyor': data.survey.metadata.surveyor, 'surveyorFirm': data.survey.metadata.firm,
    'type': data.survey.metadata.type, 'format': data.survey.metadata.format,
    'county': data.survey.metadata.jurisdiction, 'jurisdiction': data.survey.metadata.jurisdiction,
    'surveyorReference': data.survey.metadata.personnel, 'personnel': data.survey.metadata.personnel,
    'hadOfPower': data.survey.metadata.headOfPower, 
  }) # Child nodes?
  for monument in data.survey.monuments:
    ET.SubElement(surveyL, "Monument", {
        # I haven't captured monument IDs!
        'mntRef': id(monument), 'purpose': monument.type, 'state': monument.state
    })
  for instrument in data.survey.instruments:
    instrumentL = ET.SubElement(surveyL, "InstrumentSetup", {
        'id': instrument.id, 'stationName': instrument.stationName, 'instrumentHeight': instrument.height
    })
    ET.SubElement(instrumentL, "InstrumentPoint", {'pntRef': instrument.point.objID, 'pointGeometry': 'Point'})
  for name, group in data.observationGroups.items():
    groupL = ET.SubElement(surveyL, "ObservationGroup", {'id': name}):
    for observation in group:
      if isinstance(observation, ReducedObservation):
        obsL = ET.SubElement(groupL, "ReducedObservation", {
            'setupID': observation.setup.id, 'azimuth': observation.azimuth, 'horizDistance': observation.horizDist,
            'equipmentUsed': observation.equipment, 'azimuthType': observation.azimuthType,
            'distanceType': observation.distType, 'date': observation.date, 'azimuthAccuracy': observation.angleAccuracy,
            'distanceAccuracy': observation.distanceAccuracy
        })
        ET.SubElement(obsL, "TargetPoint" {'pntRef': observation.targetPoint.objID, 'pointGeometry': 'point'})
      elif isinstance(observation, ReducedArcObservation):
        # FIXME: Find a reference!
        obsL = ET.SubElement(groupL, "ReducedArcObservation", {
            'purpose': observation.purpose, 'setupID': observation.setup.id, 'targetSetupID': observation.targetSetup.id,
            'chordAzimuth': observation.chordAzimuth, 'radius': observation.radius, 'length': observation.length,
            'rot': "cw" if observation.is_clockwise else "ccw", 'equipmentUsed': observation.equipmentUsed,
            'arcLengthAccuracy': observation.arcLengthAccuracy, 'arcType': observation.arcType,
            'name': observation.name, 'date': observation.date
        })
      elif isinstance(observation, RedHorizPos):
        ... # Most LandXML variants don't use this.
      else:
        raise Exception("Unsupported observation type! " + str(type(observation)))

  monumentsL = ET.SubElement(root, "Monuments")
  for monument in data.survey.monuments:
    ET.SubElement(monumentsL, "Monument", {
        'name': monument.name, 'pntRef': monument.point.objID, 'state': monument.state, 'type': monument.type,
        'condition': monument.condition, 'oID': id(monument)
    })

  cgpointsL = ET.SubElement(root, "CgPoints")
  for point in data.points:
    ET.SubElement(cgpointsL, "CgPoint", {
        # surveyOrder?
        'name': point.name, 'pntSurv': point.survey, 'oID': point.objID
    }).text = str(point.northing) + " " + str(point.easting)

  parcelsL = ET.SubElement(root, "Parcels")
  for parcel in data.parcels:
    parcelL = ET.SubElement(parcelsL, "Parcel", {
        # Store oIDs?
        'name': parcel.name, 'oID': id(parcel), 'area': parcel.area, 'parcelType': parcel.type, 'class': parcel.klass
    })
    ET.SubElement(parcelL, "Center").text = str(parcel.center.northing) + " " str(parcel.center.easting)
    for geom in parcel.geom:
      geomL = ET.SubElement(parcelL, "CoordGeom", {'name': geom.name})
      for seg in geom.segments:
        if isinstance(seg, Line):
          lineL = ET.SubElement(geomL, "Line", {'state': seg.state, 'oID': id(seg)})
          ET.SubElement(lineL, "Start", {'pntRef': seg.start.objID, 'pointGeometry': 'point'})
          ET.SubElement(lineL, "End", {'pntRef': seg.end.objID, 'pointGeometry': 'point'})
        elif isinstance(seg, Curve):
          # Is this right?
          curveL = ET.SubElement(geomL, "Curve", {'state': seg.state, 'oID': id(seg),
              'rot': 'cw' if seg.is_clockwise else 'ccw', 'radius': seg.radius})
          ET.SubElement(curveL, "Start", {'pntRef': seg.start.objID, 'pointGeometry': 'point'})
          ET.SubElement(curveL, "Mid", {'pntRef': seg.mid.objID, 'pointGeometry': 'point'})
          ET.SubElement(curveL, "End", {'pntRef': seg.end.objID, 'pointGeometry': 'point'})
        else:
          raise Exception("Unsupported geometry-segment type! " + str(type(seg)))
    ET.SubElement(parcelL, "Title", {'name': parcel.titleDoc})

  root.write(file)
  

# TODO: Tweak according to surveyfeatures.json, circular_arc.json
def exportJSONfg(data, file = None):
  import json
  def simplifyPoly(trans, lines):
    print(trans)
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
        'geometry': {
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
        'geometry': None,
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
        'geometry': None,
        'properties': obs.properties
      } for i, obs in enumerate(group)]
    } for groupName, group in data.survey.observationGroups.items()],
    'features': exportJSONfg(data)["features"] # I (Adrian) question of the value of this...
  }
  json.dump(ret, file, indent=4)

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
      rets.append("CURVE( " + ptr(seg.start) + ", " + pt(seg.mid) + ", " pt(seg.end) + " )")
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
    ret += "\t" + monument.name + "\t" + str(monument.point) + "\t" + monument.state + "\t" + monument.condition + "\n"
  ret += "Points:\t" + str(len(data.points)) + "\n"
  for point in data.points or []:
    ret += "\t" + str(point.name) + "\t" + str(point.state) + "\t" + str(point) + "\n"
  ret += "Parcels:\t" + str(len(data.parcels)) + "\n"
  for parcel in data.parcels or []:
    ret += "\t" + str(parcel.name) + "\t" + str(parcel.type) + "\t" + str(parcel.area) + "\t" + str(parcel.center) + "\t" + \
        str(parcel.state) + "\t" + str(parcel.klass) + "\t" + str(parcel.comment) + "\n"
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

## Commandline interface

if __name__ == "__main__":
  import argparse

  argparser = argparse.ArgumentParser()
  argparser.add_argument('-X', '--LANDXML', help="LandXML input file", nargs='+')
  argparser.add_argument('-V', '--vocab', help="Jurisdictional vocabulary to align to")
  argparser.add_argument('-C', '--CSDM', help="JSON-Topology input file", nargs='+')
  argparser.add_argument('--interpolate', help="Interpolate curves, possibly specifying length in meters of each segment (default 1m)",
      const=1, type=int, nargs='?')
  argparser.add_argument('--epsg', help="Overwrite the projection being used.")
  argparser.add_argument('-j', '--jsonfg', help="JSONfg output file", const="?.json", nargs='?')
  argparser.add_argument('-r', '--rdf', help="RDF syntax to output", const="ttl", nargs='?')
  argparser.add_argument('-o', '--output', help="RDF file to output to", const="?.rdf", nargs='?')
  argparser.add_argument('-c', '--csdm', help="CSDM JSON output file", const="?.json", nargs='?')
  argparser.add_argument('-s', '--summary', help="Textual summary of parsed data", const="?.txt", nargs='?')
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
    if args.summary:
      with open(args.summary.replace("?", filename), "w") as f: f.write(exportSummary(data))

  read_input = False
  for source in args.LANDXML or []:
    with open(source) as f: export(source, importLandXML(f))
    read_input = True
  for source in args.CSDM or []:
    with open(source) as f: export(source, importCSDM(f))
    read_input = True

  if not read_input:
    print("Please specify at least one input file!")
    argparser.print_help()
