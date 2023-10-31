from data import *
from xml.etree import ElementTree as ET

def importLandXML(file):
  # Look at both CRSs, if survey header has different projection apply to places.
  file = ET.parse(file)

  # NOTE: App data has been discarded from intermediate model...
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
    if center is None:
      centerPt = None # Accept absence
    elif not center.get('pntRef'):
      centerPt = pointsIndex.get(center.get('pntRef'))
    elif center.text.strip() != "":
      centerPt = Point(None, None, None, *[float(x) for x in center.text.strip().split()])
    else:
      centerPt = None
      print("Could not parse center point!")

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
                    centerPt, geoms, titles, address,
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
                                surveyHeaderEl.get('surveyFormat', surveyHeaderEl.get('format')),
                                # TODO: Capture county, desc, endTime, surveyor, surveyorFirm, surveyorReference, class
                                # LINZ has it as SurveyPurpose
                                i(surveyHeaderEl.find('{http://www.landxml.org/schema/LandXML-1.2}PurposeOfSurvey') or surveyHeaderEl.find('{http://www.landxml.org/schema/LandXML-1.2}SurveyPurpose')).get('name', surveyHeaderL.get('surveyPurpose')),
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

def exportLandXML(data, file):
  root = ET.Element("{http://www.landxml.org/schema/LandXML-1.2}LandXML") # FIXME: I forgot to capture the attributes on the root element
  for units in data.units:
    unitsL = ET.SubElement(root, "{http://www.landxml.org/schema/LandXML-1.2}Units")
    ET.SubElement(unitsL, "{http://www.landxml.org/schema/LandXML-1.2}Metric", {
        'areaUnit': units.areaUnit, 'linearUnit': units.linearUnit, 'volumeUnit': units.volumeUnit,
        'temperatureUnit': units.temparatureUnit, 'pressureUnit': units.pressureUnit,
        'angularUnit': units.angularUnit, 'directionUnit': units.directionUnit
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
    'type': data.survey.metadata.type, 'format': data.survey.metadata.format, 'surveyFormat': data.survey.metadata.format,
    'county': data.survey.metadata.jurisdiction, 'jurisdiction': data.survey.metadata.jurisdiction,
    'surveyorReference': data.survey.metadata.personnel, 'personnel': data.survey.metadata.personnel,
    'hadOfPower': data.survey.metadata.headOfPower,
  }) # FIXME: Child nodes?
  ET.SubElement(surveyHeaderL, 'PurposeOfSurvey', {'name': data.survey.metadata.purpose})
  ET.SubElement(surveyHeaderL, 'SurveyPurpose', {'name': data.survey.metadata.purpose})
  ET.SubElement(surveyHeaderL, 'HeadOfPower', {'name': data.survey.metadata.headOfPower})
  for adminArea in data.survey.metadata.adminAreas:
    ET.SubElement(surveyHeaderL, 'AdministrativeArea', {'adminAreaType': adminArea.type, 'adminAreaName': adminArea.name, 'adminAreaCode': adminArea.code})
  for personnel in data.survey.metadata.personnel:
    ET.SubElement(surveyHeaderL, 'Personnel', {'name': personnel.name, 'role': personnel.role, 'regType': personnel.registerType, 'regNumber': personnel.registerID})
  for annotation in data.survey.metadata.annotations:
    ET.SubElement(surveyHeaderL, 'Annotation', {'type': annotation.type, 'name': annotation.name, 'desc': annotation.desc, 'pclRef': annotation.parcel.name})

  for monument in data.survey.monuments:
    ET.SubElement(surveyL, "Monument", {
        # I haven't captured monument IDs!
        'mntRef': id(monument), 'purpose': monument.type, 'state': monument.state, 'name': monument.name, 'pntRef': monument.point.objID, 'type': monument.type, 'condition': monument.condition
    })
  for instrument in data.survey.instruments:
    instrumentL = ET.SubElement(surveyL, "InstrumentSetup", {
        'id': instrument.id, 'stationName': instrument.stationName, 'instrumentHeight': instrument.height
    })
    ET.SubElement(instrumentL, "InstrumentPoint", {'pntRef': instrument.point.objID, 'pointGeometry': 'Point'})
  for name, group in data.observationGroups.items():
    groupL = ET.SubElement(surveyL, "ObservationGroup", {'id': name})
    for observation in group:
      if isinstance(observation, ReducedObservation):
        obsL = ET.SubElement(groupL, "ReducedObservation", {
            'setupID': observation.setup.id, 'azimuth': observation.azimuth, 'horizDistance': observation.horizDist,
            'targetSetupId': observation.targetSetup.id, 'horizDistance': observation.horizDist,
            'equipmentUsed': observation.equipment, 'azimuthType': observation.azimuthType,
            'distanceType': observation.distType, 'date': observation.date, 'azimuthAccuracy': observation.angleAccuracy,
            'distanceAccuracy': observation.distanceAccuracy, 'purpose': observation.purpose, 'name': observation.name
        })
        ET.SubElement(obsL, "TargetPoint", {'pntRef': observation.targetPoint.objID, 'pointGeometry': 'point'})
      elif isinstance(observation, ReducedArcObservation):
        # FIXME: Find a reference!
        ET.SubElement(groupL, "ReducedArcObservation", {
            'purpose': observation.purpose, 'setupID': observation.setup.id, 'targetSetupID': observation.targetSetup.id,
            'chordAzimuth': observation.chordAzimuth, 'radius': observation.radius, 'length': observation.length,
            'rot': "cw" if observation.is_clockwise else "ccw", 'equipmentUsed': observation.equipmentUsed,
            'arcLengthAccuracy': observation.arcLengthAccuracy, 'arcType': observation.arcType,
            'name': observation.name, 'date': observation.date
        })
      elif isinstance(observation, RedHorizPos):
        # FIXME: Find a reference!
        ET.SubElement(groupL, "RedHorizontalPosition", {'desc': observation.desc, 'name': observation.name, 'oID': observation.objID, 'setupID': observation.setup.id, 'date': observation.date, 'horizontalDatum': observation.horizDatum, 'latitude': observation.northing, 'longitude': observation.longitude, 'horizontalFix': observation.horizFix, 'order': observation.order})
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
        'name': point.name, 'pntSurv': point.survey, 'state': point.state, 'oID': point.objID
    }).text = str(point.northing) + " " + str(point.easting)

  def outputGeom(geoms):
    for geom in geoms:
      geomL = ET.SubElement(parcelL, "CoordGeom", {'name': geom.name})
      for seg in geom.segments:
        if isinstance(seg, Line):
          lineL = ET.SubElement(geomL, "Line", {'state': seg.state, 'oID': seg.id})
          ET.SubElement(lineL, "Start", {'pntRef': seg.start.objID, 'pointGeometry': 'point'})
          ET.SubElement(lineL, "End", {'pntRef': seg.end.objID, 'pointGeometry': 'point'})
        elif isinstance(seg, Curve):
          # Is this right?
          curveL = ET.SubElement(geomL, "Curve", {'state': seg.state, 'oID': seg.id,
              'rot': 'cw' if seg.is_clockwise else 'ccw', 'radius': seg.radius})
          ET.SubElement(curveL, "Start", {'pntRef': seg.start.objID, 'pointGeometry': 'point'})
          ET.SubElement(curveL, "Mid", {'pntRef': seg.mid.objID, 'pointGeometry': 'point'})
          ET.SubElement(curveL, "End", {'pntRef': seg.end.objID, 'pointGeometry': 'point'})
        else:
          raise Exception("Unsupported geometry-segment type! " + str(type(seg)))

  parcelsL = ET.SubElement(root, "Parcels")
  for parcel in data.parcels:
    parcelL = ET.SubElement(parcelsL, "Parcel", {
        # Store oIDs?
        'name': parcel.name, 'oID': id(parcel), 'area': parcel.area, 'parcelType': parcel.type, 'class': parcel.klass, 'state': parcel.state, 'parcelFormat': parcel.format, 'desc': parcel.desc
    })
    ET.SubElement(parcelL, "Center").text = str(parcel.center.northing) + " " + str(parcel.center.easting)
    outputGeom(parcel.geom)
    if isinstance(parcel.titleDoc, dict):
      for titleName, titleType in parcel.titleDoc.items():
        ET.SubElement(parcelL, "Title", {'name': titleName, 'titleType': titleType})
    else:
      ET.SubElement(parcelL, "Title", {'name': titleName})

    if parcel.address is not None:
      addressL = ET.SubElement(parcelL, 'LocationAddress', {'addressType': parcel.address.type, 'numberFirst': parcel.address.num})
      ET.SubElement(addressL, 'RoadName', {'roadNameType': parcel.address.roadNameType, 'roadName': parcel.address.roadName, 'roadType': parcel.address.roadType})
      adminArea = parcel.address.adminArea
      if adminArea is not None:
        ET.SubElement(addressL, 'AdministrativeArea', {'adminAreaType': adminArea.type, 'adminAreaName': adminArea.name, 'adminAreaCode': adminArea.code})

  for name, features in data.features.items():
    featuresL = ET.SubElement(root, 'PlanFeatures')
    for feature in features:
      # Is this element name correct?
      ET.SubElement(featuresL, 'PlanFeature', {'name': feature.name, 'desc': feature.desc})
      outputGeom(feature.geom)

  root.write(file)
