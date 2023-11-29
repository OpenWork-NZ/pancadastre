from data import *

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
    ret += "\t" + str(monument.oid) + "\t" + str(monument.name) + "\t" + str(monument.point) + "\t" + str(monument.state) + "\t" + str(monument.condition) + "\n"
  # Since CSDM schema doesn't differentiate points from lines,
  # This is a reflection of LandXML.
#  ret += "Points:\t" + str(len(data.points)) + "\n"
#  for point in data.points or []:
#    ret += "\t" + str(point.name or point.objID) + "\t" + str(point.state) + "\t" + str(point) + "\n"
  ret += "Parcels:\t" + str(len(data.parcels)) + "\n"
  for parcel in data.parcels or []:
    ret += "\t" + str(parcel.oid) + "\t" + str(parcel.name) + "\t" + str(parcel.type) + "\t" + str(parcel.area) + "\t" + str(parcel.center) + "\t" + \
        str(parcel.state) + "\t" + str(parcel.klass) + "\t" + str(parcel.desc) + "\n"
  ret += "Observations:\t" + str(len([y for x in data.survey.observationGroups.values() for y in x])) + "\n"
  for label, observations in data.survey.observationGroups.items():
    ret += "\t\t" + label + "\n"
    for observation in observations:
      # TODO: Capture input's object IDs to echo here...
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
