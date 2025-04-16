from data import *

def loadWA_CSD(file):
  # NOTE: This code is dramatically incomplete...
  # We're struggling to understand this cadastral format!
  typeHandlers = {}
  def t(i):
    def inner(func):
      typeHandlers[i] = func
      return func
    return inner
  def unsupported(i):
    def inner(*fields):
      print("Ignoring unsupported row of type", i, "with fields", fields)
    return inner

  points = []
  pointsIndex = {}
  pointsQuality = {} # Where do these go? In the observations?
  @t(10)
  def point(id, coord0, coord1, whatsit, quality, class0, class1, name, bool0):
    pt = Point("", id.strip(), "", float(coord0.strip()), float(coord1.strip()), id.strip())
    points.append(pt)
    pointsIndex[id] = id
    pointsQuality[id] = quality
  projection = None
  @t(4)
  def proj(whatsit, proj, *whatsem):
    projection = Projection(proj.strip())
  # What are the other options?

  for line in file:
    line = line.split(",")
    rowtype = int(line[0])
    typeHandlers.get(rowtype, unsupported(rowtype))(*line[1:])

  return Cadastre(projection, {}, [], [], points, [], Survey(None, [], {}))
