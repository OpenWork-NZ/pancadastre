from data import *
from xml.etree import ElementTree as ET

def importCityGML(file):
    NS = '{http://www.opengis.net/citygml/2.0}'
    data = ET.parse(file).getroot()

    parcels = []
    for member in data.findall(NS + 'cityObjectMember') or []:
        BLD = '{http://www.opengis.net/citygml/building/2.0}'
        for building in member.findall(BLD + 'BuildingPart') or []:
            area = None
            for measure in building.findall('{http://www.opengis.net/citygml/generics/2.0}measureAttribute') or []:
                if measure.get(NS + 'name') != "GrossPlannedArea": continue
                area = measure.find('{http://www.opengis.net/citygml/generics/2.0}value').text
            geom = []
            for solid in building.findall(BLD + 'lod2Solid') or []:
                for solidb in solid.findall('{http://www.opengis.net/gml}Solid') or []:
                    segments = []
                    for exterior in solidb.findall('{http://www.opengis.net/gml}exterior'):
                        for composite in exterior.findall('{http://www.opengis.net/gml}CompositeSurface'):
                            for surface in composite.findall('{http://www.opengis.net/gml}surfaceMember'):
                                for polygon in surface.findall('{http://www.opengis.net/gml}Polygon'):
                                    start = None
                                    for exteriorb in polygon.findall('{http://www.opengis.net/gml}exterior'):
                                        for ring in exteriorb.findall('{http://www.opengis.net/gml}LinearRing'):
                                            for pointProp in ring.findall('{http://www.opengis.net/gml}pointProperty'):
                                                point = pointProp.find('{http://www.opengis.net/gml}Point').find('{http://www.opengis.net/gml}pos').text.split()
                                                point = Point(None, None, None, *map(float, point))
                                                if start is not None:
                                                    segments.append(Line(None, start, point))
                                                start = point
                    geom.append(Geom(solidb.get('name'), segments))
            parcel = Parcel(building.find('{http://www.opengis.net/gml}name').text, area, building.find(NS + 'parcelType').text,
                    None, building.find(BLD + 'class').text, None, None, geom, None, 
                    desc = building.find('{http://www.opengis.net/gml}description').text, oid = building.find(NS + 'parcelID').text)
            parcels.append(parcel)

    return Cadastre(None, None, None, None, None, parcels, None)
