# pancadastre
Convert between different cadastre fileformats. Can be used as either a Python library or a commandline program.

## Example data
`example.json` & `example.txt` were generated from `example.json.input` using the command:

		./pancadastre.py -C example.json.input -j -s

`example.json.input` was taken from `../3d-csdm-schema/_sources/csdm/features/CSD/example.json`.

## Dependencies
This program requires Python3, PyProj, & RDFLib. PyProj & RDFLib can be installed using: `pip install rdflib pyproj`.

## Usage
Pancadastre can be run as a commandline program with inputs specified using uppercase flags & outputs specified using lowercase flags. Each input will be written to all outputs, substituting the input's basename in place of questionmarks (`?`) in the output filepaths. If no outputs are given then it validates the input files without writing them out.

Some flags like `--interpolate` & `--epsg` mutates the data before writing it out.

### Commandline Options
usage: pancadastre.py [-h] [-X LANDXML [LANDXML ...]] [-V VOCAB] [-C CSDM [CSDM ...]] [--interpolate [INTERPOLATE]] [--epsg EPSG] [-j [JSONFG]] [-r [RDF]]
                      [-o [OUTPUT]] [-c [CSDM]]

options:
  -h, --help            show this help message and exit
  -X LANDXML [LANDXML ...], --LANDXML LANDXML [LANDXML ...]
                        LandXML input file
  -V VOCAB, --vocab VOCAB
                        Jurisdictional vocabulary to align to
  -C CSDM [CSDM ...], --CSDM CSDM [CSDM ...]
                        JSON-Topology input file
  --interpolate [INTERPOLATE]
                        Interpolate curves, possibly specifying length in meters of each segment (default 1m)
  --epsg EPSG           Overwrite the projection being used.
  -j [JSONFG], --jsonfg [JSONFG]
                        JSONfg output file
  -r [RDF], --rdf [RDF]
                        RDF syntax to output
  -o [OUTPUT], --output [OUTPUT]
                        RDF file to output to
  -c [CSDM], --csdm [CSDM]
                        CSDM JSON output file
