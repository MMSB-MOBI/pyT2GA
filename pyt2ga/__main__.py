"""T^2 Genome Analysis.

Usage:
    py2tga <pathToFile> <pathToGO> <pathToString> <pathToSTU> [--purb=purb] [--intg=intg] [--alpha=alpha] [--ncore=ncore] [--sizelim=sizelim]
  
Options:
  -h --help     Show this screen.

"""


from . import analysis

from docopt import docopt

import json

arguments = docopt(__doc__)

ptF, ptG, ptS, ptSTU, purb, intg, alpha, ncore, sizelim = arguments.values()
purb = 1.5 if purb is None else float(purb)
intg = True if intg is None else eval(intg)
alpha = 0.01 if alpha is None else float(alpha)
ncore = 7 if ncore is None else int(ncore)
sizelim = 100 if sizelim is None else int(sizelim)

res = analysis(ptF, ptG, ptS, ptSTU, purb, intg, alpha, ncore, sizelim)

with open('T2GA-analysis.json', 'w') as outfile: 
    json.dump(res, outfile)
# Example:
# python T2GAscript.py './Data/incomplete-wt1.tsv' './Databases/GO_pws.json' './Databases/STRING_v110.tsv' './Databases/STRING_to_Uniprot.tsv' --sizelim=50 --alpha=0.001

