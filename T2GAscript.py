# #!/usr/bin/env python
# #coding: utf-8

"""T^2 Genome Analysis.

Usage:
  T2GA.py <pathToFile> <pathToGO> <pathToString> <pathToSTU> [--purb=purb] [--intg=intg] [--alpha=alpha] [--ncore=ncore] [--sizelim=sizelim]
  
Options:
  -h --help     Show this screen.

"""

from T2GAlib import *
from docopt import docopt


def flatten(l):
    return [item for sublist in l for item in sublist]

def toDict(df, pws):
    theDict = {}
    for pathway in df.values:
        theDict[pathway[1]] = {'pathway_title' : pathway[0],
                               'uniprot_ids'   : pathway[2],
                               'members'       : pathway[3],
                               'matrix_rank'   : pathway[4],
                               't_square'      : pathway[5],
                               'p_value'       : pathway[6],
                               'lineage'       : pws.loc[pathway[1], 'Lineage']}
    return theDict

def analysis(pathToFile, pathToGO, pathToSTRING, pathToSTU, purb=1.5, intg=True, alpha=0.01, ncore=7, sizelim=100):
    
    file = pd.read_csv(pathToFile, sep='\t')
    STRING = pd.read_csv(pathToSTRING, sep='\t')
    STU = pd.read_csv(pathToSTU, sep='\t')
    with open(pathToGO,'r') as json_file: GO = pd.DataFrame.from_dict(json.load(json_file), orient="index")

    GO_vex = pd.DataFrame(flatten([ [[j, prot] for prot in GO.loc[j,'Proteins']] for j in GO.index ]),
                          columns=['pathway_id','prot_id'])
    GO_pid = pd.DataFrame([ [i, GO.loc[i,'Name']] for i in GO.index ],
                          columns=['pathway_id','pathway_name'])
    
    dat = importdata(file)
    print("=================================================")
    print(" Using pathway database: " + pathToGO)
    print(" Using ppi database:     " + pathToSTRING)
    print("-------------------------------------------------")
    result = computeT2(dat, GO_vex, GO_pid, STRING, STU, intg=intg, alpha=alpha, ncore=ncore, sizelim=sizelim)
    result = toDict(result, GO)
    
    with open('./Out/T2GA-analysis.json', 'w') as outfile: 
        json.dump(result, outfile)
    
    
if __name__ == '__main__':
    arguments = docopt(__doc__)
    print(arguments)
    ptF, ptG, ptS, ptSTU, purb, intg, alpha, ncore, sizelim = arguments.values()
    purb = 1.5 if purb is None else float(purb)
    intg = True if intg is None else eval(intg)
    alpha = 0.01 if alpha is None else float(alpha)
    ncore = 7 if ncore is None else int(ncore)
    sizelim = 100 if sizelim is None else int(sizelim)
    
    analysis(ptF, ptG, ptS, ptSTU, purb, intg, alpha, ncore, sizelim)

# Example:
# python T2GAscript.py './Data/incomplete-wt1.tsv' './Databases/GO_pws.json' './Databases/STRING_v110.tsv' './Databases/STRING_to_Uniprot.tsv' --sizelim=50 --alpha=0.001

