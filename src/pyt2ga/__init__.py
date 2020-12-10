import pandas as pd
import json
from .core import computeT2, importdata

## Accepting stream like object
# as pathToFile, pathToGO, pathToSTRING, pathToSTU

def analysis(pathToFile, 
             pathToGO, 
             pathToSTRING, 
             pathToSTU, 
             purb=1.5, intg=True, alpha=0.01, ncore=7, sizelim=100):
    
    try :
        file = pd.read_csv(pathToFile, sep='\t')
    except Exception as e:
        print (f"Can't open proteomic file {e}")
    
    try:
        STRING = pd.read_csv(pathToSTRING, sep='\t')
    except Exception as e:
        print(f"Can't open STRING {e}")
    
    try:
        STU = pd.read_csv(pathToSTU, sep='\t')
    except Exception as e:
        print(f"Can't open STRING to Uniprot mapper {e}")

    try:
        with open(pathToGO,'r') as json_file: 
            GO = pd.DataFrame.from_dict(json.load(json_file), orient="index")
    except Exception as e:
        print(f"Can't open GO pathway {pathToGO} {e}")

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
    
    return result





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