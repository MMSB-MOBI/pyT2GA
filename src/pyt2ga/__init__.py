import pandas as pd
import json
from .core import computeT2, importdata
from .io import *

import sys

## Accepting stream like object
# as pathToFile, pathToGO, pathToSTRING, pathToSTU

def analysis(proteomicRessource, 
             GoRessource, 
             STRINGRessource, 
             MapperUniprotSTRING, 
             abnd_label = "Corrected_Abundance_Ratio",
             purb=1.5, intg=True, alpha=0.01, ncore=1, sizelim=100):
    """analysis
        Run a T2 analysis 
    """

    defineProteomicRessourcHeader(abnd_label)

    try :
        _                  = pd.read_csv(proteomicRessource, sep='\t')
        proteomicRessource = _
    except Exception as e:
        try :
            assertValidproteomicRessource(proteomicRessource)
        except Exception as e: 
            raise TypeError(f"Can't open proteomic file {e}")
    


    try:
        _               = pd.read_csv(STRINGRessource, sep='\t')
        STRINGRessource = _
    except Exception as e:
        try :
            assertValidSTRINGRessource(STRINGRessource)
        except Exception as e: 
            raise TypeError(f"Can't open STRING {e}")
    
    try:
        _                    = pd.read_csv(MapperUniprotSTRING, sep='\t')
        MapperUniprotSTRING  = _
    except Exception as e:
        try :
            assertValidUniprotStringMapper(MapperUniprotSTRING)
        except Exception as e: 
            raise TypeError(f"Can't open STRING to Uniprot mapper {e}")

    try:
        with open(GoRessource,'r') as json_file: 
            _           = pd.DataFrame.from_dict(json.load(json_file), orient="index")
            GoRessource = _
    except Exception as e:
        try :
            assertValidGoRessource(GoRessource)
            GoRessource = pd.DataFrame.from_dict(GoRessource, orient="index")
        except Exception as e: 
            raise TypeError(f"Can't open GO pathway {e}")

    GO = GoRessource
    GO_vex = pd.DataFrame(
                flatten([ [[j, prot] for prot in GO.loc[j,'Proteins']] for j in GO.index ]),
                          columns=['pathway_id','prot_id'])
    GO_pid = pd.DataFrame([ [i, GO.loc[i,'Name']] for i in GO.index ],
                          columns=['pathway_id','pathway_name'])
    
    dat = importdata(proteomicRessource, abnd_label=abnd_label) # proteomicRessource
    #print("=================================================")
    #print(" Using pathway database: " + pathToGO)
    #print(" Using ppi database:     " + pathToSTRING)
    #print("-------------------------------------------------")
    result = computeT2(dat, GO_vex, GO_pid, 
                       STRINGRessource, 
                       MapperUniprotSTRING, 
                       intg=intg, alpha=alpha, ncore=ncore, sizelim=sizelim)
     
    return toDict(result, GO)

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
                               'p_value'       : pathway[6]#,
                               #'lineage'       : pws.loc[pathway[1], 'Lineage']
                               }
    return theDict