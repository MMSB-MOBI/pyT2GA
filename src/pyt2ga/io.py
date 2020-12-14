
PROTEOMIC_RESSOURCE_HEADER = ['Accession']#, 'Corrected_Abundance_Ratio',
       #'LOG2(Corrected_Abundance_Ratio)', 'LOG10(Adj.P-val)']
UNIPROT_STRING_MAPPER_RESSOURCE_HEADER = [ 'String_id', 'Uniprot_id' ]
STRING_RESSOURCE_HEADER = ['protein1', 'protein2', 'experimental']


def defineProteomicRessourcHeader(columnLabel):
    global PROTEOMIC_RESSOURCE_HEADER

    PROTEOMIC_RESSOURCE_HEADER = list( set(PROTEOMIC_RESSOURCE_HEADER) | set([columnLabel]) )

def assertValidproteomicRessource(data):
    try:
        data.columns
    except AttributeError as e:
        raise ValueError("proteomic Ressource is not a panda frame")
    
    missingHead = set(PROTEOMIC_RESSOURCE_HEADER) - set(data.columns)
    
    if len(missingHead) > 0:
        raise ValueError(f"Following columns are missing in proteomic Ressource {missingHead}")
    
    return True

def assertValidGoRessource(data):
    try :
        for k,v in data.items():
            if not k.startswith("GO:"):
                raise ValueError(f"Invalid key in GO ressource{k}")
        if (not 'Proteins' in v) or (not 'Name' in v):
            raise ValueError(f"Invalid pathway key in GO ressource {v}")

    except AttributeError as e:
        raise ValueError("Go Ressource is not a dictionary")
    return True

def assertValidSTRINGRessource(data):
    try:
        data.columns
    except AttributeError as e:
        raise ValueError("STRING Ressource is not a panda frame")
    
    missingHead = set(STRING_RESSOURCE_HEADER) - set(data.columns)
    if len(missingHead) > 0:
        raise ValueError(f"Following columns are missing in STRING ressource {missingHead}")
    
    return True

def assertValidUniprotStringMapper(data):
    try:
        data.columns
    except AttributeError as e:
        raise ValueError("UniprotStringMapper Ressource is not a panda frame")
    
    missingHead = set(UNIPROT_STRING_MAPPER_RESSOURCE_HEADER) - set(data.columns)
    if len(missingHead) > 0:
        raise ValueError(f"Following columns are missing in UniprotStringMapper {missingHead}")
    
    return True