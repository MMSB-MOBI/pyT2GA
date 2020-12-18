import pandas as pd
import numpy as np
import json
from scipy.stats import chi2
from numpy import linalg as la

import multiprocessing

def nearestPD(A):
    """Find the positive-definite matrix nearest to input.
    
    Written by Ahmed Fasih: https://stackoverflow.com/users/500207/ahmed-fasih.

    A Python/Numpy port of John D'Errico's `nearestSPD` MATLAB code [1], which
    credits [2].

    [1] https://www.mathworks.com/matlabcentral/fileexchange/42885-nearestspd.

    [2] N.J. Higham, "Computing a nearest symmetric positive semidefinite
    matrix" (1988): https://doi.org/10.1016/0024-3795(88)90223-6.
    """

    B = (A + A.T) / 2
    _, s, V = la.svd(B)

    H = np.dot(V.T, np.dot(np.diag(s), V))

    A2 = (B + H) / 2

    A3 = (A2 + A2.T) / 2

    if isPD(A3):
        return A3

    spacing = np.spacing(la.norm(A))
    # The above is different from [1]. It appears that MATLAB's `chol` Cholesky
    # decomposition will accept matrixes with exactly 0-eigenvalue, whereas
    # Numpy's will not. So where [1] uses `eps(mineig)` (where `eps` is Matlab
    # for `np.spacing`), we use the above definition. CAVEAT: our `spacing`
    # will be much larger than [1]'s `eps(mineig)`, since `mineig` is usually on
    # the order of 1e-16, and `eps(1e-16)` is on the order of 1e-34, whereas
    # `spacing` will, for Gaussian random matrixes of small dimension, be on
    # the order of 1e-16. In practice, both ways converge, as the unit test
    # below suggests.
    I = np.eye(A.shape[0])
    k = 1
    while not isPD(A3):
        mineig = np.min(np.real(la.eigvals(A3)))
        A3 += I * (-mineig * k**2 + spacing)
        k += 1

    return A3

def isPD(B):
    """Returns True when input is positive-definite, via Cholesky"""
    try:
        _ = la.cholesky(B)
        return True
    except la.LinAlgError:
        return False

if __name__ == '__main__':
    for i in range(10):
        for j in range(2, 100):
            A = np.random.randn(j, j)
            B = nearestPD(A)
            assert(isPD(B))
    print('unit test passed!')

def TV(z, S):
    """
        This function calculates the T^2 value of the given pathway
        and the corresponding interaction matrix.
    """
    kik = np.where(z!=0)[0]
    if len(kik)==0:
        return 0
    else:
        return np.dot(np.dot(np.transpose(z[kik]), np.linalg.pinv(np.diag(S[kik,kik]))), z[kik])
    
def TS(pathway, ppi, stu, purb, dgv=0.4):
    """ T-square.
        For the given pathway, this function creates the corresponding interaction matrix.
        Returns the associated T^2, p-value, and other information.
    """
    # - pathway is a pandas dataframe containing the id of a pathway, its included proteins, and their abundance ratios.
    # - z contains only the abundance ratios of the proteins, to be used later to calculate the T^2 value.
    # - m contains the indexes of the pathway's proteins that can be translated form Uniprot to STRING.
    # - S is the interaction matrix to be built from STRING interaction scores.
    pathway = pathway.sort_values(by='prot_id').reset_index(drop = True)
    z = np.vectorize(float)(pathway['exp'])
    m = np.where(np.isin(stu, pathway))[0]
    S = dgv * np.identity(len(z))
    nrow, ncol = S.shape
    if nrow != 1:
        # Each possible pair of proteins in the pathway will be looked at.
        for i in range(1, nrow):
            for j in range(i):
                # x1 is the Uniprot accession to one protein i in the pathway.
                x1 = pathway.iat[i, 1]
                # x2 is the Uniprot accession to another protein j in the pathway.
                x2 = pathway.iat[j, 1]
                # s1 is the STRING id to protein i, translated using m.
                s1 = stu.iloc[m]['String_id'].to_numpy()[np.where(np.isin(stu.iloc[m]['Uniprot_id'], x1))[0]]
                # s2 is the STRING id to protein j.
                s2 = stu.iloc[m]['String_id'].to_numpy()[np.where(np.isin(stu.iloc[m]['Uniprot_id'], x2))[0]]   
                if len(s1)*len(s2) !=0:
                    # If there is one, p will contain the experimental value of interaction between the two proteins.
                    # If there are more, the mean will be used.
                    # Get all protein 1 partners in STRING
                    p = ppi.iloc[np.where(np.isin(ppi['protein1'],s1))[0]]
                    # Check for protein among partners
                    p = p.iloc[np.where(np.isin(p['protein2'], s2))]['experimental']
                    # p holds interaction score(s) from experimental evidences only
                    if len(p)>0:   
                        # Modify S to include that value at the corresponding index.   
                        if z[pathway['prot_id']==x1]*z[pathway['prot_id']==x2] <0:
                            S[i,j] = -np.mean(p)
                            S[j,i] = -np.mean(p)
                        else:
                            S[i,j] = np.mean(p)
                            S[j,i] = np.mean(p)
    # Transform STRING scores, S matrix to be positive-definite, and therefore useable for the T^2 method.
    S = nearestPD(S)
    r = np.linalg.matrix_rank(S, tol=1e-10)
    # T2 score matrix effectively used
    T2 = TV(z, S)
    I = dgv * np.identity(len(z))
    T2I = TV(z, I)
    return np.array([pathway.iat[0,0],
                    ','.join(pathway['prot_id']),
                    len(pathway),
                    r,
                    T2,
                    chi2.sf(T2, r),
                    T2I,
                    chi2.sf(T2I, r)],
                    dtype=object)

def PS(pi, cov=0):
    """
        This function returns a list of lists of pathways.
        Each sublist contains an 'delegate' pathway and all the ones included in it.
    """
    pi[1] = np.vectorize(int)(pi[1])
    a = np.array([s.split(",") for s in pi[2]])
    tag = np.array(a[pi[1].idxmax()])
    ll = np.array([len(np.setdiff1d(x, tag)) for x in a])
    g = pd.DataFrame(pi[ll<=cov].values)
    if np.ndim(g)!=1:
        g.sort_values(1, ascending=False, inplace=True)
        g.reset_index(drop = True, inplace=True)
    newpi = pd.DataFrame(pi[ll>cov].values)
    nrow, ncol = newpi.shape
    if ncol==1:
        g = [g]+[newpi]
    elif nrow==0:
        return [g]
    else:
        g = [g]+PS(newpi)
    return g

def extract_accession_abundance(df, abnd_label):
    acc_col = df.columns.get_loc("Accession")
    abnd_col = df.columns.get_loc(abnd_label)
    return df.values[:,[acc_col,abnd_col]]

def predata(data, outth=10):
    """
        This function gives input data the appropriate format for importdata to use it. 
        Namely, groups duplicates by their median value.
    """
    nrow, ncol = data.shape
    print("    #(input site/probe): {}".format(nrow))
    
    ### Remove missing
    tokp = np.logical_not([(k[1]=="NAN") or (k[1]=="NaN") or (k[1]=="NA") or (k[1]=="na") or (k[1]=="-") or (k[0]=="") or (k[1]=="") or pd.isnull(k[0]) or pd.isnull(k[1]) for k in data])
    data = data[tokp]
        
    ### Multiple ids one value [GL: Not understood, unitprot -X fragment notation patching?]
    data1 = np.array([s for s in data if len(s[0])<11])
    data2 = np.array([s for s in data if len(s[0])>10])
    if data2.shape[0] != 0:
        cop = np.empty((0,2), dtype=object)
        for protein in data2:
            for realid in protein[0].split("|"):
                cop = np.append(cop, [[realid, protein[1]]], axis=0)
        data2 = np.copy(cop)
    else:
        data2 = data2.reshape(0,2)
    data = np.concatenate((data1,data2))
    data = pd.DataFrame(data, columns=["id","exp"])
    data['id']= np.vectorize(lambda t: t[0:6])(data['id'])
    
    ### one id multiple values
    #data['exp'] = np.vectorize(float)(data['exp'])
    data['exp'] = np.vectorize(float)(
        data.applymap(lambda x: x if type(x) == 'float' else str( x.replace(',','.') ))['exp'] 
    )
    med = data.groupby(['id']).median()
    data = pd.DataFrame({'id':med.index.values, 'M':med.values.flatten()})
    
    ### normalization
    if min(data['M'])>=0:
        nul = data[data['M']==0].index
        data.iloc[nul].M = min(data[data['M']!=0]['M'])
        data['M']=np.log2(data['M'])
    
    ### outliers replacement
    for i in range(len(data)):       
        if not (-outth <= data.iloc[i]['M'] <= outth):
            data.iloc[i].M = np.sign(data.iloc[i]['M'])*max(abs(data['M']))   
            
    ### standardization
    data['M'] = data['M'] / np.std(data['M'])
    
    return data


# __Main functions__

def importdata(inputDf1, inputDf2=None, outth=100, abnd_label = "Corrected_Abundance_Ratio"):
    """ Data import and pre-process.
    
        This function removes NA, replaces extreme values, and
            standardizes the data. It also maps the data with
            Uniprot identifiers to ensp identifiers.
         
        file1: Expression data with Uniprot identifiers.
        file2: Expression data with Uniprot identifiers, optional for time-course data.
        outth: Outlier threshold, default is 10.

        We extract the uniprot ID and the "Corrected_Abundance_Ratio" foreach DF row
    """
    print("=================================================")
    print(" Dataset summary:")
    print("-------------------------------------------------")
    
    if inputDf2 is None:
        #x = file1.values[:,0:2]
        x = extract_accession_abundance(inputDf1, abnd_label)
        data = predata(x, outth)
    else:
        #x = file1.values[:,0:2]
        #y = file2.values[:,0:2]
        x = extract_accession_abundance(inputDf1, abnd_label)
        y = extract_accession_abundance(inputDf2, abnd_label)
        data1, data2 = predata(x, outth), predata(y, outth)
        #time series, divide time 1 and time 2
        data12 = np.empty((0,2), dtype=object) 
        m2 = min(data2['M'], key=abs)
        for i in range(len(data1)):
            proti = data1.iloc[i]['id']
            if proti not in data2['id'].values:
                data12 = np.concatenate((data12, [[proti, data1.iloc[i]['M']-m2]]))
                
        data21 = np.empty((0,2), dtype=object)
        m1 = min(data1['M'], key=abs)
        for i in range(len(data2)):
            proti = data2.iloc[i]['id']
            if proti not in data1['id'].values:
                data21 = np.concatenate((data21, [[proti, data2.iloc[i]['M']-m1]]))
        
        data12 = pd.DataFrame(data12, columns=['id','M'])
        data21 = pd.DataFrame(data21, columns=['id','M'])
        
        datai, d1i, d2i = np.intersect1d(data1['id'],data2['id'], return_indices=True)
        data = np.empty((0,2), dtype=object)
        for i, inter in enumerate(datai):
            data = np.concatenate((data, [[inter, data1.iloc[d1i[i]]['M'] - data2.iloc[d2i[i]]['M']]]), axis=0)
        data = pd.DataFrame(data, columns=['id','M'])
        
        data = pd.concat((data,data12,data21), ignore_index=True)
    
    data.columns = ['id', 'exp']
    data['exp'] = np.vectorize(float)(data['exp'])
    print("    #(input proteins):   {}".format(data.shape[0]))
    print("=================================================")
    return data



def computeT2(data, vex, pid, ppi, stu, purb=1.5, intg=True, alpha=0.05, ncore=1, sizelim = 100):
    """ Computes T**2 and its p-value for each pathway.
        This function computes the T**2 score and its significance level.
        
        data: Proccessed data using importdata function.
        vex: A dataframe displaying pathways and the proteins they contain.
        pid: A dataframe displaying pathways and their full-length name.
        ppi: Protein-protein interaction dataframe: obtained from STRING, for example.
        stu: STRING to Uniprot translation dataframe.
        purb: Perturbance threshold, default is 1.5 (after normalization).
        intg: Apply pathway integration or not. Default is True.
        alpha: Significance level. Default is 0.05.
        ncore: Number of parallel computing cores. Default is 7.
        sizelim: Maximum size of the pathways to return.
    """
    
    # We get rid of values under the perturbance threshold.
    data.loc[data['exp'].between(-purb,purb,inclusive=False), 'exp']=0
    
    ### data mapping
    # n is the indexes of the proteins from our data that can be found in our dataset of pathways. We print its length.
    n = np.where(np.isin(vex['prot_id'], data['id']))[0]
    print("    #(mapped entries):    {}".format(len(np.intersect1d(vex['prot_id'], data['id']))))
    # Conversely, vexData is the subset of pathways that contain proteins from our data. We print its length and ignore the rest.
    vexData = vex.iloc[n].copy()
    vexData['exp'] = np.array([[dv[1] for dv in data.values if dv[0]==vv[1]] for vv in vexData.values])
    pathwayList = np.unique(vexData['pathway_id'])
    print("    #(mapped pathways):   {}".format(len(pathwayList))) 
    
    ### pathway integration

    # This function and the following loop are only to reshape our set of pathways to something easier to use.
    def desc(cp):
        pathway = vexData[vexData['pathway_id']==cp]
        pathwaygenes = ",".join(np.vectorize(str)(pathway['prot_id']))
        return (pathway.iat[0,0], pathway.shape[0], pathwaygenes)
    pi = pd.DataFrame([desc(cp) for cp in pathwayList]).sort_values(by=1,ascending=False).reset_index(drop = True)

    # The PS function will apply the 'pathway integration' process.
    # The inpt variable will contain a reshaped result, depending on whether or not the user wants to aknowledge 'pathway integration'.
    ps = PS(pi)
    print("    #(summary pathways):  {}".format(len(ps)))
    inpt = np.array([p[0][0] for p in ps]) if intg else np.concatenate([list(p[0]) for p in ps])
    
   


#    print(inpt[:10])
#    jobs =[]
#    for pthwyID in inpt:
#        jobs.append( (pthwyID, len(vexData[vexData['pathway_id']==pthwyID]) ) )
#    print( sorted(jobs, key=lambda x:x[1]) )
    

    ### compute T2

    # This function and the following loop will allow to create a T^2 value for each pathway in inpt. See TS.
    # The result is stored in the pandas dataframe r.
    def desc2(cp): # cpis a pathway
        pathway = vexData[vexData['pathway_id']==cp]
        size = len(pathway)
        if size==0:
            return np.array([cp, '', size, 0, 0, 1, 0, 1], dtype=object) # dummy result
        if size<=sizelim:
            return TS(pathway, ppi, stu, purb)# returns the vector of T-scores 
        else:
            return np.array([pathway.iat[0,0], ','.join(pathway['prot_id']), size, 0,0,1,0,1], dtype=object) # dummy result

    def parallel_desc2(inpt, vexData, ncore):

        print(f"Total input len {len(inpt)}")
            ### We Create pathway pool by distributing same number of pathway based on their sizes.
        def createPool(inputList):
            tasks = []
            for pthwyID in inputList:
                tasks.append( (pthwyID, len(vexData[vexData['pathway_id']==pthwyID]) ) )
            jobQueue = [ [] for _ in range(ncore) ]
            
            tasks = sorted(tasks, key=lambda x:x[1], reverse=True)
            
            i=0
            while (i + len(jobQueue)) <= len(tasks):
                for j,q in enumerate(jobQueue):               
                    q.append(tasks[ i + j ][0]) # We only store patchway name, getting rid of its weight
                i = i + len(jobQueue)

            for _, jobleft in enumerate(tasks[i:]):
                jobQueue[_].append(jobleft[0]) # We only store patchway name, getting rid of its weight
            return jobQueue

        jobQueue = createPool(inpt)

        def worker(pthwList, procnum, resultStore):
            print(f"Starting:{multiprocessing.current_process().name} len={len(pthwList)}")
            _r = [desc2(cp) for cp in pthwList]
            resultStore[procnum] = _r

        for j in jobQueue:
            print(f"#{len(j)}")
        manager = multiprocessing.Manager()
        resultStore = manager.dict()
        procList = [ multiprocessing.Process( target=worker, args=(j, i, resultStore) ) for i,j in enumerate(jobQueue) ]
        for proc in procList:
            proc.start()
        for proc in procList:
            proc.join()
        print(resultStore.keys())
        r = []
        for k in sorted(list(resultStore.keys()), key=lambda x:float(x)):
            print(k)
            r += resultStore[k]
    
        return pd.DataFrame(r)


    
# original single thread call
    if ncore > 1
        r = pd.DataFrame([desc2(cp) for cp in inpt]) # stores vectors of score foreach pathway
    else:
        r = parallel_desc2(inpt, vexData, ncore)


    # This function and the following loop reshape r to provide a more satisfactory display.
    # The result is stored in the pandas dataframe rrr.
    def desc3(l):
        ttl = pid[pid['pathway_id']==l[0]]
        ttl = ttl.iat[0,1]
        return np.insert(l[:6],0,ttl) # drop trail 6> column and insert pathway name
    # reshape pathway rows
    rrr = pd.DataFrame([desc3(l) for l in r.values])
    # reshape finale DF by adding headers
    rrr.columns = ["Pathway title","Pathway ID","Uniprot IDs","#Mapped","df","T-square","p-value"]
    # Finally, the pathways with too high a p-value are excluded, and we print how many pathways are left.
    rrr = rrr[rrr["p-value"]<=alpha].reset_index(drop=True)
    print("    #(enriched pathways): {}".format(rrr.shape[0]))
    print("=================================================")
    return(rrr)

    

