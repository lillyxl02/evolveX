import numpy as np

def get_alternative_flux_bases(binSolp, notOnlyOptimal):
    nSols = binSolp.shape[1]

    #Count number of non-zero targets per solution
    tol = 1e-6
    nTargets = np.sum(binSolp > tol, axis=0)
    indeces = np.arange(nSols)

    # Filter optimal solutions (minimal number of targets) unless allowed otherwise
    if not notOnlyOptimal: #når den er false
        minvalue = np.min(nTargets) #finner minste antall aktive targets
        mask = (nTargets == minvalue) #boolean array for de som har minvalue
        optSolMatrix = binSolp[:, mask] #kun de med min antall targets
        indeces = indeces[mask]
    else:
        optSolMatrix = binSolp #hvis True, tar alle

    #find unique solutions (transponeres for å få løsninger som rader)
    altoptsols, unique_indices = np.unique(optSolMatrix.T, axis=0, return_index=True)

    #Keep columns in original orientation (reactions × solutions)
    altoptsols = altoptsols.T
    indeces = indeces[unique_indices].tolist()

    return indeces
