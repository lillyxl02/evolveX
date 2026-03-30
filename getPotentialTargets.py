import numpy as np

def getPotentialTargets(model, wtminmax, redVector, tol=1e-7):

    rxns = len(model.reactions)

    #Hvis ingen redVector gitt, start med alle tillatt
    if redVector is None or len(redVector) == 0:
        redVector = np.ones(rxns, dtype=int)
    else:
        redVector = np.asarray(redVector).astype(int)

    #Ekskluder alle exchange reactions automatisk
    for rxn in model.reactions:  
        if rxn.id.startswith("EX_"):
            redVector[model.reactions.index(rxn)] = 0

    #Hent bounds fra modellen
    lb = np.array([r.lower_bound for r in model.reactions], dtype=float)
    ub = np.array([r.upper_bound for r in model.reactions], dtype=float)

    #Finn hvilke reaksjoner som faktisk kan bli targets
    binRxns = []
    for r in range(rxns):
        if redVector[r] != 1:
            continue

        wt_lb = float(wtminmax[r, 0])
        wt_ub = float(wtminmax[r, 1])

        if lb[r] < tol:
            #reversible
            if not (wt_lb == -1000 and wt_ub == 1000):
                binRxns.append(r)

        else:
            #irreversible
            if not (wt_lb == 0 and wt_ub == 1000):
                binRxns.append(r)


    return np.array(binRxns, dtype=int)

