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

        #Sjekker om wt intervallet for reaksjonen er identisk med modellens bounds. Hvis ja, lite poeng å gi binærevariabler
        #Hvis nei som man sjekker nedenfor, interessant som targets og derfor legges det til i binRxns. 
        same_as_bounds = (
            abs(wt_lb - lb[r]) < tol and
            abs(wt_ub - ub[r]) < tol
        )

        if not same_as_bounds:
            binRxns.append(r)

    return np.array(binRxns, dtype=int)

