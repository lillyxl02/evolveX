import numpy as np

def prune_flux_basis_solutions(wtminmax, binSol, fullSol, binRxns, delta, eps, solution):
    #fjerner ikke-optimale løsninger
    
    binSol = np.round(binSol) #sikre at binsol faktisk er 0 eller 1
    n_rxns, n_sols = binSol.shape[0], binSol.shape[1] #n_rxns = antall binære reaksjoner, n_sols = antall løsninger
    dirBinSol = np.zeros_like(binSol) #matrise fylt med 0
    prune = np.ones(n_sols, dtype=bool)

    tol = 1e-8

    for j in range(n_sols):
        sol = solution[j]
        if sol.get("origStat", 0) != 1:
            prune[j] = False
            continue
        
        for i_idx in range(n_rxns):
            if binSol[i_idx, j] != 1:
                continue  # ignorer binær=0
            rxn_idx = binRxns[i_idx]
            wt_min, wt_max = wtminmax[rxn_idx, :]
            min_thres = wt_min - delta * abs(wt_min) - eps
            max_thres = wt_max + delta * abs(wt_max) + eps
            flux_val = fullSol[rxn_idx, j]

            #DOWN
            if flux_val <= min_thres + tol:
                dirBinSol[i_idx, j] = -1

            #UP
            elif flux_val >= max_thres - tol:
                dirBinSol[i_idx, j] = 1

            else:
                prune[j] = False #kaster løsningen hvis den ikke er innenfor WT-intervall
                dirBinSol[i_idx, j] = 0
                break  #markér som uendret, men pruner ikke løsningen

    binSolp = binSol[:, prune] #beholder kun kolonner der prune[i]==true
    fullSolp = fullSol[:, prune]
    dirBinSolp = dirBinSol[:, prune]
    solutionp = [solution[j] for j in range(n_sols) if prune[j]]


    return binSolp, fullSolp, dirBinSolp, solutionp