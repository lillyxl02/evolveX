
import numpy as np
from trait_flux_basis import trait_flux_basis
from prune_flux_basis_solutions import prune_flux_basis_solutions
from get_alternative_flux_bases import get_alternative_flux_bases
from getPotentialTargets import getPotentialTargets 
from get_flux_basis_directions import get_flux_basis_directions

def determine_trait_flux_basis(model, osense, wtminmax, alpha, delta, eps, noSols):

    limit = 0
    fullSol = []
    solutionA = [] #liste av flere solutions som kommer fra trait_flux_basis. 
    binSol = []

    targets, dirs, signs, dirSolp, fullSolOpt, solutionp = {}, [], [], [], [], []
    
    # Round wtminmax for numerical stability
    wtminmax = np.round(wtminmax * 1e9) / 1e9

    # Get reactions evaluated for potential flux change
    binRxns = getPotentialTargets(model, wtminmax, [], tol=1e-7)
    # Enumerate alternative solutions
    ultimate_limit = 3 * noSols
    while limit < noSols:
        solutionA, binSol = trait_flux_basis(model, osense, wtminmax, binRxns, alpha, binSol, solutionA, delta, eps=1e-5)
        #Om den finner en løsning vil den legge til solutionA[cfullSol]
        cfullSol = len(fullSol)
        if len(solutionA) > cfullSol and solutionA[cfullSol]["full"] is not None: 
            fullSol.append(solutionA[cfullSol]["full"]) #[alle fluxer v] + [alle binære y]
            limit += 1
        else:
            print("Ingen fullsol")
        ultimate_limit -= 1
        if ultimate_limit == 0:
            limit = noSols


    if binSol:
        binSolp, fullSolp, dirThresp, solutionp = prune_flux_basis_solutions(wtminmax, np.column_stack(binSol), np.column_stack(fullSol), binRxns, delta, eps, solutionA)
        binSolp = np.round(binSolp)
        indeces = get_alternative_flux_bases(binSolp, False) #True (1) betyr at alle løsninger tillates. 0 hvis man kun vil ha de med minste antall targets. 
        #indeces er en liste for hvilke kolonner i binsolp som beholdes. 
        if indeces:
            binSol_mat = np.array(binSolp)
            fullSol_mat = np.array(fullSolp)
            dirThresp_mat = np.array(dirThresp)
            targets, dirs, signs, dirSolp = get_flux_basis_directions(wtminmax, binSol_mat[:, indeces], fullSol_mat[:, indeces], binRxns, dirThresp_mat[:, indeces], delta, eps)
            fullSolOpt = [fullSolp[i] for i in indeces]
            solutionp = [solutionp[i] for i in indeces]
        else:
            print("ingen indeces")
    else:
        print("ingen binsol")
    

    return targets, dirs, signs, dirSolp, binRxns, fullSolOpt, solutionp


