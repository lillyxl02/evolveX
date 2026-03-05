import numpy as np

def get_flux_basis_directions(wtminmax, binSol, fullSol, binRxns, dirThresp, delta, eps):

    binSol = np.round(binSol) #skal være 0 eller 1, runder opp små feil. binsol består av binæreverdier
    fullSol = np.round(1e9 * fullSol) / 1e9 #fullsol er tilsvarende fluxer. Runder til 9 desimaler. 

    n_bin, n_solutions = binSol.shape #hver kolonne i binsol er en løsning fra MILP. Antall kolonner = antall løsninger

    #lager matrise med samme form som dirThresp
    dirSolp = np.zeros_like(dirThresp) 
    signChange = np.ones_like(dirThresp)

    #Tom liste for hver løsning i n_solutions som er brukervennlig å lese av. 
    dirs = [[] for _ in range(n_solutions)]
    targets = [[] for _ in range(n_solutions)]
    signs = [[] for _ in range(n_solutions)]

    for j in range(n_solutions):          
        k = 0

        for i in range(n_bin):      #for each binary reaction
            if binSol[i, j] != 1: #ignorer hvis reaksjonen ikke er en target
                continue
            rxn=binRxns[i]
            #min og max terskel for fluxendring
            minThres = wtminmax[rxn, 0] - delta * abs(wtminmax[rxn, 0]) - eps 
            maxThres = wtminmax[rxn, 1] + delta * abs(wtminmax[rxn, 1]) + eps

            if maxThres < minThres:
                print("threshold error")

            #detect sign change (har den motsatt fortegn enn WT?) 1000 hvis ja, 1 hvis nei. 
            if dirThresp[i, j] == -1 and np.sign(minThres) * np.sign(fullSol[rxn, j]) < 0: #hvis WT var negativ og fluxen er nå positiv. 
                signChange[i, j] = 1000
            elif dirThresp[i, j] == 1 and np.sign(maxThres) * np.sign(fullSol[rxn, j]) < 0: #hvis WT var positiv og fluxen er nå negativ. 
                signChange[i, j] = 1000

            #-1000/1000 betyr at det har skjedd en sign change.
            #DOWN 
            if fullSol[rxn, j] <= minThres:
                dirSolp[i, j] = -1 * signChange[i, j]
                dirs[j].append("DOWN")
                targets[j].append(rxn)
                signs[j].append(np.sign(fullSol[rxn, j]))
                k += 1

            #UP
            elif fullSol[rxn, j] >= maxThres:
                dirSolp[i, j] = 1 * signChange[i, j]
                dirs[j].append("UP")
                targets[j].append(rxn)
                signs[j].append(np.sign(fullSol[rxn, j]))
                k += 1

            #neither up nor down
            else:
                dirSolp[i, j] = 5000

    return targets, dirs, signs, dirSolp
