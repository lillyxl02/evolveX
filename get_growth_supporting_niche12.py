import itertools
import numpy as np
from run_FBA import run_FBA

def get_growth_supporting_niche12(model, condMap):

    carbons = condMap["carbon"]
    nitrogens = condMap["nitrogen"]
    growth_value=10.0 #vekstnivå som skal oppnås
    uptake_lb=-10000.0 #Må forbruke, derfor ligger ub og lb på negativ side. 
    uptake_ub=0


    #Kombinasjoner av 1 karbon og 2 nitrogenkilder
    combos = [(c, n1, n2) for c in carbons for (n1, n2) in itertools.combinations(nitrogens, 2)]

    growth = [0] * len(combos) #printer 1 hvis den kombinasjonen gir vekst. 
    uptakes = [None] * len(combos)

    for k, (c_rxn, n1_rxn, n2_rxn) in enumerate(combos, start=1):
        print(k, c_rxn, n1_rxn, n2_rxn)
        m = model.copy()
        #låser ub/lb til vekst
        biomass_rxn = m.reactions.get_by_id("BIOMASS_Ec_iML1515_core_75p37M")
        biomass_rxn.lower_bound = float(growth_value)
        biomass_rxn.upper_bound = float(growth_value)

        lb = float(uptake_lb)
        ub = float(uptake_ub)
        for r_id in (c_rxn, n1_rxn, n2_rxn):
            r = m.reactions.get_by_id(r_id)
            r.lower_bound = lb
            r.upper_bound = ub

        #objektivfunksjonen
        for r in m.reactions:
            r.objective_coefficient = 0.0
        m.reactions.get_by_id(c_rxn).objective_coefficient = 1.0
        m.reactions.get_by_id(n1_rxn).objective_coefficient = 1.0
        m.reactions.get_by_id(n2_rxn).objective_coefficient = 1.0

        #sjekker om vekst er optimal med FBA, hvis ikke settes growth til 0 og uptakes til none. 
        sol = run_FBA(m, "max")
        if sol is not None and sol.status == "optimal":
            u = np.array([float(sol.fluxes[c_rxn]),
                          float(sol.fluxes[n1_rxn]),
                          float(sol.fluxes[n2_rxn])], dtype=float)
            growth[k-1] = 1
            uptakes[k-1] = u

    return growth, uptakes, combos

