import numpy as np
import cobra
import cplex

def trait_flux_basis(model, osense, wtminmax, binRxns, alpha, binSol, solutionA, delta, eps):
    """
OUTPUT: solution   solution structure from the solver
        binSol     binary solutions defining the identified targets

1.Først løser den LP som gir obj.value og vopt (nivå vi krever at cellen fortsatt skal oppnå)
2.Finner reaksjoner som har ulik lb og ub i forhold til wtminmax. Disse får binære variabler
3.MILP løser det automatisk, hvilke reaksjoner som havner utenfor WT eller ikke. 
totalt sett skal den finne fluxløsninger som oppfyller massebalanse, produksjonskrav og endre færrest mulig reaksjoner slik at de havner utenfor WT (ROOM)
    """
    S = np.array(cobra.util.array.create_stoichiometric_matrix(model)) #støkiometrimatrise
    mets, rxns = S.shape

    lb = np.array([rxn.lower_bound for rxn in model.reactions], dtype=float)
    ub = np.array([rxn.upper_bound for rxn in model.reactions], dtype=float)
    c = np.array([rxn.objective_coefficient for rxn in model.reactions], dtype=float) #objektivkoeffisient 

    #LP: maksimer/minimer for å få vopt (hvor mye produkt som kan kreves). Finner "optimal" produksjon (FBA)
    lp = cplex.Cplex()
    lp.set_log_stream(None)
    lp.set_results_stream(None)
    lp.variables.add(obj=c.tolist(), lb=lb.tolist(), ub=ub.tolist(),
                     types=[lp.variables.type.continuous]*rxns,
                     names=[f"v_{i}" for i in range(rxns)])
    #S*v = 0 constraints
    for j in range(mets):
        lp.linear_constraints.add(
            lin_expr=[[list(range(rxns)), S[j, :].tolist()]],
            senses=["E"],
            rhs=[0.0]
        )
    #finner vopt (produksjonsnivå)
    if osense.lower() == "max":
        lp.objective.set_sense(lp.objective.sense.maximize) #maksimer hvis osense er max
    else:
        lp.objective.set_sense(lp.objective.sense.minimize)
    lp.solve()
    if not lp.solution.is_primal_feasible():
        return None, binSol
    
    objValue = lp.solution.get_objective_value()
    thr = 1e9
    vopt = np.floor(thr * objValue * alpha) / thr

    #MILP
    milp = cplex.Cplex()
    milp.set_log_stream(None)
    milp.set_results_stream(None)
    milp.objective.set_sense(milp.objective.sense.minimize)

    # MILP har alle flukser med alle binære verdiene 
    wtrxns = len(binRxns)
    milp.variables.add(
        obj=[0.0]*rxns + [1.0]*wtrxns, #de første rxns (flux) coefficient=0, de neste wtrxns (binære) =1, summen av de binære er obj. 
        lb=lb.tolist() + [0.0]*wtrxns, #grensene lb og ub. 
        ub=ub.tolist() + [1.0]*wtrxns,
        types=[milp.variables.type.continuous]*rxns +
            [milp.variables.type.binary]*wtrxns #de første rxns er continuous, de neste er binary. 
    )

    # Massebalanse
    for j in range(mets):
        milp.linear_constraints.add(lin_expr=[[[i for i in range(rxns)], S[j, :].tolist()]], senses=["E"],rhs=[0.0])

    #G betyr ≥: cᵀv >= vopt, L betyr ≤: cᵀv <= vopt. Låser produksjonen til å være like god eller lik vopt. Men avhenger om det er max eller min. 
    #lenger ned brukes produksjonskravet til å minimere antall reaksjoner som må avvike fra WT for å nå kravet. Reaksjonene er enten y=0 eller y=1 (tillate avvik)
    prod_sense = "G" if osense.lower() == "max" else "L" 

    milp.linear_constraints.add(
    lin_expr=[[[i for i in range(rxns)], c.tolist()]],
    senses=[prod_sense],
    rhs=[vopt]
    )

    # WT-terskler (trickle flow)
    wiu = wtminmax[:, 1] + delta*np.abs(wtminmax[:, 1]) + eps
    wil = wtminmax[:, 0] - delta*np.abs(wtminmax[:, 0]) - eps

    # Flux-endringsconstraints BigM constraint
    #Setter inn y_i=0, får: wilr≤vr≤wiur
    #Setter y_i=1, lbr≤vr≤ubr
    # hvis y_i=0 tving vr inn i WT intervallet , 
    # Hvis y_i=1, da kan vr være hvor som helst innenfor lb og ub
    #big M defineres her som wiu[r] - ub[r] og wil[r] - lb[r]
    for i, r in enumerate(binRxns):
        r = int(r) #binRxns er float, må gjøre om til int. 
        yi = int(rxns + i)

        # v_r + y_i*(wiu - ub) ≤ wiu
        milp.linear_constraints.add(
            lin_expr=[[[r, yi],[1.0, wiu[r] - ub[r]]]], senses=["L"], rhs=[wiu[r]])

        # v_r + y_i*(wil - lb) ≥ wil
        milp.linear_constraints.add(
            lin_expr=[[[r, yi],[1.0, wil[r] - lb[r]]]],senses=["G"], rhs=[wil[r]])

    milp.parameters.mip.tolerances.integrality.set(1e-9)
    milp.solve()
    if not milp.solution.is_primal_feasible():
        return None, binSol

    sol = milp.solution.get_values()
    flux_sol = np.array(sol[:rxns])
    binvals = np.array(sol[rxns:rxns+wtrxns])
    targets = [binRxns[i] for i in range(wtrxns) if binvals[i] > 0.5]


    solution = {
        "obj": milp.solution.get_objective_value(),
        "flux": flux_sol,
        "targets": targets,
        "full": sol,
        "origStat": 1
    }

    binSol.append(binvals) #binære verdier
    solutionA.append(solution)

    return solutionA, binSol

