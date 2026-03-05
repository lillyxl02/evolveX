from run_FBA import run_FBA

def apply_bounds(m, env):
    for rxn_id, typ in env.items():
        rxn = m.reactions.get_by_id(rxn_id)
        if typ == "inh":
            rxn.lower_bound = 0.0
            rxn.upper_bound = 0.0
        elif typ == "comp":
            rxn.lower_bound = -1000
            rxn.upper_bound = -1 

def compute_opt(m, env):
    comp_rxns = [m.reactions.get_by_id(rid) for rid,t in env.items() if t=="comp"]
    m.objective = {rxn: 1.0 for rxn in comp_rxns}
    sol = run_FBA(m, csense="max")
    return sol
