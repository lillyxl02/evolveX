import numpy as np
from EVOLVEX_SCORE import evolveX
from optimality import compute_opt, apply_bounds

def enumerate_evolvex(model, envs, targets):
    strengths = []
    comp = []

    for i, env in enumerate(envs, start=1):
        #hvert miljø har sin egen verdi på hva optimality er

        m = model.copy()
        apply_bounds(m, env)
        sol_upt = compute_opt(m, env)

        optimality = float(sol_upt.objective_value)

        strength = evolveX(m, optimality, env, targets)
        strengths.append(strength)
        comp_rxns = [rid for rid, typ in env.items() if typ == "comp"]
        comp.append(comp_rxns)

        print(f"{i}/{len(envs)} strength={strength} comp={comp_rxns}")

    return strengths, comp
