#FBA
import cobra

def run_FBA(model, csense=None):
    
    if csense is None or csense.lower() == "max":
        model.solver.objective.direction = 'max'
    else:
        model.solver.objective.direction = 'min'
    
    sol = model.optimize()
    return sol