import cobra
from cobra import Reaction
from math import floor, ceil
import cplex

#this script determines which fluxes can be exposed to positive selection using combinations of a limited set of most common nutrients
def setBounds(cbmodel): #gjør grensene bredere
	for rxn in cbmodel.reactions:
		rxn.lower_bound = rxn.lower_bound * 10
		rxn.upper_bound = rxn.upper_bound * 10


def createSecretion(model, rxn): #lager en kunstig reaksjon for sekresjon av metabolitter. 
	irxn = model.reactions.get_by_id(rxn)
	reaction = Reaction("secretion_" + rxn)
	reaction.name="secretion of rxn"
	reaction.lower_bound = 0
	reaction.upper_bound = 10000
	reaction.add_metabolites(irxn.metabolites)
	model.add_reactions([reaction])


def evolveX(model, optimality, env, targets):
	    #"""**evolvex** Predicts the relative selection on the given target fluxes in evolution conditions. Compatible also with gpr-transformed models.
    #Authors: Paula Jouhten, Sandra Castillo

    #Arguments:
        #modelin (CBModel): a constraint-based model
        #optimality (double): optimal value of the objective function of the model
        #envinronment (dict): reaction ids as keys and the type ("comp" or "inh") as values
        #targets (dict): rxn ids as keys, desired change ("UP" or "DOWN") as values
    
    #Returns:
        #double: 
    #"""

	#optimization parameters
	#Lager et cplex modell
	relax = 1
	tol = 1e-9
	problem = cplex.Cplex()
	problem.set_log_stream(None)
	problem.set_error_stream(None)
	problem.set_warning_stream(None)
	problem.set_results_stream(None)
	problem.parameters.simplex.tolerances.optimality = 1e-9
	problem.parameters.randomseed.set(1)


	for rxn in model.reactions:
		#the fluxes that are targets and can have only positive values get non-zero objective coefficients
		if rxn.id in targets:
			if rxn.lower_bound >= 0: #hvis reaksjonen bare har positive flukser
				if targets[rxn.id] == 'UP':
					#if the correpsonding targetDirs is UP the flux variable gets coefficient 1 (minimize) in the objective
					problem.variables.add(obj=[1], lb=[rxn.lower_bound], ub=[rxn.upper_bound], names=[rxn.id])
				else:
					#if the correpsonding targetDirs is DOWN the flux variable gets coefficient -1 (maximize) in the objective
					problem.variables.add(obj=[-1], lb=[rxn.lower_bound], ub=[rxn.upper_bound], names=[rxn.id])
			
			#the fluxes that are not targets or have lower bound < 0 get 0 in the objective.
			#hvis reaksjonen kan gå i begge retninger, setter den ingen påvirkning på objektivfunksjonen obj[0]
			elif rxn.lower_bound < 0:
				problem.variables.add(obj=[0], lb=[rxn.lower_bound], ub=[rxn.upper_bound], names=[rxn.id])
				
				#hvis reaksjonen skal økes, oprettes en ny variable up_rxn.id, og legger til to begresninger som samsvarer med reaksjonsflux const1,2
				#elif: samme hvis en reaksjon skal nedreguleres. legger til 4 begrensninger
				if targets[rxn.id] == 'UP':
			    	#one new continuous variable is created
					ubs = max(abs(rxn.upper_bound),abs(rxn.lower_bound))

					problem.variables.add(obj=[1], lb=[0], ub=[ubs], names=["up_"+rxn.id])
					
					#creating constraints on how these variables depend on the original flux variables
					const1 = [[rxn.id, "up_"+rxn.id],[1,-1]]
					const2 = [[rxn.id, "up_"+rxn.id],[-1,-1]]
					problem.linear_constraints.add(lin_expr=[const1, const2], senses=['L','L'], rhs=[0,0])

				elif targets[rxn.id] =='DOWN':
					#this case needs two new variables (one continuous and one binary)
					ubs = max(abs(rxn.upper_bound),abs(rxn.lower_bound))
					
					problem.variables.add(obj=[-1], lb=[0], ub=[ubs], names=["down_"+rxn.id])
					
					problem.variables.add(obj=[0], lb=[0], ub=[1], names=["down2_"+rxn.id])
					problem.variables.set_types("down2_" + rxn.id, problem.variables.type.binary)

					#creating constraints on how these variables depend on the original flux variables
					const1 = [[rxn.id, "down_"+rxn.id],[1, -1]] #v-d<=0, d>=v
					const2 = [[rxn.id, "down_"+rxn.id],[-1, -1]] #d>=-v, begge blir d>=|v|
					#here the big coefficient should be twice as big/small as the values that the original variables may have
					const3 = [[rxn.id, "down_"+rxn.id, "down2_"+rxn.id],[-1, 1, -20000]] 
					const4 = [[rxn.id, "down_"+rxn.id,"down2_"+rxn.id],[1, 1, 20000]]
					problem.linear_constraints.add(lin_expr=[const1, const2, const3, const4], senses=['L', 'L', 'L', 'L'], rhs=[0,0,0,20000])	
		else:
			problem.variables.add(obj=[0], lb=[rxn.lower_bound], ub=[rxn.upper_bound], names=[rxn.id])
		

	#adding constrains to the problem (the stoichimetries of each metabolite in each reaction are the constrains)
	#lager stoikiometrimatrise, for hver metabolitt S*v=0
	A = cobra.util.array.create_stoichiometric_matrix(model, array_type="DataFrame")
	names = A.index
	A = A.transpose()
	for name in names:
		r = A[A[name] != 0]
		problem.linear_constraints.add(lin_expr=[[r.index, r[name]]],
			senses=['E'],
			rhs=[0])

	#setting a constraint on optimality of nutrient conversion to biomass
	#optimality er optimal verdi av vekst fra forrige simulering. Tvinger bakterien til å fortsatt kunne vokse
	if optimality < 0: #maximized
		ub=(1+(1-relax))*floor(optimality/tol)*tol #rund nedover
		#ub=floor(optimality/tol)*tol
	else: #minimized
		ub=(1+(1-relax))*ceil(optimality/tol)*tol #rund oppover
		#ub=ceil(optimality/tol)*tol
	#relax er definert over, bestemmer hvor mye man tillater at veksten er lavere enn den maksimale veksten. 

#hvis en reaksjon er inhibert, settes flux til 0, hvis en reaksjon er kompetitiv, bidrar den til en sum constraints som må være større enn ub. 
	cn = []
	coeff = []
	for i, rxn_id in enumerate(env):
		rxn = model.reactions.get_by_id(rxn_id)
		if env[rxn_id] == "inh":
			inh_const = [[rxn.id],[1]]
			problem.linear_constraints.add(lin_expr=[inh_const],senses=['E'],rhs=[0]) #Sense: L;G;E. L er less than. G greater than. E equals
		if env[rxn_id] == "comp":
			coeff.append(1)
			cn.append(rxn_id)
	

#løser optimaliseringsproblemet, og hvis den feiler returnerer den -1000000. Hva er tallet som returneres her...
	const1 = [cn, coeff]
	problem.linear_constraints.add(lin_expr=[const1],senses=['G'],rhs=[ub]) #summen av comp må være større enn ub. Metabolsk aktiv. 
	problem.objective.set_sense(problem.objective.sense.minimize)
	strength = -1000000

	try:
		problem.solve()	
		strength = problem.solution.get_objective_value()
		return strength
	except:
		return strength



####
####
"""
def evolveX_gene(model_gene, model_gene_opt, env, target_genes, target_dirs):
	#input:
	#model_gene
	#model_gene_opt
	#env is a dictionary with reaction ids as keys and either "comp" or "inh as values" 090821 PJ
	#target_genes -> vector of target gene activity variable indeces in the model
	#targetDirs -> a vector of target gene activity desired directions of change (UP/DOWN)
	#output:
	#strength -> the worst case of selection pressure on the target genes created by the individual chemical environment

	#setting the tolerances
	tol = 1e-9
	problem = cplex.Cplex()
	problem.parameters.simplex.tolerances.optimality = 1e-9
	problem.parameters.randomseed.set(1)

	#adding variables to the problem, for gene activity a variable is created
	target_variables = []
	target_variable_dirs = []
	for u_gene in model_gene.u_reactions:
		if 'UP' in dir:
			#if the correpsonding targetDirs is UP the flux variable gets coefficient 1 (minimize) in the objective
			problem.variables.add(obj=[1], lb=[0], ub=[u_gene.upper_bound], names=[u_gene.id])
			target_variables.append(u_gene.id)
			target_variable_dirs.append("UP")
		else:
			#if the correpsonding targetDirs is DOWN the flux variable gets coefficient -1 (maximize) in the objective
			problem.variables.add(obj=[-1], lb=[0], ub=[u_gene.upper_bound], names=[u_gene.id])
			target_variables.append(u_gene.id)
			target_variable_dirs.append("DOWN")
			#the fluxes that are not targets_genes or have lower bound < 0 get 0 in the objective
		#else:
		#	problem.variables.add(obj=[0], lb=[rxn.lower_bound], ub=[rxn.upper_bound], names=[rxn.id])

	#adding constrains to the problem (the stoichimetries of each metabolite in each reaction are the constrains)
	A = cobra.util.array.create_stoichiometric_matrix(modelin, array_type="DataFrame")
	names = A.index
	A = A.transpose()
	for name in names:
		r = A[A[name] != 0]
		problem.linear_constraints.add(lin_expr=[[r.index, r[name]]],
			senses=['E'],
			rhs=[0])

	#setting a constraint on optimality of nutrient conversion to biomass
	ub=1.01*floor(modelin_opt/tol)*tol
	#print(ub)
	cn = []
	coeff = []
	#turned env into dictionary to handle inhibitors 090821 PJ
	for index,rxn in enumerate(env):
		rxn_reaction = modelin.reactions.get_by_id(rxn)
		if env.get(rxn) == "inh":
			inh_const = [[rxn_reaction.id],[1]]
			problem.linear_constraints.add(lin_expr=[inh_const],senses=['E'],rhs=[0])
		if env.get(rxn) == "comp":
			coeff.append(1)
			cn.append(rxn_reaction.id)
			#print(rxn_reaction.id)
	const1 = [cn, coeff]

	problem.linear_constraints.add(lin_expr=[const1],senses=['G'],rhs=[ub])

	problem.objective.set_sense(problem.objective.sense.minimize)
	
	#print(problem.linear_constraints.get_num())
	try:
		problem.solve()
		val= -1000000
#		print(problem.solution.get_status())
		#print("Status: "+ str(problem.solution.get_status))
		#if problem.solution.get_status() == 1 or problem.solution.get_status() == 101:
		values = problem.solution.get_values(target_variables)
		values = np.absolute(values)
		val = np.sum(values)
#			print(val)
		return val
	except:
		val = -1000000
#		print(val)
		return val
"""
