import cobra
from cobra.flux_analysis import flux_variability_analysis
from determine_trait_flux_basis import determine_trait_flux_basis
from enumerate_evolvex import enumerate_evolvex
from get_growth_supporting_niche12 import get_growth_supporting_niche12
from EVOLVEX_SCORE import setBounds, createSecretion

model = cobra.io.read_sbml_model("iML1515.xml")
#1.Nullstiller objektivfunksjon og setter biomasse som objektiv
for rxn in model.reactions:
    rxn.objective_coefficient = 0
model.reactions.get_by_id("BIOMASS_Ec_iML1515_core_75p37M").objective_coefficient = 1


#2.Create a copy for flux-basis determination and run FVA to define wild-type flux ranges
model_wt = model.copy()

fva_result = flux_variability_analysis(
    model_wt,
    reaction_list=model_wt.reactions,
    fraction_of_optimum=1.0)
wtminmax = fva_result[['minimum', 'maximum']].to_numpy() #benyttes i determine_trait_flux_basis etterpå
#WT intervall er når modellen optimaliserer vekst. 
#lb og ub er for trait-modellen med ønsket objektivfunksjon. 


#3.Definerer den nye objektivfunksjonen
trait_model = model.copy() #M9 medium
for rxn in trait_model.reactions:
    rxn.objective_coefficient = 0

trait=trait_model.reactions.get_by_id("EX_lac__D_e")
trait.lower_bound=0 #blokkerer opptak. Hindre at den ikke øke fluxen ved å ta opp direkte fra medium. 
trait.upper_bound=1000 
trait.objective_coefficient = 1 #setter ønsket reaksjon som objektivfunksjon


#4.Pro env
lb_open = ["EX_ala__L_e","EX_glu__L_e","EX_gln__L_e","EX_asp__L_e","EX_ser__L_e","EX_gly_e","EX_leu__L_e","EX_lys__L_e", "EX_thm_e", "EX_btn_e"]
for e in lb_open:
    if e in trait_model.reactions:
        trait_model.reactions.get_by_id(e).lower_bound = -10


#5.Determine_trait_flux_basis
targets, dirs, signs, dirSolp, binRxns, fullSolOpt, solutionp = determine_trait_flux_basis(trait_model, 'max', wtminmax, 1, 1, 0.00001, 5)


#6.parse-funksjon som fikser evolvex input
targets_flat = targets[0]
dirs_flat = dirs[0]
evolveX_targets = [int(x) for x in targets_flat]
evolveX_target_dirs = list(dirs_flat)
# lag dict til evolveX_score
targets_dict = {model.reactions[i].id: d for i, d in zip(evolveX_targets, evolveX_target_dirs)}


#7.Prepare EvolveX-model som stenger glukose for å finne alternative ruter
model_evolvex = model.copy()
model_evolvex.reactions.get_by_id("EX_glc__D_e").lower_bound = 0


#8.lager kombinasjoner som støtter vekst med alternative kilder
condMap = {
    "carbon":  ["EX_gal_e", "EX_fru_e", "EX_ac_e", "EX_pyr_e", "EX_succ_e", "EX_glyc_e", "EX_xyl__D_e"],
    "nitrogen": ["EX_gln__L_e", "EX_glu__L_e","EX_ala__L_e", "EX_ser__L_e", "EX_nh4_e"]
}
growth, uptakes, combos = get_growth_supporting_niche12(model_evolvex, condMap)


#9.lager envs fra get_growth_supporting_niche1 og finner optimality
envs = []
for (c, n1, n2) in combos:   
    env = {
        c: "comp",
        n1: "comp",
        n2: "comp",
        "EX_glc__D_e": "inh",
    }
    envs.append(env)

components = sorted(set(rxn for combo in combos for rxn in combo))
setBounds(model_evolvex)
for ex in components:
    createSecretion(model_evolvex, ex)

for ex in components:
    r = model_evolvex.reactions.get_by_id(ex)
    r.lower_bound = 0
    r.upper_bound = 0


#10.sammenligner targets med alle mulige envs for å gi en evolvex score på beste miljø
strength, comp = enumerate_evolvex(model_evolvex, envs, targets_dict)
print(targets_dict)
