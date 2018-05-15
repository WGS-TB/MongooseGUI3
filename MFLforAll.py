import os
from ModelParsing import *
import shelve
import sys
import time
import pandas as pd
import json
from MinimalFreeLunch import *
import copy
from collections import defaultdict
from FeaturePreparation import convertFormula
#import elementallyBalancedCheck

output_dict = defaultdict(list)
indexes_output_dict = defaultdict(list)
ModelCounter = 0
for modelNameWext in os.listdir('/home/hzabeti/qsopt-ex/build/python-qsoptex/MongooseGUI3/Models/DB/'):
    if modelNameWext.endswith(".db"):
        # model's name with out ext.
        modelName = os.path.splitext(modelNameWext)[0]
        print("Model Name: ", modelName)
        s = shelve.open('/home/hzabeti/qsopt-ex/build/python-qsoptex/MongooseGUI3/Models/DB/'+modelNameWext)
        model= s[modelName]
        # find sloppy reactions
        N, removed_column_index, p, ne = remove_unusful_reaction(model.fullMatrix)
        # find biomass reaction
        biomass_reaction = model.findBiomassReaction()
        # make a copy of sloppy reactions - blocked reactions are the reactions that would cause false free lunches
        # because of their representation in stoichiometric matrix
        blocked_reactions = copy.copy(removed_column_index)
        # add biomass reaction to blocked reactions - we are not interested in biomass reaction for now
        blocked_reactions.append(biomass_reaction)
        # find additional biomass reactions
        biomass_reaction_candidates = []
        for r in model.reactions:
            reaction_name_biomass = r.name.lower()
            if reaction_name_biomass.find('biomass') != -1:
                biomass_reaction_candidates.append({'name': r.name, 'index': r.index})
                blocked_reactions.append(r.index)
        # size of matrix
        m,n = getSize(model.fullMatrix)

        constraint_rhs = np.zeros(m)
        # Remove desired blocked reactions
        deleted_reactions_bd = [None]*n
        for i in blocked_reactions:
            deleted_reactions_bd[i]=0
        
        p, status, n = minimum_free_lunches(model.fullMatrix,constraint_rhs,None,deleted_reactions_bd,deleted_reactions_bd)
        solution_non_zero_elements, solution_non_zero_indexes_w, ans_vector_x_w, ans_vector_z = l1_equiv_lp_ans(p, status, n)
        output_dict[modelName] = copy.deepcopy(ans_vector_x_w)
        indexes_output_dict[modelName] = copy.deepcopy(solution_non_zero_indexes_w)
        print('---------------------------------------------------OUTPUT----------------------------------------------')
        print('solution_non_zero_indexes', solution_non_zero_indexes_w)
        print('solution_non_zero_elements',solution_non_zero_elements)
# save dicts
        output_shelved_dict = shelve.open("/home/hzabeti/Dropbox/SFU/MFLforALLOutPut.db")
        output_shelved_dict.update(output_dict)
        output_shelved_dict.close()
    ####
        indexes_output_shelved_dict = shelve.open("/home/hzabeti/Dropbox/SFU/MFLforALLIndexOutPut.db")
        indexes_output_shelved_dict.update(indexes_output_dict)
        indexes_output_shelved_dict.close()
        ModelCounter += 1
        f = open('/home/hzabeti/Dropbox/SFU/MFL.txt','a')
        f.write('Model counter:'+str(ModelCounter)+'\n')
        f.write('Model name:'+modelName+'\n')
        f.close()
f = open('/home/hzabeti/Dropbox/SFU/MFLforAll.txt','w')
f.close()
db = shelve.open('/home/hzabeti/Dropbox/SFU/MFLforALLIndexOutPut.db')
models_list = list(db.keys())
f = open('/home/hzabeti/Dropbox/SFU/MFLforAll.txt','a')
for i in models_list:
    f.write('Model Name:'+i+'\n')
    if len(db[i])!= 0:
        f.write('Free Lunches: True\n')
        f.write('Reactions that cause free lunch with minimum l1 norm:'+str(db[i])+'\n')
    else:
        f.write('Free Lunches: False\n')
    f.write('------------------------------------------------------------\n')
f.close()