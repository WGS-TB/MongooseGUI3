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
from AllFL import *

# Following script would find all free lunch metabolites for a given set of models with db format
# and save the results on db file. Please set up the path to folder contains models that parsed and saved by Mongoose
# Please set up the path to save your results


output_dict = defaultdict(list)
indexes_output_dict = defaultdict(list)
ModelCounter = 0
# Setup the correct path
for modelNameWext in os.listdir('/path/to/your/Models/DB/files/'):
    if modelNameWext.endswith(".db"):
        # model's name with out ext.
        modelName = os.path.splitext(modelNameWext)[0]
        print("Model Name: ", modelName)
        s = shelve.open('/path/to/your/Models/DB/files/' + modelNameWext)
        model = s[modelName]
        # find sloppy reactions
        N, removed_column_index, p, ne = remove_unusful_reaction(model.fullMatrix)
        # find biomass reaction
        biomass_reaction = model.findBiomassReaction()
        # make a copy of sloppy reactions - blocked reactions are the reactions that would cause false free lunches
        #  because of their representation in stoichiometry matrix
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
        m, n = getSize(model.fullMatrix)
        constraint_rhs = np.zeros(m)
        # Remove desired blocked reactions
        deleted_reactions_bd = [None] * n
        for i in blocked_reactions:
            deleted_reactions_bd[i] = 0
        p, status, n = all_free_lunches(model.fullMatrix, constraint_rhs, None, deleted_reactions_bd,
                                        deleted_reactions_bd)
        ans_vector_x, solution_non_zero_elements, solution_non_zero_indexes, ans_vector_y, y_non_zero_elements, y_non_zero_indexes = all_FL_ans(
            p, status, m, n)
        output_dict[modelName] = copy.deepcopy(ans_vector_y)
        indexes_output_dict[modelName] = copy.deepcopy(y_non_zero_indexes)
        output_shelved_dict = shelve.open("/path/to/save/the/results/AllFLOutPut.db")
        output_shelved_dict.update(output_dict)
        output_shelved_dict.close()
        ####
        indexes_output_shelved_dict = shelve.open("/path/to/save/the/results/AllFLIndexOutPut.db")
        indexes_output_shelved_dict.update(indexes_output_dict)
        indexes_output_shelved_dict.close()
        ModelCounter += 1
        print(ModelCounter)
