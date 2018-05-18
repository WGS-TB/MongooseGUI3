import copy
import math
import cplex
import os
from MinimalFreeLunch import *

# Following functions would help us find a set of free lunch witnesses for specific model
# For further analysis, first, we save the problem in mps format
# and we solve the problem via Cplex 12.8.0 API for Python afterward.
# To find FLWs:
# 1) import FreeLunchWitness
# 2) use function model_mps_writer(model,'/path/you/want/to/save/mps/file/.mps') to save the problem in mps format
# 3) use function flw_finder_via_mps('/path/of/saved/mps/file/.mps')


def model_to_mps(N, y_lower_bound, y_upper_bound, deleted_columns, mps_file_path):
    m, n = getSize(N)
    mps_file = os.path.abspath(mps_file_path)
    opt_file = open(mps_file, 'w')
    opt_file.write('NAME\tProblem\n')
    opt_file.write('OBJSENSE\n')
    opt_file.write('\tMIN\n')
    opt_file.write('OBJNAME\n')
    opt_file.write('\tOBJ\n')
    # ---------------------- Constraints -------------------------
    '''
    rows_dict = {}
    for j in range(n):
        current_dict = {'y' + str(i): N[i][j] for i in range(m)}
        current_dict.update({'x' + str(j): 1})
        rows_dict.update({'row' + str(j): copy.copy(current_dict)})
    for j in range(n):
        current_dict = {'y' + str(i): N[i][j] for i in range(m)}
        current_dict.update({'x' + str(j): -1})
        rows_dict.update({'row' + str(n + j): copy.copy(current_dict)})
    '''
    # ---------------------- Print ROWS -------------------------
    opt_file.write('ROWS\n')
    opt_file.write(' N   OBJ\n')
    for i in range(n):
        if i not in deleted_columns:
            opt_file.write(' G   c' + str(i) + '\n')
            opt_file.write(' G   c' + str(n + i) + '\n')
    # ---------------------- Print COLUMNS -------------------------
    opt_file.write('COLUMNS\n')
    for i in range(n):
        if i not in deleted_columns:
            opt_file.write('   x' + str(i) + '    ' + 'OBJ\t1\n')
            opt_file.write('   x' + str(i) + '    ' + 'c' + str(i) + '\t1\n')
            opt_file.write('   x' + str(i) + '    ' + 'c' + str(n + i) + '\t1\n')
    for i in range(m):
        for j in range(n):
            if j not in deleted_columns:
                opt_file.write('   y' + str(i) + '    ' + 'c' + str(j) + '\t' + str(float(N[i][j])) + '\n')
                opt_file.write('   y' + str(i) + '    ' + 'c' + str(n + j) + '\t' + str(-float(N[i][j])) + '\n')
    # ----------------------- RHS ------------------------------
    # ALL of them are zero so we don't need to mention them
    # ----------------------- BOUNDS ------------------------------
    opt_file.write('BOUNDS\n')
    for i in range(m):    
        opt_file.write('   LO   BND1\ty' + str(i) + '\t' + str(0.0001) + '\n')
    for j in range(n):
        if j not in deleted_columns:
            opt_file.write('   BV   BND1\tx' + str(j) + '\n')
    opt_file.write('ENDATA')
    opt_file.close()


def matrix_norm_two(N):
    m, n = getSize(N)
    norm_array = []
    for j in range(n):
        current_norm = 0
        for i in range(m):
            current_norm += math.pow(N[i][j],2)
        norm_array.append(current_norm)
    return max(norm_array), norm_array

def y_bounds(N,norm_max):
    m, n = getSize(N)
    max_meta_atoms = 1000
    y_lower_bound = 1/(max_meta_atoms*(int(math.sqrt(m))+1)*(int(norm_max)+1))
    y_upper_bound = 1/((int(math.sqrt(m))+1)*(int(norm_max)+1))
    return y_lower_bound,y_upper_bound

def mps_deleted_columns(model):
    N, removed_column_index, p, ne = remove_unusful_reaction(model.fullMatrix)
    biomass_reaction = model.findBiomassReaction()
    blocked_reactions = copy.copy(removed_column_index)
    blocked_reactions.append(biomass_reaction)
    biomass_reaction_candidates = []
    for r in model.reactions:
        reaction_name_biomass = r.name.lower()
        if reaction_name_biomass.find('biomass') != -1:
            biomass_reaction_candidates.append({'name': r.name, 'index': r.index})
            blocked_reactions.append(r.index)
    return blocked_reactions


def model_mps_writer(model, mps_file_path):
    N = model.fullMatrix
    deleted_columns = mps_deleted_columns(model)
    max_norm, norm_array = matrix_norm_two(N)
    y_lower_bound,y_upper_bound = y_bounds(N,max_norm)
    model_to_mps(N,y_lower_bound,y_upper_bound,deleted_columns,mps_file_path)


def flw_finder_via_mps(file_path):
    mps_file =os.path.abspath(file_path)
    flw_problem = cplex.Cplex(mps_file)
    flw_problem.solve()
    flw_objective_value = flw_problem.solution.get_objective_value()
    flw_list = [v for v in flw_problem.variables.get_names() if v[0] == 'x' and flw_problem.solution.get_values(v)!=0]
    return flw_objective_value,flw_list
