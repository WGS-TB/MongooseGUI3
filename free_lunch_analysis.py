import qsoptex
from Utilities import *
from fractions import Fraction
import copy
import shelve
import numpy as np
from ModelParsing import *
#import cplex

# Please see the associated paper for further details:
# 'A duality-based method for identifying elemental balance violations in metabolic network models' \
# Hooman Zabeti, Tamon Stephen, Bonnie Berger and Leonid Chindelevitch

#TO DO: add ref for cplex and qsopt_ex and Boyd's paper


def find_import_export_reactions(N, reactions_type ='all'):
    # This function returns list of import and export reactions
    # i.e. those reactions that have only positive or only negative 
    # stoichiometric coefficients. N here is the stoichiometric matrix
    # of the model (fullMatrix) and reaction_type represent type of reactions
    # (i.e. you can select either 'import' for import reactions, 'export' for export reactions
    # or 'all' for both import and export. Default value of reactions_type is 'all').
    # Here we assumed that all reactions are irreversible.
    m,n = getSize(N)
    im_reactions_index=[]
    ex_reactions_index=[]
    for i in range(n):
        ith_reaction = [row[i] for row in N if row[i]!=0]
        if all(x> 0 for x in ith_reaction): ex_reactions_index.append(i)
        elif all(x< 0 for x in ith_reaction): im_reactions_index.append(i)
    if reactions_type=='all':
        return sorted(im_reactions_index+ex_reactions_index)
    elif reactions_type =='import':
        return im_reactions_index
    elif reactions_type == 'export':
        return ex_reactions_index


def find_biomass_reaction_candidates(model_obj):
    # This function returns list of biomass reaction candidates
    # i.e. those containing the word 'biomass' or 'growth' in their name.
    #biomass_reactions = [r.index for r in model_obj.reactions if r.name.lower().find('biomass')!=-1\
        #or r.name.lower().find('growth')!=-1]
    biomass_reactions = [r.index for r in model_obj.reactions if r.name.lower().find('biomass')!=-1 or r.name.lower().find('growth')!=-1]
    if model_obj.findBiomassReaction()!=-1:
        biomass_reactions.append(model_obj.findBiomassReaction())
    return biomass_reactions


def excluded_reactions(model_obj):
    # This function returns list of reactions which we want to exclude from free lunch process
    # i.e. pseudo reactions (biomass and growth) and import/export reactions.
    m,n = getSize(model_obj.fullMatrix)
    excluded_reactions_index = find_biomass_reaction_candidates(model_obj)+find_import_export_reactions(model_obj.fullMatrix)
    excluded_reactions_list = [None]*n
    for i in excluded_reactions_index:
        excluded_reactions_list[i]=0
    return excluded_reactions_list


def free_lunch_metabolites(N, constraint_rhs=None, objective_weight_vector=None, variable_lower_bound=None, variable_upper_bound=None):
    # This function returns the solution of the linear program that we use to verify No Free Lunch condition and find all free lunch metabolites
    # i.e. eq.(2) in Step 1 in the paper.
    # max(1^Ty) subject to Sx-y =>0, 0<= y<= 1 where S is the stoichiometric matrix. 
    m, n = getSize(N)
    if constraint_rhs is None:
        constraint_rhs = [0]*m
    if objective_weight_vector is None:
        objective_weight_vector = [1]*m
    if variable_lower_bound is None:
        variable_lower_bound = [None]*n
    if variable_upper_bound is None:
        variable_upper_bound = [None]*n
    p = qsoptex.ExactProblem()
    for i in range(n):
        p.add_variable(name='x' + str(i), objective=0, lower=variable_lower_bound[i], upper=variable_upper_bound[i])
    for i in range(m):
        p.add_variable(name='y' + str(i), objective=objective_weight_vector[i], lower=0, upper=1)
    for j in range(m):
        constraints_dict = {'x' + str(i): N[j][i] for i in range(n)}
        constraints_dict.update({'y' + str(j): -1})
        p.add_linear_constraint(qsoptex.ConstraintSense.GREATER, constraints_dict, rhs=constraint_rhs[j])
        constraints_dict.clear()
    p.set_objective_sense(qsoptex.ObjectiveSense.MAXIMIZE)
    p.set_param(qsoptex.Parameter.SIMPLEX_DISPLAY, 1)
    status = p.solve()
    return p, status

# ---------------------------------------------- Step 1------------------------------------------------
# No Free Lunch condition and Free Lunch Metabolites

def find_free_lunch_metabolites(model_obj):
    # This function returns list of indexes of all free lunch metabolites.i.e. Step 1. 
    # Note that if optimum solution is zero or if returned list is empty, then No Free Lunch Condition is satisfied.
    #optimum_solution = -1
    free_lunch_meta_index = []
    m,n  = getSize(model_obj.fullMatrix)
    excluded_reactions_list = excluded_reactions(model_obj)
    p, status = free_lunch_metabolites(model_obj.fullMatrix,None,None,excluded_reactions_list,excluded_reactions_list)
    if status == qsoptex.SolutionStatus.OPTIMAL:
        optimum_solution = p.get_objective_value()
        free_lunch_meta_index = [j for j in range(m) if p.get_value('y'+str(j))!=0]
    # ---------------------------------------- Print section ------------------------------------
    # Print section can be removed if needed.
        if optimum_solution!=0:
            print('No Free Lunch condition: Not satisfied!')
            print('Optimum solution and the number of all possible free lunch metabolites',optimum_solution)
        else:
            print('No Free Lunch condition: Satisfied!')
    # -------------------------------------------------------------------------------------------
    flm_list = [model_obj.metabolites[int(i)].species.name for i in free_lunch_meta_index]
    return free_lunch_meta_index,flm_list

# ---------------------------------------------- Step 2------------------------------------------------
# Minimal Free Lunch

def free_lunch_finder(N, constraint_rhs, objective_weight_vector=None, variable_lower_bound=None, variable_upper_bound=None):
    # This function returns an instance of Free Lunch for a given list of free lunch metabolits
    # (list of FL metabolites respresented here as constraints_rhs). This function along with 
    # re_weighted_linear_program function (below), have been used as first part of Step 2 in the paper, which
    # provide an approximation for eq 5. min1^Tz subject to Sx=> constraint_rhs, z-x=>0, z+x=>0.
    m, n = getSize(N)
    if objective_weight_vector is None:
        objective_weight_vector = [1]*n
    if variable_lower_bound is None:
        variable_lower_bound = [None]*n
    if variable_upper_bound is None:
        variable_upper_bound = [None]*n
    p = qsoptex.ExactProblem()
    for i in range(n):
        p.add_variable(name='z' + str(i), objective=objective_weight_vector[i], lower=None, upper=None)
        p.add_variable(name='x' + str(i), objective=0, lower=variable_lower_bound[i], upper=variable_upper_bound[i])
    for j in range(m):
        p.add_linear_constraint(qsoptex.ConstraintSense.GREATER,
                                {'x' + str(i): N[j][i] for i in range(n)}, rhs=constraint_rhs[j])
    for i in range(n):
        p.add_linear_constraint(qsoptex.ConstraintSense.GREATER,
                                {'z' + str(i): 1, 'x' + str(i): 1}, rhs=0)
        p.add_linear_constraint(qsoptex.ConstraintSense.GREATER,
                                {'z' + str(i): 1, 'x' + str(i): -1}, rhs=0)
    p.set_objective_sense(qsoptex.ObjectiveSense.MINIMIZE)
    p.set_param(qsoptex.Parameter.SIMPLEX_DISPLAY, 1)
    status = p.solve()
    return p, status

def re_weighted_linear_program(N, constraint_rhs, variable_lower_bound=None, variable_upper_bound=None):
    m, n = getSize(N)
    # Initial values
    flag = True
    optimal_solution_length = n+1
    iteration = 0
    w = [1]*n
    epsilon = Fraction(1, 10)
    stop_num = 2
    solution_non_zero =[]
    if variable_lower_bound is None:
        variable_lower_bound = [None]*n
    if variable_upper_bound is None:
        variable_upper_bound = [None]*n

    while iteration < stop_num and flag:
        p, status = free_lunch_finder(N, constraint_rhs, w, variable_lower_bound, variable_upper_bound)
        if status == qsoptex.SolutionStatus.OPTIMAL:
            solution_non_zero = [(i,p.get_value('x' + str(i))) for i in range(n) if p.get_value('x' + str(i)) != 0]
            if len(solution_non_zero) >= optimal_solution_length:
                flag = False
            else:
                for i in range(n):
                    w[i] = Fraction(1, p.get_value('z' + str(i)) + epsilon)
                optimal_solution_length = len(solution_non_zero)
                #------------------------Print-------------------------------------
                # Print section can be removed if needed
                print('Number of iteration: ',iteration)
                print('Number of reactions involved in FL:', len(solution_non_zero))
                print('------------------------------------------------------------')
                #------------------------ end of print -----------------------------
        else:
            print('LP does not have a solution for these constraints')
            iteration = stop_num
        iteration += 1

    return solution_non_zero


def free_lunch_check(N, temp_vector, constraint_rhs):
    free_lunch_status = False
    p, status = free_lunch_finder(N, constraint_rhs, None, temp_vector, temp_vector)
    if status == qsoptex.SolutionStatus.OPTIMAL:
        free_lunch_status = True
    return free_lunch_status,p,status


def extract_minimal_set(N, reactions_list, constraint_rhs):
    m,n = getSize(N)
    final_minimal = []
    minimal_set_indexes = [i[0] for i in reactions_list]
    minimal_vector = [None if i in minimal_set_indexes else 0 for i in range(n)]
    final_p = None
    for i in minimal_set_indexes:
        temp_vector = copy.copy(minimal_vector)
        temp_vector[i] = 0
        print('Checking index: ', i)
        free_lunch_status,p,status = free_lunch_check(N, temp_vector, constraint_rhs)
        if free_lunch_status:
            print('This index is not in the minimal', i)
            minimal_vector = copy.copy(temp_vector)
        else:
            print('This index is in the minimal set', i)
            final_minimal.append(i)
    return final_minimal

def minimum_free_lunches(N, objective_weight_vector=None, variable_lower_bound=None, variable_upper_bound=None):
    m, n = getSize(N)
    if objective_weight_vector is None:
        objective_weight_vector = [1]*n
    if variable_lower_bound is None:
        variable_lower_bound = [None]*n
    if variable_upper_bound is None:
        variable_upper_bound = [None]*n
    p = qsoptex.ExactProblem()
    for i in range(n):
        p.add_variable(name='z' + str(i), objective=objective_weight_vector[i], lower=None, upper=None)
        p.add_variable(name='x' + str(i), objective=0, lower=variable_lower_bound[i], upper=variable_upper_bound[i])
    # The following part will add non-negative slack variable y
    for i in range(m):
        p.add_variable(name='y' + str(i), objective=0, lower=0, upper=None)
    for j in range(m):
        constraints_dict = {'x' + str(i): N[j][i] for i in range(n)}
        constraints_dict.update({'y' + str(j): -1})
        p.add_linear_constraint(qsoptex.ConstraintSense.EQUAL, constraints_dict, rhs=0)
        constraints_dict.clear()
    for i in range(n):
        p.add_linear_constraint(qsoptex.ConstraintSense.GREATER,
                                {'z' + str(i): 1, 'x' + str(i): 1}, rhs=0)
        p.add_linear_constraint(qsoptex.ConstraintSense.GREATER,
                                {'z' + str(i): 1, 'x' + str(i): -1}, rhs=0)
    # Adding tolerance to whole system
    p.add_linear_constraint(qsoptex.ConstraintSense.GREATER, {'y' + str(i): 1 for i in range(m)},
                            rhs=1)
    # ------------
    p.set_objective_sense(qsoptex.ObjectiveSense.MINIMIZE)
    p.set_param(qsoptex.Parameter.SIMPLEX_DISPLAY, 1)
    status = p.solve()
    return p, status


def min_free_lunch_finder_ans(p, status, m,n):
    if status == qsoptex.SolutionStatus.OPTIMAL:
        constraints_non_zero = [(i,p.get_value('y' + str(i))) for i in range(m) if p.get_value('y' + str(i)) != 0]
        solution_non_zero = [(i, p.get_value('x' + str(i))) for i in range(n) if p.get_value('x' + str(i)) != 0]
    elif status == qsoptex.SolutionStatus.INFEASIBLE:
        print('Problem is infeasible!')
    return solution_non_zero, constraints_non_zero


def minimal_free_lunch_finder(model_obj,FL_meta_index_list = None,print_mess = True):
    m, n = getSize(model_obj.fullMatrix)
    excluded_reactions_list = excluded_reactions(model_obj)
    minimal_free_lunch = []
    free_lunch_list = []
    cons= []
    if FL_meta_index_list is None:
        p, status =  minimum_free_lunches(model_obj.fullMatrix,None,excluded_reactions_list,excluded_reactions_list)
        if status == qsoptex.SolutionStatus.OPTIMAL:
            FL_meta_index_list = [i for i in range(m) if p.get_value('y' + str(i)) != 0]
            FL_meta_list = [0]*m
            FL_meta_list[FL_meta_index_list[0]]=1
            reactions_list = re_weighted_linear_program(model_obj.fullMatrix, FL_meta_list, excluded_reactions_list,excluded_reactions_list)
            minimal_free_lunch = extract_minimal_set(model_obj.fullMatrix, reactions_list,FL_meta_list)
            print("Minimal free lunch indexes:", minimal_free_lunch)
    elif FL_meta_index_list is not None:
        FL_meta_list = [1 if i in FL_meta_index_list else 0 for i in range(m)]
        reactions_list = re_weighted_linear_program(model_obj.fullMatrix, FL_meta_list, excluded_reactions_list,excluded_reactions_list)
        minimal_free_lunch = extract_minimal_set(model_obj.fullMatrix, reactions_list,FL_meta_list)
        print("Minimal free lunch indexes:", minimal_free_lunch)
    #------------------------Print message-------------------------------------
    if print_mess==True and len(minimal_free_lunch)!=0 :
        print('------------------ minimal free lunch equation----------------')
        minimal_vector = [None if i in minimal_free_lunch else 0 for i in range(n)]
        a,b,c = free_lunch_check(model_obj.fullMatrix,minimal_vector,FL_meta_list)
        if c== qsoptex.SolutionStatus.OPTIMAL:
            free_lunch_list = [(i,b.get_value('x'+str(i))) for i in range(n) if b.get_value('x' + str(i)) != 0]
        ans_x = [0]*n
        for i in free_lunch_list:
            ans_x[i[0]] = i[1]
        cons_list = np.dot(model_obj.fullMatrix,ans_x)
        cons = [(i,cons_list[i]) for i in range(m) if cons_list[i]!=0]
        for i in free_lunch_list:
            print(str(i[1])+'*('+model_obj.printReactionFormula(i[0])+')')
        print('-----------')
        for j in cons:
            print(str(j[1])+'*('+model_obj.metabolites[j[0]].species.name +')')
        print('-------------------------------------------------------------')
    #-------------------------------------------------------------------
    return minimal_free_lunch,free_lunch_list,cons

def free_lunch_finder_ans(p, status, n):
    if status == qsoptex.SolutionStatus.OPTIMAL:
        solution_non_zero = [(i,p.get_value('x' + str(i))) for i in range(n) if p.get_value('x' + str(i)) != 0]
    elif status == qsoptex.SolutionStatus.INFEASIBLE:
        print('Problem is infeasible!')
    return solution_non_zero



def standard_basis(basis_dimention, basis_index):
    # This method is faster than numpy.eye
    e = np.zeros(basis_dimention)
    e[basis_index] = 1
    return e


def ans_check(N, ans_vector_x, e, n):
    check_array = np.dot(N, ans_vector_x)
    check_array_non_zero = {}
    flag = True
    for i in range(n):
        if check_array[i]!=0:
            check_array_non_zero.update({'index:'+str(i): check_array[i]})
        if check_array[i] - e[i] < 0:
            flag = False
            print(i)
    print('This solution is:')
    print(flag)
    return check_array, check_array_non_zero, flag

# ---------------------------------------------- Step 3------------------------------------------------
# Individually producible free lunch metabolites

def find_ind_flm(model_obj):
    flm_index_list = []
    flm_list = []
    excluded_reactions_index = find_biomass_reaction_candidates(model_obj)+find_import_export_reactions(model_obj.fullMatrix)
    updated_matrix = copy.copy(model_obj.fullMatrix)
    for i in updated_matrix:
        for j in range(len(i)):
            if j in excluded_reactions_index:
                i[j]=0
    MT = transpose(updated_matrix)
    RM = GaussJordan(MT)
    for i in range(len(RM[0])):
        current_row = []
        current_row = [k for k in range(len(RM[0][i])) if RM[0][i][k]!=0]
        if len(current_row)==1:
            flm_index_list.append(current_row[0])
    flm_list = [model_obj.metabolites[int(i)].species.name for i in flm_index_list]
    #----------------------------------Print----------------------------------------------------------
    print('Number of individually producible free lunch metabolites',len(flm_index_list))
    #-------------------------------------------------------------------------------------------------
    return flm_index_list,flm_list
# ---------------------------------------------- Additional--------------------------------------------
# No free Lunch Condition - Faster approach

def no_free_lunch_condition(model_obj):
    N = model_obj.fullMatrix
    m, n = getSize(N)
    excluded_reactions_list = excluded_reactions(model_obj)
    p = qsoptex.ExactProblem()
    for i in range(n):
        p.add_variable(name='x' + str(i), objective=0, lower=excluded_reactions_list[i], upper=excluded_reactions_list[i])
    # The following part will add non-negative slack variable y
    for i in range(m):
        p.add_variable(name='y' + str(i), objective=0, lower=0, upper=None)
    for j in range(m):
        constraints_dict = {'x' + str(i): N[j][i] for i in range(n)}
        constraints_dict.update({'y' + str(j): -1})
        p.add_linear_constraint(qsoptex.ConstraintSense.EQUAL, constraints_dict, rhs=0)
        constraints_dict.clear()
    # Adding tolerance to whole system
    p.add_linear_constraint(qsoptex.ConstraintSense.GREATER, {'y' + str(i): 1 for i in range(m)},
                            rhs=1)
    # ------------
    p.set_objective_sense(qsoptex.ObjectiveSense.MINIMIZE)
    p.set_param(qsoptex.Parameter.SIMPLEX_DISPLAY, 1)
    status = p.solve()
    return p, status