import qsoptex
import numpy as np
from Utilities import *
from fractions import Fraction
import copy
import shelve

def l1_equiv_lp(N, constraint_rhs, objective_weight_vector=None, variable_lower_bound=None, variable_upper_bound=None):
    m, n = getSize(N)
    # alpha = Fraction(1, 10)
    # following part set a default weight vector
    if objective_weight_vector is None:
        objective_weight_vector = np.ones(n)
    if variable_lower_bound is None:
        variable_lower_bound = [None]*n
    if variable_upper_bound is None:
        variable_upper_bound = [None]*n
    p = qsoptex.ExactProblem()
    # In the following part we define the objective of the linear programming
    for i in range(n):
        p.add_variable(name='z' + str(i), objective=objective_weight_vector[i], lower=None, upper=None)
        p.add_variable(name='x' + str(i), objective=0, lower=variable_lower_bound[i], upper=variable_upper_bound[i])
    # In the following part we define constrains of the linear Programming
    for j in range(m):
        p.add_linear_constraint(qsoptex.ConstraintSense.GREATER,
                                {'x' + str(i): N[j][i] for i in range(n)}, rhs=constraint_rhs[j])
    for i in range(n):
        p.add_linear_constraint(qsoptex.ConstraintSense.GREATER,
                                {'z' + str(i): 1, 'x' + str(i): 1}, rhs=0)
        p.add_linear_constraint(qsoptex.ConstraintSense.GREATER,
                                {'z' + str(i): 1, 'x' + str(i): -1}, rhs=0)
    # Following part would solve the LP with qsoptex
    p.set_objective_sense(qsoptex.ObjectiveSense.MINIMIZE)
    p.set_param(qsoptex.Parameter.SIMPLEX_DISPLAY, 1)
    status = p.solve()
    return p, status, n


def minimum_free_lunches(N, constraint_rhs, objective_weight_vector=None, variable_lower_bound=None, variable_upper_bound=None):
    m, n = getSize(N)
    # alpha = Fraction(1, 10)
    # Set default weight vector
    if objective_weight_vector is None:
        objective_weight_vector = np.ones(n)
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
        p.add_linear_constraint(qsoptex.ConstraintSense.EQUAL, constraints_dict, rhs=constraint_rhs[j])
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
    return p, status, n


def l1_equiv_lp_ans(p, status, reactions_num):
    solution_non_zero_elements = []
    solution_non_zero_indexes = []
    ans_vector_x = []
    ans_vector_z = []

    if status == qsoptex.SolutionStatus.OPTIMAL:
        print('Optimal solution')
        print(p.get_objective_value())
        for i in range(reactions_num):
            ans_vector_x.append(p.get_value('x' + str(i)))
            ans_vector_z.append(p.get_value('z' + str(i)))
            if p.get_value('x' + str(i)) != 0:
                # print('x' + str(i), p.get_value('x' + str(i)))
                # print('z' + str(i), p.get_value('z' + str(i)))
                solution_non_zero_elements.append(p.get_value('x' + str(i)))
                solution_non_zero_indexes.append(i)
    print('l0 approximation')
    print(len(solution_non_zero_elements))
    # print('Non-zero elements of the solution')
    # print(solution_non_zero_elements)
    # print('Indexes of non-zero elements of the solution ')
    # print(solution_non_zero_indexes)
    return solution_non_zero_elements, solution_non_zero_indexes, ans_vector_x, ans_vector_z


def re_weighted_linear_program(N, e, problem_case,variable_lower_bound, variable_upper_bound):
    flag = True
    optimial_solution_lenght = 0
    m, n = getSize(N)
    iteration = 0
    w = np.ones(n)
    epsilon = Fraction(1, 10)

    while iteration <= 5 and flag:
        if int(problem_case) == 1:
            p, status, reactions_num = l1_equiv_lp(N, e, w,variable_lower_bound, variable_upper_bound)
        else:
            p, status, reactions_num = minimum_free_lunches(N, e, w,variable_lower_bound, variable_upper_bound)
        print('%d iteration' % iteration)
        solution_non_zero_elements, solution_non_zero_indexes, ans_vector_x, ans_vector_z \
            = l1_equiv_lp_ans(p, status, n)
        if status == qsoptex.SolutionStatus.OPTIMAL:
            for i in range(n):
                w[i] = Fraction(1, p.get_value('z' + str(i)) + epsilon)
            check_array, check_array_non_zero, flag = ans_check(N, ans_vector_x, e, m)
            print('check_array_non_zero:', check_array_non_zero)
            print('solution non-zero element', solution_non_zero_elements)
            print('solutions non-zero indexes', solution_non_zero_indexes)
        else:
            print('LP does not have a solution for these constraints')
            iteration = 6
        iteration += 1
        if len(solution_non_zero_elements) == optimial_solution_lenght:
            flag = False
        else:
            optimial_solution_lenght = len(solution_non_zero_elements)
    return p, status, reactions_num, solution_non_zero_elements, solution_non_zero_indexes, ans_vector_x, ans_vector_z


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


def free_lunch_check(N, reaction_vector, constraint_rhs):
    #non_zero_indexes = list(np.nonzero(reaction_vector)[0])
    non_zero_indexes = []
    for j in range(len(reaction_vector)):
        if reaction_vector[j] != 0:
            non_zero_indexes.append(j)
    lower_upper_bound = [0]*(len(reaction_vector))
    # print("TRACE 1st lower_upper_bound:", lower_upper_bound)
    free_lunch_status = False
    for i in non_zero_indexes:
        lower_upper_bound[i] = None
    # print("TRACE 2st lower_upper_bound:", lower_upper_bound)
    p, status, n = l1_equiv_lp(N, constraint_rhs, None, lower_upper_bound, lower_upper_bound)
    if status == qsoptex.SolutionStatus.OPTIMAL:
        free_lunch_status = True
    print("TRACE free lunch statuse:",free_lunch_status)
    return free_lunch_status, p, status, n


def extract_minimal_set(N, reaction_vector, constraint_rhs):
    non_zero_indexes = []
    minimal_set_indexes = []
    for j in range(len(reaction_vector)):
        if reaction_vector[j] != 0:
            non_zero_indexes.append(j)
    # non_zero_indexes = list(np.nonzero(reaction_vector)[0])
    for i in non_zero_indexes:
        temp_vector = copy.deepcopy(reaction_vector)
        # print("TRACE 1st Temp-vector:", temp_vector)
        # print("TRACE reaction vector:", reaction_vector)
        temp_vector[i] = 0
        # print("TRACE 2st Temp-vector:", temp_vector)
        free_lunch_status, p, status, n = free_lunch_check(N,temp_vector,constraint_rhs)
        if free_lunch_status:
            print('This index is not in the minimal', i)
            reaction_vector = copy.deepcopy(temp_vector)
        else:
            print('This index is in the miniaml set', i)
            minimal_set_indexes.append(i)
    return minimal_set_indexes, reaction_vector


def remove_unusful_reaction(N):
        m,n = getSize(N)
        removed_column_index = []
        p =[]
        ne=[]
        for i in range(n):
            # change the name to non_zero_elements
            non_zero_index = []
            column_ith = [row[i] for row in N]
            for j in column_ith:
                if j != 0:
                    non_zero_index.append(j)
            num_positive_sign = list(np.sign(non_zero_index)).count(1)
            p.append(num_positive_sign)
            num_negative_sign = list(np.sign(non_zero_index)).count(-1)
            ne.append(num_negative_sign)
            if num_positive_sign == 0 or num_negative_sign == 0:
                removed_column_index.append(i)
        #for i in range(n-1,-1,-1):
            #if i in removed_column_index:
                #for row in N:
                    #row.pop(i)
        return N, removed_column_index,p,ne


def scripts_mfl(N,irrev_reactions):
    from collections import defaultdict

    m,n = getSize(N)
    output_dict = defaultdict(list)
    indexes_output_dict = defaultdict(list)
    irrev_reactions_bd = [None]*n
    for i in irrev_reactions:
        irrev_reactions_bd[i]=0
    print(irrev_reactions_bd)
    for i in range(1,10,1):
        e = standard_basis(m,i)
        print('-------------------------------------->>>>   '+str(i)+'   <<<<------------------------------------------------------')
        print('-----------------------------------------------WEAK---------------------------------------------------')
        p, status, reactions_num, solution_non_zero_elements, solution_non_zero_indexes_w, ans_vector_x_w, ans_vector_z=\
            re_weighted_linear_program(N,e,1,None,None)
        output_dict['weak'+str(i)] = copy.deepcopy(ans_vector_x_w)
        indexes_output_dict['weak' + str(i)] = copy.deepcopy(solution_non_zero_indexes_w)
        print('---------------------------------------------------OUTPUT----------------------------------------------')
        print('solution_non_zero_indexes', solution_non_zero_indexes_w)
        print('solution_non_zero_elements',solution_non_zero_elements)
        print('-----------------------------------------------Minimal-WEAK---------------------------------------------------')
        minimal_set_indexes_w, reaction_vector_w = extract_minimal_set(N,ans_vector_x_w,e)
        output_dict['minimalWeak' + str(i)] = copy.deepcopy(reaction_vector_w)
        indexes_output_dict['minimalWeak' + str(i)] = copy.deepcopy(minimal_set_indexes_w)
        print('---------------------------------------------------OUTPUT----------------------------------------------')
        print(output_dict)
        print('-----------------------------------------------STRONG---------------------------------------------------')
        p, status, reactions_num, solution_non_zero_elements, solution_non_zero_indexes_s, ans_vector_x_s, ans_vector_z = \
            re_weighted_linear_program(N, e, 1, irrev_reactions_bd, None)
        output_dict['strong' + str(i)] = copy.deepcopy(ans_vector_x_s)
        indexes_output_dict['strong' + str(i)] = copy.deepcopy(solution_non_zero_indexes_s)
        print('---------------------------------------------------OUTPUT----------------------------------------------')
        print('solution_non_zero_indexes', solution_non_zero_indexes_s)
        print('solution_non_zero_elements',solution_non_zero_elements)
        print('-----------------------------------------------Minimal-STRONG---------------------------------------------------')
        minimal_set_indexes_s, reaction_vector_s = extract_minimal_set(N, ans_vector_x_s, e)
        output_dict['minimalStrong' + str(i)] = copy.deepcopy(reaction_vector_s)
        indexes_output_dict['minimalStrong' + str(i)] = copy.deepcopy(minimal_set_indexes_s)
        print('---------------------------------------------------OUTPUT----------------------------------------------')
        print(output_dict)
    # save dicts
    output_shelved_dict = shelve.open("/home/hzabeti/Dropbox/SFU/output.db")
    output_shelved_dict.update(output_dict)
    output_shelved_dict.close()
    ####
    indexes_output_shelved_dict = shelve.open("/home/hzabeti/Dropbox/SFU/indexesOutput.db")
    indexes_output_shelved_dict.update(indexes_output_dict)
    indexes_output_shelved_dict.close()
    return output_dict, indexes_output_dict


def scripts_mfl_after_matrix_reduction(M, irrev_reactions):
    from collections import defaultdict
    # update Matrix
    N, removed_column_index, p, ne = remove_unusful_reaction(M)
    print("TRACE removed column indexes length:", len(removed_column_index))
    # update irreversible reactions vector
    #for j in irrev_reactions:
        #if j in removed_column_index:
            #irrev_reactions.remove(j)
    # re-ordering ir irrev_reactions indexes
    #irrev_reactions_updated=copy.deepcopy(irrev_reactions)
    #for k in removed_column_index:
        #for l in irrev_reactions:
            #if k<l:
                #irrev_reactions_updated[:] = irrev_reactions_updated[:irrev_reactions.index(l)] + \
                                             #[x - 1 for x in irrev_reactions_updated[irrev_reactions.index(l):]]
                #break
    #irrev_reactions = copy.deepcopy(irrev_reactions_updated)
    m,n = getSize(N)
    print("TRACE m and n",m,n)
    output_dict = defaultdict(list)
    indexes_output_dict = defaultdict(list)
    irrev_reactions_bd = [None]*n
    deleted_reactions_bd = [None]*n
    finial_reactions_lower_bd = [None]*n
    for i in irrev_reactions:
        irrev_reactions_bd[i]=0
    print(irrev_reactions_bd)
    for i in removed_column_index:
       deleted_reactions_bd[i]=0
    print(deleted_reactions_bd)
    for i in irrev_reactions:
        finial_reactions_lower_bd[i]=0
    for i in removed_column_index:
        finial_reactions_lower_bd[i]=0
    deleted_reactions_bd[479]=0
    for i in range(m):
        e = standard_basis(m,i)
        print('-------------------------------------->>>>   '+str(i)+'   <<<<------------------------------------------------------')
        print('-----------------------------------------------WEAK---------------------------------------------------')
        p, status, reactions_num, solution_non_zero_elements, solution_non_zero_indexes_w, ans_vector_x_w, ans_vector_z=\
            re_weighted_linear_program(N,e,1,deleted_reactions_bd,deleted_reactions_bd)
        output_dict['weak'+str(i)] = copy.deepcopy(ans_vector_x_w)
        indexes_output_dict['weak' + str(i)] = copy.deepcopy(solution_non_zero_indexes_w)
        print('---------------------------------------------------OUTPUT----------------------------------------------')
        print('solution_non_zero_indexes', solution_non_zero_indexes_w)
        print('solution_non_zero_elements',solution_non_zero_elements)
        print('-----------------------------------------------Minimal-WEAK---------------------------------------------------')
        minimal_set_indexes_w, reaction_vector_w = extract_minimal_set(N,ans_vector_x_w,e)
        output_dict['minimalWeak' + str(i)] = copy.deepcopy(reaction_vector_w)
        indexes_output_dict['minimalWeak' + str(i)] = copy.deepcopy(minimal_set_indexes_w)
        print('---------------------------------------------------OUTPUT----------------------------------------------')
        print(output_dict)
        '''
        print('-----------------------------------------------STRONG---------------------------------------------------')
        p, status, reactions_num, solution_non_zero_elements, solution_non_zero_indexes_s, ans_vector_x_s, ans_vector_z = \
            re_weighted_linear_program(N, e, 1, finial_reactions_lower_bd, deleted_reactions_bd)
        output_dict['strong' + str(i)] = copy.deepcopy(ans_vector_x_s)
        indexes_output_dict['strong' + str(i)] = copy.deepcopy(solution_non_zero_indexes_s)
        print('---------------------------------------------------OUTPUT----------------------------------------------')
        print('solution_non_zero_indexes', solution_non_zero_indexes_s)
        print('solution_non_zero_elements',solution_non_zero_elements)
        print('-----------------------------------------------Minimal-STRONG---------------------------------------------------')
        minimal_set_indexes_s, reaction_vector_s= extract_minimal_set(N, ans_vector_x_s, e)
        output_dict['minimalStrong' + str(i)] = copy.deepcopy(reaction_vector_s)
        indexes_output_dict['minimalStrong' + str(i)] = copy.deepcopy(minimal_set_indexes_s)
        print('---------------------------------------------------OUTPUT----------------------------------------------')
        print(output_dict)
        '''
    # save dicts
    output_shelved_dict = shelve.open("/home/hzabeti/Dropbox/SFU/afterReductionOutput2.db")
    output_shelved_dict.update(output_dict)
    output_shelved_dict.close()
    ####
    indexes_output_shelved_dict = shelve.open("/home/hzabeti/Dropbox/SFU/afterReductionIndexesOutput2.db")
    indexes_output_shelved_dict.update(indexes_output_dict)
    indexes_output_shelved_dict.close()
    return output_dict, indexes_output_dict, removed_column_index
