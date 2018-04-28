from MinimalFreeLunch import *
import numpy as np
import copy


def all_free_lunches(N, constraint_rhs=None, objective_weight_vector=None, variable_lower_bound=None,
                     variable_upper_bound=None):
    m, n = getSize(N)
    # Set default weight vector
    if objective_weight_vector is None:
        objective_weight_vector = np.ones(n)
    if variable_lower_bound is None:
        variable_lower_bound = [None] * n
    if variable_upper_bound is None:
        variable_upper_bound = [None] * n
    if constraint_rhs is None:
        constraint_rhs = [0] * m
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
    # ------------
    p.set_objective_sense(qsoptex.ObjectiveSense.MAXIMIZE)
    p.set_param(qsoptex.Parameter.SIMPLEX_DISPLAY, 1)
    status = p.solve()
    return p, status, n


def all_FL_ans(p, status, meta_num, reactions_num):
    solution_non_zero_elements = []
    solution_non_zero_indexes = []
    y_non_zero_elements = []
    y_non_zero_indexes = []
    ans_vector_x = []
    ans_vector_y = []
    if status == qsoptex.SolutionStatus.OPTIMAL:
        print('Optimal solution')
        print(p.get_objective_value())
        for i in range(reactions_num):
            ans_vector_x.append(p.get_value('x' + str(i)))
            if p.get_value('x' + str(i)) != 0:
                solution_non_zero_elements.append(p.get_value('x' + str(i)))
                solution_non_zero_indexes.append(i)
        for j in range(meta_num):
            ans_vector_y.append(['y' + str(j), p.get_value('y' + str(j))])
            if p.get_value('y' + str(j)) != 0:
                y_non_zero_elements.append(p.get_value('y' + str(j)))
                y_non_zero_indexes.append(j)
    print('l0 approximation')
    print(len(y_non_zero_elements))
    return ans_vector_x, solution_non_zero_elements, solution_non_zero_indexes, ans_vector_y, y_non_zero_elements, y_non_zero_indexes


def all_fl_metabolites(model):
    # find sloppy reactions
    N, removed_column_index, p, ne = remove_unusful_reaction(model.fullMatrix)
    # find biomass reaction
    biomass_reaction = model.findBiomassReaction()
    # make a copy of sloppy reactions - blocked reactions are the reactions that would cause false free lunches
    # because of their representation in stoichiomatric matrix
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
    # blocked reactions
    deleted_reactions_bd = [None] * n
    for i in blocked_reactions:
        deleted_reactions_bd[i] = 0
    p, status, n = all_free_lunches(model.fullMatrix, None, None, deleted_reactions_bd, deleted_reactions_bd)
    ans_vector_x, solution_non_zero_elements, solution_non_zero_indexes, ans_vector_y, y_non_zero_elements, y_non_zero_indexes = all_FL_ans(p, status, m, n)
    return y_non_zero_indexes, ans_vector_y, ans_vector_x
