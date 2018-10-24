import qsoptex
import numpy as np
from Utilities import *
from fractions import Fraction
import copy
import shelve

def conservation_relations(N, constraint_rhs, blocked_reactions, objective_weight_vector=None, variable_lower_bound=None, variable_upper_bound=None):
    m, n = getSize(N)
    # alpha = Fraction(1, 10)
    # following part set a default weight vector
    if objective_weight_vector is None:
        objective_weight_vector = np.zeros(m)
    if variable_lower_bound is None:
        variable_lower_bound = [None]*m
    if variable_upper_bound is None:
        variable_upper_bound = [None]*m
    p = qsoptex.ExactProblem()
    # In the following part we define the objective of the linear programming
    for i in range(m):
        p.add_variable(name='y' + str(i), objective=objective_weight_vector[i], lower=variable_lower_bound[i], upper=variable_upper_bound[i])
    # In the following part we define constrains of the linear Programming
    for j in range(n):
        if j not in blocked_reactions:
            p.add_linear_constraint(qsoptex.ConstraintSense.EQUAL,{'y' + str(i): N[i][j] for i in range(m)}, rhs=constraint_rhs[j])
    # Following part would solve the LP with qsoptex
    p.set_objective_sense(qsoptex.ObjectiveSense.MINIMIZE)
    p.set_param(qsoptex.Parameter.SIMPLEX_DISPLAY, 1)
    status = p.solve()
    return p, status, m


def conservation_relations_ans(p, status, meta_num):
    ans_vector_y = []
    ans_dict_y = {}
    if status == qsoptex.SolutionStatus.OPTIMAL:
        print('Optimal solution')
        print(p.get_objective_value())
        for i in range(meta_num):
            ans_vector_y.append(p.get_value('y' + str(i)))
            ans_dict_y.update({i:p.get_value('y' + str(i))})
    return ans_vector_y,ans_dict_y
