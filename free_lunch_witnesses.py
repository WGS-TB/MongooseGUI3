import copy
from free_lunch_analysis import *
import math

def free_lunch_witness_finder(model_obj,y_lower_bound = None,max_num_atoms=1e3):
    import cplex
    cpx_feasibility_tol = 1e-9
    excluded_reactions_index = find_biomass_reaction_candidates(model_obj)+find_import_export_reactions(model_obj.fullMatrix)
    m,r_num = getSize(model_obj.fullMatrix)
    n = r_num - len(set(excluded_reactions_index))
    if y_lower_bound == None:
        #y_lower_bound = 1e-4
        norm_list = []
        Max_norm = 0
        m_sq = int(math.sqrt(m))+1
        S = model_obj.fullMatrix
        for i in range(r_num):
            S_c = [float(c[i]) for c in S]
            current_column = np.linalg.norm(S_c)
            norm_list.append(current_column)
        Max_norm = int(max(norm_list))+1
        y_lower_bound = 1/(Max_norm*m_sq*max_num_atoms)
        if y_lower_bound <= cpx_feasibility_tol:
            y_lower_bound = 1/(1/cpx_feasibility_tol-1)
# -------------------------- Problem initials ----------------------------
    p_obj = [1.0]*n+[0.0]*m
    p_rhs = [0]*(2*n)
    p_lb =[0]*n+[y_lower_bound]*m
    p_ub = [1]*n+[cplex.infinity]*m
    p_ctype ="B"*n+"C"*m
    p_sense = "G"*(2*n)
    p_colnames =["x"+str(i) for i in range(r_num) if i not in excluded_reactions_index ]+["y"+str(j) for j in range(m)] 
    p_rownames =["r"+str(i) for i in range(2*n)]
    rows = []
    prob = cplex.Cplex()
    # --------------------------- Param tolerances -----------------------
    prob.parameters.mip.tolerances.integrality.set(1e-15)
    prob.parameters.simplex.tolerances.feasibility.set(cpx_feasibility_tol)
    print(prob.parameters.mip.tolerances.integrality)
    prob.objective.set_sense(prob.objective.sense.minimize)
    prob.variables.add(obj = p_obj, lb = p_lb, ub = p_ub, types = p_ctype,names = p_colnames)
    
    for i in range(r_num):
        if i not in excluded_reactions_index:
            r = model_obj.reactions[i]
            row_inner_f = ["x"+str(i)]+["y"+str(j[0]) for j in r.pairs ]
            row_inner_s = [1]+[float(j[1]) for j in r.pairs]
            rows.append(copy.copy([row_inner_f,row_inner_s]))
    for i in range(r_num):
        if i not in excluded_reactions_index:
            r = model_obj.reactions[i]
            row_inner_f = ["x"+str(i)]+["y"+str(j[0]) for j in r.pairs ]
            row_inner_s = [1]+[-float(j[1]) for j in r.pairs]
            rows.append(copy.copy([row_inner_f,row_inner_s]))
        
    prob.linear_constraints.add(lin_expr = rows, senses = p_sense,rhs = p_rhs, names = p_rownames)

    #prob.set_log_stream(None)
    #prob.set_error_stream(None)
    #prob.set_warning_stream(None)
    #prob.set_results_stream(None)
    prob.solve()
    flw_objective_value = prob.solution.get_objective_value()
    #print([(v[1:],prob.solution.get_values(v))  for v in prob.variables.get_names() if v[0] == 'x' and prob.solution.get_values(v)!=0])
    flw_list = [int(v[1:]) for v in prob.variables.get_names() if v[0] == 'x' and prob.solution.get_values(v)!=0]
    #------------------------------ Print----------------------------------------
    print('Number of free lunch witnesses:',len(flw_list))
    print('List of free lunch witnesses:',flw_list)
    print('Computed y lower bound', y_lower_bound)
    #----------------------------------------------------------------------------
    return flw_list,y_lower_bound
