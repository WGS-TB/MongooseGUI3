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


def check_emlemental_balance_formula(model_obj, reaction_num, metabolits_data):
    r = model_obj.reactions[reaction_num]
    r_pairs = r.pairs
    r_formula_list_rhs = []
    r_formula_list_lhs = []
    r_formula_dict_list=[]
    for i in r_pairs:
        cur_met = model_obj.metabolites[i[0]]
        cur_met_name = cur_met.species.name
        for j in metabolits_data:
            if cur_met_name[2:]== j['id']:
                if i[1]>0:
                    r_formula_list_rhs.append((i[1],j['formula'],j['id']))
                elif i[1]<0:
                    r_formula_list_lhs.append((i[1],j['formula'],j['id']))
                f =convertFormula(j['formula'])
                for k in f.keys():
                    f[k] *=i[1]
                r_formula_dict_list.append(f)
    return r_formula_dict_list,r_formula_list_lhs,r_formula_list_rhs


def final_check_EB(r_formula_dict_list):
    chem_set= set()
    final_check = {}
    for l in r_formula_dict_list:
        for k in l.keys():
            chem_set.add(k)
    for i in chem_set:
        total_temp = 0
        for j in r_formula_dict_list:
            if i in j.keys():
                total_temp +=j[i]
        final_check.update({i:total_temp})
    return final_check,chem_set


def reactions_correspond_to_metabolite(model,meta_id_list):
    reations_list=[]
    meta_index = []
    for i in model.metabolites:
        if i.species.name[2:] in meta_id_list:
            meta_index.append(i.index)
    for m_index in meta_index:
        row = model.fullMatrix[m_index]
        for j in range(len(row)):
            if row[j] != 0:
                reations_list.append(j)
    return reations_list,meta_index
    
    

def script_code(model,blocked_reactions,n):
    final_list_dict = {}
    for r_num in range(n):
        if r_num not in blocked_reactions:
            r = model.reactions[r_num]
            r_name = r.name
            r_formula_dict_list,r_formula_list_lhs,r_formula_list_rhs = check_emlemental_balance_formula(model,r_num,metabolits_data)
            final_check,chem_set = final_check_EB(r_formula_dict_list)
            final_list_dict.update({r_num:{'name':r.name,'balance':final_check,'reaction_rhs':r_formula_list_rhs, 'reaction_lhs': r_formula_list_lhs}})
    fl = set()
    fl_dict = {}
    for i in final_list_dict:
        test = final_list_dict[i]
        test2= test['balance']
        for j in test2.keys():
            if test2[j]!=0:
                fl.add(i)
                fl_dict.update({i:final_list_dict[i]})
    return fl, fl_dict

output_dict = defaultdict(list)
f = open('/home/hzabeti/Dropbox/SFU/EB_check_text.txt','w')
f.write('Elementally balanced check with chemical formula text version\n')
for modelNameWext in os.listdir('/home/hzabeti/qsopt-ex/build/python-qsoptex/MongooseGUI3/Models/DB/'):
    if modelNameWext.endswith(".db"):
        # model's name with out ext.
        modelName = os.path.splitext(modelNameWext)[0]
        f.write('\n'+'Model Name:\n')
        f.write(modelName+'\n')
        print("Model Name: ",modelName)
        s = shelve.open('/home/hzabeti/qsopt-ex/build/python-qsoptex/MongooseGUI3/Models/DB/'+modelNameWext)
        model= s[modelName]
        #find sloppy reactions
        N, removed_column_index,p,ne = remove_unusful_reaction(model.fullMatrix)
        #find biomass reaction
        biomass_reaction = model.findBiomassReaction()
        # make a copy of sloppy reactions - blocked reactions are the reactions that would cause false free lunches because of their representation in stoichomatric matrix
        blocked_reactions = copy.copy(removed_column_index)
        # add biomass reaction to blocked reactions - we are not intrested in biomass reaction for now
        blocked_reactions.append(biomass_reaction)
        # find additional biomass reactions
        biomass_reaction_candidates = []
        for r in model.reactions:
            reaction_name_biomass = r.name.lower()
            if reaction_name_biomass.find('biomass')!=-1:
                biomass_reaction_candidates.append({'name':r.name,'index':r.index})
                blocked_reactions.append(r.index)
        f = open('/home/hzabeti/Dropbox/SFU/EB_check_text.txt','a')
        f.write('Biomass reaction index:'+str(biomass_reaction)+'\n')
        f.write('All biomass reaction candidates:'+str(biomass_reaction_candidates)+'\n')
        with open('/home/hzabeti/qsopt-ex/build/python-qsoptex/MongooseGUI3/Models/json/'+modelName+'.json', 'r') as f:
            model_jason = json.load(f)
            metabolits_data = model_jason['metabolites']
        # size of matrix
        m,n = getSize(model.fullMatrix)
        # metabolites with out formula
        meta_without_formula = []
        
        for i in range(m):
            try:
                meta_formula = metabolits_data[i]['formula']
            except KeyError:
                meta_without_formula.append(metabolits_data[i]['id'])
        print("metabolites without formula: ",meta_without_formula)
        if len(meta_without_formula)==0:
            f = open('/home/hzabeti/Dropbox/SFU/EB_check_text.txt','a')
            f.write('Number of metabolites without formula:'+str(len(meta_without_formula))+'\n')
            fl,fl_dict  = script_code(model,blocked_reactions,n)
            f.write('Number of imbalance reactions:'+str(len(fl))+'\n')
            '''
            f.write('+++++++++++++++++ Details ++++++++++++++++++++\n\n')
            f.write('metabolites without formula:\n None \n')
            for i in fl:
                f.write('Imbalanced reactions:\n'+'Index:\n'+str(i)+'\n'+'Name:\n'+model.reactions[i].name+'\n')
                cur_reaction=fl_dict[i]
                f.write('Balance:\n'+str(cur_reaction['balance'])+'\n\n')
                f.write('Reaction LHS:\n'+str(cur_reaction['reaction_lhs'])+'\n\n')
                f.write('Reaction RHS:\n'+str(cur_reaction['reaction_rhs'])+'\n\n')
                '''
            print("Imbalance reactions",fl)
            if len(fl) == len(fl_dict):
                output_dict[modelName] = {'metas without formula':None,'correspond reactions':None,'imbalanced reactions':copy.copy(fl_dict)}
            else:
                print('Problem fl does not agree with fl_dict',modelName)
        elif len(meta_without_formula)!=m and len(meta_without_formula)>0:
            f = open('/home/hzabeti/Dropbox/SFU/EB_check_text.txt','a')
            f.write('Number of metabolites without formula:'+str(len(meta_without_formula))+'\n')
            reations_list,meta_index = reactions_correspond_to_metabolite(model,meta_without_formula)
            new_blocked_reactions =blocked_reactions+reations_list 
            fl,fl_dict  = script_code(model,new_blocked_reactions,n)
            f.write('Number of imbalance reactions:'+str(len(fl))+'\n')
            '''
            f.write('+++++++++++++++++ Details ++++++++++++++++++++\n\n')
            f.write('metabolites without formula:\n'+str(meta_without_formula)+'\n')
            f.write('correspond reactions index:\n'+str(reations_list)+'\n\n\n')
            for i in fl:
                f.write('Imbalanced reactions:\n'+'Index:\n'+str(i)+'\n'+'Name:\n'+model.reactions[i].name+'\n')
                cur_reaction=fl_dict[i]
                f.write('Balance:\n'+str(cur_reaction['balance'])+'\n\n')
                f.write('Reaction LHS:\n'+str(cur_reaction['reaction_lhs'])+'\n\n')
                f.write('Reaction RHS:\n'+str(cur_reaction['reaction_rhs'])+'\n\n')
                '''
            print("Imbalance reactions",fl)
            if len(fl) == len(fl_dict):
                output_dict[modelName] = {'metas without formula':copy.copy(meta_without_formula),'correspond reactions':copy.copy(reations_list),'imbalanced reactions':copy.copy(fl_dict)}
            else:
                print('Problem fl does not agree with fl_dict',modelName)
        else:
            f = open('/home/hzabeti/Dropbox/SFU/EB_check_text.txt','a')
            f.write('Number of metabolites without formula:'+str(len(meta_without_formula))+'\n')
            f.write('Number of imbalance reactions: Not available\n')
            '''
            f.write('+++++++++++++++++ Details ++++++++++++++++++++\n\n')
            f.write('metabolites without formula:\n ALL')            
            output_dict[modelName] = {'metas without formula':'ALL','correspond reactions':'ALL','imbalanced reactions':None}
            '''
        f = open('/home/hzabeti/Dropbox/SFU/EB_check_text.txt','a')    
        f.write('\n\n\n-------------------------------------------------------------\n\n\n')
f.close()
output_shelved_dict = shelve.open("/home/hzabeti/Dropbox/SFU/EB_check_output.db")
output_shelved_dict.update(output_dict)
output_shelved_dict.close()
