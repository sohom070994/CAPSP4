# -*- coding: utf-8 -*-
"""
Created on Sun Mar  5 16:41:05 2023

@author: hebaish
"""
#%%
import numpy as np
import math
import pandas as pd
from tqdm import tqdm
import time
from multiprocessing import Pool

#%%
def memoize(fn):
    cache = {}
    def wrapped(*args):
        if args in cache:
            return cache[args]
        else:
            result = fn(*args)
            cache[args] = result
            return result
    return wrapped

@memoize
def Pr(n,k,rho):
    if n == 1:
        P = 1
    if n >= 2:
        if k == 1:
            P = 1 - sum(Pr(n,i,rho) for i in range(2,n)) - (rho**(n-1))/((rho+1)**(n-1))
        if k >= 2 and k <= n-1:
            term = 0
            for j in range(k-1, n):
                term = term + ((1/(rho+1))**(j-k+1))*Pr(n-1,j, rho)
            P = (rho/(rho+1))*term
        if k == n:
            P = (rho**(n-1))/((rho+1)**(n-1)) 
    return P

def weight(lam,U,Node_1,Node_2):
    wt=0
    for i in range(Node_1,Node_2):
        lamm = lam[i]
        st = Node_2 - i
        sr = 1/st
        rho = lamm/sr
        for n in range(1,U+1):
            for k in range(1, n+1):
                wt = wt + Pr(n,k,rho)*k*st
            wt = wt - st
        wt = wt + (wt/U +1)*lamm
    return wt


def Bellman(cost,B): 
    global pop
    d=[[None for i in range(pop+1)] for i in range(B+1)]
    d[0]=[float("inf") for i in range(pop+1)]
    d[0][0]=0.0
    pathNode=[[None for i in range(pop+1)] for i in range(B+1)]
    for k in range(B):
        for i in range(pop+1):
            tempMin=float("inf")
            tempIndex=-1.0
            for j in range(i):
                if d[k][j]+cost[j][i]<tempMin:
                    tempMin=d[k][j]+cost[j][i]
                    tempIndex=j
            if tempMin<d[k][i]:
                d[k+1][i]=tempMin
                pathNode[k+1][i]=tempIndex
            else:
                d[k+1][i]=d[k][i]
                pathNode[k+1][i]=pathNode[k][i]
    
    path=[]
    node=pop
    # path.append(node)
    k=B
    while node>0:
        node=pathNode[k][node]
        path.append(node)
        k=k-1
        
    node_list=path
    of=0
    of+=cost[node_list[0]][slots]
    for ele in range(0,len(node_list)-1):
        of+=cost[node_list[ele+1]][node_list[ele]]
    return node_list,of

def calc_meas(cost,node_list):
    of=0
    of+=cost[node_list[0]][slots]
    for ele in range(0,len(node_list)-1):
        of+=cost[node_list[ele+1]][node_list[ele]]
    return of

def get_wts(lam,U):
    global slots
    wt=[[100000000 for i in range(slots+1)]for i in range(slots)]
    for i in tqdm(range(0,slots)):
       for j in range(i+1,slots+1):
           wt[i][j]=weight(lam,U,i,j)
    return wt

def chk_conf_mult(solution_dict):
    """creates a dictionary where each key is a time slot that exists in more than one shortest paths, and the value is the list of service types that have a conflict in that time slot.

    Parameters
    ----------
    solution_dict : dict
        keys are service types, values are nodes/time slots in the shortest path

    Returns
    -------
    dict
        a dictionary where each key is a time slot that exists in more than one shortest paths, and the value is the list of service types that have a conflict in that time slot.
    """    
    # Dictionary to store the conflicts (time slot -> list of service types)
    conflict_dict = {}

    #for all service types, remove the last element from the values list
    for key in solution_dict:
        solution_dict[key] = solution_dict[key][:-1]
    
    # Iterate over the solution_dict
    for service_type, time_slots in solution_dict.items():
        for time_slot in time_slots[0]:
            # If the time slot already exists in conflict_dict, add the service_type to the list
            if time_slot in conflict_dict:
                conflict_dict[time_slot].append(service_type)
            else:
                conflict_dict[time_slot] = [service_type]
    
    # Filter out time slots that do not have conflicts (i.e., only 1 service type)
    conflict_dict = {time_slot: services for time_slot, services in conflict_dict.items() if len(services) > 1}
    
    return conflict_dict


def chk_conf(a,b):
    # Remove last element (0) 
    a_set = set(a[:-1])
    b_set = set(b[:-1])
    if (a_set & b_set):
        return list(a_set & b_set)
    else:
        return []
    
def modifygraph(wt,conf):
    for i in range(len(wt)):
        wt[i][conf]=100000000
    return wt    

def Branch_Bound_NoBnB(wr,wc, wp, regular_slots, crisis_slots, psych_slots, solution_list):
    reg_sol=Bellman(wr,regular_slots)
    cris_sol=Bellman(wc,crisis_slots)
    psych_sol=Bellman(wp,psych_slots)
    solution_list=solution_list+[reg_sol,cris_sol,psych_sol]

    return solution_list    


def Branch_Bound_BnB(wr,wc, wp):
    wt_crisis = wc
    wt_reg = wr
    wt_psych = wp
    reg_sol=Bellman(wr,regular_slots)
    cris_sol=Bellman(wc,crisis_slots)
    psych_sol=Bellman(wp,psych_slots)
    sol_dict = {}
    sol_dict[1] = reg_sol
    sol_dict[2] = cris_sol
    sol_dict[3] = psych_sol
# Check Conflicts
    conflict_dict=chk_conf_mult(sol_dict)
# if no conflicts, add solution and OF to a global list and exit
    if not conflict_dict:
        global solution_list
        list_of_solutions = [sol_dict[service] for service in sol_dict]
        solution_list=solution_list+list_of_solutions
        print('Leaf node!')
        # print([reg_sol,cris_sol])
        return
# if conflicts, resolve conflicts, 
    else: 
        # for time_slot in conflict_dict:
        # Get the first key
        time_slot = next(iter(conflict_dict))
        for service in conflict_dict[time_slot]:
            # find list of services in conflict_dict[time_slot] except current service
            disallowed_services = [s for s in conflict_dict[time_slot] if s != service]
            # modify graph for disallowed service
            if 1 in disallowed_services:
                wt_r = modifygraph(wr,time_slot)
                Branch_Bound_BnB(wt_r,wt_crisis,wt_psych)
            if 2 in disallowed_services:
                wt_c = modifygraph(wc,time_slot)
                Branch_Bound_BnB(wt_reg,wt_c,wt_psych)
            if 3 in disallowed_services:
                wt_p = modifygraph(wp,time_slot)
                Branch_Bound_BnB(wt_reg,wt_crisis,wt_p)
                

#%%
def in_pr_2s(pr, pc, pp, ar, ac, ap, g_triage, g_crisis, g_psych, BnB):
    global slots 
    global pop
    global regular_slots
    global crisis_slots
    global psych_slots
    global wt_crisis
    global wt_reg
    global wt_psych
    global feas_c
    global feas_r
    global feas_p
    global solution_list
    global weeks_
    global slots_per_week
    U_triage = math.floor((((weeks_*slots_per_week)))*pr[0]/100)
    U_crisis = math.floor((weeks_*slots_per_week)*pc[0]/100) 
    U_psych = math.floor((weeks_*slots_per_week)*pp[0]/100) 
    weekly_lam_c = ac
    lambda_c=[]
    for i in weekly_lam_c:
        lambda_c=lambda_c+[i]*slots_per_week
    
    weekly_lam_r = ar
    lambda_r=[]
    for i in weekly_lam_r:
        lambda_r=lambda_r+[i]*slots_per_week
        
    weekly_lam_p = ap
    lambda_p=[]
    for i in weekly_lam_p:
        lambda_p=lambda_p+[i]*slots_per_week

    slots=len(lambda_r)
    pop=slots
    
    reg_props = pr
    crisis_props = pc  
    psych_props = pp
    df = []
    for rp in tqdm(reg_props):
        for cp in tqdm(crisis_props):
            for pp in tqdm(psych_props):
                regular_slots=math.ceil(rp/100*slots)
                crisis_slots=math.ceil(cp/100*slots)
                psych_slots = math.ceil(pp/100*slots)
                solution_list=[]
                if BnB == 'BnB':
                    Branch_Bound_BnB(g_triage, g_crisis,g_psych)
                    df.append(solution_list[0][0][:-1])
                    df.append(solution_list[1][0][:-1])
                    df.append(solution_list[2][0][:-1])
                if BnB == 'NoBnB':
                    sol_list = Branch_Bound_NoBnB(g_triage, g_crisis, g_psych, regular_slots, crisis_slots, psych_slots, solution_list)
                    df.append(sol_list[0][0][:-1])
                    df.append(sol_list[1][0][:-1])
                    df.append(sol_list[2][0][:-1])
    return df
#%%
ac = [0.0111111111111111,0.0118055555555556,0.0201388888888889,0.0125,
              0.0208333333333333,0.0229166666666667,0.0201388888888889,
              0.0194444444444444,0.0229166666666667,0.0131944444444444,
              0.0215277777777778,0.0166666666666667,0.0125,0.00138888888888889,
              0.000694444444444444,0.000694444444444444,0,0]

# # Regular Arrival Rate
# "C:\Users\Sohom\Documents\GitHub\CAPS\CAPS_Paper3\ARRIVAL_RATES.xlsx" SHEET REGULAR COLUMN E REGARRIVALS
ar = [0.0823295454463125,0.103177897995375,0.081846910096875,
            0.0801925233928875,0.084468831868575,0.078303409110675,
            0.0824852064396937,0.0711768749619375,0.0667458333100125,
            0.06103125001575,0.0570360668810625,0.0245395161167625,
            0.0070811831538,0.0170455645206,0.00590625,0.00039375,0.00039375,
            0.00118125]

ap = [0.00138888888888889,0.00138888888888889,0.00416666666666667,0.00347222222222222,0.00208333333333333,0.00208333333333333,0.00347222222222222,0.00416666666666667,0.00347222222222222,0.00208333333333333,0.00208333333333333,0.00138888888888889,0.00138888888888889,0.000694444444444444,0.00208333333333333,0,0,0]
ap =[ele *10 for ele in ap]
#%%
# #Condense to 6
# new_ac = []
# for i in range(0, 18, 3):
#     print(i)
#     print(sum(ac[i:i+3]))
#     new_ac.append(sum(ac[i:i+3])/3)
    
# new_ar = []
# for i in range(0, 18, 3):
#     new_ar.append(sum(ar[i:i+3])/3)


# ac = new_ac
# ar = new_ar
#%%

weeks_ = 18
slots_per_week = 15
weekly_lam_c = ac
lambda_c=[]
for i in weekly_lam_c:
    lambda_c=lambda_c+[i]*slots_per_week

weekly_lam_r = ar
lambda_r=[]
for i in weekly_lam_r:
    lambda_r=lambda_r+[i]*slots_per_week
    
weekly_lam_p = ap
lambda_p=[]
for i in weekly_lam_p:
    lambda_p=lambda_p+[i]*slots_per_week
#%%    
slots=len(lambda_r)
pop=slots

rp = 8
cp = 4
pp = 6
regular_slots=math.ceil(rp/100*slots)
crisis_slots=math.ceil(cp/100*slots)
psych_slots=math.ceil(pp/100*slots)
# #%%
# g_triage = get_wts(lambda_r,regular_slots)
# print('G_triage_generated')
# g_crisis = get_wts(lambda_c,crisis_slots)
# print('C_triage_generated')
# g_psych = get_wts(lambda_p,psych_slots)
# print('P_triage_generated')
