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

def Branch_Bound_NoBnB(wr,wc, regular_slots, crisis_slots, solution_list):
    reg_sol=Bellman(wr,regular_slots)
    cris_sol=Bellman(wc,crisis_slots)
    solution_list=solution_list+[reg_sol,cris_sol]

    return solution_list    

def Branch_Bound_NoBnB_seq(wr,wc, regular_slots, crisis_slots, solution_list):
    reg_sol=Bellman(wr,regular_slots)
    wt_cris_mod = wc.copy()
    for conflict_slot in reg_sol[0][:-1]:
        wt_cris_mod = modifygraph(wt_cris_mod,conflict_slot)
    cris_sol=Bellman(wt_cris_mod,crisis_slots)
    solution_list=solution_list+[reg_sol,cris_sol]

    return solution_list 


def Branch_Bound_BnB(wr,wc):
    wt_crisis = wc
    wt_reg = wr
    reg_sol=Bellman(wr,regular_slots)
    cris_sol=Bellman(wc,crisis_slots)
# Check Conflicts
    conflict_list=chk_conf(reg_sol[0],cris_sol[0])
# if no conflicts, add solution and OF to a global list and exit
    if conflict_list==[]:
        global solution_list
        solution_list=solution_list+[reg_sol,cris_sol]
        # print([reg_sol,cris_sol])
        return
# if conflicts, resolve conflicts, 
    else: 
# solve P1 Disallow 1st regular conflict
        conf=conflict_list[0]
        
        wt_reg_up=modifygraph(wr,conf)
        Branch_Bound_BnB(wt_reg_up,wt_crisis)

# solve P2 Disallow 1st crisis conflict
        wt_crisis_up=modifygraph(wc,conf)
        Branch_Bound_BnB(wt_reg,wt_crisis_up)

#%%
def in_pr_2s(pr, pc, ar, ac, g_triage, g_crisis, BnB):
    global slots 
    global pop
    global regular_slots
    global crisis_slots
    global wt_crisis
    global wt_reg
    global feas_c
    global feas_r
    global solution_list
    global weeks_
    global slots_per_week
    U_triage = math.floor((((weeks_*slots_per_week)))*pr[0]/100)
    U_crisis = math.floor((weeks_*slots_per_week)*pc[0]/100) 
    weekly_lam_c = ac
    lambda_c=[]
    for i in weekly_lam_c:
        lambda_c=lambda_c+[i]*slots_per_week
    
    weekly_lam_r = ar
    lambda_r=[]
    for i in weekly_lam_r:
        lambda_r=lambda_r+[i]*slots_per_week

    slots=len(lambda_r)
    pop=slots
    
    reg_props = pr
    crisis_props = pc  
    df = []
    for rp in tqdm(reg_props):
        for cp in tqdm(crisis_props):
            regular_slots=math.ceil(rp/100*slots)
            crisis_slots=math.ceil(cp/100*slots)
            solution_list=[]
            if BnB == 'BnB':
                Branch_Bound_BnB(g_triage, g_crisis)
                df.append(solution_list[0][0][:-1])
                df.append(solution_list[1][0][:-1])
            if BnB == 'NoBnB':
                sol_list = Branch_Bound_NoBnB(g_triage, g_crisis, regular_slots, crisis_slots, solution_list)
                df.append(sol_list[0][0][:-1])
                df.append(sol_list[1][0][:-1])
            if BnB == 'seq':
                sol_list = Branch_Bound_NoBnB_seq(g_triage, g_crisis, regular_slots, crisis_slots, solution_list)
                df.append(sol_list[0][0][:-1])
                df.append(sol_list[1][0][:-1])
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
slots_per_week = 10
weekly_lam_c = ac
lambda_c=[]
for i in weekly_lam_c:
    lambda_c=lambda_c+[i]*slots_per_week

weekly_lam_r = ar
lambda_r=[]
for i in weekly_lam_r:
    lambda_r=lambda_r+[i]*slots_per_week
    
slots=len(lambda_r)
pop=slots

rp = 8
cp = 4
regular_slots=math.ceil(rp/100*slots)
crisis_slots=math.ceil(cp/100*slots)
# #%%
# g_triage = get_wts(lambda_r,regular_slots)
# print('G_triage_generated')
# g_crisis = get_wts(lambda_c,crisis_slots)
