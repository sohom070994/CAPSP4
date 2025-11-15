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
    for i in range(0,slots):
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
    U_triage = math.floor(720*pr[0]/100)
    U_crisis = math.floor(720*pc[0]/100) 
    global slots 
    global pop
    global regular_slots
    global crisis_slots
    global wt_crisis
    global wt_reg
    global feas_c
    global feas_r
    global solution_list
    weekly_lam_c = ac
    lambda_c=[]
    for i in weekly_lam_c:
        lambda_c=lambda_c+[i]*40
    
    weekly_lam_r = ar
    lambda_r=[]
    for i in weekly_lam_r:
        lambda_r=lambda_r+[i]*40

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
    return df


