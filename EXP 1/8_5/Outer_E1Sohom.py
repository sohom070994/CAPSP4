# -*- coding: utf-8 -*-
"""
Created on Sun Mar  5 16:58:58 2023

@author: hebaish
"""
#%%
import numpy as np
import math
#import SPfuncs
from InnerProblem_create_graph import in_pr_2s
import pandas as pd
import os
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


def wtcalc_old(a,s,n1,n2,U,log_):
    global weeks_
    global slots_per_week
    if log_==1:
        print("n1="+str(n1)+"; n2="+str(n2))
    arr = a
    
    slots = s
    
    
    #to get the number of slots for the edge between n1 and n2
    n = sum([1 if ((x>=(n1-1)*slots_per_week)&(x<(n2-1)*slots_per_week)) else 0 for x in slots])
    if log_==1:
        print("n="+str(n))
    
    #to get the span of n1 and n2 (total number of slots)
    N = (n2-n1)*slots_per_week
    
    n_2=n2*slots_per_week
    n_1=n1*slots_per_week
    
    p_s = n/N
    if log_==1:
        print("p_s="+str(p_s))
    
    if n!=0:
        lam = [arr[i] for i in range(n1-1, n2-1) for j in range(slots_per_week)]
        H = [0] * N
        # for i in range(N):
        #     try:
        #         H[i] = (((1-(1-p_s)**(n_2-i-1))/p_s)-((n_2-i-1)*(1-p_s)**(n_2-i-1))+((n_2-i)*(p_s)*(1-p_s)**(n_2-i-1)))/(1-(1-p_s)**(n_2-i))
        #     except ZeroDivisionError:
        #         print("p_s: " + str(p_s)+ "; N =" + str(N)+"; i =" + str(i))
        xbar = [0] * N
        xbar[-1] = 1
        if n > 2:
            freq = math.floor((N-1)/(n-1))
            for i in range(1, n):
                idx = -1 - i*freq
                xbar[idx] = 1
                
            for i in range(N):
                count  = 1
                j = i  
                while xbar[j] != 1:
                    count += 1
                    j += 1
                H[i] = count
        elif n == 2:
            freq = math.floor((N-1)/(n))
            ind = math.floor(N/2)
            xbar[ind] = 1
            for i in range(N):
                count  = 1
                j = i  
                while xbar[j] != 1:
                    count += 1
                    j += 1
                H[i] = count
                    
        elif n == 1:
            for i in range(N):
                count  = 1
                j = i  
                while xbar[j] != 1:
                    count += 1
                    j += 1
                H[i] = count
        try:
            mu = [1/H[i] for i in range(len(H))]
        except ZeroDivisionError:
            print ("ERRORRR")
        if log_==1:
            print("H_old="+str(H))   
        
        wt=0
        for i in range(N):
            rho = lam[i]/mu[i]
            for n in range(1,U+1):
                for k in range(1, n+1):
                    wt = wt + Pr(n,k,rho)*k/mu[i]
                wt = wt - (1/mu[i])
            # wt = wt/U + 1
            wt = wt + (wt/U +1)*lam[i]
        
    else: 
        wt = float('inf')
    
    return wt

def wtcalc(a,s,n1,n2,U,log_):
    global weeks_
    global slots_per_week
    if log_==1:
        print("n1="+str(n1)+"; n2="+str(n2))
    
    arr = a
    
    slots = s
    
    
    #to get the number of slots for the edge between n1 and n2
    n = sum([1 if ((x>=(n1-1)*slots_per_week)&(x<(n2-1)*slots_per_week)) else 0 for x in slots])
    if log_==1:
        print("n="+str(n))
    
    #to get the span of n1 and n2 (total number of slots)
    N = (n2-n1)*slots_per_week
    
    n_2=n2*slots_per_week
    n_1=n1*slots_per_week
    
    p_s = n/N
    
    if log_==1:
        print("p_s="+str(p_s))
    
    if n!=0:
        lam = [arr[i] for i in range(n1-1, n2-1) for j in range(slots_per_week)]
        H = [0] * N
        for i in range(N):
            try:
                H[i] = (((1-(1-p_s)**(n_2-i-1))/p_s)-((n_2-i-1)*(1-p_s)**(n_2-i-1))+((n_2-i)*(p_s)*(1-p_s)**(n_2-i-1)))/(1-(1-p_s)**(n_2-i))
            except ZeroDivisionError:
                print("p_s: " + str(p_s)+ "; N =" + str(N)+"; i =" + str(i))
        # xbar = [0] * N
        # xbar[-1] = 1
        # if n > 2:
        #     freq = math.floor((N-1)/(n-1))
        #     for i in range(1, n):
        #         idx = -1 - i*freq
        #         xbar[idx] = 1
                
        #     for i in range(N):
        #         count  = 1
        #         j = i  
        #         while xbar[j] != 1:
        #             count += 1
        #             j += 1
        #         H[i] = count
        # elif n == 2:
        #     freq = math.floor((N-1)/(n))
        #     ind = math.floor(N/2)
        #     xbar[ind] = 1
        #     for i in range(N):
        #         count  = 1
        #         j = i  
        #         while xbar[j] != 1:
        #             count += 1
        #             j += 1
        #         H[i] = count
                    
        # elif n == 1:
        #     for i in range(N):
        #         count  = 1
        #         j = i  
        #         while xbar[j] != 1:
        #             count += 1
        #             j += 1
        #         H[i] = count
        try:
            mu = [1/H[i] for i in range(len(H))]
        except ZeroDivisionError:
            print ("ERRORRR")
        
        if log_==1:
            print("H_new="+str(H))    
        wt=0
        for i in range(N):
            rho = lam[i]/mu[i]
            for n in range(1,U+1):
                for k in range(1, n+1):
                    wt = wt + Pr(n,k,rho)*k/mu[i]
                wt = wt - (1/mu[i])
            # wt = wt/U + 1
            wt = wt + (wt/U +1)*lam[i]
        
    else: 
        wt = float('inf')
    
    return wt


def weightscs(n1, n2):
    global weeks_
    global slots_per_week
    global wtr, wtc, ar, sr, ac, sc, U_triage, U_crisis
    log_=0
    
    wtd=(wtcalc(ar,sr,n1,n2,U_triage,log_)/wtr)+(wtcalc(ac,sc,n1,n2,U_crisis,log_)/wtc)
    
    log_ = 1 
    check_NewT = wtcalc(ar,sr,n1,n2,U_triage,log_)
    check_NewC = wtcalc(ac,sc,n1,n2,U_crisis,log_)
    check_OldT = wtcalc_old(ar,sr,n1,n2,U_triage,log_)
    Check_OldC = wtcalc_old(ac,sc,n1,n2,U_crisis,log_)
    log_ = 0
    
    return wtd
    
def weightsrs(a, s, n1, n2, U):
    
    log_ = 0
    wt=wtcalc(a, s, n1, n2, U,log_)
    
    # log_ = 1 
    # check_NewT = wtcalc(ar,sr,n1,n2,U_triage,log_)
    # check_NewC = wtcalc(ac,sc,n1,n2,U_crisis,log_)
    # check_OldT = wtcalc_old(ar,sr,n1,n2,U_triage,log_)
    # Check_OldC = wtcalc_old(ac,sc,n1,n2,U_crisis,log_)
    # log_ = 0
    
    return wt  

'''
This function takes the nodes and the number of changes and creates a directed
graph.
'''
def create_graph(nodes, num_changes, ns, s, ar, sr, ac, sc):
    global U_triage, U_crisis
    #Check for a valid value for the number of changes
    if num_changes<1 or num_changes>(len(nodes)-1):
        raise ValueError('Please enter a value greater than 1 less than ' 
                         + str(len(nodes)-2))
        
    '''
    This loop creates the edges from each node to the following nodes in a 
    directed manner.
    '''
    '''
    To force the number of changes, and consequently the edges between nodes, 
    add the following to len(nodes) instead of the (-1) (-num_changes+k) in the
    for loop of the j below.
    '''

    edges = []
    for i in range(1,len(nodes)):
        for j in range(i+1,(len(nodes)+1)):
            if ns == 1:
                if s ==1:
                    edges.append([i, j, weightsrs(ar, sr, i, j, U_triage)])
                if s ==2:
                    edges.append([i, j, weightsrs(ac, sc, i, j, U_crisis)])
            if ns == 2:
                edges.append([i, j, weightscs(i,j)])
                    
    
    
    print('Graph created')        
            
    '''        
    This loop creates the graph matrix using the edges and the weights of each
    edge.
    '''    
    edges = np.array(edges)       
    ind=[]
    graph = np.array([[float('inf')]*len(nodes)]*len(nodes))

    for l in range(1,len(nodes)+1):
        for m in range(1,len(nodes)+1):
            if l == m:
                graph[l-1,m-1] = 0
            else:
                ind = np.array(np.where((edges[:,0]==l)
                                       &(edges[:,1]==m))).tolist()[0]
                if len(ind)>0:
                    graph[l-1,m-1] = edges[ind[0],2]				
									       

    return edges, graph


def shortestPath(graph, u, v, k):
    
    global V

    sp = [[None] * V for i in range(V)]
    path = [[None] * V for i in range(V)]

    for i in range(V):
        for j in range(V):
            sp[i][j] = [None] * (k + 1)
            path[i][j] = [None] * (k + 1)

    for e in range(k + 1):
        for i in range(V):
            for j in range(V):
                sp[i][j][e] = float('inf')
                path[i][j][e] = float('inf')

                if (e == 0 and i == j):
                    sp[i][j][e] = 0
                if (e == 1 and graph[i][j] != float('inf')):
                    sp[i][j][e] = graph[i][j]
                    path[i][j][e] = []

                if (e > 1):
                    for a in range(V):
                        if (graph[i][a] != float('inf') and i != a and
							j!= a and sp[a][j][e - 1] != float('inf')):
                            if graph[i][a] + sp[a][j][e - 1] < sp[i][j][e]:
                                sp[i][j][e] = graph[i][a] + sp[a][j][e - 1]
                                path[i][j][e] = path[a][j][e-1] + [a]
 
    PathsLen = []
    Paths = []
    for i in range(len(sp[u][v])):
        if sp[u][v][i] < float('inf'):
            PathsLen.append(sp[u][v][i])
            path[u][v][i] = [x+1 for x in path[u][v][i]]
            path[u][v][i] = [u+1] + path[u][v][i][::-1] + [v+1]
            Paths.append(path[u][v][i])

    kk = sp[u][v].index(min(sp[u][v]))
    
    path[u][v][kk] = [x+1 for x in path[u][v][kk]]
    path[u][v][kk] = [u+1] + path[u][v][kk][::-1] + [v+1]
    return Paths, PathsLen

#%%
slice_ = weeks_

nodes = list(range(1,slice_+2))            #Create the nodes to be used in the graph
num_changes =  slice_                  #The number of desired changes
u = 1                                #source node (usually node 1)
v = slice_+1                              #sink node 
V = len(nodes)                       #The number of vertices in the graph


#%%
#arrival rates already loaded
# ac = [0.0111111111111111,0.0118055555555556,0.0201388888888889,0.0125,
#               0.0208333333333333,0.0229166666666667,0.0201388888888889,
#               0.0194444444444444,0.0229166666666667,0.0131944444444444,
#               0.0215277777777778,0.0166666666666667,0.0125,0.00138888888888889,
#               0.000694444444444444,0.000694444444444444,0,0]

# # # Regular Arrival Rate
# # "C:\Users\Sohom\Documents\GitHub\CAPS\CAPS_Paper3\ARRIVAL_RATES.xlsx" SHEET REGULAR COLUMN E REGARRIVALS
# ar = [0.0823295454463125,0.103177897995375,0.081846910096875,
#             0.0801925233928875,0.084468831868575,0.078303409110675,
#             0.0824852064396937,0.0711768749619375,0.0667458333100125,
#             0.06103125001575,0.0570360668810625,0.0245395161167625,
#             0.0070811831538,0.0170455645206,0.00590625,0.00039375,0.00039375,
#             0.00118125]

# ac = ac[:slice_]
# ar = ar[:slice_]
#%%
pr = [9]
pc = [5]


U_triage = math.floor(slots*pr[0]/100)
U_crisis = math.floor(slots*pc[0]/100)

t_graph = g_triage
c_graph = g_crisis

s = in_pr_2s(pr, pc, ar, ac, t_graph, c_graph, "NoBnB")
sr = s[0]
sc = s[1]

proptbl = np.zeros((slice_,3))
propsr = []
propsc = []


#Solvig for each service
print('FOR REGULAR:')
edges, graph = create_graph(nodes, num_changes, 1, 1, ar, sr, ac, sc)
paths, lengths= shortestPath(graph, u-1, v-1,num_changes+1)
for i in range(len(paths)):
    print('A path from node ' + str(u)
          + ' to node ' + str(v) + ' Exist, having the following nodes: ')
    print(*paths[i], sep = ' --> ') 
    print('The length of the above path from node ' + str(u)
          + ' to node ' + str(v) + ' is ' + str(lengths[i]))
    
print('The shortest path from node ' + str(u)
      + ' to node ' + str(v) + ' has the following nodes: ')
print(*paths[lengths.index(min(lengths))], sep = ' --> ') 
print('The length of the shortest path from node ' + str(u)
      + ' to node ' + str(v) + ' is ' + str(lengths[lengths.index(min(lengths))]))

shrtspthr = paths[lengths.index(min(lengths))]
count = 0
for i in range(0,len(shrtspthr)-1):
    for j in range(shrtspthr[i+1]-shrtspthr[i]):
        propsr.append(sum((x < (shrtspthr[i+1]-1)*slots_per_week and x >= (shrtspthr[i]-1)*slots_per_week) for x in sr)/(slots_per_week*(shrtspthr[i+1]-shrtspthr[i])))
        proptbl[count][0] = count+1
        proptbl[count][1] = propsr[count]
        count += 1
        
        
print('FOR CRISIS:')
edges, graph = create_graph(nodes, num_changes,1, 2, ar, sr, ac, sc)
paths, lengths= shortestPath(graph, u-1, v-1,num_changes+1)
for i in range(len(paths)):
    print('A path from node ' + str(u)
          + ' to node ' + str(v) + ' Exist, having the following nodes: ')
    print(*paths[i], sep = ' --> ') 
    print('The length of the above path from node ' + str(u)
          + ' to node ' + str(v) + ' is ' + str(lengths[i]))
    
print('The shortest path from node ' + str(u)
      + ' to node ' + str(v) + ' has the following nodes: ')
print(*paths[lengths.index(min(lengths))], sep = ' --> ') 
print('The length of the shortest path from node ' + str(u)
      + ' to node ' + str(v) + ' is ' + str(lengths[lengths.index(min(lengths))]))

shrtspthc = paths[lengths.index(min(lengths))]
count = 0
for i in range(0,len(shrtspthc)-1):
    for j in range(shrtspthc[i+1]-shrtspthc[i]):
        propsc.append(sum((x < (shrtspthc[i+1]-1)*slots_per_week and x >= (shrtspthc[i]-1)*slots_per_week) for x in sc)/(slots_per_week*(shrtspthc[i+1]-shrtspthc[i])))
        proptbl[count][0] = count+1
        proptbl[count][2] = propsc[count]
        count += 1
        
        

#%%
policy = pd.DataFrame(proptbl, columns = ['Week', 'Triage', 'Crisis'])

current_path = os.path.dirname(os.path.abspath(__file__))
fn = f'Policy_E1_Sohom_{slice_}_{pr[0]}_{pc[0]}.csv'
current_fn = os.path.join(current_path,fn)
policy.to_csv(current_fn, index = False)

