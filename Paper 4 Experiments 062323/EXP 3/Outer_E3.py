# -*- coding: utf-8 -*-
"""
Created on Sun Mar  5 16:58:58 2023

@author: hebaish
"""
#%%
import numpy as np
import math
#import SPfuncs
from InnerProblem import in_pr_2s
import pandas as pd
from matplotlib import pyplot as plt
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

def wtcalc(a,s,n1,n2,U):
    arr = a
    
    slots = s
    
    
    #to get the number of slots for the edge between n1 and n2
    n = sum([1 if ((x>=(n1-1)*40)&(x<(n2-1)*40)) else 0 for x in slots])
    
    #to get the span of n1 and n2 (total number of slots)
    N = (n2-n1)*40
    
    if n!=0:
        lam = [arr[i] for i in range(n1-1, n2-1) for j in range(40)]
        H = [0] * N
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
    
        mu = [1/H[i] for i in range(len(H))]
        
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
    global wtr, wtc, ar, s2r, ac, s2c, U_triage, U_crisis
    
    wtd=(wtcalc(ar,s2r,n1,n2,U_triage)/wtr)+(wtcalc(ac,s2c,n1,n2,U_crisis)/wtc)
    
    return wtd
    
def weightsrs(a, s, n1, n2, U):

    wt=wtcalc(a, s, n1, n2, U)
    
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
nodes = list(range(1,20))            #Create the nodes to be used in the graph
num_changes = 18                     #The number of desired changes
u = 1                                #source node (usually node 1)
v = 19                               #sink node 
V = len(nodes)                       #The number of vertices in the graph

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

pr = [8]
pc = [4]

U_triage = math.floor(720*pr[0]/100)
U_crisis = math.floor(720*pc[0]/100)

c_graph = crisis_graph
t_graph = triage_graph

s1 = in_pr_2s(pr, pc, ar, ac, t_graph, c_graph, "NoBnB")
s1r = s1[0]
s1c = s1[1]

s2 = in_pr_2s(pr, pc, ar, ac, t_graph, c_graph, "BnB")
s2r = s2[0]
s2c = s2[1]


#getting the weights for service 1 and 2
#ns is the number of services, and s is which service to solve for.
#if ns is 1 then s has to be specified. If ns is 2, then s doesn't matter.
edges, graph = create_graph(nodes, num_changes, 1, 1, ar, s1r, ac, s1c)
paths, lengths= shortestPath(graph, u-1, v-1,num_changes+1)
mod1 = len(paths[lengths.index(min(lengths))])-2
wtr = lengths[lengths.index(min(lengths))]


edges, graph = create_graph(nodes, num_changes, 1, 2, ar, s1r, ac, s1c)
paths, lengths= shortestPath(graph, u-1, v-1,num_changes+1)
mod2 = len(paths[lengths.index(min(lengths))])-2
wtc = lengths[lengths.index(min(lengths))]
#%%
edges, graph = create_graph(nodes, num_changes, 2, 1, ar, s2r, ac, s2c)
paths, lengths= shortestPath(graph, u-1, v-1,num_changes+1)
mod3 = len(paths[lengths.index(min(lengths))])-2



#Solvig for the two services combined
ub = max(mod1, mod2, mod3)

res = []
graph_triage = [None]*ub
length_triage = [None]*ub
beta_triage = [None]*ub
path_triage = [None]*ub
graph_crisis = [None]*ub
length_crisis = [None]*ub
beta_crisis = [None]*ub
path_crisis = [None]*ub
graph_comb = [None]*ub
length_comb = [None]*ub
path_comb = [None]*ub
shrtspth = [None]*ub

#resolving the separate problem with an upper bound
for i in range(1,ub+1):
    edges, graph_triage[i-1] = create_graph(nodes, i, 1, 1, ar, s1r, ac, s1c)
    path_triage[i-1], length_triage[i-1]= shortestPath(graph_triage[i-1], u-1, v-1,i+1)
    beta_triage[i-1] = length_triage[i-1][length_triage[i-1].index(min(length_triage[i-1]))]
    
    edges, graph_crisis[i-1] = create_graph(nodes, i, 1, 2, ar, s1r, ac, s1c)
    path_crisis[i-1], length_crisis[i-1]= shortestPath(graph_crisis[i-1], u-1, v-1,i+1)
    beta_crisis[i-1] = length_crisis[i-1][length_crisis[i-1].index(min(length_crisis[i-1]))]
    
    edges, graph_comb[i-1] = create_graph(nodes, i, 2, 1, ar, s2r, ac, s2c)
    path_comb[i-1], length_comb[i-1]= shortestPath(graph_comb[i-1], u-1, v-1,i+1)
    shrtspth[i-1] = path_comb[i-1][length_comb[i-1].index(min(length_comb[i-1]))]
    for j in range(len(path_comb[i-1])):
        print('A path from node ' + str(u)
              + ' to node ' + str(v) + ' Exist, having the following nodes: ')
        print(*path_comb[i-1][j], sep = ' --> ') 
        print('The length of the above path from node ' + str(u)
              + ' to node ' + str(v) + ' is ' + str(length_comb[i-1][j]))
        
    print('The shortest path from node ' + str(u)
          + ' to node ' + str(v) + ' has the following nodes: ')
    print(*path_comb[i-1][length_comb[i-1].index(min(length_comb[i-1]))], sep = ' --> ') 
    print('The length of the shortest path from node ' + str(u)
          + ' to node ' + str(v) + ' is ' + str(length_comb[i-1][length_comb[i-1].index(min(length_comb[i-1]))]))
    
    proptbl = np.zeros((18,3))
    propsr = []
    propsc = []
    count = 0
    for k in range(0,len(shrtspth[i-1])-1):
        for l in range(shrtspth[i-1][k+1]-shrtspth[i-1][k]):
            propsr.append(sum((x < (shrtspth[i-1][k+1]-1)*40 and x >= (shrtspth[i-1][k]-1)*40) for x in s2r)/(40*(shrtspth[i-1][k+1]-shrtspth[i-1][k])))
            propsc.append(sum((x < (shrtspth[i-1][k+1]-1)*40 and x >= (shrtspth[i-1][k]-1)*40) for x in s2c)/(40*(shrtspth[i-1][k+1]-shrtspth[i-1][k])))
            proptbl[count][0] = count+1
            proptbl[count][1] = propsr[count]
            proptbl[count][2] = propsc[count]
            count += 1
    res.append(proptbl)

#%%

for i in range(len(res)):
   pd.DataFrame(res[i], columns=['Week', 'Triage', 'Crisis']).to_csv(f'E3_Policy_{pr[0]}_{pc[0]}_{i+1}mod.csv', index = False)


#%%
#%%

# Define the 7 points
x =  list(range(1, 6))

y = [min(length_comb[i]) for i in range(5)]

# Plot the curve and the points
plt.scatter(x, y)
plt.plot(x, y, '-')

# Set the plot labels and title
plt.xlabel('# of modifications')
plt.ylabel('Weighted objective function')
# plt.title('Curve Fitting')

# Show the legend
# plt.legend()

# Show the plot
plt.show()

