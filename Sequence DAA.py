# -*- coding: utf-8 -*-
"""
Created on Thu Jan 12 08:19:10 2023

@author: Sanford Jone
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Jan 12 08:00:17 2023

@author: Sanford Jone
"""

import numpy as np
import random

s2="AGCACACA"
s1="ACACACTA"
m1 = len(s1)+1
n = len(s2)+1

m = np.ones((m1,n)) * -100
parent = np.zeros((m1,n,2))
m[:,0] = 0
m[0,:] = 0

match = 2
mismatch = -1

def allign(i,j):
    
    global s1,s2,m,parent
    
    if i>0 and j>0:
        if m[j][i-1] == -100:
            allign(i-1,j)
        if m[j-1][i-1] == -100:
            allign(i-1,j-1)
        if m[j-1][i] == -100:
            allign(i,j-1)
        
        if s1[i-1] == s2[j-1]:
            if m[j-1][i-1] != -100:
                m[j][i] = m[j-1][i-1]+match
                
        else:
            l = np.array([m[j][i-1],m[j-1][i-1],m[j-1][i]])
            l = l + mismatch
            
            if l[1] >= l[0] and l[1]>=l[2]:
                m[j][i] = l[1]
                parent[j][i] = [i-1,j-1]
            elif l[0]>=l[0] and l[0]>=l[2]:
                m[j][i] = l[0]
                parent[j][i] = [i-1,j]
            else:
                m[j][i] = l[2]
                parent[j][i] = [i,j-1]
    
    
    

allign(len(s1),len(s2))    
#trace(m,trm,len(s1),len(s2))
print(m)
#print(trm)


#%%

def seq(seq1, seq2, match=1, mismatch=-1, gap=-1):
    m, n = len(seq1), len(seq2)
    dp = [[0 for _ in range(n + 1)] for _ in range(m + 1)]
    traceback = [[0 for _ in range(n + 1)] for _ in range(m + 1)]
    
    for i in range(m + 1):
        for j in range(n + 1):
            if i == 0 and j == 0:
                dp[i][j] = 0
                traceback[i][j] = "E"
            elif i == 0:
                dp[i][j] = dp[i][j - 1] + gap
                traceback[i][j] = "L"
            elif j == 0:
                dp[i][j] = dp[i - 1][j] + gap
                traceback[i][j] = "U"
            else:
                dp[i][j], traceback[i][j] = max((dp[i - 1][j - 1] + (match if seq1[i - 1] == seq2[j - 1] else mismatch), "D"), 
                                               (dp[i][j - 1] + gap, "L"), 
                                               (dp[i - 1][j] + gap, "U"), 
                                               key=lambda x: x[0])
    align1, align2 = "", ""
    i, j = m, n
    while i > 0 or j > 0:
        if traceback[i][j] == "D":
            align1 = seq1[i - 1] + align1
            align2 = seq2[j - 1] + align2
            i -= 1
            j -= 1
        elif traceback[i][j] == "L":
            align1 = "-" + align1
            align2 = seq2[j - 1] + align2
            j -= 1
        elif traceback[i][j] == "U":
            align1 = seq1[i - 1] + align1
            align2 = "-" + align2
            i -= 1
    return (dp[m][n], align1, align2)

score, align1, align2 = seq("AGCACACA", "ACACACTA")
print("Score:", score)
print("Alignment:")
print(align1)
print(align2)

#%%

i=len(s1)
j=len(s2)
trm_a=[]

#%%
while i!=0 and j!=0:
    print(i,j)
    trm_a+=[[i,j]]
    i,j=int(parent[j][i][0]),int(parent[j][i][1])
trm_a+=[[0,0]]
#%%
 trm=[]    
 def trace(mat , trm, x, y):
     if x==1 and y==1:
         return trm
     elif (mat[x,y] == (mat[x - 1,y - 1] + mismatch) or mat[x,y] == (mat[x - 1,y - 1] + match)):
         trm.append([x - 1, y - 1])
         return trace(mat, trm, x-1, y-1)
     elif (mat[x,y] == (mat[x-1,y] + mismatch) ):
         trm.append([x-1, y])
         return trace(mat, trm, x-1, y)
     elif (mat[x,y] == (mat[x,y-1] + mismatch) ):
         trm.append([x, y-1])
         return trace(mat, trm, x, y-1)    
#%%
trace(m,trm,len(s1),len(s2))
print(trm)
#%%


#%%
def seq(seq1, seq2, match_score=1, mismatch_score=-1, gap_score=-1):
    m = len(seq1)
    n = len(seq2)
    dp = [[0 for _ in range(n + 1)] for _ in range(m + 1)]
    path = [['' for _ in range(n + 1)] for _ in range(m + 1)]
    
    for i in range(m + 1):
        for j in range(n + 1):
            if i == 0:
                dp[i][j] = j * gap_score
                path[i][j] = 'L'
            elif j == 0:
                dp[i][j] = i * gap_score
                path[i][j] = 'U'
            elif seq1[i-1] == seq2[j-1]:
                dp[i][j] = dp[i-1][j-1] + match_score
                path[i][j] = 'D'
            else:
                if dp[i-1][j-1] + mismatch_score > dp[i-1][j] + gap_score and dp[i-1][j-1] + mismatch_score > dp[i][j-1] + gap_score:
                    dp[i][j] = dp[i-1][j-1] + mismatch_score
                    path[i][j] = 'D'
                elif dp[i-1][j] + gap_score > dp[i-1][j-1] + mismatch_score and dp[i-1][j] + gap_score > dp[i][j-1] + gap_score:
                    dp[i][j] = dp[i-1][j] + gap_score
                    path[i][j] = 'U'
                else:
                    dp[i][j] = dp[i][j-1] + gap_score
                    path[i][j] = 'L'

    align1 = ""
    align2 = ""
    i = m
    j = n

    while i > 0 or j > 0:
        if path[i][j] == 'D':
            align1 = seq1[i-1] + align1
            align2 = seq2[j-1] + align2
            i -= 1
            j -= 1
        elif path[i][j] == 'U':
            align1 = seq1[i-1] + align1
            align2 = '-' + align2
            i -= 1
        elif path[i][j] == 'L':
            align1 = '-' + align1
            align2 = seq2[j-1] + align2
            j -= 1

    return align1, align2        

x,y=seq(s1,s2)       
print(x)
print(y)

#%%
def scoring_model(s1, s2, match_score=1, mismatch_score=-1, gap_score=-1):
    m = len(s1)
    n = len(s2)
    dp = [[0 for _ in range(n + 1)] for _ in range(m + 1)]
    
    for i in range(m + 1):
        for j in range(n + 1):
            if i == 0:
                dp[i][j] = j * gap_score
            elif j == 0:
                dp[i][j] = i * gap_score
            elif s1[i-1] == s2[j-1]:
                dp[i][j] = dp[i-1][j-1] + match_score
            else:
                dp[i][j] = max(dp[i-1][j-1] + mismatch_score, dp[i-1][j] + gap_score, dp[i][j-1] + gap_score)
    return dp[m][n]

s=scoring_model(s1,s2)       
print(s)