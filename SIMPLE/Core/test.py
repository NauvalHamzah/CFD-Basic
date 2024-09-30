import math
import numpy as np
import matplotlib.pyplot as plt
import os

n_points=6
Ae=np.zeros((n_points-1,n_points-1))
Aw=np.zeros((n_points-1,n_points-1))
An=np.zeros((n_points-1,n_points-1))
As=np.zeros((n_points-1,n_points-1))
Ap=np.zeros((n_points-1,n_points-1))

for i in range (n_points-1):
    for j in range (n_points-1):
        Ap[i,j] = n_points+3 + j + i*(n_points+1)

Ae = Ap+1
Aw = Ap-1
An = Ap-(n_points+1)
As = Ap+(n_points+1)

M=np.zeros(((n_points-1)*(n_points-1),(n_points-1)*(n_points-1)))

e=0 #Diagonal
for i in range (n_points-1):
    for j in range (n_points-1):
        M[e,e]=Ap[i,j]
        e=e+1

e=0 #West
for i in range (n_points-1):
    for j in range (n_points-1):
        if e==0:
            pass
        elif j==0:
            M[e,e]=M[e,e]-Aw[i,j]
        else:
            M[e,e-1]=Aw[i,j]
        e=e+1

e=0 #East
for i in range (n_points-1):
    for j in range (n_points-1):
        if e==(n_points-1)*(n_points-1)-1:
            pass
        elif j==(n_points-2):
            M[e,e]=M[e,e]-Ae[i,j]
        else:
            M[e,e+1]=Ae[i,j]
        e=e+1

e=0 #North
for i in range (n_points-1):
    for j in range (n_points-1):
        if i==0:
            M[e,e] = M[e,e] - An[i,j]
        else:
            M[e,e-(n_points-1)] = An[i,j]
        e=e+1

e=0 #South
for i in range (n_points-1):
    for j in range (n_points-1):
        if i==n_points-2:
            M[e,e] = M[e,e] - As[i,j]
        else:
            M[e,e+(n_points-1)] = As[i,j] 
        e=e+1

np.savetxt('M.txt',M,delimiter=',')