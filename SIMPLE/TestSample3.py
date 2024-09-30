import math
import numpy as np
import matplotlib.pyplot as plt
import os
import Core.CoreMethod as CoreMethod
import Core.PlotMethod as PlotMethod

# Directory
# Create directory
dirName = "results"

if not os.path.exists(dirName):
    os.mkdir(dirName)
else:
    pass


#Defining the problem domain
n_inlet = 3
N = (n_inlet-1)*5+1
n_points = N
L= 1
A = L/(n_points-1)
x = np.linspace(0,L,N) #X domain span
y = np.linspace(0,L,N) #Y domain span
Error = []

rho = 1.225
miu = 1.7894e-5
alfa = 0.8
v1 = 0.001
MaxIter = 1000
Threshold = 1e-5

#Initialization
u,u_star,u_new,u_final,d_u,v,v_star,v_new,v_final,d_v,p,p_star,p_new,p_final,p_prime,b = CoreMethod.SIMPLEinit(n_points,n_inlet,v1)

for i in range (MaxIter):

    #-----------SOLVING MOMENTUM EQUATIONS------------#
    u_star,d_u=CoreMethod.SolveU_star(u,u_star,d_u,v,p,n_points,rho,A,miu)
    np.savetxt('1/u_star.txt',u_star,delimiter=',')
    v_star,d_v=CoreMethod.SolveV_star(u,v,v_star,d_v,p,n_points,n_inlet,rho,A,miu,v1)
    np.savetxt('1/v_star.txt',v_star,delimiter=',')

    #-----------SOLVING PRESSURE CORRECTION------------#
    p_prime,b = CoreMethod.SolveP_prime(u_star,d_u,v_star,d_v,n_points,p_prime,b,rho,A)

    #-----------UPDATE VALUE-------------------#
    u_new,v_new,p_new = CoreMethod.UpdateValue(u_star,u_new,d_u,v_star,v_new,d_v,p,p_prime,p_new,n_points,n_inlet,alfa,v1)
    
    u=u_new
    v=v_new
    p=p_new

    #Error
    error = CoreMethod.Error(b,n_points)
    Error.append(error)
    print(i,error)
    if error < Threshold: break



Error_all = np.zeros(len(Error))
numiter = np.zeros(len(Error))

for i in range (len(Error)):
    numiter[i] = i+1
    Error_all[i] = Error[i]

u_final,v_final,p_final,velocity = CoreMethod.FinalMapping(u,v,p,u_final,v_final,p_final,n_points)
'''
PlotMethod.Residual(numiter,Error_all)
PlotMethod.Contour(x,y,velocity,20)
PlotMethod.Vector(x,y, u_final, v_final,1)

PlotMethod.ContourP(x,y,p_final,25)
'''

np.savetxt('1/u.txt',u,delimiter=',')

np.savetxt('1/u_new.txt',u_new,delimiter=',')
np.savetxt('1/u_final.txt',u_final,delimiter=',')
np.savetxt('1/d_u.txt',d_u,delimiter=',')
np.savetxt('1/v.txt',v,delimiter=',')

np.savetxt('1/v_new.txt',v_new,delimiter=',')
np.savetxt('1/v_final.txt',v_final,delimiter=',')
np.savetxt('1/d_v.txt',d_v,delimiter=',')
np.savetxt('1/p.txt',p,delimiter=',')
np.savetxt('1/p_star.txt',p_star,delimiter=',')
np.savetxt('1/p_new.txt',p_new,delimiter=',')
np.savetxt('1/p_final.txt',p_final,delimiter=',')
np.savetxt('1/p_prime.txt',p_prime,delimiter=',')
np.savetxt('1/b.txt',b,delimiter=',')
np.savetxt('1/velocity.txt',velocity,delimiter=',')