import math
import numpy as np
import matplotlib.pyplot as plt
import os
from numpy.linalg import matrix_rank


def Gauss(M,c):
    n=len(M)
    a = np.zeros((n,n+1))
    x = np.zeros(n)
    for i in range (n):
        a[i,n]=c[i]
        for j in range (n):
            a[i,j] = M[i][j]

    for i in range (n):
        for j in range (i+1,n):
            ratio = a[j,i]/a[i,i]
            for k in range (n+1):
                a[j,k] = a[j,k]-ratio*a[i,k]

    x[n-1] = a[n-1,n]/a[n-1,n-1]

    for i in range (n-2,-1,-1):
        x[i] = a [i,n]
        for j in range (i+1,n):
            x[i]=x[i]-a[i,j]*x[j]
        
        x[i]=x[i]/a[i,i]
    return (x)

def SIMPLEinit(n_points,n_inlet,v1):
    #Initializing the variables
    #Final collocated variables
    u_final = np.zeros((n_points,n_points))
    v_final = np.zeros((n_points,n_points))
    p_final = np.zeros((n_points,n_points))
    p_final[n_points-1,n_points-1] = 1
    
    for i in range (n_inlet):
        v_final[n_points-1,2*(n_inlet-1)+i] = v1
        v_final[0,n_points-1-i] = v1
    

    #Staggered variables
    u = np.zeros((n_points+1,n_points))
    u_star = np.zeros((n_points+1,n_points))
    d_u = np.zeros((n_points+1,n_points))
    v = np.zeros((n_points,n_points+1))
    v_star = np.zeros((n_points,n_points+1))
    d_v = np.zeros((n_points,n_points+1))
    p = np.zeros((n_points+1,n_points+1))
    p[n_points,n_points] = 1
    p_star = np.zeros((n_points+1,n_points+1))
    p_star[n_points,n_points] = 1
    p_prime = np.zeros((n_points+1,n_points+1))
    b = np.zeros((n_points+1,n_points+1))
    for i in range (n_inlet-1):
        v[n_points-1,2*(n_inlet-1)+1+i] = v1
        v[0,n_points-1-i] = v1
    

    u_new = np.zeros((n_points+1,n_points))
    v_new = np.zeros((n_points,n_points+1))
    p_new = np.zeros((n_points+1,n_points+1))
    p_new[n_points,n_points] = 1
    for i in range (n_inlet-1):
        v_new[n_points-1,2*(n_inlet-1)+1+i] = v1
        v_new[0,n_points-1-i] = v1

    return(u,u_star,u_new,u_final,d_u,v,v_star,v_new,v_final,d_v,p,p_star,p_new,p_final,p_prime,b)

def SolveU_star(u,u_star,d_u,v,p,n_points,rho,A,miu):
    #Interior
    Ae=np.zeros((n_points-1,n_points-2))
    Aw=np.zeros((n_points-1,n_points-2))
    An=np.zeros((n_points-1,n_points-2))
    As=np.zeros((n_points-1,n_points-2))
    Ap=np.zeros((n_points-1,n_points-2))
    Su=np.zeros((n_points-1,n_points-2))
    for i in range (1,n_points):
        for j in range (1,n_points-1):
            Fe = 0.5*rho*A*(u[i,j]+u[i,j+1])
            Fw = 0.5*rho*A*(u[i,j]+u[i,j-1])
            Fn = 0.5*rho*A*(v[i-1,j]+v[i-1,j+1])
            Fs = 0.5*rho*A*(v[i,j]+v[i,j+1])
            
            
            Ae[i-1,j-1] = max((-0.5*Fe + (miu/A)),-Fe,0)
            Aw[i-1,j-1] = max((0.5*Fw + (miu/A)), Fw, 0)
            An[i-1,j-1] = max((-0.5*Fn + (miu/A)),-Fn,0)
            As[i-1,j-1] = max((0.5*Fs + (miu/A)), Fs, 0)
            
            Ap[i-1,j-1] = Ae[i-1,j-1]+Aw[i-1,j-1]+An[i-1,j-1]+As[i-1,j-1]+(Fe-Fw+Fn-Fs)
            
            d_u[i,j] = -A/Ap[i-1,j-1]
            Su[i-1,j-1] = d_u[i,j]*(p[i,j+1] - p[i,j])
            
    M=np.zeros(((n_points-1)*(n_points-2),(n_points-1)*(n_points-2)))

    e=0 #Diagonal
    for i in range (n_points-1):
        for j in range (n_points-2):
            M[e,e]=Ap[i,j]
            e=e+1

    e=0 #West
    for i in range (n_points-1):
        for j in range (n_points-2):
            if e==0:
                pass
            elif j==0:
                M[e,e-1]=0
            else:
                M[e,e-1]=-Aw[i,j]
            e=e+1

    e=0 #East
    for i in range (n_points-1):
        for j in range (n_points-2):
            if e==(n_points-1)*(n_points-2)-1:
                pass
            elif j==(n_points-3):
                M[e,e+1]=0
            else:
                M[e,e+1]=-Ae[i,j]
            e=e+1

    e=0 #North
    for i in range (n_points-1):
        for j in range (n_points-2):
            if i==0:
                M[e,e] = M[e,e]+An[i,j]
            else:
                M[e,e-(n_points-2)] = -An[i,j]
            e=e+1

    e=0 #South
    for i in range (n_points-1):
        for j in range (n_points-2):
            if i==n_points-2:
                M[e,e] = M[e,e]+As[i,j]
            else:
                M[e,e+(n_points-2)] = -As[i,j] 
            e=e+1

    c=np.zeros((n_points-1)*(n_points-2))

    e=0 #Constant
    for i in range (n_points-1):
        for j in range (n_points-2):
            c[e]=Su[i,j]
            e=e+1
    print('uSolver ', matrix_rank(M))
    np.savetxt('2/Mu.txt',M,delimiter=',')
    u_sol = Gauss(M,c)
    
    e=0 #Solution Array
    for i in range (1,n_points):
        for j in range (1,n_points-1):
            u_star[i,j] = u_sol[e]
            e=e+1

    #Boundary
    u_star[1:n_points-1,0] = 0
    u_star[1:n_points-1,n_points-1] = 0
    u_star[0,:] = -u_star[1,:]
    u_star[n_points,:] = -u_star[n_points-1,:]
    
    return(u_star,d_u)

def SolveV_star(u,v,v_star,d_v,p,n_points,n_inlet,rho,A,miu,v1):
    #Interior
    Ae=np.zeros((n_points-2,n_points-1))
    Aw=np.zeros((n_points-2,n_points-1))
    An=np.zeros((n_points-2,n_points-1))
    As=np.zeros((n_points-2,n_points-1))
    Ap=np.zeros((n_points-2,n_points-1))
    Su=np.zeros((n_points-2,n_points-1))
    for i in range (1,n_points-1):
        for j in range (1,n_points):
            Fe = 0.5*rho*A*(u[i,j]+u[i+1,j])
            Fw = 0.5*rho*A*(u[i,j-1]+u[i+1,j-1])
            Fn = 0.5*rho*A*(v[i-1,j]+v[i,j])
            Fs = 0.5*rho*A*(v[i,j]+v[i+1,j])
            
            Ae[i-1,j-1] = max((-0.5*Fe + (miu/A)),-Fe,0)
            Aw[i-1,j-1] = max((0.5*Fw + (miu/A)), Fw, 0)
            An[i-1,j-1] = max((-0.5*Fn + (miu/A)),-Fn,0)
            As[i-1,j-1] = max((0.5*Fs + (miu/A)), Fs, 0)
            
            
            Ap[i-1,j-1] = Ae[i-1,j-1]+Aw[i-1,j-1]+An[i-1,j-1]+As[i-1,j-1]+(Fe-Fw+Fn-Fs)
            
            d_v[i,j] = -A/Ap[i-1,j-1]
            Su[i-1,j-1] = d_v[i,j]*(p[i,j] - p[i+1,j])
    
    M=np.zeros(((n_points-2)*(n_points-1),(n_points-2)*(n_points-1)))

    e=0 #Diagonal
    for i in range (n_points-2):
        for j in range (n_points-1):
            M[e,e]=Ap[i,j]
            e=e+1

    e=0 #West
    for i in range (n_points-2):
        for j in range (n_points-1):
            if e==0:
                pass
            elif j==0:
                M[e,e]=M[e,e]+Aw[i,j]
            else:
                M[e,e-1]=-Aw[i,j]
            e=e+1

    e=0 #East
    for i in range (n_points-2):
        for j in range (n_points-1):
            if e==(n_points-2)*(n_points-1)-1:
                pass
            elif j==(n_points-2):
                M[e,e]=M[e,e]+Ae[i,j]
            else:
                M[e,e+1]=-Ae[i,j]
            e=e+1

    e=0 #North
    for i in range (n_points-2):
        for j in range (n_points-1):
            if i==0:
                pass
            else:
                M[e,e-(n_points-1)] = -An[i,j]
            e=e+1

    e=0 #South
    for i in range (n_points-2):
        for j in range (n_points-1):
            if i==n_points-3:
                pass
            else:
                M[e,e+(n_points-1)] = -As[i,j] 
            e=e+1
    
    c=np.zeros((n_points-1)*(n_points-2))

    e=0 #Constant
    for i in range (n_points-2):
        for j in range (n_points-1):
            c[e]=Su[i,j]
            e=e+1
    np.savetxt('2/Mv.txt',M,delimiter=',')
    print('vSolver ', matrix_rank(M))
    v_sol = Gauss(M,c)
    
    e=0 #Solution Array
    for i in range (1,n_points-1):
        for j in range (1,n_points):
            v_star[i,j] = v_sol[e]
            e=e+1

    #Boundary
    v_star[0,1:n_points-1] = 0
    v_star[n_points-1,1:n_points-1] = 0
    v_star[:,0] = -v_star[:,1]
    v_star[:,n_points] = -v_star[:,n_points-1]
    
    for i in range (n_inlet-1):
        v_star[n_points-1,2*(n_inlet-1)+1+i] = v1
        v_star[0,n_points-1-i] = v1

    return (v_star,d_v)

def SolveP_prime(u_star,d_u,v_star,d_v,n_points,p_prime,b,rho,A):
    #Zeroing the corrections to begin with
    p_prime[:,:]=0
    
    #Continuity equation a.k.a. pressure correction - Interior
    for i in range (1,n_points):
        for j in range (1,n_points):
            Ae = -rho*A*d_u[i,j]
            Aw = -rho*A*d_u[i,j-1]
            An = -rho*A*d_v[i-1,j]
            As = -rho*A*d_v[i,j]
            Ap = Ae + Aw + An + As
            b[i,j] = -(u_star[i,j] - u_star[i,j-1])*A + (v_star[i,j] - v_star[i-1,j])*A
            
            p_prime[i,j] = (Ae*p_prime[i,j+1] + Aw*p_prime[i,j-1] + An*p_prime[i-1,j] + As*p_prime[i+1,j] + b[i,j])/Ap
    '''
    #Zeroing the corrections to begin with
    p_prime[:,:]=0

    #Continuity equation a.k.a. pressure correction - Interior
    #Zeroing the corrections to begin with
    p_prime[:,:]=0
    Ae=np.zeros((n_points-1,n_points-1))
    Aw=np.zeros((n_points-1,n_points-1))
    An=np.zeros((n_points-1,n_points-1))
    As=np.zeros((n_points-1,n_points-1))
    Ap=np.zeros((n_points-1,n_points-1))
    Su=np.zeros((n_points-1,n_points-1))
    #Continuity equation a.k.a. pressure correction - Interior
    for i in range (1,n_points):
        for j in range (1,n_points):
            Ae[i-1,j-1] = -rho*A*d_u[i,j]
            Aw[i-1,j-1] = -rho*A*d_u[i,j-1]
            An[i-1,j-1] = -rho*A*d_v[i-1,j]
            As[i-1,j-1] = -rho*A*d_v[i,j]
            Ap[i-1,j-1] = Ae[i-1,j-1] + Aw[i-1,j-1] + An[i-1,j-1] + As[i-1,j-1]
            b[i,j] = -(u_star[i,j] - u_star[i,j-1])*A + (v_star[i,j] - v_star[i-1,j])*A
            Su[i-1,j-1] = b[i,j]

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

    c=np.zeros((n_points-1)*(n_points-1))

    e=0 #Constant
    for i in range (n_points-1):
        for j in range (n_points-1):
            c[e]=Su[i,j]
            e=e+1
    np.savetxt('2/Mp.txt',M,delimiter=',')
    np.savetxt('2/cp.txt',c,delimiter=',')
    print('pSolver ', matrix_rank(M))
    p_sol = Gauss(M,c)
    
    e=0 #Solution Array
    for i in range (1,n_points):
        for j in range (1,n_points):
            p_prime[i,j] = p_sol[e]
            e=e+1
    '''
    return (p_prime,b)

def UpdateValue(u_star,u_new,d_u,v_star,v_new,d_v,p,p_prime,p_new,n_points,n_inlet,alfa,v1):
    #Correcting the pressure field
    for i in range (1,n_points):
        for j in range (1,n_points):
            p_new[i,j] = p[i,j] + alfa*p_prime[i,j]
        
    #Continuity eq. - Boundary
    p_new[0,:] = p_new[1,:]
    p_new[n_points,:] = p_new[n_points-1,:]
    p_new[:,0] = p_new[:,1]
    p_new[:,n_points] = p_new[:,n_points-1]
    
    #Correcting the velocities
    for i in range (1,n_points):
        for j in range (1,n_points-1):
            u_new[i,j] = u_star[i,j] + alfa*d_u[i,j]*(p_prime[i,j+1] - p_prime[i,j])
        
    
    #x-momentum eq. - Boundary
    u_new[1:n_points-1,0] = 0
    u_new[1:n_points-1,n_points-1] = 0
    u_new[0,:] = -u_new[1,:]
    u_new[n_points,:] = -u_new[n_points-1,:]
    
    
    for i in range (1,n_points-1):
        for j in range (1,n_points):
            v_new[i,j] = v_star[i,j] + alfa*d_v[i,j]*(p_prime[i,j] - p_prime[i+1,j])
       
    
    #y-momentum eq. - Boundary
    v_new[0,1:n_points-1] = 0
    v_new[n_points-1,1:n_points-1] = 0
    v_new[:,0] = -v_new[:,1]
    v_new[:,n_points] = -v_new[:,n_points-1]
    

    for i in range (n_inlet-1):
        v_new[n_points-1,2*(n_inlet-1)+1+i] =v1
        v_new[0,n_points-1-i] = v1
    

    return (u_new,v_new,p_new)

def Error(b,n_points):
    # Continuity residual as error measure
    error = 0
    for i in range (1,n_points):
        for j in range (1,n_points):
            error = error + abs(b[i,j])
    return (error)
        
def FinalMapping(u,v,p,u_final,v_final,p_final,n_points):
    velocity=np.zeros((len(u_final),len(u_final[0])))
    for i in range (n_points):
        for j in range (n_points):
            u_final[i,j] = 0.5*(u[i,j] + u[i+1,j])
            v_final[i,j] = 0.5*(v[i,j] + v[i,j+1])
            p_final[i,j] = 0.25*(p[i,j] + p[i,j+1] + p[i+1,j] + p[i+1,j+1])
            velocity[i,j] = np.sqrt((u_final[i,j])**2 + (v_final[i,j])**2)
    return(u_final,v_final,p_final,velocity)
        
