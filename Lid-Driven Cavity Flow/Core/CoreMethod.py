import math
import numpy as np
import matplotlib.pyplot as plt
import os


def SIMPLEinit(n_points,n_inlet,v1):
    #Initializing the variables
    #Final collocated variables
    u_final = np.zeros((n_points,n_points))
    v_final = np.zeros((n_points,n_points))
    p_final = np.zeros((n_points,n_points))
    p_final[n_points-1,n_points-1] = 1
    u_final[0,:] = 1
    

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
    u[0,:] = 2
    

    u_new = np.zeros((n_points+1,n_points))
    v_new = np.zeros((n_points,n_points+1))
    p_new = np.zeros((n_points+1,n_points+1))
    p_new[n_points,n_points] = 1
    u_new[0,:] = 2

    return(u,u_star,u_new,u_final,d_u,v,v_star,v_new,v_final,d_v,p,p_star,p_new,p_final,p_prime,b)

def SolveU_star(u,u_star,d_u,v,p,n_points,rho,A,miu):
    #Interior
    for i in range (1,n_points):
        for j in range (1,n_points-1):
            Fe = 0.5*rho*A*(u[i,j]+u[i,j+1])
            Fw = 0.5*rho*A*(u[i,j]+u[i,j-1])
            Fn = 0.5*rho*A*(v[i-1,j]+v[i-1,j+1])
            Fs = 0.5*rho*A*(v[i,j]+v[i,j+1])
            
            Ae = -0.5*Fe + (miu/A)
            Aw = 0.5*Fw + (miu/A)
            An = -0.5*Fn + (miu/A)
            As = 0.5*Fs + (miu/A)
            
            Ap = 0.5*Fe - 0.5*Fw + 0.5*Fn - 0.5*Fs + 4*(miu/A)
            
            d_u[i,j] = -A/Ap
            
            u_star[i,j] = (Ae*u[i,j+1] + Aw*u[i,j-1] + An*u[i-1,j] + As*u[i+1,j])/Ap + d_u[i,j]*(p[i,j+1] - p[i,j])
    
    #Boundary
    u_star[0,:] = 2 - u_star[1,:]
    u_star[n_points,:] = -u_star[n_points-1,:]
    u_star[1:n_points-1,0] = 0
    u_star[1:n_points-1,n_points-1] = 0

    return(u_star,d_u)

def SolveV_star(u,v,v_star,d_v,p,n_points,rho,A,miu):
    #Interior
    for i in range (1,n_points-1):
        for j in range (1,n_points):
            Fe = 0.5*rho*A*(u[i,j]+u[i+1,j])
            Fw = 0.5*rho*A*(u[i,j-1]+u[i+1,j-1])
            Fn = 0.5*rho*A*(v[i-1,j]+v[i,j])
            Fs = 0.5*rho*A*(v[i,j]+v[i+1,j])
            
            Ae = -0.5*Fe + (miu/A)
            Aw = 0.5*Fw + (miu/A)
            An = -0.5*Fn + (miu/A)
            As = 0.5*Fs + (miu/A)
            
            Ap = 0.5*Fe - 0.5*Fw + 0.5*Fn - 0.5*Fs + 4*(miu/A)
            
            d_v[i,j] = -A/Ap
            
            v_star[i,j] = (Ae*v[i,j+1] + Aw*v[i,j-1] + An*v[i-1,j] + As*v[i+1,j])/Ap + d_v[i,j]*(p[i,j] - p[i+1,j])
    
    #Boundary
    v_star[:,0] = -v_star[:,1]
    v_star[:,n_points] = -v_star[:,n_points-1]
    v_star[0,1:n_points-1] = 0
    v_star[n_points-1,1:n_points-1] = 0

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
    u_new[0,:] = 2 - u_new[1,:]
    u_new[n_points,:] = -u_new[n_points-1,:]
    u_new[1:n_points-1,0] = 0
    u_new[1:n_points-1,n_points-1] = 0
    
    for i in range (1,n_points-1):
        for j in range (1,n_points):
            v_new[i,j] = v_star[i,j] + alfa*d_v[i,j]*(p_prime[i,j] - p_prime[i+1,j])
       
    
    #y-momentum eq. - Boundary
    v_new[:,0] = -v_new[:,1]
    v_new[:,n_points] = -v_new[:,n_points-1]
    v_new[0,1:n_points-1] = 0
    v_new[n_points-1,1:n_points-1] = 0
    

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
        
