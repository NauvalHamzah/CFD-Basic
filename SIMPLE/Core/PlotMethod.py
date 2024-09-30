import math
import numpy as np
import matplotlib.pyplot as plt
import os

def Residual(numiter,Error):
    fig = plt.figure(figsize=(9, 6))
    plt.semilogy(numiter, Error)
    plt.xlabel('Iterations')
    plt.ylabel("Residual Error")
    plt.grid(True,which='both')
    plt.show()

def Contour(x,y,z,n):
    x_dom = x
    y_dom = 1-y
    X, Y = np.meshgrid(x_dom, y_dom)
    plt.figure(2)
    plt.contourf(X, Y, z, n, alpha=0.5, cmap='jet')  
    plt.colorbar()

def ContourP(x,y,z,n):
    x_dom = x
    y_dom = 1-y
    X, Y = np.meshgrid(x_dom, y_dom)
    plt.figure(2)
    plt.contourf(X, Y, z, n, alpha=0.5, cmap='jet')  
    plt.colorbar()
    plt.show()
    

def Vector(x,y,u,v,s):
    x_dom = x
    y_dom = 1-y
    X, Y = np.meshgrid(x_dom, y_dom)
    plt.quiver(X[::s, ::s], Y[::s, ::s], u[::s, ::s], v[::s, ::s]) 
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.show()


