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

def Contour(x,y,z):
    x_dom = x
    y_dom = 1-y
    X, Y = np.meshgrid(x_dom, y_dom)
    plt.figure(2)
    plt.contourf(X, Y, z, 50, alpha=0.5, cmap='jet')  
    plt.colorbar()
    

def Vector(x,y,u,v):
    x_dom = x
    y_dom = 1-y
    X, Y = np.meshgrid(x_dom, y_dom)
    plt.quiver(X[::2, ::2], Y[::2, ::2], u[::2, ::2], v[::2, ::2]) 
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.show()

