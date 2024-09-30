import math
import matplotlib.pyplot as plt
import numpy as np
import csv


# ---------------------------------------------------------------------------- #


def PlotArea(X, A, nMax, L, title="Chart Area"):
    # generate plot
    r = np.zeros(nMax)
    for i in range(nMax):
        r[i] = np.sqrt(A[i]*0.0025/math.pi)

    fig = plt.figure(figsize=(13.5,4.5))
    plt.plot(X, r, 'r', label=r'$A(x) = 0.0025 + 0.01(x-1)^{2} ,0<x<1$ and $A(x) = 0.0025 + 0.0025(x-1)^{2} ,1<x<2$')
    plt.plot(X, -r, 'r')
    plt.xlabel(r'$X^{*}$')
    plt.ylabel(r"$r(x)$")
    plt.legend(prop={"size":10}, loc='upper center')
    plt.xlim(0, L)
    plt.minorticks_on()
    plt.grid(b=True, which='major', color='0.65', linestyle='-')
    plt.grid(b=True, which='minor', color='0.65', linestyle='--')

    plt.show()      # show plot

    # save plot
    fName = "results/" + title + ".png"
    fig.savefig(fName, dpi = 150)

def PlotResidue(numiter, resRho, resT, resV, resU, iters, title="Residual"):
    fig = plt.figure(figsize=(9, 6))
    plt.subplot(211)
    plt.plot(numiter, resRho, color='r', linestyle='-', linewidth='1', label=r'$\rho$')
    plt.plot(numiter, resT, color='g', linestyle='--', linewidth='1', label=r'T')
    plt.plot(numiter, resV, color='b', linestyle='-.', linewidth='1', label=r'V')
    plt.legend(prop={"size":10}, loc='upper right')
    plt.xlim(1, iters )
    plt.minorticks_on()
    plt.grid(b=True, which='major', color='0.65', linestyle='-')
    plt.grid(b=True, which='minor', color='0.65', linestyle='--')
    plt.ylabel("Residuals")
    #plt.yscale("log")

    plt.subplot(212)
    plt.plot(numiter, resU[:,0], color='r', linestyle='-', linewidth='1', label=r'$U_{1}$')
    plt.plot(numiter, resU[:,1], color='g', linestyle='--', linewidth='1', label=r'$U_{2}$')
    plt.plot(numiter, resU[:,2], color='b', linestyle='-.', linewidth='1', label=r'$U_{3}$')
    plt.legend(prop={"size":10}, loc='upper right')
    plt.xlim(1, iters)
    plt.minorticks_on()
    plt.grid(b=True, which='major', color='0.65', linestyle='-')
    plt.grid(b=True, which='minor', color='0.65', linestyle='--')
    plt.ylabel("Residuals")
    plt.xlabel("Timesteps")
    #plt.yscale("log")

    # save plot
    fName = "results/" + title + ".png"
    fig.savefig(fName, dpi = 150)

def PlotPrimitiveVariables(X, A, rho, T, V, L, nMax, title="Primitive Variables"):
    r = np.zeros(nMax)
    for i in range(nMax):
        r[i] = np.sqrt(A[i]/math.pi)

    fig = plt.figure(figsize=(9,8))

    plt.subplot(411)
    plt.plot(X, r, color='b', linewidth='2', label=r'$A(x) = 1 + 2.2(x-1.5)^{2}$')
    plt.plot(X, -r, color='b', linewidth='2')
    plt.legend(prop={"size":10}, loc='center right')
    plt.xlim(0, L)
    plt.ylabel(r"$r(x)$")
    plt.minorticks_on()
    plt.grid(b=True, which='major', color='0.65', linestyle='-')
    plt.grid(b=True, which='minor', color='0.65', linestyle='--')

    plt.subplot(412)
    plt.plot(X, rho, color='r', linestyle='-', linewidth='1', label="Numerical Solution")
    plt.ylabel(r'$\rho/\rho_{0}$')   
    plt.legend(prop={"size":10}, loc='upper right')
    plt.minorticks_on()
    plt.grid(b=True, which='major', color='0.65', linestyle='-')
    plt.grid(b=True, which='minor', color='0.65', linestyle='--')
    plt.xlim(0.0, L)

    plt.subplot(413)
    plt.plot(X, T, color='r', linestyle='-', linewidth='1')
    plt.ylabel(r'$T/T_{0}$')   
    plt.minorticks_on()
    plt.grid(b=True, which='major', color='0.65', linestyle='-')
    plt.grid(b=True, which='minor', color='0.65', linestyle='--')
    plt.xlim(0.0, L)

    plt.subplot(414)
    plt.plot(X, V, color='r', linestyle='-', linewidth='1')
    plt.ylabel(r'$V/a_{0}$')   
    plt.xlabel(r'$X$')
    plt.minorticks_on()
    plt.grid(b=True, which='major', color='0.65', linestyle='-')
    plt.grid(b=True, which='minor', color='0.65', linestyle='--')
    plt.xlim(0.0, L)

    plt.show()    

    # save plot
    fName = "results/" + title + ".png"
    fig.savefig(fName, dpi = 150)

def PlotMDot(X, mDot, L, nMax, title="Mass Flow Comparison"):
    fig = plt.figure(figsize=(13.5,4.5))

    plt.plot(X, mDot, color='r', linestyle='-', linewidth='1', label="Numerical solution")
    plt.plot(X, 0.579*np.ones(nMax), color='k', linestyle='--', linewidth='2', label='exact solution')
    plt.ylabel(r'$(\rho A V)^{*}$')
    plt.xlabel(r'$X^{*}$')   
    plt.legend(prop={"size":10}, loc='upper right')
    plt.minorticks_on()
    plt.grid(b=True, which='major', color='0.65', linestyle='-')
    plt.grid(b=True, which='minor', color='0.65', linestyle='--')
    plt.xlim(0.0, L)

    plt.show()    

    # save plot
    fName = "results/" + title + ".png"
    fig.savefig(fName, dpi = 150)   

def CSVdat(Name):
    with open(Name, newline='') as csvfile:
        data = list(csv.reader(csvfile))
        A=np.zeros(len(data[0]))
        for i in range (len(data[0])):
            A[i] = float(data[0][i])
    return A

def Plotdata3(X1, X2, X3, data1, data2, data3, L, title, ylab):
    fig = plt.figure(figsize=(13.5,4.5))

    plt.plot(X1, data1, color='b', linestyle='--', linewidth='2', label='Analytical')
    plt.plot(X2, data2, color='r', linestyle='-.', linewidth='2', label='N=21')
    plt.plot(X3, data3, color='k', linestyle='-', linewidth='2', label='N=201')
    plt.ylabel(ylab)
    plt.xlabel(r'$X^{*}$')   
    plt.legend(prop={"size":10}, loc='upper right')
    plt.minorticks_on()
    plt.grid(b=True, which='major', color='0.65', linestyle='-')
    plt.grid(b=True, which='minor', color='0.65', linestyle='--')
    plt.xlim(0.0, L)  

    # save plot
    fName = "results/" + title + ".png"
    fig.savefig(fName, dpi = 150)   

def Plotdata2(X1, X2, data1, data2, L, title, ylab):
    fig = plt.figure(figsize=(13.5,4.5))

    plt.plot(X1, data1, color='r', linestyle='-.', linewidth='2', label='N=21')
    plt.plot(X2, data2, color='k', linestyle='-', linewidth='2', label='N=201')
    plt.ylabel(ylab)
    plt.xlabel(r'$X^{*}$')   
    plt.legend(prop={"size":10}, loc='upper right')
    plt.minorticks_on()
    plt.grid(b=True, which='major', color='0.65', linestyle='-')
    plt.grid(b=True, which='minor', color='0.65', linestyle='--')
    plt.xlim(0.0, L)   

    # save plot
    fName = "results/" + title + ".png"
    fig.savefig(fName, dpi = 150)   
# ---------------------------------------------------------------------------- #