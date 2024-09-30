import copy
import numpy as np
import os

import cores.CoreMethod as CoreMethod
import cores.PlotMethod as PlotMethod


# --------------------------------------------------------------------------------- #

# Directory
# Create directory
dirName = "results"

if not os.path.exists(dirName):
    os.mkdir(dirName)
else:
    pass

# --------------------------------------------------------------------------------- #

# --- INPUT: Discretization --- #
nMax = 201

# --- INITIALIZE VARIABLES --- #
X, A, p, rho, V, T, mDot, M, U, Up, F, J, S, dUp, dUc, dUav = CoreMethod.InitializeVariables(nMax)

# --- INPUT: Physical Condition --- #
# Geometry
L = 2         # Length of Nozzles

# Define boundary condition
# All the total properties is calculated in inlet condition
rhoi = 1.0      # non-dimensional inlet static density (rho/rho_t)
Ti = 1.0        # non-dimensional inlet static temperature (T/T_t)
pi = rhoi*Ti    # non-dimensional inlet static pressure (p inlet/p_t)

print ("Pilih Kasus\n")
print ("1. Isentropik\n")
print ("2. Non-Isentropik\n")
kasus = 0
pilihan = True


#Pilihan kasus isentropik atau non isentropik
while pilihan is True:
    kasus=input ("(1/2)?")
    if kasus == '1':
        pe = 0
        pilihan = False
        print ("Isentropik")
    elif kasus == '2':
        p_exit = False
        while p_exit is False:
            print ("Masukkan tekanan keluar antara 0.1-0.9\n")
            a=input ()
            pe_val=float(a)
            if pe_val > 0.1 and pe_val<0.9:
                pe=pe_val
                p_exit = True
                pilihan = False
                print ("Non-Isentropik")
            else:
                p_exit=False
    else:
        pilihan = True
        print("Ulangi pilihan\n")


# Constant
gamma = 1.4     # ratio of specific heat capacity

# --- INITIALIZATION --- #
dx, X, A, rho, T, p, V, M, mDot, U = CoreMethod.GenerateInitialCondition(X, A, rho, T, p, V, M, mDot, U,
                                                                        nMax, L, rhoi, Ti, pi, pe, gamma)

# if there is flow (pe != 1), start calculating
if (pe != 1):
    # --- INPUT: MacCormack Iteration & Artificial Viscosity --- #
    # Artifical viscosity input
    C_x = 0.2               # Artificial viscosity constant

    # Time step input
    C = 0.5                 # Courant number

    # Iteration input
    maxIter = 10000          # number of iteration

    # residual variables
    resRho = np.ones(maxIter)
    resT = np.ones(maxIter)
    resV = np.ones(maxIter)
    resU = np.ones((maxIter, 3))
    run = True
    iters = 0
    error = 0.0001

    # Start iteration
    while run is True:
        # resiudal temporary storage
        rho_old = copy.deepcopy(rho)
        V_old = copy.deepcopy(V)
        T_old = copy.deepcopy(T)
        U_old = copy.deepcopy(U)

        # Generate timestep
        dt = CoreMethod.GenerateTimeStep(C, dx, T, V, nMax)
        
        # Predictor loop
        U, Up, dUp, F, J, S, A, rho, T, p, V = CoreMethod.PredictorLoop(U, Up, dUp, F, J, S, A, rho, 
                                                                        T, p, V, gamma, dt, dx, C_x, 
                                                                        pe, nMax)
        
        # Corrector loop
        U, dUav, rho, T, p, V = CoreMethod.CorrectorLoop(U, Up, dUp, dUc, dUav, F, J, S, A, rho, T, 
                                                            p, V, gamma, dt, dx, C_x, pe, nMax)

        # Generate Residual
        resRho_temp = abs(rho-rho_old)/rho_old
        resT_temp = abs(T-T_old)/T_old
        resV_temp = abs(V-V_old)/V_old
        resU1_temp = abs(U[:,0]-U_old[:,0])/U_old[:,0]
        resU2_temp = abs(U[:,1]-U_old[:,1])/U_old[:,1]
        resU3_temp = abs(U[:,2]-U_old[:,2])/U_old[:,2]

        resRho[iters] = np.sum(resRho_temp)
        resT[iters] = np.sum(resT_temp)
        resV[iters] = np.sum(resV_temp)
        resU[iters,0] = np.sum(resU1_temp)
        resU[iters,1] = np.sum(resU2_temp)
        resU[iters,2] = np.sum(resU3_temp)

        for i in range(nMax):
            mDot[i] = rho[i]*A[i]*V[i]
            M[i] = V[i]/(np.sqrt(T[i]))
    
        if iters == maxIter-2:
            run = False
        elif resRho[iters] < error and resT[iters] <error and resV[iters] <error and resU[iters,0] <error and resU[iters,1] <error and resU[iters,2] <error:
            run = False

        iters = iters + 1

numiter = np.zeros(iters)
resRho_final = np.ones(iters)
resT_final = np.ones(iters)
resV_final = np.ones(iters)
resU_final = np.ones((iters, 3))
for n in range (iters):
    numiter[n]=n+1
    resRho_final[n]=resRho[n]
    resT_final[n]=resT[n]
    resV_final[n]=resV[n]
    resU_final[n,0]=resU[n,0]
    resU_final[n,1]=resU[n,1]
    resU_final[n,2]=resU[n,2]

# Post processing
# Plotting
PlotMethod.PlotArea(X, A, nMax, L, title="Chart Area")
PlotMethod.PlotResidue(numiter, resRho_final, resT_final, resV_final, resU_final, iters, title="Residual")
PlotMethod.PlotPrimitiveVariables(X, A, rho, T, V, L, nMax, title="Primitive Variables")
PlotMethod.PlotMDot(X, mDot, L, nMax, title="Mass Flow Comparison")

# Write results
CoreMethod.WriteResults(rho, T, V, p, mDot, M, X, A, nMax)
