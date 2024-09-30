import math
import numpy as np
import matplotlib.pyplot as plt

import core.Airfoil as Foil
import core.Plot as Plot

def Initialization(chord, dInlet, dGap, dFlap, dOutflow, dFarfield, dWake, thickness, nInlet,
                   nAirfoil, nGap, nFlap, nWake, nOutflow, iMax, jMax, camber, camberPos):
    
    # Status
    print("")
    print("Initializing...")
    
    # Initialize
    X = np.zeros(shape=(iMax,2*jMax))
    Y = np.zeros(shape=(iMax,2*jMax))

    # Initialize coordinates
    # point A
    X[0, jMax-1] = -dInlet
    Y[0, jMax-1] = 0

    # point B
    X[nInlet-1, jMax-1] = 0
    Y[nInlet-1, jMax-1] = 0

    # point C
    X[nInlet-1 + nAirfoil-1, jMax-1] = chord
    Y[nInlet-1 + nAirfoil-1, jMax-1] = 0

    # point D
    X[nInlet-1 + nAirfoil-1 + nGap-1, jMax-1] = chord + dGap
    Y[nInlet-1 + nAirfoil-1 + nGap-1, jMax-1] = 0

    # point E
    X[nInlet-1 + nAirfoil-1 + nGap-1 + nFlap-1, jMax-1] = chord + dGap + dFlap
    Y[nInlet-1 + nAirfoil-1 + nGap-1 + nFlap-1, jMax-1] = 0

    # point F
    X[nInlet-1 + nAirfoil-1 + nGap-1 + nFlap-1 + nWake-1, jMax-1] = chord + dGap + dFlap + dWake
    Y[nInlet-1 + nAirfoil-1 + nGap-1 + nFlap-1 + nWake-1, jMax-1] = 0

    # point G
    X[iMax-1, jMax-1] = chord + dGap + dFlap + dWake + dOutflow
    Y[iMax-1, jMax-1] = 0

    # point H
    X[iMax-1, 0] = chord + dGap + dFlap + dWake + dOutflow
    Y[iMax-1, 0] = -dFarfield

    # point I
    X[nInlet-1 + nAirfoil-1 + nGap-1 + nFlap-1 + nWake-1, 0] = chord + dGap + dFlap + dWake
    Y[nInlet-1 + nAirfoil-1 + nGap-1 + nFlap-1 + nWake-1, 0] = -dFarfield

    # point J
    X[nInlet-1 + nAirfoil-1 + nGap-1 + nFlap-1, 0] = chord + dGap + dFlap
    Y[nInlet-1 + nAirfoil-1 + nGap-1 + nFlap-1, 0] = -dFarfield    
    
    # point K
    X[nInlet-1 + nAirfoil-1 + nGap-1, 0] = chord + dGap
    Y[nInlet-1 + nAirfoil-1 + nGap-1, 0] = -dFarfield

    # point L
    X[nInlet-1 + nAirfoil-1, 0] = chord
    Y[nInlet-1 + nAirfoil-1, 0] = -dFarfield

    # point M
    X[nInlet-1, 0] = 0
    Y[nInlet-1, 0] = -dFarfield

    # point N
    X[0,0] = -dInlet
    Y[0,0] = -dFarfield

    # point A'
    X[0, jMax] = -dInlet
    Y[0, jMax] = 0

    # point B'
    X[nInlet-1, jMax] = 0
    Y[nInlet-1, jMax] = 0

    # point C'
    X[nInlet-1 + nAirfoil-1, jMax] = chord
    Y[nInlet-1 + nAirfoil-1, jMax] = 0

    # point D'
    X[nInlet-1 + nAirfoil-1 + nGap-1, jMax] = chord + dGap
    Y[nInlet-1 + nAirfoil-1 + nGap-1, jMax] = 0

    # point E'
    X[nInlet-1 + nAirfoil-1 + nGap-1 + nFlap-1, jMax] = chord + dGap + dFlap
    Y[nInlet-1 + nAirfoil-1 + nGap-1 + nFlap-1, jMax] = 0

    # point F'
    X[nInlet-1 + nAirfoil-1 + nGap-1 + nFlap-1 + nWake-1, jMax] = chord + dGap + dFlap + dWake
    Y[nInlet-1 + nAirfoil-1 + nGap-1 + nFlap-1 + nWake-1, jMax] = 0

    # point G'
    X[iMax-1, jMax] = chord + dGap + dFlap + dWake + dOutflow
    Y[iMax-1, jMax] = 0

    # point H'
    X[iMax-1, 2*jMax-1] = chord + dGap + dFlap + dWake + dOutflow
    Y[iMax-1, 2*jMax-1] = dFarfield

    # point I'
    X[nInlet-1 + nAirfoil-1 + nGap-1 + nFlap-1 + nWake-1, 2*jMax-1] = chord + dGap + dFlap + dWake
    Y[nInlet-1 + nAirfoil-1 + nGap-1 + nFlap-1 + nWake-1, 2*jMax-1] = dFarfield

    # point J'
    X[nInlet-1 + nAirfoil-1 + nGap-1 + nFlap-1, 2*jMax-1] = chord + dGap + dFlap
    Y[nInlet-1 + nAirfoil-1 + nGap-1 + nFlap-1, 2*jMax-1] = dFarfield    
    
    # point K'
    X[nInlet-1 + nAirfoil-1 + nGap-1, 2*jMax-1] = chord + dGap
    Y[nInlet-1 + nAirfoil-1 + nGap-1, 2*jMax-1] = dFarfield

    # point L'
    X[nInlet-1 + nAirfoil-1, 2*jMax-1] = chord
    Y[nInlet-1 + nAirfoil-1, 2*jMax-1] = dFarfield

    # point M'
    X[nInlet-1, 2*jMax-1] = 0
    Y[nInlet-1, 2*jMax-1] = dFarfield

    # point N'
    X[0,2*jMax-1] = -dInlet
    Y[0,2*jMax-1] = dFarfield

    # Initialize Boundary coordinates
    # lower domain
    for j in range(1, jMax-1):
        eta = j/(jMax-1)
        m = dFarfield
        startIndex = jMax-1
        
        # Distribution
        A = 1.2   # exponential
        
        # N->A
        X[0,startIndex-j] = -dInlet
        #Y[0,startIndex-j] = -m*eta # linear
        Y[0,startIndex-j] = -m*(math.exp(A*eta)-1)/(math.exp(A)-1) # exponential

        # H->G
        X[iMax-1,startIndex-j] = chord + dGap + dFlap + dWake + dOutflow
        #Y[iMax-1,startIndex-j] = -m*eta # linear
        Y[iMax-1,startIndex-j] = -m*(math.exp(A*eta)-1)/(math.exp(A)-1) # exponential

    # upper domain
    for j in range(1, jMax-1):
        eta = j/(jMax-1)
        m = dFarfield
        startIndex = jMax
        
        # Distribution
        A = 1.2   # exponential

        # A'->N'
        X[0,startIndex+j] = -dInlet
        #Y[0,startIndex+j] = m*eta # linear
        Y[0,startIndex+j] = m*(math.exp(A*eta)-1)/(math.exp(A)-1) # exponential
        
        # G'->H'
        X[iMax-1,startIndex+j] = chord + dGap + dFlap + dWake + dOutflow
        #Y[iMax-1,startIndex+j] = m*eta # linear
        Y[iMax-1,startIndex+j] = m*(math.exp(A*eta)-1)/(math.exp(A)-1) # exponential


    # A->B, A'->B', N->M, and N'->M'
    for i in range(1, nInlet-1):
        xi = i/(nInlet-1)
        m = dInlet
        startIndex = nInlet-1

        # Distribution
        A = 1.2   # exponential
        
        # A'->B'
        #X[startIndex-i,jMax] = -m*xi # linear
        X[startIndex-i,jMax] = -m*(math.exp(A*xi)-1)/(math.exp(A)-1) # exponential
        Y[startIndex-i,jMax] = 0

        # A->B
        X[startIndex-i,jMax-1] = X[startIndex-i,jMax]
        Y[startIndex-i,jMax-1] = 0

        # N'->M'
        X[startIndex-i,2*jMax-1] = X[startIndex-i,jMax]
        Y[startIndex-i,2*jMax-1] = dFarfield

        # N->M
        X[startIndex-i,0] = X[startIndex-i,jMax]
        Y[startIndex-i,0] = -dFarfield

    # B->C, B'->C', M->L, and M'->L'
    for i in range(1, nAirfoil-1):
        xi = i/(nAirfoil-1)
        m = chord
        startIndex = nInlet-1 + nAirfoil-1
                
        # ---------------------------------------------------------------------------- #
        # SYMMETRIC
        # ---------------------------------------------------------------------------- #
        # B'->C'
        X[startIndex-i,jMax] = m - m*xi # linear
        #X[startIndex-i,jMax] = m - m*math.sin(0.5*math.pi*xi) # cosinus
        Xrel = X[startIndex-i,jMax]/m
        Y[startIndex-i,jMax] = Foil.NACA00TT(Xrel, thickness, m)

        # B->C
        X[startIndex-i,jMax-1] = X[startIndex-i,jMax]
        Y[startIndex-i,jMax-1] = -Y[startIndex-i,jMax]

        # M'->L'
        X[startIndex-i,2*jMax-1] = X[startIndex-i,jMax]
        Y[startIndex-i,2*jMax-1] = dFarfield               
        
        # M->L
        X[startIndex-i,0] = X[startIndex-i,jMax]
        Y[startIndex-i,0] = -dFarfield     

        # ---------------------------------------------------------------------------- #
        # ASYMMETRIC
        # ---------------------------------------------------------------------------- #
        #X[startIndex-i,jMax] = m - m*xi # linear
        '''
        X[startIndex-i,jMax] = m - m*math.sin(0.5*math.pi*xi) # cosinus
        Xrel = X[startIndex-i,jMax]/m
        Xu, Xl, Yu, Yl = Foil.NACAMPTT(Xrel, camber, camberPos, thickness, m)
        
        # B'->C'
        X[startIndex-i,jMax] = Xu
        Y[startIndex-i,jMax] = Yu

        # B->C
        X[startIndex-i,jMax-1] = Xl
        Y[startIndex-i,jMax-1] = Yl
       
        # M'->L'
        X[startIndex-i,2*jMax-1] = X[startIndex-i,jMax]
        Y[startIndex-i,2*jMax-1] = dFarfield               
        
        # M->L
        X[startIndex-i,0] = X[startIndex-i,jMax-1]
        Y[startIndex-i,0] = -dFarfield
        '''       

    # C->D, C'->D', L->K, and L'->K'
    for i in range(1, nGap-1):
        xi = i/(nGap-1)
        m = dGap
        startIndex = nInlet-1 + nAirfoil-1 + nGap-1

        # C'->D'
        #X[startIndex-i,jMax] = m - m*xi + chord # linear
        X[startIndex-i,jMax] = m - m*(math.exp(A*xi)-1)/(math.exp(A)-1) + chord # exponential
        Y[startIndex-i,jMax] = 0

        # C->D
        X[startIndex-i,jMax-1] = X[startIndex-i,jMax]
        Y[startIndex-i,jMax-1] = 0

        # L'->K'
        X[startIndex-i,2*jMax-1] = X[startIndex-i,jMax]
        Y[startIndex-i,2*jMax-1] = dFarfield

        # L->K
        X[startIndex-i,0] = X[startIndex-i,jMax]
        Y[startIndex-i,0] = -dFarfield

    # D->E, D'->E', K->J, and K'->J'
    for i in range(1, nFlap-1):
        xi = i/(nFlap-1)
        m = dFlap
        startIndex = nInlet-1 + nAirfoil-1 + nGap-1 + nFlap-1

        # ---------------------------------------------------------------------------- #
        # SYMMETRIC
        # ---------------------------------------------------------------------------- #
        # D'->E'
        X[startIndex-i,jMax] = m - m*xi + chord + dGap # linear
        #X[startIndex-i,jMax] = m - m*math.sin(0.5*math.pi*xi) + chord + dGap # cosinus
        Xrel = (X[startIndex-i,jMax] - chord - dGap)/m
        Y[startIndex-i,jMax] = Foil.NACA00TT(Xrel, thickness, m)

        # D->E
        X[startIndex-i,jMax-1] = X[startIndex-i,jMax]
        Y[startIndex-i,jMax-1] = -Y[startIndex-i,jMax]

        # K'->J'
        X[startIndex-i,2*jMax-1] = X[startIndex-i,jMax]
        Y[startIndex-i,2*jMax-1] = dFarfield               
        
        # K->J
        X[startIndex-i,0] = X[startIndex-i,jMax]
        Y[startIndex-i,0] = -dFarfield     

        # ---------------------------------------------------------------------------- #
        # ASYMMETRIC
        # ---------------------------------------------------------------------------- #
        #X[startIndex-i,jMax] = m - m*xi + chord + dGap # linear
        '''
        X[startIndex-i,jMax] = m - m*math.sin(0.5*math.pi*xi) # cosinus
        Xrel = X[startIndex-i,jMax]/m
        Xu, Xl, Yu, Yl = Foil.NACAMPTT(Xrel, camber, camberPos, thickness, m)
        
        # D'->E'
        X[startIndex-i,jMax] = Xu + chord + dGap
        Y[startIndex-i,jMax] = Yu

        # D->E
        X[startIndex-i,jMax-1] = Xl + chord + dGap
        Y[startIndex-i,jMax-1] = Yl
       
        # K'->J'
        X[startIndex-i,2*jMax-1] = X[startIndex-i,jMax]
        Y[startIndex-i,2*jMax-1] = dFarfield               
        
        # K->J
        X[startIndex-i,0] = X[startIndex-i,jMax-1]
        Y[startIndex-i,0] = -dFarfield
        '''

    # E->F, E'->F', J->I, and J'->I'
    for i in range(1, nWake-1):
        xi = i/(nWake-1)
        m = dWake
        startIndex = nInlet-1 + nAirfoil-1 + nGap-1 + nFlap-1
        
        # E'->F'
        X[startIndex+i,jMax] = m*xi + chord + dGap + dFlap # linear
        Y[startIndex+i,jMax] = 0

        # E->F
        X[startIndex+i,jMax-1] = X[startIndex+i,jMax]
        Y[startIndex+i,jMax-1] = 0

        # J'->I'
        X[startIndex+i,2*jMax-1] = X[startIndex+i,jMax]
        Y[startIndex+i,2*jMax-1] = dFarfield

        # J->I
        X[startIndex+i,0] = X[startIndex+i,jMax]
        Y[startIndex+i,0] = -dFarfield

    # F->G, F'->G', I->H, and I'->H'
    for i in range(1, nOutflow-1):
        xi = i/(nOutflow-1)
        m = dOutflow
        startIndex = nInlet-1 + nAirfoil-1 + nGap-1 + nFlap-1 + nWake-1

        # Distribution
        A = 1.2   # exponential
       
        # F'->G'
        #X[startIndex+i,jMax] = m*xi + chord + dGap + dFlap + dWake # linear
        X[startIndex+i,jMax] = m*(math.exp(A*xi)-1)/(math.exp(A)-1) + chord + dGap + dFlap + dWake #exponential
        Y[startIndex+i,jMax] = 0

        # E->F
        X[startIndex+i,jMax-1] = X[startIndex+i,jMax]
        Y[startIndex+i,jMax-1] = 0

        # J'->I'
        X[startIndex+i,2*jMax-1] = X[startIndex+i,jMax]
        Y[startIndex+i,2*jMax-1] = dFarfield

        # J->I
        X[startIndex+i,0] = X[startIndex+i,jMax]
        Y[startIndex+i,0] = -dFarfield

    return (X, Y)

def LaplaceSmoothing(X, Y, iMax, jMax, r, omega, targetError):
    
    # Status
    print("")
    print("Smoothing...")   
    
    # Laplace Smoothing
    Rx = np.zeros(shape=(iMax-1, jMax-1))
    Ry = np.zeros(shape=(iMax-1, jMax-1))
    
    iteration = 0
    error = 1
    lastRValue = 1
    residual = []

    jMax = int(jMax/2)

    while (error>targetError):
        for i in range(1,iMax-1):
            for j in range(1,2*jMax-1):
                if j == jMax-1 or j == jMax: continue
                
                xXi = (X[i+1,j]-X[i-1,j])/2
                yXi = (Y[i+1,j]-Y[i-1,j])/2
                xEta = (X[i,j+1]-X[i,j-1])/2
                yEta = (Y[i,j+1]-Y[i,j-1])/2
                #J = xXi*yEta - xEta*yXi

                alpha = xEta**2 + yEta**2
                beta = xXi*xEta + yXi*yEta
                gamma = xXi**2 + yXi**2

                # Finding Y
                #Ry1 = alpha*(Y[i+1,j] - 2*Y[i,j] + Y[i-1,j])
                #Ry2 = (-0.5)*beta*(Y[i+1,j+1] - Y[i-1,j+1] - Y[i+1,j-1] + Y[i-1,j-1])
                #Ry3 = gamma*(Y[i,j+1] - 2*Y[i,j] + Y[i,j-1])
                #Ry[i,j] = Ry1 + Ry2 + Ry3
                
                #Y[i,j] = Y[i,j] + omega*((Ry[i,j])/(2*(alpha + gamma)))

                # Finding X
                if i >= r: continue

                Rx1 = alpha*(X[i+1,j] - 2*X[i,j] + X[i-1,j])
                Rx2 = (-0.5)*beta*(X[i+1,j+1] - X[i-1,j+1] - X[i+1,j-1] + X[i-1,j-1])
                Rx3 = gamma*(X[i,j+1] - 2*X[i,j] + X[i,j-1])
                Rx[i,j] = Rx1 + Rx2 + Rx3

                X[i,j] = X[i,j] + omega*((Rx[i,j])/(2*(alpha + gamma)))
        
        # Boundary Shifting 
        for i in range(1,r):
            X[i,0] = X[i,1] #N->J
            X[i,2*jMax-1] = X[i,2*jMax-2] #N'->J'
        
        #for j in range(1,2*jMax-1):
            #if j == jMax-1 or j == jMax: continue
            #Y[0,j] = Y[1,j] #N->A and A'->N'
            #Y[iMax-1,j] = Y[iMax-2,j]  #H->G and G'->H'

        # Find residual
        currentRValue = np.sqrt(np.sum(Rx)**2 + np.sum(Ry)**2)
        error = abs(lastRValue - currentRValue)
        
        # Store residual
        iteration = iteration + 1
        
        # Other escape routes
        if (iteration>1000):
            break
        
        # Status
        print("iteration #%4d residual = %.3e" % (iteration, error))
        
        residual.append(error*100)
        
        # Update value
        lastRValue = currentRValue

    return (X, Y, residual)

def MeshQuality(X, Y, iMax, jMax, title):
    
    # Status
    print("")
    print("Checking Mesh Quality...")
    
    # Mesh Quality
    meshArea = np.zeros(shape=(iMax-1, jMax-1))
    meshSkewness = np.zeros(shape=(iMax-1, jMax-1))

    for i in range(iMax-1):
        for j in range(jMax-1):
            p = np.array([X[i,j+1] - X[i,j], Y[i,j+1] - Y[i,j]])
            q = np.array([X[i+1,j+1] - X[i,j+1], Y[i+1,j+1] - Y[i,j+1]])
            r = np.array([X[i+1,j+1] - X[i+1,j], Y[i+1,j+1] - Y[i+1,j]])
            s = np.array([X[i+1,j] - X[i,j], Y[i+1,j] - Y[i,j]])

            # Mesh Area
            area1 = -np.cross(p,s)
            area2 = np.cross(-q,-r)
            meshArea[i,j] = 0.5*(area1 + area2)

            # Skewness
            teta1 = math.degrees(math.atan2(-np.cross(p,s),np.dot(p,s)))
            teta2 = math.degrees(math.atan2(-np.cross(-s,r),np.dot(-s,r)))
            teta3 = math.degrees(math.atan2(-np.cross(-r,-q),np.dot(-r,-q)))
            teta4 = 360 - (teta1 + teta2 + teta3)
            teta = np.array([teta1, teta2, teta3, teta4])        
            tetaMaxOnMesh = (np.max(teta)-90)/90
            tetaMinOnMesh = (90-np.min(teta))/90
            skewness = max(np.array([tetaMaxOnMesh, tetaMinOnMesh]))
            meshSkewness[i,j] = skewness
            
    row = int((jMax-2)/2)
    meshArea = np.delete(meshArea,row,1)
    meshSkewness = np.delete(meshSkewness,row,1)
            
    # Create bar chart of skewness
    data1 = []      # storing cell skewness < 0.2
    data2 = []      # storing cell skewness < 0.4
    data3 = []      # storing cell skewness < 0.6
    data4 = []      # storing cell skewness < 0.8
    data5 = []      # storing cell skewness > 0.8

    for i in range(iMax-1):
        for j in range(jMax-1 - 1):
            if (meshSkewness[i,j] < 0.2):
                data1.append(meshSkewness[i,j])
            elif (meshSkewness[i,j] < 0.4):
                data2.append(meshSkewness[i,j])
            elif (meshSkewness[i,j] < 0.6):
                data3.append(meshSkewness[i,j])
            elif (meshSkewness[i,j] < 0.8):
                data4.append(meshSkewness[i,j])
            else:
                data5.append(meshSkewness[i,j])
    
    skewnessData = [len(data1), len(data2), len(data3), 
                    len(data4), len(data5)]
    Plot.plotBarQuality(skewnessData, title)

    print("")
    print("# --------------------------------------- #")
    print("             MESH QUALITY CHECK")
    print("# --------------------------------------- #")
    print("minimum mesh area = %.2e" % np.min(meshArea))
    print("maximum mesh area = %.4f" % np.max(meshArea))
    print("ortoghonality scale. 0 = Orthogonal")
    print("minimum value = %.2e" % np.min(meshSkewness))
    print("maximum value = %.4f" % np.max(meshSkewness))
    print("")

    return (skewnessData)

def plotGrid(X, Y, jMax, chartTitle ='', lineColor='b', lineWidth=1, 
             activatePoint=True, pointColor='r', pointSize=10):
    
    """
    plot: To plot structured grid

        plot(X, Y, lineColor, lineWidth, activatePoint, pointColor, pointSize)

        INPUT:
            X (matrix)      - matrix with x-coordinates of gridpoints
            Y (matrix)      - matrix with y-coordinates of gridpoints
            lineColor       - color of mesh lines. Default blue
            lineWidth       - width of mesh lines.
            activatePoint   - to activate mesh points
            pointColor      - mesh points color
            pointSize       - size of mesh points

    """

    xdiv, ydiv = X.shape    # extracting size, X and Y should be the same

    # Generating plot
    fig = plt.figure(figsize=(12, 9))
    ax = plt.subplot(111)

    # point scatter
    if (activatePoint==True):
        ax.scatter(X, Y, c=pointColor, s=pointSize)
    else:
        pass        
    
    # line plot
    for i in range(xdiv):
        ax.plot(X[i,:jMax], Y[i,:jMax], lineColor, lineWidth)
        ax.plot(X[i,jMax:], Y[i,jMax:], lineColor, lineWidth)
    
    for j in range(ydiv):
        ax.plot(X[:,j], Y[:,j], lineColor, lineWidth)

    ax.minorticks_on()
    ax.grid(b=True, which='major', color='0.65', linestyle='-')
    ax.grid(b=True, which='minor', color='0.65', linestyle='--')
    plt.title("2D Structured Grid Generation")
    plt.show()

    # Save figure
    fName = "results/" + chartTitle + ".png"
    fig.savefig(fName, dpi=150)

def TFI(X, Y, u, v, iMax, jMin, jMax):
    
    # Status
    print("")
    print("Meshing...")     
    
    # Transfinite Interpolation
    for i in range(1,iMax-1):
        for j in range(jMin+1,jMax-1):
            U = (1-u[i,j])*X[0,j] + u[i,j]*X[iMax-1,j]
            V = (1-v[i,j])*X[i,jMin] + v[i,j]*X[i,jMax-1]
            UV = u[i,j]*v[i,j]*X[iMax-1,jMax-1] + u[i,j]*(1-v[i,j])*X[iMax-1,jMin] +\
                (1-u[i,j])*v[i,j]*X[0,jMax-1] + (1-u[i,j])*(1-v[i,j])*X[0,jMin]
            X[i,j] = U + V - UV

            U = (1-u[i,j])*Y[0,j] + u[i,j]*Y[iMax-1,j]
            V = (1-v[i,j])*Y[i,jMin] + v[i,j]*Y[i,jMax-1]
            UV = u[i,j]*v[i,j]*Y[iMax-1,jMax-1] + u[i,j]*(1-v[i,j])*Y[iMax-1,jMin] +\
                (1-u[i,j])*v[i,j]*Y[0,jMax-1] + (1-u[i,j])*(1-v[i,j])*Y[0,jMin]
            Y[i,j] = U + V - UV

    return (X, Y)

def NodesCoordinates(X, Y, iMax, jMax, nAirfoil, nInlet, nFlap, nGap, thickness):
    # Create Basic Points
    nInnerPoints = iMax*(2*jMax)
    nGhostPoints = 2*(iMax + (2*jMax)) + 2*(nAirfoil+nFlap)
    nTotalPoints = nInnerPoints + nGhostPoints

    basicCoor = np.zeros(shape=(3, nTotalPoints))

    # Point internal mesh
    for j in range(2*jMax):
        for i in range(iMax):
            index = i + j*iMax

            # Store
            basicCoor[0, index] = index + 1
            basicCoor[1, index] = X[i,j]
            basicCoor[2, index] = Y[i,j]

    # Ghost Points
    #Bottom
    Index_bot = int(nInnerPoints/1000)*1000 + 1000
    for i in range(iMax):
        index = i+nInnerPoints
        basicCoor[0, index] = Index_bot + 1 + i
        basicCoor[1, index] = X[i,0]
        basicCoor[2, index] = Y[i,0] - abs(Y[i,0] - Y[i,1])

    #Right
    Index_right = int(Index_bot/1000)*1000 + 1000
    for j in range((2*jMax)):
        index = j + nInnerPoints + iMax
        basicCoor[0, index] = Index_right + 1 + j
        basicCoor[1, index] = X[iMax-1, j] + abs(X[iMax-1, j] - X[iMax-2, j])
        basicCoor[2, index] = Y[iMax-1, j]

    #Top
    Index_top = int(Index_right/1000)*1000 + 1000
    for i in range(iMax):
        index = i + nInnerPoints + iMax + 2*jMax
        basicCoor[0, index] = Index_top + iMax - i
        basicCoor[1, index] = X[iMax-1-i,(2*jMax)-1]
        basicCoor[2, index] = Y[iMax-1-i,(2*jMax)-1] + abs(Y[iMax-1-i,(2*jMax)-1] - Y[iMax-1-i,(2*jMax)-2])

    #Left
    Index_left = int(Index_top/1000)*1000 + 1000
    for j in range((2*jMax)):
        index = j + nInnerPoints + 2*(iMax) + 2*jMax 
        basicCoor[0, index] = Index_left + (2*jMax) - j
        basicCoor[1, index] = X[0, j] - abs(X[0, j] - X[1, j])
        basicCoor[2, index] = Y[0, j]

    #Lower - Upper Airfoil 
    Index_lower = int(Index_left/1000)*1000 + 1000
    Index_upper = int(Index_lower/1000)*1000 + 1000
    for k in range(nAirfoil):
        # Lower 
        index = nTotalPoints-(nAirfoil+2*nFlap)-1-k
        basicCoor[0, index] = Index_lower + (nAirfoil) - k
        basicCoor[1, index] = X[nInlet+nAirfoil-2-k,jMax-1]
        basicCoor[2, index] = Y[nInlet+nAirfoil-2-k,jMax-1] + 0.05*thickness

        # Upper 
        index = nTotalPoints-(nFlap)-1-k
        basicCoor[0, index] = Index_upper + (nAirfoil) - k
        basicCoor[1, index] = X[nInlet+nAirfoil-2-k,jMax]
        basicCoor[2, index] = Y[nInlet+nAirfoil-2-k,jMax] - 0.05*thickness

    #Lower - Upper Flap 
    for k in range(nFlap):
        # Lower 
        index = nTotalPoints-(nAirfoil + nFlap)-1-k
        basicCoor[0, index] = Index_lower + (nAirfoil+nFlap) - k
        basicCoor[1, index] = X[nInlet+(nAirfoil+nGap+nFlap-4)-k,jMax-1]
        basicCoor[2, index] = Y[nInlet+(nAirfoil+nGap+nFlap-4)-k,jMax-1] + 0.2*0.05*thickness

        # Upper 
        index = nTotalPoints-1-k
        basicCoor[0, index] = Index_upper + (nAirfoil+nFlap) - k
        basicCoor[1, index] = X[nInlet+(nAirfoil+nGap+nFlap-4)-k,jMax]
        basicCoor[2, index] = Y[nInlet+(nAirfoil+nGap+nFlap-4)-k,jMax] - 0.2*0.05*thickness

    #Leading Edge and Trailing Edge Correction
    #Airfoil
    basicCoor[1, nTotalPoints-2*(nAirfoil+nFlap)] = 0.5*(X[nInlet,jMax-1]+X[nInlet-1,jMax-1])
    basicCoor[1, nTotalPoints-(nAirfoil+2*nFlap)-1] = 0.5*(X[nInlet+nAirfoil-2,jMax-1]+X[nInlet+nAirfoil-3,jMax-1])
    basicCoor[1, nTotalPoints-(nAirfoil+nFlap)] = 0.5*(X[nInlet,jMax]+X[nInlet-1,jMax])
    basicCoor[1, nTotalPoints-(nFlap)-1] = 0.5*(X[nInlet+nAirfoil-2,jMax]+X[nInlet+nAirfoil-3,jMax])

    basicCoor[2, nTotalPoints-2*(nAirfoil+nFlap)] = 0           #Lower LE
    basicCoor[2, nTotalPoints-(nAirfoil+2*nFlap)-1] = 0         #Lower TE
    basicCoor[2, nTotalPoints-(nAirfoil+nFlap)] = 0             #upper LE
    basicCoor[2, nTotalPoints-(nFlap)-1] = 0                    #upper TE
    
    #Flap
    basicCoor[1, nTotalPoints-(nAirfoil+2*nFlap)] = 0.5*(X[nInlet+nAirfoil+nGap-3,jMax-1]+X[nInlet+nAirfoil+nGap-2,jMax-1])
    basicCoor[1, nTotalPoints-(nAirfoil+nFlap)-1] = 0.5*(X[nInlet+(nAirfoil+nGap+nFlap-5),jMax-1]+X[nInlet+(nAirfoil+nGap+nFlap-4),jMax-1])
    basicCoor[1, nTotalPoints-nFlap] = 0.5*(X[nInlet+nAirfoil+nGap-3,jMax]+X[nInlet+nAirfoil+nGap-2,jMax])
    basicCoor[1, nTotalPoints-1] = 0.5*(X[nInlet+(nAirfoil+nGap+nFlap-5),jMax]+X[nInlet+(nAirfoil+nGap+nFlap-4),jMax])

    basicCoor[2, nTotalPoints-(nAirfoil+2*nFlap)] = 0       #Lower LE
    basicCoor[2, nTotalPoints-(nAirfoil+nFlap)-1] = 0       #Lower TE
    basicCoor[2, nTotalPoints-nFlap] = 0                    #upper LE
    basicCoor[2, nTotalPoints-1] = 0                        #upper TE
    
    return (basicCoor)

def CellNumber(X, Y, iMax, jMax, nAirfoil, nFlap):
    # Cell number
    nInnerCell = (iMax-1)*(2*(jMax)-2)
    nGhostCell = 2*((iMax-1) + (2*(jMax)-2)) + 2*(nAirfoil+nFlap-2)
    nTotalCell = nInnerCell+nGhostCell

    cellNumber = np.zeros(shape=(nTotalCell,1))
    n=0
    # inner cell
    for j in range(0,2*(jMax)):
        for i in range(0,iMax-1):
            index = i+j*(iMax) + 1
            if j == jMax-1:
                continue
            else:
                cellNumber[n,0] = index
                n = n+1

    # ghost cell: Bottom
    Index_bot = int(nInnerCell/1000)*1000 + 1000
    for i in range(1,iMax):
        # Bottom
        index = i + Index_bot
        cellNumber[nInnerCell-1+i,0] = index

    # ghost cell: Right
    Index_right = int(Index_bot/1000)*1000 + 1000
    for i in range(1,2*(jMax)-1):
        index = i + Index_right 
        cellNumber[nInnerCell+(iMax-2)+i,0] = index

    # ghost cell: Top
    Index_top = int(Index_right/1000)*1000 + 1000
    for i in range(1,iMax):
        index = i + Index_top
        cellNumber[nInnerCell+2*(iMax-1)+2*(jMax)-2-i,0] = index

    # ghost cell: Left
    Index_left = int(Index_top/1000)*1000 + 1000
    for i in range(1,2*(jMax)-1):
        # left
        index = i + Index_left 
        cellNumber[nInnerCell+2*(iMax-1+2*(jMax)-2)-i,0] = index
    
    # ghost cell: Airfoil-Flap
    Index_lower = int(Index_left/1000)*1000 + 1000
    Index_upper = int(Index_lower/1000)*1000 + 1000
    for k in range(1,nAirfoil+nFlap-1):
        # Lower
        index = k + Index_lower 
        cellNumber[nTotalCell-2*(nAirfoil+nFlap)+3+k,0] = index

        # Upper
        index = k + Index_upper 
        cellNumber[nTotalCell-(nAirfoil+nFlap)+1+k,0] = index

    return (cellNumber)

def CellNeighbor(X, Y, cellNumber, iMax, jMax, nAirfoil, nFlap, nInlet, nGap):
    # Cell neighbor
    nInnerPoints = iMax*(2*jMax)
    nGhostPoints = 2*(iMax + (2*jMax)) + 2*(nAirfoil+nFlap)
    nTotalPoints = nInnerPoints + nGhostPoints
    nInnerCell = (iMax-1)*(2*(jMax)-2)
    nGhostCell = 2*((iMax-1) + (2*(jMax)-2)) + 2*(nAirfoil+nFlap-2)
    nTotalCell = nInnerCell+nGhostCell

    cellNeighboring = []

    # Cell neighbor
    # Internal: general
    for j in range((2*jMax)-2):
        for i in range(iMax-1):
            center = int(cellNumber[i+j*(iMax-1),0])

            left = cellNumber[i+j*(iMax-1)-1,0]
            bottom = cellNumber[i+j*(iMax-1)-(iMax-1),0]
            right = cellNumber[i+j*(iMax-1)+1,0]
            top = cellNumber[i+j*(iMax-1)+(iMax-1),0]

            cellNeighboring.append([left, bottom, right, top])

    # Internal: correction: bottom-top
    for i in range(iMax-1):
        # bottom
        cellNeighboring[i][1] = cellNumber[nInnerCell+i,0]
        # top
        cellNeighboring[nInnerCell-(iMax-1)+i][3] = cellNumber[nInnerCell + 2*(iMax-1) + ((2*jMax)-2) - 1 - i,0]

    # Internal: correction: left-right
    for j in range((2*jMax)-2):
        # left
        cellNeighboring[j*(iMax-1)][0] = cellNumber[nTotalCell-2*(nAirfoil+nFlap-2)-1-j,0]
        
        # right
        cellNeighboring[(j+1)*(iMax-1) - 1][2] = cellNumber[nInnerCell+iMax-1+j,0]

    # Internal: correction: Airfoil
    for k in range(nAirfoil-1):
        # Lower
        cellNeighboring[(jMax-2)*(iMax-1)+(nInlet-1)+k][3] = cellNumber[nTotalCell-2*(nAirfoil+nFlap-2)+k,0]
        
        # Upper
        cellNeighboring[(jMax-1)*(iMax-1)+(nInlet-1)+k][1] = cellNumber[nTotalCell-(nAirfoil+nFlap-2)+k,0]
    
    # Internal: correction: Flap
    for k in range(nFlap-1):
        # Lower
        cellNeighboring[(jMax-2)*(iMax-1)+(nInlet+nAirfoil+nGap-3)+k][3] = cellNumber[nTotalCell-(nAirfoil-1)-2*(nFlap-1)+k,0]
        
        # Upper
        cellNeighboring[(jMax-1)*(iMax-1)+(nInlet+nAirfoil+nGap-3)+k][1] = cellNumber[nTotalCell-(nFlap-1)+k,0]

    # Ghost: bottom
    for i in range(iMax-1):
        center = int(cellNumber[i+nInnerCell,0])
        
        left = center - 1
        bottom = -1
        right = center + 1
        top = i+1
        
        if (i==0):
            left = -1
        elif (i==iMax-2):
            right = -1
            
        cellNeighboring.append([left, bottom, right, top])

    # Ghost: right
    for j in range((2*jMax)-2):
        center = int(cellNumber[j+nInnerCell+(iMax-1),0])

        left = cellNumber[(j+1)*(iMax-1)-1,0]
        bottom = cellNumber[j+nInnerCell+(iMax-1)-1,0]
        right = -1
        top = cellNumber[j+nInnerCell+(iMax-1)+1,0]

        if (j==0):
            bottom = -1
        elif (j==((2*jMax)-3)):
            top = -1

        cellNeighboring.append([left, bottom, right, top])

    # Ghost: top
    for i in range(iMax-1):
        center = int(cellNumber[nInnerCell+(iMax-1)+((2*jMax)-2)+i,0])
        
        left = cellNumber[nInnerCell+(iMax-1)+((2*jMax)-2)+i+1,0]
        bottom = cellNumber[nInnerCell - 1 - i,0]
        right = cellNumber[nInnerCell+(iMax-1)+((2*jMax)-2)+i-1,0]
        top = -1
        
        if (i==0):
            right = -1
        elif (i==iMax-2):
            left = -1
            
        cellNeighboring.append([left, bottom, right, top])

    # Ghost: left
    for j in range((2*jMax)-2):
        center = int(cellNumber[nInnerCell+2*(iMax-1)+((2*jMax)-2)+j,0])

        left = -1
        bottom = cellNumber[nInnerCell+2*(iMax-1)+((2*jMax)-2)+j+1,0]
        right = cellNumber[((2*jMax)-3-j)*(iMax-1),0]
        top = cellNumber[nInnerCell+2*(iMax-1)+((2*jMax)-2)+j-1,0]

        if (j==0):
            top = -1
        elif (j==((2*jMax)-3)):
            bottom = -1

        cellNeighboring.append([left, bottom, right, top])

    # Ghost: lower airfoil
    for k in range(nAirfoil-1):
        center = int(cellNumber[nTotalCell-2*(nAirfoil+nFlap-2)+k,0])

        left = cellNumber[nTotalCell-2*(nAirfoil+nFlap-2)+k-1,0]
        bottom = cellNumber[(jMax-2)*(iMax-1)+(nInlet-1)+k,0]
        right = cellNumber[nTotalCell-2*(nAirfoil+nFlap-2)+k+1,0]
        top = -1

        if (k==0):
            left = -1
        elif (k==(nAirfoil-2)):
            right = -1

        cellNeighboring.append([left, bottom, right, top])

    # Ghost: lower Flap
    for k in range(nFlap-1):
        center = int(cellNumber[nTotalCell-(nAirfoil-1)-2*(nFlap-1)+k,0])

        left = cellNumber[nTotalCell-(nAirfoil-1)-2*(nFlap-1)+k-1,0]
        bottom = cellNumber[(jMax-2)*(iMax-1)+(nInlet+nAirfoil+nGap-3)+k,0]
        right = cellNumber[nTotalCell-(nAirfoil-1)-2*(nFlap-1)+k+1,0]
        top = -1

        if (k==0):
            left = -1
        elif (k==(nFlap-2)):
            right = -1

        cellNeighboring.append([left, bottom, right, top])

    # Ghost: Upper airfoil
    for k in range(nAirfoil-1):
        center = int(cellNumber[nTotalCell-(nAirfoil+nFlap-2)+k,0])

        left = cellNumber[nTotalCell-(nAirfoil+nFlap-2)+k-1,0]
        bottom = -1
        right = cellNumber[nTotalCell-(nAirfoil+nFlap-2)+k+1,0]
        top = cellNumber[(jMax-1)*(iMax-1)+(nInlet-1)+k,0]

        if (k==0):
            left = -1
        elif (k==(nAirfoil-2)):
            right = -1

        cellNeighboring.append([left, bottom, right, top]) 
 
  
    # Ghost: upper Flap
    for k in range(nFlap-1):
        center = int(cellNumber[nTotalCell-(nFlap-1)+k,0])

        if (k==0):
            left = -1
            bottom = -1
            right = cellNumber[nTotalCell-(nFlap-1)+k+1,0]
            top = cellNumber[(jMax-1)*(iMax-1)+(nInlet+nAirfoil+nGap-3)+k,0]
        elif (k==(nFlap-2)):
            left = cellNumber[nTotalCell-(nFlap-1)+k-1,0]
            bottom = -1
            right = -1
            top = cellNumber[(jMax-1)*(iMax-1)+(nInlet+nAirfoil+nGap-3)+k,0]
        else:
            left = cellNumber[nTotalCell-(nFlap-1)+k-1,0]
            bottom = -1
            right = cellNumber[nTotalCell-(nFlap-1)+k+1,0]
            top = cellNumber[(jMax-1)*(iMax-1)+(nInlet+nAirfoil+nGap-3)+k,0]

        cellNeighboring.append([left, bottom, right, top])

    for i in range (len(cellNeighboring)):
        for j in range (4):
            cellNeighboring [i][j] = int (cellNeighboring [i][j])

    return (cellNeighboring)

def CellNodalNumber(X, Y, Nodes, cellNumber, iMax, jMax, nAirfoil, nFlap, nInlet, nGap):
    # Cell nodal number
    nInnerPoints = iMax*(2*jMax)
    nGhostPoints = 2*(iMax + (2*jMax)) + 2*(nAirfoil+nFlap)
    nTotalPoints = nInnerPoints + nGhostPoints
    nInnerCell = (iMax-1)*(2*(jMax)-2)
    nGhostCell = 2*((iMax-1) + (2*(jMax)-2)) + 2*(nAirfoil+nFlap-2)
    nTotalCell = nInnerCell+nGhostCell

    cellNodalNumber = []

    # Cell Nodal Number
    # inner cell
    for j in range(0,(2*jMax)-1):
        for i in range(0,iMax-1):
            if j == jMax-1:
                continue
            else:
                bottomLeft = Nodes[0,i+j*iMax]
                bottomRight = Nodes[0,i+j*iMax+1]
                topRight = Nodes[0,i+j*iMax+1+iMax]
                topLeft = Nodes[0,i+j*iMax+iMax]
            
                cellNodalNumber.append([bottomLeft, bottomRight, topRight, topLeft])
    
    # ghost cell: bottom
    for i in range(0,iMax-1):
        index = nInnerPoints + i
        
        bottomLeft = Nodes[0,index]
        bottomRight = Nodes[0,index+1]
        topRight = Nodes[0,index+1 - nInnerPoints]
        topLeft = Nodes[0,index - nInnerPoints]

        cellNodalNumber.append([bottomLeft, bottomRight, topRight, topLeft])

    # ghost cell: right
    for j in range(0,(2*jMax)-1):
        if j == jMax-1:
            continue
        else:    
            bottomLeft = Nodes[0,iMax*(j+1)-1]
            bottomRight = Nodes[0,nInnerPoints + iMax + j]
            topRight = Nodes[0,nInnerPoints + iMax + j + 1]
            topLeft = Nodes[0,iMax*(j+2)-1]

            cellNodalNumber.append([bottomLeft, bottomRight, topRight, topLeft])

    # ghost cell: top
    for i in range(0,iMax-1):
        index = nInnerPoints + iMax + (2*jMax) + i
        
        bottomLeft = Nodes[0,iMax*(2*jMax) - i - 2]
        bottomRight = Nodes[0,iMax*(2*jMax) - i - 1]
        topRight = Nodes[0,index]
        topLeft = Nodes[0,index+1]

        cellNodalNumber.append([bottomLeft, bottomRight, topRight, topLeft])

    # ghost cell: left
    for j in range(0,(2*jMax)-1):
        if j==jMax-1:
            continue
        else:    
            bottomLeft = Nodes[0,nInnerPoints + 2*iMax + (2*jMax) + j + 1]
            bottomRight = Nodes[0,iMax*((2*jMax)-2-j)]  #nInnerPoints + iMax + j + 1
            topRight = Nodes[0,iMax*((2*jMax)-2-j)+ iMax]
            topLeft = Nodes[0,nInnerPoints + 2*iMax + (2*jMax) + j]

            cellNodalNumber.append([bottomLeft, bottomRight, topRight, topLeft])

    # ghost cell: lower airfoil
    for k in range(0,nAirfoil-1):    
        bottomLeft = Nodes[0,iMax*(jMax-1)+nInlet+k-1]
        bottomRight = Nodes[0,iMax*(jMax-1)+nInlet+k]  
        topRight = Nodes[0,nInnerPoints + 2*(iMax + (2*jMax)) + k +1]
        topLeft = Nodes[0,nInnerPoints + 2*(iMax + (2*jMax)) + k ]

        cellNodalNumber.append([bottomLeft, bottomRight, topRight, topLeft])

    # ghost cell: lower Flap
    for k in range(0,nFlap-1):    
        bottomLeft = Nodes[0,iMax*(jMax-1)+(nInlet+nAirfoil+nGap-3)+k]
        bottomRight = Nodes[0,iMax*(jMax-1)+(nInlet+nAirfoil+nGap-2)+k]  
        topRight = Nodes[0,nInnerPoints + 2*(iMax + (2*jMax)) + (nAirfoil) + k +1]
        topLeft = Nodes[0,nInnerPoints + 2*(iMax + (2*jMax)) + (nAirfoil) + k ]

        cellNodalNumber.append([bottomLeft, bottomRight, topRight, topLeft])
    
    # ghost cell: upper airfoil
    for k in range(0,nAirfoil-1):    
        bottomLeft = Nodes[0,nInnerPoints + 2*(iMax + (2*jMax)) + (nAirfoil+nFlap) + k ]
        bottomRight = Nodes[0,nInnerPoints + 2*(iMax + (2*jMax)) + (nAirfoil+nFlap) + k + 1]  #nInnerPoints + iMax + j + 1
        topRight = Nodes[0,iMax*(jMax)+nInlet+k]
        topLeft = Nodes[0,iMax*(jMax)+nInlet+k-1]

        cellNodalNumber.append([bottomLeft, bottomRight, topRight, topLeft])
    
    # ghost cell: upper flap
    for k in range(0,nFlap-1):    
        bottomLeft = Nodes[0,nTotalPoints - nFlap + k ]
        bottomRight = Nodes[0,nTotalPoints - nFlap + k + 1]  #nInnerPoints + iMax + j + 1
        topRight = Nodes[0,iMax*(jMax)+(nInlet+nAirfoil+nGap-2)+k]
        topLeft = Nodes[0,iMax*(jMax)+(nInlet+nAirfoil+nGap-3)+k]

        cellNodalNumber.append([bottomLeft, bottomRight, topRight, topLeft])

    for i in range (len(cellNodalNumber)):
        for j in range (4):
            cellNodalNumber [i][j] = int (cellNodalNumber [i][j])

    return (cellNodalNumber)

def CellTypes(cellNodalNumber, iMax, jMax, nAirfoil, nFlap):
    nInnerCell = (iMax-1)*(2*(jMax)-2)
    nGhostCell = 2*((iMax-1) + (2*(jMax)-2)) + 2*(nAirfoil+nFlap-2)
    nTotalCell = nInnerCell+nGhostCell

    cellType = np.zeros(shape=(nTotalCell,1))

    # Cell types
    for i in range(len(cellNodalNumber)):
        cellType[i,0] = int(len(cellNodalNumber[i]))

    return (cellType)

def BoundaryFlags(iMax, jMax, cellNumber, nAirfoil, nFlap, nInlet):
    nInnerCell = (iMax-1)*(2*(jMax)-2)
    nGhostCell = 2*((iMax-1) + (2*(jMax)-2)) + 2*(nAirfoil+nFlap-2)
    nTotalCell = nInnerCell+nGhostCell
    dataSolid = []
    dataInlet = []
    dataOutlet = []
    dataFarfield = []
    dataSymmetric = []
    dataPeriodic = []

    # for top | the tops are farfield
    # farfield: bottom
    for i in range(iMax-1):
        index = nInnerCell + i
        dataFarfield.append(int(cellNumber[index,0]))

    # farfield: top
    for i in range(iMax-1):
        index = nInnerCell + 2*(iMax-1) + ((2*jMax)-2) - 1 - i
        dataFarfield.append(int(cellNumber[index,0]))

    # Left is Inlet
    for j in range((2*jMax)-2):
        index = nTotalCell - 2*(nAirfoil+nFlap-2) - 1 - j
        dataInlet.append(int(cellNumber[index,0]))

    # Right is Outflow
    for j in range((2*jMax)-2):
        index = nInnerCell + iMax-1 + j
        dataOutlet.append(int(cellNumber[index,0]))

    # For Lower and upper are Solid
    #Lower
    for k in range(nAirfoil+nFlap-2):
        index = nTotalCell - 2*(nAirfoil+nFlap-2) + k
        dataSolid.append(int(cellNumber[index,0]))
    
    #Upper
    for k in range(nAirfoil+nFlap-2):
        index = nTotalCell - (nAirfoil+nFlap-2) + k
        dataSolid.append(int(cellNumber[index,0]))

    dataFlags = [dataSolid, dataInlet, dataOutlet,
                    dataFarfield, dataSymmetric, dataPeriodic]

    return dataFlags

def WriteDataStructures(basicCoor, cellNumber, cellNeighboring, cellNodalNumber, cellType,boundaryFlags, iMax, jMax, nAirfoil, nFlap):

    nInnerPoints = iMax*(2*jMax)
    nGhostPoints = 2*(iMax + (2*jMax)) + 2*(nAirfoil+nFlap)
    nTotalPoints = nInnerPoints + nGhostPoints

    nInnerCell = (iMax-1)*(2*(jMax)-2)
    nGhostCell = 2*((iMax-1) + (2*(jMax)-2)) + 2*(nAirfoil+nFlap-2)
    nTotalCell = nInnerCell+nGhostCell

    # WRITE POINT COORDINATES
    f = open("results/nodes-coordinates", "w+")
    f.write("/*---------------------------------------------------------------------------*\\")
    f.write("\nDescription\n")
    f.write("{\n")
    f.write("\tobject\t\t\t\tnodes-coordinates;\n")
    f.write("\tnumber-of-nodes\t\t%d;\n" % nTotalPoints)
    f.write("\tinternal-nodes\t\t%d;\n" %nInnerPoints)
    f.write("\tghost-nodes\t\t\t%d;\n" %nGhostPoints)
    f.write("}\n")
    f.write("\*---------------------------------------------------------------------------*/")
    f.write("\n\n")
    f.write("%d\n" % nTotalPoints)
    f.write("(\n")
    for i in range(nTotalPoints):
        f.write("\t(%d, " % basicCoor[0, i])
        coor = np.array([round(basicCoor[1,i],4), round(basicCoor[2,i],4)])
        f.write("%s)\n " % coor)
    f.write(")")
    f.close()

    # WRITE CELL CONNECTIVITY
    f = open("results/cells-connectivity", "w+")
    f.write("/*---------------------------------------------------------------------------*\\")
    f.write("\nDescription\n")
    f.write("{\n")
    f.write("\tobject\t\t\t\tcells-connectivity;\n")
    f.write("\tnumber-of-cells\t\t%d;\n" %nTotalCell)
    f.write("\tinternal-cells\t\t%d;\n" %nInnerCell)
    f.write("\tghost-cells\t\t\t%d;\n" %nGhostCell)
    f.write("}\n")
    f.write("\*---------------------------------------------------------------------------*/")
    f.write("\n\n")
    f.write("%d\n" % nInnerCell)
    f.write("(\n")
    for i in range(nInnerCell):
        f.write("\t(%d, " % (cellNumber[i]))
        f.write("%s, " % cellNeighboring[i])
        f.write("%s, " % cellNodalNumber[i])
        f.write("%d)\n" % cellType[i])
    f.write(")")
    f.close()

    # WRITE BOUNDARY FLAGS
    f = open("results/boundary-flags", "w+")
    f.write("/*---------------------------------------------------------------------------*\\")
    f.write("\nDescription\n")
    f.write("{\n")
    f.write("\tobject\t\tboundary-flags;\n")
    f.write("}\n")
    f.write("\*---------------------------------------------------------------------------*/")
    f.write("\n\n")
    f.write("%d\n" % len(boundaryFlags))
    f.write("(\n")
    f.write("\tsolid\t\t1\t: %s\n"% boundaryFlags[0])
    f.write("\tinlet\t\t2\t: %s\n"% boundaryFlags[1])
    f.write("\toutlet\t\t3\t: %s\n"% boundaryFlags[2])
    f.write("\tfarfield\t4\t: %s\n"% boundaryFlags[3])
    f.write("\tsymmetry\t5\t: %s\n"% boundaryFlags[4])
    f.write("\tperiodic\t6\t: %s\n"% boundaryFlags[5])
    f.write(")")
    f.close()

