import math
import numpy as np

import Core.Airfoil as Foil
import Core.PlotMethod as PlotMethod

def Initialization(chord, gap, flap, rOutflow, rFarfield, wakeChord, thickness, nAirfoil, ngap, nflap, 
                        nAirfoilWake, nOutflow, nWake, iMax, jMax):
    
    # Status
    print("")
    print("Initializing...")

    # Initialize
    X = np.zeros(shape=(iMax,jMax))
    Y = np.zeros(shape=(iMax,jMax))

    # Initialize coordinates
    # point 1
    X[0,0] = rOutflow
    Y[0,0] = 0

    # point 2
    X[nWake-1,0] = chord+gap+flap
    Y[nWake-1,0] = 0

    # point 3
    X[nWake-1 + nflap-1,0] = chord+gap
    Y[nWake-1 + nflap-1,0] = 0

    # point 4
    X[nWake-1 + nflap-1 + ngap-1,0] = chord
    Y[nWake-1 + nflap-1 + ngap-1,0] = 0

    # point 5
    X[nWake-1 + nflap-1 + ngap-1 + nAirfoil-1 ,0] = 0
    Y[nWake-1 + nflap-1 + ngap-1 + nAirfoil-1 ,0] = 0

    # point 9
    X[iMax-1,0] = rOutflow
    Y[iMax-1,0] = 0

    # point 10
    X[iMax-1,jMax-1] = rOutflow
    Y[iMax-1,jMax-1] = rFarfield

    # point 11
    X[((iMax-1)-(nWake+nflap+ngap-3)),jMax-1] = chord
    Y[((iMax-1)-(nWake+nflap+ngap-3)),jMax-1] = rFarfield

    # point 12
    X[nWake+nflap+ngap-3,jMax-1] = chord
    Y[nWake+nflap+ngap-3,jMax-1] = -rFarfield

    # point 13
    X[0,jMax-1] = rOutflow
    Y[0,jMax-1] = -rFarfield

    # Initialize Boundary coordinates
    # Vertical left and right
    for j in range(1, jMax-1):
        eta = j/(jMax-1)
        m = rFarfield

        # Distribution
        A = 4   # exponential

        # Left
        X[0,j] = rOutflow
        #Y[0,j] = -m*eta     # linear
        Y[0,j] = -m*(math.exp(A*eta)-1)/(math.exp(A)-1) #exponential

        # Right
        X[iMax-1,j] = rOutflow
        #Y[iMax-1,j] = m*eta # linear
        Y[iMax-1,j] = m*(math.exp(A*eta)-1)/(math.exp(A)-1) # exponential

    # Top Boundary: Outflow
    for i in range(1, (nOutflow)):
        xi = i/(nOutflow-1)
        m = rOutflow-(chord+gap+flap+wakeChord)

        #Lower
        X[i,jMax-1] = rOutflow - m*xi   # linear
        #X[i,jMax-1] = rOutflow - m*math.sin(0.5*math.pi*xi)    # cosinus
        Y[i,jMax-1] = -rFarfield
        
        #Upper
        startIndex = iMax-1 
        X[startIndex-i,jMax-1] = rOutflow - m*xi
        Y[startIndex-i,jMax-1] = rFarfield

    # Top Boundary: wake
    for i in range(1, (nAirfoilWake)):
        xi = i/(nAirfoilWake-1)
        m = wakeChord

        #Lower
        X[nOutflow-1+i,jMax-1] = (chord+gap+flap+wakeChord) - m*xi   # linear
        #X[i,jMax-1] = rOutflow - m*math.sin(0.5*math.pi*xi)    # cosinus
        Y[nOutflow-1+i,jMax-1] = -rFarfield
        
        #Upper
        startIndex = (iMax-nOutflow) 
        X[startIndex-i,jMax-1] = (chord+gap+flap+wakeChord) - m*xi
        Y[startIndex-i,jMax-1] = rFarfield

    # Top Boundary: flap
    for i in range(1, (nflap)):
        xi = i/(nflap-1)
        m = flap

        # Lower
        X[nWake-1+i,jMax-1] = (chord+gap+flap) - m*xi   # linear
        #X[i,jMax-1] = rOutflow - m*math.sin(0.5*math.pi*xi)    # cosinus
        Y[nWake-1+i,jMax-1] = -rFarfield
        
        # Upper
        startIndex = (iMax-nWake) 
        X[startIndex-i,jMax-1] = (chord+gap+flap) - m*xi
        Y[startIndex-i,jMax-1] = rFarfield

    # Top Boundary: gap
    for i in range(1, (ngap)):
        xi = i/(ngap-1)
        m = gap

        # Lower
        X[nWake+nflap-2+i,jMax-1] = (chord+gap) - m*xi   # linear
        #X[i,jMax-1] = rOutflow - m*math.sin(0.5*math.pi*xi)    # cosinus
        Y[nWake+nflap-2+i,jMax-1] = -rFarfield
        
        # Upper
        startIndex = (iMax-(nWake+nflap-1))
        X[startIndex-i,jMax-1] = (chord+gap) - m*xi
        Y[startIndex-i,jMax-1] = rFarfield

    # Top Boundary: C-shaped
    for i in range(1, 2*(nAirfoil-1)):
        m = rFarfield
        xi = i/(2*(nAirfoil-1))   # linear
        #xi = 1-math.cos(0.5*math.pi*xi)     # cosinus

        # linear distribution
        startIndex = nWake+nflap+ngap-3
        X[startIndex+i,jMax-1] = m*math.sin(-math.pi*xi) + (chord)
        Y[startIndex+i,jMax-1] = -m*math.cos(-math.pi*xi)

    # Bottom Boundary: Outflow
    for i in range(1,nOutflow):
        xi = i/(nOutflow-1)
        m = rOutflow-(chord+gap+flap+wakeChord)

        X[i,0] = rOutflow - m*xi                            # linear
        #X[i,0] = rOutflow - m*math.sin(0.5*math.pi*xi)     # cosinus
        Y[i,0] = 0

        X[(iMax-1)-i,0] = X[i,0]                                    # linear
        #X[iMax-1)-i,0] = rOutflow - m*math.sin(0.5*math.pi*xi)     # cosinus
        Y[(iMax-1)-i,0] = Y[i,0]

    # Bottom Boundary: Airfoil wake
    for i in range(1,nAirfoilWake):
        xi = i/(nAirfoilWake-1)
        m = wakeChord 

        # Lower
        startIndex1 = nOutflow-1
        #X[startIndex+i,0] = chord + m*xi   # linear
        X[startIndex1+i,0] = chord + gap + flap + wakeChord - m*math.sin(0.5*math.pi*xi)    # cosinus
        Y[startIndex1+i,0] = 0
        
        # Upper
        startIndex2 = (iMax-nOutflow)
        X[startIndex2-i,0] = X[startIndex1+i,0]
        Y[startIndex2-i,0] = 0

    # Bottom Boundary: flap
    for i in range(1,nflap):
        xi = i/(nflap-1)
        m = flap

        # Lower
        startIndex1 = nWake-1
        #X[startIndex1+i,0] = chord + gap + flap - m*xi   # linear
        X[startIndex1+i,0] = chord + gap + flap - m*math.sin(0.5*math.pi*xi)    # cosinus
        Xrel = (X[startIndex1+i,0]- (chord + gap))/m
        Y[startIndex1+i,0] = Foil.NACA0012(Xrel, thickness, flap)

        # Upper
        startIndex2 = (iMax-nWake)
        X[startIndex2-i,0] = X[startIndex1+i,0]
        Y[startIndex2-i,0] = -Y[startIndex1+i,0]

    # Bottom Boundary: gap
    for i in range(1,ngap):
        xi = i/(ngap-1)
        m = gap

        # Lower
        startIndex1 = nWake + nflap - 2
        X[startIndex1+i,0] = chord + gap - m*xi   # linear
        #X[startIndex1+i,0] = chord + gap - m*math.sin(0.5*math.pi*xi)    # cosinus
        Y[startIndex1+i,0] = 0

        # Upper
        startIndex2 = (iMax- (nWake + nflap -1))
        X[startIndex2-i,0] = X[startIndex1+i,0]
        Y[startIndex2-i,0] = 0

    # Bottom Boundary: Airfoil
    for i in range(1,nAirfoil):
        xi = i/(nAirfoil-1)
        m = chord

        # Lower
        startIndex1 = nWake + nflap + ngap - 3
        #X[startIndex+i,0] = m - m*xi   # linear
        X[startIndex1+i,0] = chord - m*math.sin(0.5*math.pi*xi)    # cosinus
        Xrel = (X[startIndex1+i,0]-0)/chord
        Y[startIndex1+i,0] = Foil.NACA0012(Xrel, thickness, chord)

    for i in range(1,nAirfoil-1):
        # Upper
        startIndex1 = nWake + nflap + ngap - 3
        startIndex2 = (iMax-(nWake + nflap + ngap -2))
        X[startIndex2-i,0] = X[startIndex1+i,0]
        Y[startIndex2-i,0] = -Y[startIndex1+i,0]
    return (X, Y)

def BoundaryNormalization(X, Y, iMax, jMax):
    print("")
    print("Normalizing...") 

    # Normalization at boundary
    meshLength = np.zeros(shape=(iMax-1, jMax-1))
    u = np.zeros(shape=(iMax, jMax))
    v = np.zeros(shape=(iMax, jMax))

    # Bottom
    totalLength = 0
    for i in range(iMax-1):
        dx = X[i+1,0] - X[i,0]
        dy = Y[i+1,0] - Y[i,0]
        dLength = math.sqrt(dx**2 + dy**2)
        totalLength = totalLength + dLength
        meshLength[i,0] = totalLength

    for i in range(iMax-1):
        u[i+1,0] = meshLength[i,0]/totalLength
        v[i+1,0] = 0

    # Top
    totalLength = 0
    for i in range(iMax-1):
        dx = X[i+1,jMax-1] - X[i,jMax-1]
        dy = Y[i+1,jMax-1] - Y[i,jMax-1]
        dLength = math.sqrt(dx**2 + dy**2)
        totalLength = totalLength + dLength
        meshLength[i,jMax-2] = totalLength
        
    for i in range(iMax-1):
        u[i+1,jMax-1] = meshLength[i,jMax-2]/totalLength
        v[i+1,jMax-1] = 1

    # reset
    meshLength = np.zeros(shape=(iMax-1, jMax-1))

    # Left
    totalLength = 0
    for i in range(jMax-1):
        dx = X[0,i+1] - X[0,i]
        dy = Y[0,i+1] - Y[0,i]
        dLength = math.sqrt(dx**2 + dy**2)
        totalLength = totalLength + dLength
        meshLength[0,i] = totalLength

    for i in range(jMax-1):
        u[0,i+1] = 0
        v[0,i+1] = meshLength[0,i]/totalLength
        
    # Right
    totalLength = 0
    for i in range(jMax-1):
        dx = X[iMax-1,i+1] - X[iMax-1,i]
        dy = Y[iMax-1,i+1] - Y[iMax-1,i]
        dLength = math.sqrt(dx**2 + dy**2)
        totalLength = totalLength + dLength
        meshLength[iMax-2,i] = totalLength

    for i in range(jMax-1):
        u[iMax-1,i+1] = 1
        v[iMax-1,i+1] = meshLength[iMax-2,i]/totalLength

    return (u, v)

def BoundaryBlendedControlFunction(u, v, iMax, jMax):
    # Status
    print("")
    print("Blending...")
    
    # Boundary-Blended Control Functions
    for i in range(iMax-1):
        for j in range(jMax-1):
            u[i,j] = u[i,0] + v[0,j]*(u[i,jMax-1]-u[i,0])
            v[i,j] = v[0,j] + u[i,0]*(v[iMax-1,j]-v[0,j])
    return (u, v)

def TFI(X, Y, xi, eta, iMax, jMax):
    # Status
    print("")
    print("Meshing...")  
    
    # Transfinite Interpolation
    for i in range(1,iMax-1):
        for j in range(1,jMax-1):
            U = (1-xi[i,j])*X[0,j] + xi[i,j]*X[iMax-1,j]
            V = (1-eta[i,j])*X[i,0] + eta[i,j]*X[i,jMax-1]
            UV = xi[i,j]*eta[i,j]*X[iMax-1,jMax-1] + xi[i,j]*(1-eta[i,j])*X[iMax-1,0] +\
                (1-xi[i,j])*eta[i,j]*X[0,jMax-1] + (1-xi[i,j])*(1-eta[i,j])*X[0,0]
            X[i,j] = U + V - UV

            U = (1-xi[i,j])*Y[0,j] + xi[i,j]*Y[iMax-1,j]
            V = (1-eta[i,j])*Y[i,0] + eta[i,j]*Y[i,jMax-1]
            UV = xi[i,j]*eta[i,j]*Y[iMax-1,jMax-1] + xi[i,j]*(1-eta[i,j])*Y[iMax-1,0] +\
                (1-xi[i,j])*eta[i,j]*Y[0,jMax-1] + (1-xi[i,j])*(1-eta[i,j])*Y[0,0]
            Y[i,j] = U + V - UV

    return (X, Y)

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
    
    # Create bar chart of skewness
    data1 = []      # storing cell skewness < 0.2
    data2 = []      # storing cell skewness < 0.4
    data3 = []      # storing cell skewness < 0.6
    data4 = []      # storing cell skewness < 0.8
    data5 = []      # storing cell skewness > 0.8

    for i in range(iMax-1):
        for j in range(jMax-1):
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
    PlotMethod.plotBarQuality(skewnessData, title)

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

    while (error>targetError):
        for i in range(1,iMax-1):
            for j in range(1,jMax-1):
                xXi = (X[i+1,j]-X[i-1,j])/2
                yXi = (Y[i+1,j]-Y[i-1,j])/2
                xEta = (X[i,j+1]-X[i,j-1])/2
                yEta = (Y[i,j+1]-Y[i,j-1])/2
                #J = xXi*yEta - xEta*yXi

                alpha = xEta**2 + yEta**2
                beta = xXi*xEta + yXi*yEta
                gamma = xXi**2 + yXi**2

                # Finding Y
                Ry1 = alpha*(Y[i+1,j] - 2*Y[i,j] + Y[i-1,j])
                Ry2 = (-0.5)*beta*(Y[i+1,j+1] - Y[i-1,j+1] - Y[i+1,j-1] + Y[i-1,j-1])
                Ry3 = gamma*(Y[i,j+1] - 2*Y[i,j] + Y[i,j-1])
                Ry[i,j] = Ry1 + Ry2 + Ry3

                Y[i,j] = Y[i,j] + omega*((Ry[i,j])/(2*(alpha + gamma)))

                # Finding X
                if i == r or i == iMax-1 - r: continue

                Rx1 = alpha*(X[i+1,j] - 2*X[i,j] + X[i-1,j])
                Rx2 = (-0.5)*beta*(X[i+1,j+1] - X[i-1,j+1] - X[i+1,j-1] + X[i-1,j-1])
                Rx3 = gamma*(X[i,j+1] - 2*X[i,j] + X[i,j-1])
                Rx[i,j] = Rx1 + Rx2 + Rx3

                X[i,j] = X[i,j] + omega*((Rx[i,j])/(2*(alpha + gamma)))

        
        # Boundary Shifting 
        for i in range(1,iMax-1):
            if r <= i <= iMax-1 - r: continue 
            X[i,jMax-1] = X[i,jMax-2]
        
        for j in range(1,jMax-1):
            Y[0,j] = Y[1,j]
            Y[iMax-1,j] = Y[iMax-2,j]

        # Find residual
        currentRValue = np.sqrt(np.sum(Rx)**2 + np.sum(Ry)**2)
        error = abs(lastRValue - currentRValue)
        
        # Store residual
        iteration = iteration + 1
        
        # Other escape routes
        if (iteration>1000):
            break
        
        print("iteration #%4d residual = %.3e" % (iteration, error))

        residual.append(error*100)
        
        # Update value
        lastRValue = currentRValue

    return (X, Y, residual)

def NodesCoordinates(X, Y, iMax, jMax, nWake, nFlap, nGap, nAirfoil, thickness):
    # Create Basic Points
    nInnerPoints = iMax*jMax
    nGhostPoints = iMax + 2*jMax + 2*nFlap + 2*nAirfoil-1 
    nTotalPoints = nInnerPoints + nGhostPoints

    basicCoor = np.zeros(shape=(3, nTotalPoints))

    # Point internal mesh
    for j in range(jMax):
        for i in range(iMax):
            index = i + j*iMax

            # Store
            basicCoor[0, index] = index + 1
            basicCoor[1, index] = X[i,j]
            basicCoor[2, index] = Y[i,j]

    # Ghost Points:Left Boundary (Computational)
    Index_left = int(nInnerPoints/1000)*1000 + 1000
    for j in range(jMax):
        index = nInnerPoints + j
        basicCoor[0, index] = Index_left + 1 + j
        basicCoor[1, index] = X[0,j] + abs(X[0,j] - X[1,j])
        basicCoor[2, index] = Y[0,j]

    # Ghost Points : Outer Boundary
    #Bottom
    Index_out = int((Index_left+jMax)/1000)*1000 + 1000
    for i in range(nWake+nFlap+nGap-2):
        index = nInnerPoints + jMax + i
        basicCoor[0, index] = Index_out + 1 + i
        basicCoor[1, index] = X[i,jMax-1] 
        basicCoor[2, index] = Y[i,jMax-1] - abs(Y[i,jMax-1] - Y[i,jMax-2])
    
    # Ghost Points : Outer Boundary
    #C-Shaped
    for j in range(2*nAirfoil-1):
        index = nInnerPoints + jMax + i + j
        basicCoor[0, index] = Index_out + 1 + i + j
        basicCoor[1, index] = X[i+j,jMax-1] - abs(X[i+j,jMax-1]-X[i+j,jMax-2])
        if Y[i+j,jMax-1] < 0 :
            basicCoor[2, index] = Y[i+j,jMax-1] - abs(Y[i+j,jMax-1] - Y[i+j,jMax-2])
        else:
            basicCoor[2, index] = Y[i+j,jMax-1] + abs(Y[i+j,jMax-1] - Y[i+j,jMax-2])

    # Ghost Points : Outer Boundary
    #Top
    for k in range(nWake+nFlap+nGap-2):
        index = nInnerPoints + jMax + i + j + k
        basicCoor[0, index] = Index_out + 1 + i + j + k
        basicCoor[1, index] = X[(iMax-(nWake+nFlap+nGap-2))+k,jMax-1] 
        basicCoor[2, index] = Y[(iMax-(nWake+nFlap+nGap-2))+k,jMax-1] + abs(Y[(iMax-(nWake+nFlap+nGap-2))+k,jMax-1] - Y[(iMax-(nWake+nFlap+nGap-2))+k,jMax-2])

    # Ghost Points : Right Boundary (Computational)
    Index_right = int((Index_out+iMax)/1000)*1000 + 1000
    for j in range(1,jMax+1):
        index = nInnerPoints + jMax + iMax - 1 + j
        basicCoor[0, index] = Index_right + jMax + 1 - j
        basicCoor[1, index] = X[iMax-1,jMax-j] + abs(X[iMax-1,jMax-j] - X[iMax-2,jMax-j])
        basicCoor[2, index] = Y[iMax-1,jMax-j]

    # Ghost Points : AIrfoil-Flap Boundary
    # Flap 
    Index_solid = int((Index_right+jMax)/1000)*1000 + 1000
    for k in range(nFlap):
        index_1 = nInnerPoints + iMax + 2*(jMax) + k
        index_2 = nTotalPoints - nFlap + k
        basicCoor[0, index_1] = Index_solid + 1 + k
        basicCoor[0, index_2] = Index_solid + nFlap + 2*nAirfoil + k
        if k == 0 :
            basicCoor[1, index_1] = 0.5*(X[nWake,0]+X[nWake-1,0])
            basicCoor[2, index_1] = 0
            basicCoor[1, index_2] = 0.5*(X[iMax-(nWake+nFlap-2),0]+X[iMax-(nWake+nFlap-1),0])
            basicCoor[2, index_2] = 0
        elif k == nFlap-1:
            basicCoor[1, index_1] = 0.5*(X[nWake+nFlap-2,0]+X[nWake+nFlap-3,0])
            basicCoor[2, index_1] = 0
            basicCoor[1, index_2] = 0.5*(X[iMax-(nWake+1),0]+X[iMax-(nWake),0])
            basicCoor[2, index_2] = 0
        else :
            basicCoor[1, index_1] = X[nWake-1+k,0]
            basicCoor[2, index_1] = Y[nWake-1+k,0] + 0.2*0.05*thickness
            basicCoor[1, index_2] = X[iMax-(nWake+nFlap-1)+k,0]
            basicCoor[2, index_2] = Y[iMax-(nWake+nFlap-1)+k,0] - 0.2*0.05*thickness

    # Ghost Points : AIrfoil-Flap Boundary
    # Airfoil
    for k in range(2*nAirfoil-1):
        index = nInnerPoints + iMax + 2*(jMax) + nFlap + k
        basicCoor[0, index] = Index_solid + nFlap + 1 + k
        if k == 0 :
            basicCoor[1, index] = 0.5*(X[nWake+nFlap+nGap-3,0]+X[nWake+nFlap+nGap-2,0])
            basicCoor[2, index] = 0
        elif k == nAirfoil-1:
            basicCoor[1, index] = 0.5*(X[nWake+nFlap+nGap+nAirfoil-3,0]+X[nWake+nFlap+nGap+nAirfoil-4,0])
            basicCoor[2, index] = 0
        elif k == 2*nAirfoil-2:
            basicCoor[1, index] = 0.5*(X[nWake+nFlap+nGap+2*nAirfoil-5,0]+X[nWake+nFlap+nGap+2*nAirfoil-6,0])
            basicCoor[2, index] = 0
        else :
            if k < nAirfoil-1:
                basicCoor[1, index] = X[nWake+nFlap+nGap-3+k,0]
                basicCoor[2, index] = Y[nWake+nFlap+nGap-3+k,0] + 0.05*thickness
            else :
                basicCoor[1, index] = X[nWake+nFlap+nGap-3+k,0]
                basicCoor[2, index] = Y[nWake+nFlap+nGap-3+k,0] - 0.05*thickness

    return (basicCoor)

def CellNumber(X, Y, iMax, jMax, nAirfoil, nFlap):
    # Cell number
    nInnerCell = (iMax-1)*(jMax-1)
    nGhostCell = (iMax-1)+2*(jMax-1)+2*(nAirfoil+nFlap-2)
    nTotalCell = nInnerCell+nGhostCell

    cellNumber = np.zeros(shape=(nTotalCell,1))

    # inner cell
    for j in range(jMax-1):
        for i in range(iMax-1):
            index = i+j*(iMax-1) + 1
            cell = i+j*(iMax) + 1
            cellNumber[index-1,0] = cell

    # ghost cell: left
    Index_left = int(nInnerCell/1000)*1000 + 1000
    for j in range(jMax-1):
        # Bottom
        index = (nInnerCell) + j
        cellNumber[index,0] = Index_left+1+j

    # ghost cell: Outer
    Index_out = int((Index_left+jMax-1)/1000)*1000 + 1000
    for i in range(iMax-1):
        # left
        index = nInnerCell + (jMax-1) + i
        cellNumber[index,0] = Index_out + 1 + i

    # ghost cell: right
    Index_right = int((Index_out+iMax-1)/1000)*1000 + 1000
    for j in range(jMax-1):
        # Bottom
        index = (nInnerCell) + (jMax-1) + (iMax-1) + j
        cellNumber[index,0] = Index_right+(jMax-1)-j
    
    # ghost cell: Airfoil-Flap
    Index_solid = int((Index_right+jMax-1)/1000)*1000 + 1000
    for k in range(2*(nAirfoil+nFlap-2)):
        # Bottom
        index = nTotalCell - 2*(nAirfoil+nFlap-2) + k
        cellNumber[index,0] = Index_solid+1+k

    return (cellNumber)

def CellNeighbor(X, Y, cellNumber, iMax, jMax, nAirfoil, nFlap, nWake, nGap):
    # Cell neighbor
    nInnerPoints = iMax*jMax
    nGhostPoints = iMax + 2*jMax + 2*nFlap + 2*nAirfoil-1 
    nTotalPoints = nInnerPoints + nGhostPoints
    nInnerCell = (iMax-1)*(jMax-1)
    nGhostCell = (iMax-1)+2*(jMax-1)+2*(nAirfoil+nFlap-2)
    nTotalCell = nInnerCell+nGhostCell

    cellNeighboring = []

    # Cell neighbor
    # Internal: general
    n=0
    for j in range(jMax-1):
        for i in range(iMax-1):
            center = i+j*(iMax-1)
            
            if i == 0:
                left = int(cellNumber[nInnerCell+j])
                right = int(cellNumber[center+1])
            elif i == iMax-2:
                left = int(cellNumber[center-1])
                right = int(cellNumber[nInnerCell+2*(jMax-1)+(iMax-1)-1-j])
            else :
                left = int(cellNumber[center-1])
                right = int(cellNumber[center+1])
            
            if j == 0:
                top = int(cellNumber[center+(iMax-1)])
                if nWake-2<i<nWake+nFlap-2:
                    bottom = int(cellNumber[nTotalCell-(2*(nAirfoil+nFlap-2))+n])
                    n=n+1
                elif nWake+nFlap+nGap-4<i<nWake+nFlap+nGap+2*nAirfoil-5:
                    bottom = int(cellNumber[nTotalCell-(2*(nAirfoil+nFlap-2))+n])
                    n=n+1
                elif iMax-(nWake+nFlap)<i<iMax-nWake:
                    bottom = int(cellNumber[nTotalCell-(2*(nAirfoil+nFlap-2))+n])
                    n=n+1
                else:
                    bottom = int(cellNumber[(iMax-2)-center])
            elif j == jMax-2:
                bottom = int(cellNumber[center-(iMax-1)])
                top = int(cellNumber[center+(iMax-1)+(jMax-1)])
            else :
                bottom = int(cellNumber[center-(iMax-1)])
                top = int(cellNumber[center+(iMax-1)])

            cellNeighboring.append([left, bottom, right, top])
    
    #Ghost Cell : Left
    for j in range (jMax-1):
        center = nInnerCell + j

        left = -1
        right = int(cellNumber[j*(iMax-1)])
        bottom = int(cellNumber[center-1])
        top =  int(cellNumber[center+1])

        #correction
        if j==0:
            bottom = -1
        elif j==jMax-2:
            top = -1

        cellNeighboring.append([left, bottom, right, top])

    #Ghost Cell : Outer
    for i in range (iMax-1):
        center = nInnerCell + jMax-1 + i

        left = int(cellNumber[center-1])
        right = int(cellNumber[center+1])
        bottom = int(cellNumber[nInnerCell-(iMax-1) + i])
        top =  -1

        #correction
        if i==0:
            left = -1
        elif i==iMax-2:
            right = -1

        cellNeighboring.append([left, bottom, right, top])

    #Ghost Cell : Right
    for j in range (jMax-1):
        center = nInnerCell + (jMax-1) + (iMax-1) + j

        left = int(cellNumber[(jMax-1-j)*(iMax-1)-1])
        right = - 1
        bottom = int(cellNumber[center+1])
        top =  int(cellNumber[center-1])

        #correction
        if j==jMax-2:
            bottom = -1
        elif j==0:
            top = -1

        cellNeighboring.append([left, bottom, right, top])
    
    # Ghost Cell : Lower Flap
    for k in range (nFlap-1):
        center = nTotalCell-2*(nAirfoil+nFlap-2) + k

        left = int(cellNumber[center-1])
        right = int(cellNumber[center+1])
        bottom = -1
        top = int(cellNumber[nWake - 1 + k])

        #correction
        if k==0:
            left = -1
        elif k==nFlap-2:
            right = -1
        
        cellNeighboring.append([left, bottom, right, top])

    #Ghost Cell: Airfoil
    for k in range (2*(nAirfoil-1)):
        center = nTotalCell-2*(nAirfoil-1)-(nFlap-1) + k

        left = int(cellNumber[center-1])
        right = int(cellNumber[center+1])
        bottom = -1
        top = int(cellNumber[(nWake+nFlap+nGap-3) + k])

        #correction
        if k==0:
            left = -1
        elif k==2*(nAirfoil-1)-1:
            right = -1
        
        cellNeighboring.append([left, bottom, right, top])

    # Ghost Cell : Upper Flap
    for k in range (nFlap-1):
        center = nTotalCell - (nFlap-1) + k

        bottom = -1
        top = int(cellNumber[iMax-(nWake+nFlap-1) + k])

        if k==0:
            left = -1
            right = int(cellNumber[center+1])
        elif k==nFlap-2:
            left = int(cellNumber[center-1])
            right = -1
        else:
            left = int(cellNumber[center-1])
            right = int(cellNumber[center+1])
        
        cellNeighboring.append([left, bottom, right, top])

    return (cellNeighboring)

def CellNodalNumber(X, Y, Nodes, CellNumber, iMax, jMax, nAirfoil, nFlap, nWake, nGap):
    # Cell nodal number
    nInnerPoints = iMax*jMax
    nGhostPoints = iMax + 2*jMax + 2*nFlap + 2*nAirfoil-1 
    nTotalPoints = nInnerPoints + nGhostPoints
    nInnerCell = (iMax-1)*(jMax-1)
    nGhostCell = (iMax-1)+2*(jMax-1)+2*(nAirfoil+nFlap-2)
    nTotalCell = nInnerCell+nGhostCell

    cellNodalNumber = []

    # Cell Nodal Number
    # inner cell
    for j in range(jMax-1):
        for i in range(iMax-1):
            center = i+j*iMax

            bottomLeft = int(Nodes[0,center])
            bottomRight = int(Nodes[0,center+1])
            upperRight = int(Nodes[0,iMax+center+1])
            upperLeft = int(Nodes[0,iMax+center])
            
            cellNodalNumber.append([bottomLeft, bottomRight, upperRight, upperLeft])
    
    #Ghost Cell : Left
    for j in range(jMax-1):
        center = nInnerPoints + j

        bottomLeft = int(Nodes[0,center])
        bottomRight = int(Nodes[0,j*(iMax)])
        upperRight = int(Nodes[0,(j+1)*(iMax)])
        upperLeft = int(Nodes[0,center+1])
            
        cellNodalNumber.append([bottomLeft, bottomRight, upperRight, upperLeft])

    #Ghost Cell : Outer
    for i in range(iMax-1):
        upper = nInnerPoints + jMax + i
        bottom = nInnerPoints-iMax+i

        bottomLeft = int(Nodes[0,bottom])
        bottomRight = int(Nodes[0,bottom+1])
        upperRight = int(Nodes[0,upper+1])
        upperLeft = int(Nodes[0,upper])
            
        cellNodalNumber.append([bottomLeft, bottomRight, upperRight, upperLeft])

    #Ghost Cell : Right
    for j in range(jMax-1):
        center = nInnerPoints + (iMax+jMax) + j

        bottomLeft = int(Nodes[0,(jMax-1-j)*(iMax)-1])
        bottomRight = int(Nodes[0,center+1])
        upperRight = int(Nodes[0,center])
        upperLeft = int(Nodes[0,(jMax-j)*(iMax)-1])
            
        cellNodalNumber.append([bottomLeft, bottomRight, upperRight, upperLeft])

    # Ghost Cell : Lower Flap
    for k in range (nFlap-1):
        center = nTotalPoints-2*(nAirfoil+nFlap) + k

        bottomLeft = int(Nodes[0,center+1])
        bottomRight = int(Nodes[0,center+2])
        upperRight = int(Nodes[0,nWake+k])
        upperLeft = int(Nodes[0,nWake-1+k])
            
        cellNodalNumber.append([bottomLeft, bottomRight, upperRight, upperLeft])

    # Ghost Cell : Lower Airfoil
    for k in range (nAirfoil-1):
        center = nTotalPoints-2*(nAirfoil)-nFlap + k

        bottomLeft = int(Nodes[0,center+1])
        bottomRight = int(Nodes[0,center+2])
        upperRight = int(Nodes[0,(nWake+nFlap+nGap-2)+k])
        upperLeft = int(Nodes[0,(nWake+nFlap+nGap-3)+k])
            
        cellNodalNumber.append([bottomLeft, bottomRight, upperRight, upperLeft])
    
    # Ghost Cell : Upper Airfoil
    for k in range (nAirfoil-1):
        center = nTotalPoints-nAirfoil-nFlap-1 + k

        bottomLeft = int(Nodes[0,center+1])
        bottomRight = int(Nodes[0,center+2])
        upperRight = int(Nodes[0,(nWake+nFlap+nGap+nAirfoil-3)+k])
        upperLeft = int(Nodes[0,(nWake+nFlap+nGap+nAirfoil-4)+k])
            
        cellNodalNumber.append([bottomLeft, bottomRight, upperRight, upperLeft])

    # Ghost Cell : Upper Flap
    for k in range (nFlap-1):
        center = nTotalPoints-nFlap + k

        bottomLeft = int(Nodes[0,center])
        bottomRight = int(Nodes[0,center+1])
        upperRight = int(Nodes[0,iMax-(nWake+nFlap-2)+k])
        upperLeft = int(Nodes[0,iMax-(nWake+nFlap-1)+k])
            
        cellNodalNumber.append([bottomLeft, bottomRight, upperRight, upperLeft])

    return (cellNodalNumber)

def CellTypes(cellNodalNumber, iMax, jMax, nAirfoil, nFlap):
    nInnerCell = (iMax-1)*(jMax-1)
    nGhostCell = (iMax-1)+2*(jMax-1)+2*(nAirfoil+nFlap-2)
    nTotalCell = nInnerCell+nGhostCell

    cellType = np.zeros(shape=(nTotalCell,1))

    # Cell types
    for i in range(len(cellNodalNumber)):
        cellType[i,0] = int(len(cellNodalNumber[i]))

    return (cellType)

def BoundaryFlags(cellNumber, iMax, jMax, nAirfoil, nFlap, nWake, nGap):
    dataSolid = []
    dataInlet = []
    dataOutlet = []
    dataFarfield = []
    dataSymmetric = []
    dataPeriodic = []

    nInnerCell = (iMax-1)*(jMax-1)
    nGhostCell = (iMax-1)+2*(jMax-1)+2*(nAirfoil+nFlap-2)
    nTotalCell = nInnerCell+nGhostCell

    # for Outer are farfield and inlet
    # farfield: bottom
    for i in range(nWake+nFlap+nGap-3):
        index = nInnerCell + (jMax-1) + i
        dataFarfield.append(int(cellNumber[index]))

    # inlet: c-shaped
    for i in range(2*(nAirfoil-1)):
        index = nInnerCell + (jMax-1) + (nWake+nFlap+nGap-3) + i
        dataInlet.append(int(cellNumber[index]))

    # farfield: top
    for i in range(nWake+nFlap+nGap-3):
        index = nInnerCell + iMax - (nWake+nFlap+nGap-5) + i
        dataFarfield.append(int(cellNumber[index]))
            
    # for left and right | the whole left and right is outflow
    # Left
    for j in range(jMax-1):
        index = (iMax-1)*(jMax-1) + j
        dataOutlet.append(int(cellNumber[index]))

    # Right
    for j in range(jMax-1):
        index = (iMax)*(jMax) + 1 - j
        dataOutlet.append(int(cellNumber[index]))
    
    # for Airfoil and Flap
    for k in range(2*(nAirfoil+nFlap-2)):
        index = nTotalCell-2*(nAirfoil+nFlap-2) + k
        dataSolid.append(int(cellNumber[index]))


    dataFlags = [dataSolid, dataInlet, dataOutlet,
                    dataFarfield, dataSymmetric, dataPeriodic]

    return dataFlags

def WriteDataStructures(basicCoor, cellNumber, cellNeighboring, cellNodalNumber, cellType,boundaryFlags, iMax, jMax, nAirfoil, nFlap):

    nInnerPoints = iMax*jMax
    nGhostPoints = iMax + 2*jMax + 2*nFlap + 2*nAirfoil-1 
    nTotalPoints = nInnerPoints + nGhostPoints

    nInnerCell = (iMax-1)*(jMax-1)
    nGhostCell = (iMax-1)+2*(jMax-1)+2*(nAirfoil+nFlap-2)
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

