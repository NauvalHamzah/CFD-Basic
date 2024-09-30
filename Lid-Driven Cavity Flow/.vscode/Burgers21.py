import numpy as np # Import numpy
import matplotlib.pyplot as plt # Import the plot module matplotlib
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
plt.style.use('classic')

N=21
xlist = np.linspace(0, 1, N)
ylist = np.linspace(0, 1, N)

#pembulatan linspace
for i in range (N):
    xlist[i] = round (xlist[i],2)
    ylist[i] = round (ylist[i],2)

delta_t=0.01
delta_x=0.05
X, Y = np.meshgrid(xlist, ylist)

#loop boundary condition
A = np.zeros ((N,N))
A[N-1,:]=1.5
A_new = np.zeros ((N,N))
for i in range (N):
    A[i,0]=1.5
    A[i,N-1]=-0.5
    A[0,i]=1.5-2*xlist[i]

A_new=A

#loop isi matriks A
for t in range (1):
    for i in range (N-2):
        for j in range (N-2):
            A_new[i+1,j+1]= A[i+1,j+1]-delta_t*((0.5*(((A[i+1,j+2])**2)-((A[i+1,j])**2))/(2*delta_x))+((A[i+2,j+1]-A[i,j+1])/(2*delta_x)))
    A=A_new

print (A)
print (A_new)

np.savetxt('Burgers21.out', A_new, delimiter= ' ')

"""fig,ax=plt.subplots(1,1)
cp = ax.contourf(X, Y, A, 50)
fig.colorbar(cp) # Add a colorbar to a plot
ax.set_title('Velocity u Contour')
ax.set_xlabel('x')
ax.set_ylabel('y')
plt.show()"""

"""fig = plt.figure()
ax = plt.axes(projection='3d')
surf = ax.plot_surface(X, Y, A, cmap='jet',linewidth=0, antialiased=False)
fig.colorbar(surf)
ax.set_title('Velocity u Contour')
plt.show()"""