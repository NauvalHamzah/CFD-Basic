import numpy
from matplotlib import pyplot, cm
from mpl_toolkits.mplot3d import Axes3D

##variable declarations
nx = 11
ny = 11
dx = 2 / (nx - 1)
dy = 2 / (ny - 1)
target = 0.01

##initial conditions
p = numpy.zeros((ny, nx))  # create a XxY vector of 0's
delta = numpy.zeros((ny, nx))

##plotting aids
x = numpy.linspace(0, 2, nx)
y = numpy.linspace(0, 1, ny)

##boundary conditions
p[:, 0] = 0  # p = 0 @ x = 0
p[:, -1] = y  # p = y @ x = 2
p[0, :] = p[1, :]  # dp/dy = 0 @ y = 0
p[-1, :] = p[-2, :]  # dp/dy = 0 @ y = 1

##iteration
n = 0
run = True

while run is True:
    pn=p.copy()
    error = 0
    p[1:-1, 1:-1] = ((dy**2 * (pn[1:-1, 2:] + pn[1:-1, 0:-2]) +
                         dx**2 * (pn[2:, 1:-1] + pn[0:-2, 1:-1])) /
                        (2 * (dx**2 + dy**2)))
              
    p[:, 0] = 0  # p = 0 @ x = 0
    p[:, -1] = y  # p = y @ x = 2
    p[0, :] = p[1, :]  # dp/dy = 0 @ y = 0
    p[-1, :] = p[-2, :]  # dp/dy = 0 @ y = 1
    
    for j in range (ny):
        for i in range (nx):
            delta[j,i]=abs(p[j,i]-pn[j,i])
    
    for j in range (ny):
        for i in range (nx):
            error=error+delta[j,i]
    if error<target:
        run = False
    n=n+1


fig = pyplot.figure(figsize=(11, 7), dpi=100)
ax = fig.gca(projection='3d')
X, Y = numpy.meshgrid(x, y)
print(Y)
'''
surf = ax.plot_surface(X, Y, p[:], rstride=1, cstride=1, cmap='jet',linewidth=0, antialiased=False)
ax.set_xlim(0, 2)
ax.set_ylim(0, 1)
ax.view_init(30, 225)
ax.set_xlabel('$x$')
ax.set_ylabel('$y$')
pyplot.show()



numpy.savetxt("X.txt",X)
numpy.savetxt("Y.txt",Y)
numpy.savetxt("Z.txt",p)
'''