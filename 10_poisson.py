import numpy
from matplotlib import pyplot, cm
from mpl_toolkits.mplot3d import Axes3D

##variable declarations
nx = 50
ny = 50
nt  = 20
xmin = 0
xmax = 2
ymin = 0
ymax = 1
target = 0.01

dx = (xmax - xmin) / (nx - 1)
dy = (ymax - ymin) / (ny - 1)


##initial conditions
p = numpy.zeros((ny, nx))  # create a XxY vector of 0's
delta = numpy.zeros((ny, nx))
b  = numpy.zeros((ny, nx))

##plotting aids
x = numpy.linspace(0, 2, nx)
y = numpy.linspace(0, 1, ny)

# Source
b[int(ny / 4), int(nx / 4)]  = 100
b[int(3 * ny / 4), int(3 * nx / 4)] = -100

##iteration
n = 0
run = True

while run is True:
    pn=p.copy()
    error = 0
    p[1:-1, 1:-1] = ((dy**2 * (p[1:-1, 2:] + p[1:-1, 0:-2]) +dx**2 * (p[2:, 1:-1] + p[0:-2, 1:-1])\
                             - b[1:-1,1:-1]*dx**2*dy**2) /(2 * (dx**2 + dy**2)))
              
    p[:, 0] = 0  # p = 0 @ x = 0
    p[:, -1] = 0  # p = y @ x = 2
    p[0, :] = 0  # dp/dy = 0 @ y = 0
    p[-1, :] = 0  # dp/dy = 0 @ y = 1
    
    for j in range (ny):
        for i in range (nx):
            delta[j,i]=abs(p[j,i]-pn[j,i])
    
    for j in range (ny):
        for i in range (nx):
            error=error+delta[j,i]
    if error<target:
        run = False
    n=n+1

print(n)
fig = pyplot.figure(figsize=(11, 7), dpi=100)
ax = fig.gca(projection='3d')
X, Y = numpy.meshgrid(x, y)
surf = ax.plot_surface(X, Y, p[:], rstride=1, cstride=1, cmap='jet',linewidth=0, antialiased=False)
ax.set_xlim(0, 2)
ax.set_ylim(0, 1)
ax.view_init(30, 225)
ax.set_xlabel('$x$')
ax.set_ylabel('$y$')
pyplot.show()