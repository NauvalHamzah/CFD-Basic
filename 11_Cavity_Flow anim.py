import numpy
from matplotlib import pyplot, cm
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation

def build_up_b(b, rho, dt, u, v, dx, dy, i):
    b[1:-1, 1:-1] = (rho * (1 / dt * 
                    ((u[1:-1, 2:, i] - u[1:-1, 0:-2, i]) / 
                     (2 * dx) + (v[2:, 1:-1, i] - v[0:-2, 1:-1, i]) / (2 * dy)) -
                    ((u[1:-1, 2:, i] - u[1:-1, 0:-2, i]) / (2 * dx))**2 -
                      2 * ((u[2:, 1:-1, i] - u[0:-2, 1:-1, i]) / (2 * dy) *
                           (v[1:-1, 2:, i] - v[1:-1, 0:-2, i]) / (2 * dx))-
                          ((v[2:, 1:-1, i] - v[0:-2, 1:-1, i]) / (2 * dy))**2))
    
    return b

def pressure_poisson(p,dx,dy,b,j):
    for i in range (nit-1):
        p[1:-1, 1:-1, j] = (((p[1:-1, 2:, j-1] + p[1:-1, 0:-2, j-1]) * dy**2 + 
                          (p[2:, 1:-1, j-1] + p[0:-2, 1:-1, j-1]) * dx**2) /
                          (2 * (dx**2 + dy**2)) -
                          dx**2 * dy**2 / (2 * (dx**2 + dy**2)) * 
                          b[1:-1,1:-1])

        p[:, -1, j] = p[:, -2, j] # dp/dx = 0 at x = 2
        p[0, :, j] = p[1, :,j]   # dp/dy = 0 at y = 0
        p[:, 0, j] = p[:, 1, j]   # dp/dx = 0 at x = 0
        p[-1, :, j] = 0        # p = 0 at y = 2
    
    return p

def cavity_flow(nt, u, v, dt, dx, dy, p, rho, nu):
    b = numpy.zeros((ny, nx))    
    for i in range (nt-1):
        b = build_up_b(b, rho, dt, u, v, dx, dy, i)
        p = pressure_poisson(p, dx, dy, b, i)
        u[1:-1, 1:-1, i+1] = (u[1:-1, 1:-1, i]-
                         u[1:-1, 1:-1, i] * dt / dx *
                        (u[1:-1, 1:-1, i] - u[1:-1, 0:-2, i]) -
                         v[1:-1, 1:-1, i] * dt / dy *
                        (u[1:-1, 1:-1, i] - u[0:-2, 1:-1, i]) -
                         dt / (2 * rho * dx) * (p[1:-1, 2:, i] - p[1:-1, 0:-2, i]) +
                         nu * (dt / dx**2 *
                        (u[1:-1, 2:, i] - 2 * u[1:-1, 1:-1, i] + u[1:-1, 0:-2, i]) +
                         dt / dy**2 *
                        (u[2:, 1:-1, i] - 2 * u[1:-1, 1:-1, i] + u[0:-2, 1:-1, i])))
                
        v[1:-1,1:-1, i+1] = (v[1:-1, 1:-1, i] -
                        u[1:-1, 1:-1, i] * dt / dx *
                       (v[1:-1, 1:-1, i] - v[1:-1, 0:-2, i]) -
                        v[1:-1, 1:-1, i] * dt / dy *
                       (v[1:-1, 1:-1, i] - v[0:-2, 1:-1, i]) -
                        dt / (2 * rho * dy) * (p[2:, 1:-1, i] - p[0:-2, 1:-1, i]) +
                        nu * (dt / dx**2 *
                       (v[1:-1, 2:, i] - 2 * v[1:-1, 1:-1, i] + v[1:-1, 0:-2, i]) +
                        dt / dy**2 *
                        (v[2:, 1:-1, i] - 2 * v[1:-1, 1:-1, i] + v[0:-2, 1:-1, i])))
        u[0, :, i+1]  = 0
        u[:, 0, i+1]  = 0
        u[:, nx-1, i+1] = 0
        u[ny-1, :, i+1] = 1    # set velocity on cavity lid equal to 1
        v[0, :, i+1]  = 0
        v[ny-1, :, i+1] = 0
        v[:, 0, i+1]  = 0
        v[:, nx-1, i+1] = 0
    
    return u, v, p

nx = 41
ny = 41
nt = 201
nit = 50
c = 1
dx = 2 / (nx - 1)
dy = 2 / (ny - 1)
x = numpy.linspace(0, 2, nx)
y = numpy.linspace(0, 2, ny)
X, Y = numpy.meshgrid(x, y)

rho = 1
nu = .1
dt = .001

u = numpy.zeros((ny, nx, nt))
v = numpy.zeros((ny, nx, nt))
p = numpy.zeros((ny, nx, nt)) 
b = numpy.zeros((ny, nx))

u, v, p = cavity_flow(nt, u, v, dt, dx, dy, p, rho, nu)

'''fig = pyplot.figure(figsize=(11,7), dpi=100)
# plotting the pressure field as a contour
pyplot.contourf(X, Y, p[:,:,nt-2], alpha=0.5, cmap='jet')  
pyplot.colorbar()
# plotting the pressure field outlines
#pyplot.contour(X, Y, p[:,:,nt-2], cmap='jet')  
# plotting velocity field
pyplot.quiver(X[::2, ::2], Y[::2, ::2], u[::2, ::2, nt-1], v[::2, ::2, nt-1]) 
pyplot.xlabel('X')
pyplot.ylabel('Y')
pyplot.show()'''


def animate (i):
    cont = pyplot.contourf(X, Y, p[:,:,i], alpha=0.5, cmap='jet')
    quiv = pyplot.streamplot(X, Y, u[:,:,i], v[:,:,i])
    return tuple([cont]) + tuple([quiv])

fig = pyplot.figure(figsize=(11,7), dpi=100)
plot = [pyplot.contourf(X, Y, p[:,:,0], cmap='jet')]
anim = animation.FuncAnimation(fig, animate, nt-2, interval=20)
pyplot.xlabel('X')
pyplot.ylabel('Y')
pyplot.show()

f = r"g:/Kuliah/Python/Python CFD/Animation/11. Cavity.mp4" 
writervideo = animation.FFMpegWriter(fps=30) 
anim.save(f, writer=writervideo)


