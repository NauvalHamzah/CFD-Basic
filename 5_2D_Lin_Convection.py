from mpl_toolkits.mplot3d import Axes3D    ##New Library required for projected 3d plots
import numpy
from matplotlib import pyplot, cm
from matplotlib import animation

###variable declarations
nx = 81
ny = 81
nt = 100
cx = 2
cy = 1
dx = 2 / (nx - 1)
dy = 2 / (ny - 1)
sigma = .2
dt = sigma * dx

x = numpy.linspace(0, 2, nx)
y = numpy.linspace(0, 2, ny)

u = numpy.ones((ny, nx, nt)) ##create a 1xn vector of 1's
un = numpy.ones((ny, nx, nt)) ##

###Assign initial conditions

##set hat function I.C. : u(.5<=x<=1 && .5<=y<=1 ) is 2
u[int(.5 / dy):int(1 / dy + 1),int(.5 / dx):int(1 / dx + 1),0] = 2 

for n in range(nt-1): ##loop across number of time steps
    un = u.copy()
    u[1:, 1:, n+1] = (un[1:, 1:, n] - (cx * dt / dx * (un[1:, 1:, n] - un[1:, :-1, n])) -
                                  (cy * dt / dy * (un[1:, 1:, n] - un[:-1, 1:, n])))
    u[0, :, n] = 1
    u[-1, :, n] = 1
    u[:, 0, n] = 1
    u[:, -1, n] = 1

X, Y = numpy.meshgrid(x, y) 

zarray=u.copy()

def update_plot(frame_number, zarray, plot):
    #plot[0].remove()
    plot[0] = ax.plot_surface(X, Y, zarray[:,:,frame_number], cmap="jet", edgecolors='k')
    return plot[0]
    

fig = pyplot.figure()
ax = fig.add_subplot(111, projection='3d')
plot = [ax.plot_surface(X, Y, zarray[:,:,0], cmap='jet', rstride=1, cstride=1)]
anim = animation.FuncAnimation(fig, update_plot, nt, fargs=(zarray, plot), interval=20)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('u')
ax.set_zlim(0.9, 2.1)
fig.colorbar(plot[0], shrink=0.5, aspect=5)
pyplot.show()

f = r"g:/Kuliah/Python/Python CFD/Animation/5. 2D Linear Convection.mp4" 
writervideo = animation.FFMpegWriter(fps=10) 
anim.save(f, writer=writervideo)