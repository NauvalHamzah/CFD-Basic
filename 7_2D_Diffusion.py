from mpl_toolkits.mplot3d import Axes3D    ##New Library required for projected 3d plots
import numpy
from matplotlib import pyplot, cm
from matplotlib import animation

nx = 31
ny = 31
nt = 100
nu = .05
dx = 2 / (nx - 1)
dy = 2 / (ny - 1)
sigma = .25
dt = sigma * dx * dy / nu

x = numpy.linspace(0, 2, nx)
y = numpy.linspace(0, 2, ny)

u = numpy.ones((ny, nx, nt))  # create a 1xn vector of 1's
un = numpy.ones((ny, nx, nt))

###Assign initial conditions
# set hat function I.C. : u(.5<=x<=1 && .5<=y<=1 ) is 2
u[int(.75 / dy):int(1.25 / dy + 1),int(.75 / dx):int(1.25 / dx + 1),0] = 2 

for n in range(nt-1): ##loop across number of time steps
    un = u.copy()
    u[1:-1, 1:-1, n+1] = (un[1:-1, 1:-1, n] + nu*dt/dx**2*(un[1:-1, 2:, n]-2*un[1:-1, 1:-1, n]+un[1:-1, 0:-2, n])\
                + nu*dt/dy**2*(un[2:, 1:-1, n]-2*un[1:-1, 1:-1, n]+un[0:-2, 1:-1, n]))
    u[0, :, n] = 1
    u[-1, :, n] = 1
    u[:, 0, n] = 1
    u[:, -1, n] = 1

X, Y = numpy.meshgrid(x, y) 

zarray=u.copy()

def update_plot(frame_number, zarray, plot):
    plot[0].remove()
    plot[0] = ax.plot_surface(X, Y, zarray[:,:,frame_number], cmap='jet', edgecolors='k')

fig = pyplot.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_zlim(1, 6)
plot = [ax.plot_surface(X, Y, zarray[:,:,0], cmap='jet', color='0.75', rstride=1, cstride=1)]
anim = animation.FuncAnimation(fig, update_plot, nt, fargs=(zarray, plot), interval=20)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('u')
ax.set_zlim(0.9, 2.1)
pyplot.show()

'''f = r"g:/Kuliah/Python/Python CFD/Animation/7. 2D Diffusion.mp4" 
writervideo = animation.FFMpegWriter(fps=10) 
anim.save(f, writer=writervideo)'''
