from mpl_toolkits.mplot3d import Axes3D    ##New Library required for projected 3d plots
import numpy
from matplotlib import pyplot, cm
from matplotlib import animation

###variable declarations
nx = 81
ny = 81
nt = 100
dx = 2 / (nx - 1)
dy = 2 / (ny - 1)
sigma = .2
dt = sigma * dx

x = numpy.linspace(0, 2, nx)
y = numpy.linspace(0, 2, ny)

u = numpy.ones((ny, nx, nt)) ##create a 1xn vector of 1's
v = numpy.ones((ny, nx, nt))
un = numpy.ones((ny, nx, nt))
vn = numpy.ones((ny, nx, nt))


##set hat function I.C. : u(.5<=x<=1 && .5<=y<=1 ) is 2
u[int(.5 / dy):int(1 / dy + 1),int(.5 / dx):int(1 / dx + 1),0] = 2
v[int(.5 / dy):int(1 / dy + 1),int(.5 / dx):int(1 / dx + 1),0] = 2 

for n in range(nt-1): ##loop across number of time steps
    un = u.copy()
    vn = v.copy()
    u[1:, 1:, n+1] = (un[1:, 1:, n] - (un[1:, 1:, n] * dt / dx * (un[1:, 1:, n] - un[1:, :-1, n])) -
                                  (vn[1:, 1:, n] * dt / dy * (un[1:, 1:, n] - un[:-1, 1:, n])))
    v[1:, 1:, n+1] = (vn[1:, 1:, n] - (un[1:, 1:, n] * dt / dx * (vn[1:, 1:, n] - vn[1:, :-1, n])) -
                                  (vn[1:, 1:, n] * dt / dy * (vn[1:, 1:, n] - vn[:-1, 1:, n])))
    
    u[0, :, n] = 1
    u[-1, :, n] = 1
    u[:, 0, n] = 1
    u[:, -1, n] = 1

    v[0, :, n] = 1
    v[-1, :, n] = 1
    v[:, 0, n] = 1
    v[:, -1, n] = 1

X, Y = numpy.meshgrid(x, y) 

zarray=v.copy()

def update_plot(frame_number, zarray, plot):
    plot[0].remove()
    plot[0] = ax.plot_surface(X, Y, zarray[:,:,frame_number], cmap='jet', edgecolors='k')

fig = pyplot.figure()
ax = fig.add_subplot(111, projection='3d')

plot = [ax.plot_surface(X, Y, zarray[:,:,0], cmap='jet', color='0.75', rstride=1, cstride=1)]
anim = animation.FuncAnimation(fig, update_plot, nt, fargs=(zarray, plot), interval=20)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('u')
ax.set_zlim(0.9, 2.1)
fig.colorbar(plot[0], shrink=0.5, aspect=5)
pyplot.show()

'''f = r"g:/Kuliah/Python/Python CFD/Animation/6. 2D non-Linear Convection.mp4" 
writervideo = animation.FFMpegWriter(fps=10) 
anim.save(f, writer=writervideo)'''
