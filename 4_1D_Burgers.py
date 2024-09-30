import numpy
import sympy
from matplotlib import pyplot as plt
from sympy import init_printing
from sympy.utilities.lambdify import lambdify
from matplotlib import animation

init_printing(use_latex=True)
x, nu, t = sympy.symbols('x nu t')
phi = (sympy.exp(-(x - 4 * t)**2 / (4 * nu * (t + 1))) + sympy.exp(-(x - 4 * t - 2 * sympy.pi)**2 / (4 * nu * (t + 1))))
phiprime = phi.diff(x)
u = -2 * nu * (phiprime / phi) + 4
ufunc = lambdify((t, x, nu), u)

###variable declarations
nx = 101
nt = 200
dx = 2 * numpy.pi / (nx - 1)
nu = .07
dt = dx * nu

x = numpy.linspace(0, 2 * numpy.pi, nx)
un = numpy.empty((nt,nx))
u = numpy.empty((nt,nx))

#initial condition
for i in range (nx):
    u[0,i]=ufunc(0 * dt, x[i], nu)

#calculation
for n in range(1,nt):
    un = u.copy()
    u[n,1:-1] = un[n-1,1:-1] - un[n-1,1:-1] * dt / dx *(un[n-1,1:-1] - un[n-1,0:-2]) + nu * dt / dx**2 *\
                (un[n-1,2:] - 2 * un[n-1,1:-1] + un[n-1,0:-2])
    u[n,0] = un[n-1,0] - un[n-1,0] * dt / dx * (un[n-1,0] - un[n-1,-2]) + nu * dt / dx**2 *\
                (un[n-1,1] - 2 * un[n-1,0] + un[n-1,-2])
    u[n,-1] = u[n,0]

#analytical solution
u_analytical = numpy.empty((nt,nx))        
for n in range (1,nt):
    for i in range (nx):
        u_analytical[n,i]=ufunc(n*dt, x[i], nu)
    
#ANIMATE
fig, ax = plt.subplots()
ax = plt.axes(xlim=(0, 6.5), ylim=(0, 10))
line1, = ax.plot([], [], lw=1, color='blue', marker='x')
line2, = ax.plot([], [], lw=1, color='red')
line, = ax.plot([], [], lw=1)

# initialization function: plot the background of each frame
def init():
    line.set_data([], [])
    return line,

def animate(i):
    line1.set_data(x, u[i,:])
    line2.set_data(x, u_analytical[i,:])
    return [line1,line2]

anim = animation.FuncAnimation(fig, animate, init_func=init,frames=nt, interval=10, blit=True)
plt.grid()
plt.show()

'''f = r"g:/Kuliah/Python/Python CFD/Animation/4. 1D Burgers.mp4" 
writervideo = animation.FFMpegWriter(fps=60) 
anim.save(f, writer=writervideo)'''