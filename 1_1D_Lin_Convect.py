import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation

#1. define grid and time step
nx=81
dx=3/(nx-1)
nt=400
dt=0.01
c=2

#2. initial condition
u = np.ones((nt, nx))
u[0, int(.5 / dx):int(1 / dx + 1)] = 3

#3. Calculation
for n in range(1,nt):  
    un = u.copy() 
    u[n,1:] = un[n-1,1:] - c * dt / dx * (un[n-1,1:] - un[n-1,:-1])
    u[n,0] = un[n-1,0] - c * dt / dx * (un[n-1,0] - un[n-1,-2])
    u[n,-1] = u[n,0]


#4. Animate
fig = plt.figure()
ax = plt.axes(xlim=(0, 3.1), ylim=(0.9, 3.1))
line, = ax.plot([], [], lw=2)


# initialization function: plot the background of each frame
def init():
    line.set_data([], [])
    return line,


# animation function.  This is called sequentially
def animate(i):
    x = np.linspace(0, 3, nx)
    y = u[i,:]
    line.set_data(x, y)
    return line,

anim = animation.FuncAnimation(fig, animate, init_func=init,frames=nt, interval=20, blit=True)
plt.grid()
plt.show()

'''f = r"g:/Kuliah/Python/Python CFD/Animation/1. 1D Linear Convection.mp4" 
writervideo = animation.FFMpegWriter(fps=60) 
anim.save(f, writer=writervideo)'''
