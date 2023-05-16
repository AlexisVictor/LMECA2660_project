import sys
import os
import glob
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.animation as animation

m        = 60
dt       = 1/30
SimTime  = 1000
usemixer = 1

n = int(2/3*m)
H = 1
L = H/(1.5)
x = np.linspace(0,L,n)
y = np.linspace(0,H,m)
X, Y = np.meshgrid(x,y)
order = 5e-2
encod = 50
colors = np.linspace(0,order,num=100)

def fillarray(string, i, xlength, ylength):
    file = "results\{}-{}.txt".format(string, i)
    with open(file, "r") as f: 
        content = []
        for line in f:
            content.append(line.strip().split(";"))
    array = np.empty((xlength, ylength))
    for elem in content:
        xpos = int(elem[0])
        ypos = int(elem[1])
        value = elem[2]
        array[xpos][ypos] = value
    if xlength!=n:
        if xlength-n==2:
            array = np.delete(array, [0,n+1], 0)
        elif xlength-n==1:
            array = np.delete(array, [n], 0)

    if ylength!=n:
        if ylength-m == 2:
            array = np.delete(array, [0,m+1], 1 )
        elif ylength-m==1:
            array = np.delete(array, [m], 1)

    array = np.transpose(array)
    
    return array 

cmap = plt.cm.get_cmap('copper').copy()
#cmap = plt.cm.get_cmap('copper').copy()
N = cmap.N
cmap.set_under(cmap(1))
cmap.set_over(cmap(N-1))

def plotmixer(bool, iter,ax):
    if (bool == 1):
        theta = np.linspace(0,2*np.pi,1000)
        omega = 0.1
        r = 0.2*np.cos(3*(theta + iter*encod*omega*dt))
        x = r*np.cos(theta) + L/2
        y = r*np.sin(theta) + L/2
        color = "gainsboro"
        ax.fill_between(x,0,y,facecolor=color)
        x = 0.04*np.cos(theta) + L/2
        y = 0.04*np.sin(theta) + L/2
        ax.fill_between(x,0,y,facecolor=color)


def plot(string, name, xname, yname):
    file = "results\{}.txt".format(string)
    with open(file, "r") as f: 
        content = []
        for line in f:
            content.append(line.strip().split(";"))
    content = np.array(content)
    array1 = np.zeros(len(content))
    time = np.zeros(len(content))
    for i in range(len(content)):
        array1[i] = np.float64(content[i][0])
        time[i]  = np.float64(content[i][1])*dt
    plt.plot(time,array1)
    plt.xlabel(xname, fontsize = 15)
    plt.ylabel(yname, fontsize = 15)
    plt.title(name)
    plt.show()


def initTvw(ax):
    array = fillarray("Velocity", 0, n, m)
    # array = np.abs(array)
    CS = ax.contourf(X,Y,array,colors,cmap=cmap,extend="both")
    plt.colorbar(CS,ax=ax)
    plotmixer(usemixer,0,ax)
        
	
	
def animatearray(iter,mf,ax):
    array = fillarray("Velocity", iter, n, m)
    # array = np.abs(array)
    ax.clear()
    ax.set_xlim(0,L)
    ax.set_ylim(0,H)
    ax.set_title(r"$\frac{{{}U}}{{H}}$ = {}".format('t', np.round(encod*iter * dt,2)))
    ax.contourf(X,Y,array,colors,cmap=cmap,extend="both")
    plotmixer(usemixer,iter,ax)
    return ax
    

f, ax = plt.subplots(1,1)
ax.set_xlim(0,L)
ax.set_ylim(0,H)
ax.set_aspect('equal')

initTvw(ax)
maxframe = (int) (SimTime/(dt*encod))     


anim = animation.FuncAnimation(f,animatearray,interval=1,fargs=(maxframe,ax),frames=maxframe, blit = False)   
plt.show()

plot("Fluxaverage", "Averaged free-surface heat flux density" , r"$\frac{tU}{H}$", r"$\frac{\langle q_e \rangle (t)}{q_w}$")
plot("Temperatureaverage", "Averaged fluid temperature" , r"$\frac{tU}{H}$", r"$\frac{\langle T \rangle (t) - T_{\infty}}{\Delta T}$")
plot("TcylAverage", "Averaged mixer temperature" , r"$\frac{tU}{H}$", r"$\frac{T_{cyl}(t) - T_{\infty}}{\Delta T}$")
plot("RmsTemperature", "rms of the temperature" , r"$\frac{tU}{H}$", r"$\frac{T_{rms} - T_{\infty}}{\Delta T}$")

# f = r"Velocitymixercopper300-1000.gif" 
# writergif = animation.PillowWriter(fps=10) 
# anim.save(f, writer=writergif)

