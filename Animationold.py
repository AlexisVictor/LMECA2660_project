import sys
import os
import glob
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.animation as animation

m        = 75
dt       = 1/300
SimTime  = 300
usemixer = 1

n = int(2/3*m)
H = 1
L = H/(1.5)
x = np.linspace(0,L,n)
y = np.linspace(0,H,m)
X, Y = np.meshgrid(x,y)
order = 5e-3
encod = 300
colors = np.linspace(0,order,num=101)




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


cmap = plt.cm.get_cmap("hot").copy()
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
        ax.fill_between(x,0,y,facecolor='white')
        x = 0.04*np.cos(theta) + L/2
        y = 0.04*np.sin(theta) + L/2
        ax.fill_between(x,0,y,facecolor='white')


def initTvw(ax):
    array = fillarray("T", 0, n+2, m+2)
    array = np.abs(array)
    CS = ax.contourf(X,Y,array,colors,cmap=cmap,extend="both")
    plt.colorbar(CS,ax=ax)
    plotmixer(usemixer,0,ax)
        
	
	
def animatearray(iter,mf,ax):
    array = fillarray("T", iter, n+2, m+2)
    array = np.abs(array)
    ax.clear()
    ax.set_xlim(0,L)
    ax.set_ylim(0,H)
    ax.set_title(r"$\frac{{{}U}}{{H}}$ = {}".format('t', np.round(iter * dt,2)))
    ax.contourf(X,Y,array,colors,cmap=cmap,extend="both")
    plotmixer(usemixer,iter,ax)
    return ax
    

f, ax = plt.subplots(1,1)
ax.set_xlim(0,L)
ax.set_ylim(0,H)
ax.set_aspect('equal')

initTvw(ax)
maxframe = 300 #(int) (SimTime/dt)     

anim = animation.FuncAnimation(f,animatearray,interval=1,fargs=(maxframe,ax),frames=maxframe, blit = False)   
plt.show()

# f = r"VelocityMixe.gif" 
# writergif = animation.PillowWriter(fps=15) 
# anim.save(f, writer=writergif)
