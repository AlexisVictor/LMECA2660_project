import sys
import os
import glob
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.animation as animation

n        = 50
dt       = 1/428
SimTime  = 5e1
saveIter = 1 
usemixer = 0 
m = int(1.5*n)
L = 1
H = 1.5*L
x = np.linspace(0,L,num=n)
y = np.linspace(0,H,num=m)
X, Y = np.meshgrid(x,y)
Tcolors = np.linspace(0,1e-11,num=500)
# vcolors = np.linspace(0,0.05,num=101)
# wcolors = np.linspace(-2,2,num=101)

print(np.shape(X), np.shape(Y))


def fillarray(string, i, xlength, ylength):
    file = "results/{}-{}.txt".format(string, i)
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


cmap = plt.cm.get_cmap("jet").copy()
N = cmap.N
cmap.set_under(cmap(1))
cmap.set_over(cmap(N-1))

def plotmixer(bool, iter,ax):
    if (bool == 1):
        theta = np.linspace(0,2*np.pi,1000)
        omega = 0.1
        r = 0.2*np.cos(3*(theta + iter*saveIter*omega/dt))
        x = r*np.cos(theta) + 0.5
        y = r*np.sin(theta) + + 0.5
        ax.fill_between(x,0,y,facecolor='black')
        x = 0.04*np.cos(theta) + 0.5
        y = 0.04*np.sin(theta) + 0.5
        ax.fill_between(x,0,y,facecolor='black')


def initTvw(ax):
    array = fillarray("Velocity", 1, n, m)
    #CS = ax.contourf(X,Y,array,Tcolors,cmap=cmap,extend="both")
    #plt.colorbar(CS,ax=ax)
    plotmixer(usemixer,0,ax)
        
	
	
def animatearray(iter,mf,ax):
    array = fillarray("Velocity", iter, n, m)
    ax.clear()
    ax.set_xlim(0,L)
    ax.set_ylim(0,H)
    ax.set_title("iter = {}".format(iter+1))
    ax.contourf(X,Y,array,Tcolors,cmap=cmap,extend="both")
    plotmixer(usemixer,iter,ax)
    return ax
    

f, ax = plt.subplots(1,1)
ax.set_xlim(0,L)
ax.set_ylim(0,H)
ax.set_aspect('equal')

initTvw(ax)
maxframe = 18 #(int) (SimTime/dt)
anim = animation.FuncAnimation(f,animatearray,interval=100,fargs=(maxframe,ax),frames=maxframe, blit = False)
plt.show()

# f = r"Ufoirée.gif" 
# writergif = animation.PillowWriter(fps=15) 
# anim.save(f, writer=writergif)