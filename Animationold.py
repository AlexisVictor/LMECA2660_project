import sys
import os
import glob
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.animation as animation

m = 150
n  = int(2/3*m)
dt       = 1/300
SimTime  = 1e4
saveIter = 1 
usemixer = 0 

H = 1
L = 2/3*H
x = np.linspace(0,L,num=n)
y = np.linspace(0,H,num=m)
X, Y = np.meshgrid(x,y)
# Tcolors = np.linspace(-1e-1,1e-1,num=100)
Tcolors = np.linspace(0,5e-3,num=100)
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

cmap1 = plt.get_cmap("coolwarm")
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


# def initTvw(ax):
#     array = fillarray("U",0, n+1, m+2)
#     #CS = ax.contourf(X,Y,array,Tcolors,cmap=cmap,extend="both")
#     #plt.colorbar(CS,ax=ax)
#     plotmixer(usemixer,0,ax)
        
	
	
# def animatearray(iter,mf,ax,c):
#     array = fillarray("U", iter, n+1, m+2)
#     ax.clear()
#     ax.set_xlim(0,L)
#     ax.set_ylim(0,H)
#     ax.set_title("iter = {}".format(iter))
#     ax.contourf(X,Y,array,Tcolors,cmap=c,extend="both")
#     plotmixer(usemixer,iter,ax)
#     return ax

# def initTvw(ax):
#     array = fillarray("Vorticity",0, n+1, m+1)
#     #CS = ax.contourf(X,Y,array,Tcolors,cmap=cmap,extend="both")
#     #plt.colorbar(CS,ax=ax)
#     plotmixer(usemixer,0,ax)
        
	
	
# def animatearray(iter,mf,ax,c):
#     array = fillarray("Vorticity", iter, n+1, m+1)
#     ax.clear()
#     ax.set_xlim(0,L)
#     ax.set_ylim(0,H)
#     ax.set_title("iter = {}".format(iter))
#     ax.contourf(X,Y,array,Tcolors,cmap=c,extend="both")
#     plotmixer(usemixer,iter,ax)
#     return ax

# def initTvw(ax):
#     array = fillarray("Velocity",0, n, m)
#     #CS = ax.contourf(X,Y,array,Tcolors,cmap=cmap,extend="both")
#     #plt.colorbar(CS,ax=ax)
#     plotmixer(usemixer,0,ax)
        
	
	
# def animatearray(iter,mf,ax,c):
#     array = fillarray("Velocity", iter, n, m)
#     ax.clear()
#     ax.set_xlim(0,L)
#     ax.set_ylim(0,H)
#     ax.set_title("iter = {}".format(iter))
#     ax.contourf(X,Y,array,Tcolors,cmap=c,extend="both")
#     plotmixer(usemixer,iter,ax)
#     return ax

def initTvw(ax):
    array = fillarray("T",1, n+2, m+2)
    #CS = ax.contourf(X,Y,array,Tcolors,cmap=cmap,extend="both")
    #plt.colorbar(CS,ax=ax)
    plotmixer(usemixer,0,ax)
        
	
    
def animatearray(iter,mf,ax,c):
    array = fillarray("T", iter, n+2, m+2)
    ax.clear()
    ax.set_xlim(0,L)
    ax.set_ylim(0,H)
    ax.set_title("iter = {}".format(iter))
    ax.contourf(X,Y,array,Tcolors,cmap=c,extend="both")
    plotmixer(usemixer,iter,ax)
    return ax
    

f, ax = plt.subplots(1,1)
ax.set_xlim(0,L)
ax.set_ylim(0,H)
ax.set_aspect('equal')

initTvw(ax)
maxframe = 375#(int) (SimTime/dt)
anim = animation.FuncAnimation(f,animatearray,interval=100,fargs=(maxframe,ax,cmap1),frames=maxframe, blit = False)
plt.show()

# f = r"Ufoirée.gif" 
# writergif = animation.PillowWriter(fps=15) 
# anim.save(f, writer=writergif)