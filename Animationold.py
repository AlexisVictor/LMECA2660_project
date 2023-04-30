import sys
import os
import glob
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.animation as animation

Nx       = 50 #int(sys.argv[1])
dt       = 2e-3 #int(sys.argv[2])
SimTime  = 1e-1 #int(sys.argv[3])
saveIter = 0 #int(sys.argv[4])
usemixer = 1 #int(sys.argv[5])

Ny = int(1.5*Nx)
L = 1
H = 1.5*L
x = np.linspace(0,L,num=Nx)
y = np.linspace(0,H,num=Ny)
XT, YT = np.meshgrid(x,y)
x = np.linspace(0,L,num=Nx)
y = np.linspace(0,H,num=Ny)
Xv, Yv = np.meshgrid(x,y)
Tcolors = np.linspace(-0.005,0.005,num=101)
vcolors = np.linspace(0,0.05,num=101)
wcolors = np.linspace(-2,2,num=101)

def openfile(file):
    with open(file, "r") as f:
         contenu = f.read()
    return contenu


def fillarray(string, i, array):
    file = "results\{}-{}.txt".format(string, i)
    with open(file, "r") as f: 
        content = []
        for line in f:
            content.append(line.strip().split(";"))
    return content


print(fillarray("T", 19))

exit()

cmap = plt.cm.get_cmap("jet").copy()
N = cmap.N
cmap.set_under(cmap(1))
cmap.set_over(cmap(N-1))

def plotmixer(um,i,ax):
    
    if (um == 1):
        theta = np.linspace(0,2*np.pi,1000)
        omega = 0.1
        r = 0.2*np.cos(3*(theta - i*saveIter*omega/dt))
        x = r*np.cos(theta) + 1/3
        y = r*np.sin(theta) + 1/3
        ax.fill_between(x,0,y,facecolor='black')
        x = 0.04*np.cos(theta) + 1/3
        y = 0.04*np.sin(theta) + 1/3
        ax.fill_between(x,0,y,facecolor='black')


def initTvw(ax):
	# T = 0*XT**0 #np.sin(XT)+np.cos(YT)
	# CS = ax.contourf(XT,YT,T,Tcolors,cmap=cmap,extend="both")
	plt.colorbar(CS,ax=ax)
	plotmixer(usemixer,0,ax)
        
	
	
def animateTvw(i,mf,ax):
    ax.clear()
    ax.set_xlim(0,L)
    ax.set_ylim(0,H)
    plotmixer(usemixer,i,ax)
    return ax
    

f, ax = plt.subplots(1,1)
plt.subplots_adjust(wspace=0.4,hspace=0.3)
ax.set_xlim(0,L)
ax.set_ylim(0,H)
ax.set_aspect('equal')

initTvw(ax)
maxframe = 2000
anim = animation.FuncAnimation(f,animateTvw,interval=1,fargs=(maxframe,ax),frames=maxframe, blit = False)
plt.show()
# anim.save('results/Tvw_Nx%d_dt%d_mixing%d.mp4' % (Nx,dt,usemixer),bitrate=30000)
exit()
