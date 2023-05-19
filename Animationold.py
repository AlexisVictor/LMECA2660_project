import sys
import os
import glob
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.animation as animation

m        = 300
dt       = 1/300
SimTime  = 3000
usemixer = 1

n = int(2/3*m)
H = 1
L = H/(1.5)
x = np.linspace(0,L,n)
y = np.linspace(0,H,m)
X, Y = np.meshgrid(x,y)
order = 1e-0
encod = 500
colors = np.linspace(-order,order,num=100)

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

cmap = plt.cm.get_cmap('seismic').copy()
cmap1 = plt.cm.get_cmap('copper').copy()
N = cmap.N
cmap.set_under(cmap(1))
cmap.set_over(cmap(N-1))

def plotmixer(bool, iter,ax):
    if (bool == 1):
        theta = np.linspace(0,2*np.pi,1000)
        omega = 0.1
        r = 0.2*np.cos(3*(theta - iter*encod*omega*dt))
        x = r*np.cos(theta) + L/2
        y = r*np.sin(theta) + L/2
        color = "gainsboro"
        ax.fill_between(x,0,y,facecolor=color)
        x = 0.04*np.cos(theta) + L/2
        y = 0.04*np.sin(theta) + L/2
        ax.fill_between(x,0,y,facecolor=color)


def plot(string, name, xname, yname):
    file = "results/{}.txt".format(string)
    with open(file, "r") as f: 
        content = []
        for line in f:
            content.append(line.strip().split(";"))
    content = np.array(content)
    array1 = np.zeros(len(content))
    time = np.zeros(len(content))
    for i in range(len(content)):
        # print(i)
        array1[i] = np.float64(content[i][0])
        time[i]  = np.float64(content[i][1])*dt
    plt.plot(time,array1)
    plt.xlabel(xname, fontsize = 15)
    plt.ylabel(yname, fontsize = 15)
    plt.title(name)
    plt.tight_layout()
    # plt.savefig(name+".eps")
    plt.show()


def initTvw(ax):
    array = fillarray("Vorticity", 0, n+1, m+1)
    # array = np.abs(array)
    CS = ax.contourf(X,Y,array,colors,cmap=cmap,extend="both")
    plt.colorbar(CS,ax=ax)
    plotmixer(usemixer,0,ax)
        
	
	
def animatearray(iter,mf,ax):
    array = fillarray("Vorticity", iter, n+1, m+1)
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


anim = animation.FuncAnimation(f,animatearray,interval=5,fargs=(maxframe,ax),frames=maxframe, blit = False)   
plt.show()

plot("Fluxaverage", "Averaged free-surface heat flux density" , r"$\frac{tU}{H}$", r"$\frac{\langle q_e \rangle (t)}{q_w}$")
plot("Temperatureaverage", "Averaged fluid temperature" , r"$\frac{tU}{H}$", r"$\frac{\langle T \rangle (t) - T_{\infty}}{\Delta T}$")
plot("TcylAverage", "Averaged mixer temperature" , r"$\frac{tU}{H}$", r"$\frac{T_{cyl}(t) - T_{\infty}}{\Delta T}$")
plot("RmsTemperature", "rms of the temperature" , r"$\frac{tU}{H}$", r"$\frac{T_{rms} - T_{\infty}}{\Delta T}$")



# cmap = plt.cm.get_cmap('seismic').copy()
# cmapcopper = plt.cm.get_cmap('binary').copy()
# cmaphot = plt.cm.get_cmap('hot').copy()
# N = cmap.N
# cmap.set_under(cmap(1))
# cmap.set_over(cmap(N-1))
# iter=1
# array = fillarray("Velocity", 600, n, m)
# print(np.amax(array))
# print(np.shape(array))

# time =np.array([0,20,100])//5*3

# fig, axs = plt.subplots(1, 3,  layout='constrained',figsize=(12.9, 6))
# for i, ax in enumerate(axs.flat):
#     array = fillarray("Vorticity", time[i], n+1, m+1)
#     ax.set_xlim(0,L)
#     ax.set_ylim(0,H)
#     ax.set_title(r"$\frac{{{}U}}{{H}}$ = {}".format('t', np.round(encod*time[i] * dt,2)))
#     pcm = ax.contourf(X,Y,array,colors,cmap=cmap,extend="both")
#     plotmixer(1,time[i],ax, 'black')
# plt.colorbar(pcm, ax=axs[2], shrink=0.8, location='right')
# # plt.tight_layout()
# # plt.savefig("Vorticity_0.png")
# plt.savefig("Vorticity_mixer_0.png")
# plt.show()

# fig, axs = plt.subplots(1, 3,  layout='constrained',figsize=(12.9, 6))
# for i, ax in enumerate(axs.flat):
#     array = fillarray("T", time[i], n+2, m+2)
#     ax.set_xlim(0,L)
#     ax.set_ylim(0,H)
#     ax.set_title(r"$\frac{{{}U}}{{H}}$ = {}".format('t', np.round(encod*time[i] * dt,2)))
#     pcm = ax.contourf(X,Y,array,colorsh,cmap=cmaphot,extend="both")
#     plotmixer(1,time[i],ax,'white')
# plt.colorbar(pcm, ax=axs[2], shrink=0.8, location='right')
# # plt.tight_layout()
# plt.savefig("T_mixer_0.png")
# # plt.savefig("T_0.png")
# plt.show()

# fig, axs = plt.subplots(1, 3, layout='constrained', figsize=(12.9, 6))
# for i, ax in enumerate(axs.flat):
#     array = fillarray("Velocity", time[i], n, m)
#     ax.set_xlim(0,L)
#     ax.set_ylim(0,H)
#     ax.set_title(r"$\frac{{{}U}}{{H}}$ = {}".format('t', np.round(encod*time[i] * dt,2)))
#     pcm = ax.contourf(X,Y,array,colorsv,cmap=cmapcopper,extend="both")
#     plotmixer(1,time[i],ax,'black')
# plt.colorbar(pcm, ax=axs[2], shrink=0.8, location='right')
# # # plt.tight_layout()
# plt.savefig("Velocity_mixer_0.png")
# # plt.savefig("Velocity_0.png")
# # # plt.show()
# plt.show()


# time = np.array([200,500,1000])//5*3

# fig, axs = plt.subplots(1, 3, layout='constrained', figsize=(12.9, 6))
# for i, ax in enumerate(axs.flat):
#     array = fillarray("Vorticity", time[i], n+1, m+1)
#     ax.set_xlim(0,L)
#     ax.set_ylim(0,H)
#     ax.set_title(r"$\frac{{{}U}}{{H}}$ = {}".format('t', np.round(encod*time[i] * dt,2)))
#     pcm = ax.contourf(X,Y,array,colors,cmap=cmap,extend="both")
#     plotmixer(1,time[i],ax, 'black')
# plt.colorbar(pcm, ax=axs[2], shrink=0.8, location='right')
# # plt.tight_layout()
# plt.savefig("Vorticity_mixer_1.png")
# # plt.savefig("Vorticity_1.png")
# plt.show()

# fig, axs = plt.subplots(1, 3, layout='constrained', figsize=(12.9, 6))
# for i, ax in enumerate(axs.flat):
#     array = fillarray("T", time[i], n+2, m+2)
#     ax.set_xlim(0,L)
#     ax.set_ylim(0,H)
#     ax.set_title(r"$\frac{{{}U}}{{H}}$ = {}".format('t', np.round(encod*time[i] * dt,2)))
#     pcm = ax.contourf(X,Y,array,colorsh,cmap=cmaphot,extend="both")
#     plotmixer(1,time[i],ax, 'white')
# plt.colorbar(pcm, ax=axs[2], shrink=0.8, location='right')
# # plt.tight_layout()
# plt.savefig("T_mixer_1.png")
# # plt.savefig("T_1.png")
# plt.show()
# # plt.close()

# fig, axs = plt.subplots(1, 3, layout='constrained', figsize=(12.9, 6))
# for i, ax in enumerate(axs.flat):
#     array = fillarray("Velocity", time[i], n, m)
#     ax.set_xlim(0,L)
#     ax.set_ylim(0,H)
#     ax.set_title(r"$\frac{{{}U}}{{H}}$ = {}".format('t', np.round(encod*time[i] * dt,2)))
#     pcm = ax.contourf(X,Y,array,colorsv,cmap=cmapcopper,extend="both")
#     plotmixer(1,time[i],ax, 'black')
# plt.colorbar(pcm, ax=axs[2], shrink=0.8, location='right')
# # plt.tight_layout()
# plt.savefig("Velocity_mixer_1.png")
# # plt.savefig("Velocity_1.png")
# plt.show()

# def plot(string, name, xname, yname):
#     file = "results/{}.txt".format(string)
#     with open(file, "r") as f: 
#         content = []
#         for line in f:
#             content.append(line.strip().split(";"))
#     content = np.array(content)
#     array1 = np.zeros(len(content))
#     time = np.zeros(len(content))
#     for i in range(len(content)):
#         # print(i)
#         array1[i] = np.float64(content[i][0])
#         time[i]  = np.float64(content[i][1])*dt
#     plt.plot(time,array1,'black')
#     plt.xlabel(xname, fontsize = 15)
#     plt.ylabel(yname, fontsize = 15)
#     plt.title(name)
#     plt.tight_layout()
#     plt.savefig(name+".eps")
#     plt.show()

# plot("Fluxaverage", "Averaged free-surface heat flux density with mixer" , r"$\frac{tU}{H}$", r"$\frac{\langle q_e \rangle (t)}{q_w}$")
# plot("Temperatureaverage", "Averaged fluid temperature with mixer" , r"$\frac{tU}{H}$", r"$\frac{\langle T \rangle (t) - T_{\infty}}{\Delta T}$")
# plot("TcylAverage", "Averaged mixer temperature" , r"$\frac{tU}{H}$", r"$\frac{T_{cyl}(t) - T_{\infty}}{\Delta T}$")
# plot("RmsTemperature", "rms of the temperature with mixer" , r"$\frac{tU}{H}$", r"$\frac{T_{rms} - T_{\infty}}{\Delta T}$")


# time =np.arange(1,3000//5*3+1,1)
# Khey = np.zeros(3000//5*3)
# KheyV = np.zeros(3000//5*3)
# for i in range(3000//5*3):
#     array = fillarray("Vorticity", time[i], n+1, m+1)
#     array = abs(array)*(1/300)**2*10**5
#     Khey[i] = np.amax(array)
#     arrayv = fillarray("Velocity", time[i], n, m)
#     arrayv = abs(arrayv)*(1/300)**2*10**5
#     KheyV[i] = np.amax(arrayv)

# plt.figure(figsize=(10,6))
# plt.plot(time,Khey,"black",label="$Re_{h,\omega}$")
# plt.plot((0,time[-1]),(np.mean(Khey),np.mean(Khey)),'red')
# plt.legend()
# plt.xlim(0,time[-1])
# plt.savefig("kheyomega_withoutmixer.eps")
# plt.show()

# plt.figure(figsize=(10,6))
# plt.plot(time,KheyV,"black",label="$Re_{h}$")
# plt.plot((0,time[-1]),(np.mean(KheyV),np.mean(KheyV)),'red')
# plt.legend()
# plt.xlim(0,time[-1])
# plt.savefig("Kheyh_without_mixer.eps")
# plt.show()
# # plt.tight_layout()
# # plt.savefig("Vorticity_0.png")
# # plt.savefig("Vorticity_mixer_0.png")
# # plt.show()