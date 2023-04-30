import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Données

n = 50
m = 1.5*n
dt = 2e-3
x = np.linspace(0, 1, n)
y = np.linspace(0, 1.5, m)
X, Y = np.meshgrid(x, y)
Z = X**2+Y**2

# Création de la figure
fig, ax = plt.subplots()

# Initialisation du graphique
im = ax.imshow(Z, cmap='viridis', extent=[-np.pi, np.pi, -np.pi, np.pi])
fig.colorbar(im)


def plotmixer(um,i,ax):
    if (um == 1):
        theta = np.linspace(0,2*np.pi,1000)
        omega = 0.1
        t = i*dt
        r = 0.2*np.cos(3*(theta - omega*t))
        x = r*np.cos(theta) + 1/3
        y = r*np.sin(theta) + 1/3
        ax.fill_between(x,0,y,facecolor='black')
        x = 0.04*np.cos(theta) + 1/3
        y = 0.04*np.sin(theta) + 1/3
        ax.fill_between(x,0,y,facecolor='black')

# Fonction d'animation
def update(frame):
    # Mise à jour des données
    Z = np.sin(X + frame / 10) * np.cos(Y + frame / 10)
    # Mise à jour de l'image
    im.set_data(Z)
    # Retourne l'image mise à jour
    return im,

# Création de l'animation
ani = FuncAnimation(fig, update, frames=200, blit=True, interval=50)

# Affichage du graphique
plt.show()