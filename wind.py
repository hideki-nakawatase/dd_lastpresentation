import numpy as np
import matplotlib.pyplot as plt

filename='ty_applied'

t0=-0.1

im=20
jm=13

tx=np.zeros((im,jm))
ty=np.zeros((im,jm))

x=np.linspace(0,20000,im)
y=np.linspace(0,6000,jm)

X,Y=np.meshgrid(x,y,indexing='ij')

for i in range(1,im):
  for j in range(1,jm):
    tx[i,j] = t0 * np.cos( np.pi * (j-jm/2) / jm/2 )
    ty[i,j] = t0 * np.sin( np.pi/2 * (j-jm/2) / jm/2 )

plt.figure(figsize=(8,6))
Q = plt.quiver(X[1:-1,1:-1],Y[1:-1,1:-1],
               tx[1:-1,1:-1],ty[1:-1,1:-1],cmap='jet',pivot='mid')
plt.colorbar(Q,label='wind stress [N/m^2]')
plt.title('wind stress')
plt.xlabel('x [km]')
plt.ylabel('y [km]')
plt.savefig(filename+'/wind.png')