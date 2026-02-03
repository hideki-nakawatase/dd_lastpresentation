import numpy as np
import matplotlib.pyplot as plt

rdir='./'
expid = 'ex4-case1'
nbgn = 36
nend = 36
nskp = 1

im = 100
jm = 62
km = 12
#jm = 62
#km = 1
dx = 200.e3
dy = 100.e3
dz = 50.e0

# set coordinates
xc = np.arange(im+1) * dx / 1000.
yc = np.arange(jm+1) * dy / 1000.

XC,YC = np.meshgrid(xc, yc, indexing = 'ij') 

# Loop

dtyp = np.dtype( [('time','<d'),\
                  ('u','<'+str((im+2)*(jm+2)*(km+2))+'d'),\
                  ('v','<'+str((im+2)*(jm+2)*(km+2))+'d'),\
                  ('e','<'+str((im+2)*(jm+2))+'d'),\
                  ('w','<'+str((im+2)*(jm+2)*(km+2))+'d'),\
                  ('p','<'+str((im+2)*(jm+2)*(km+2))+'d'),\
                  ('t','<'+str((im+2)*(jm+2)*(km+2))+'d'),\
                  ('s','<'+str((im+2)*(jm+2)*(km+2))+'d'),\
                  ('r','<'+str((im+2)*(jm+2)*(km+2))+'d'),\
                  ] )
psi = np.zeros((im+1,jm+1))

for n in range(nbgn,nend+1,nskp):

    # read data
    fname = rdir + expid + f'.n{n:06}'
    fp=open(fname,'rb')
    chunk = np.fromfile(fp, dtype=dtyp)
    fp.close()
    time = chunk['time'][0]
    u = chunk['u'][0].reshape((im+2,jm+2,km+2),order="F")
    um = np.sum(u[:im+1,:jm+1,1:km+1],axis=2) * dz
    for i in range(1,im):
        for j in range(1,jm):
            psi[i,j] = psi[i,j-1] + um[i,j]*dy

    fig, ax = plt.subplots(figsize=(8,6))

    cntr = ax.contour(XC,YC,psi)
    plt.clabel(cntr)
    ax.grid()
    #figb.set_xlim(-0.05,0.1)
    #figb.set_ylim(-400.,0.)
    ax.set_title('Stream function')
    ax.set_xlabel('X [km]')
    ax.set_ylabel('Y [km]')

    ax.set_aspect('equal', adjustable='box')
    fig.tight_layout()
    plt.show()
