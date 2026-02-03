import numpy as np
import matplotlib.pyplot as plt

# set parameters
filename='basic_condition'

expid = filename+'/ex4-case1'

nbgn = 36
nend = 36
nskp = 1

im = 100
jm = 62
km = 12
dx = 200.e3
dy = 100.e3
dz = 50.e0

# set coordinates

xc = (np.arange(im+2)-0.5 ) * dx / 1000.
yc = (np.arange(jm+2)-0.5 ) * dy / 1000.
zc = (np.arange(km+2)-0.5 - km) * dz 

YC,ZC = np.meshgrid(yc,zc, indexing = 'ij') 

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

for n in range(nbgn,nend+1,nskp):

    # read data
    
    fname = expid + f'.n{n:06}'
    fp=open(fname,'rb')
    chunk = np.fromfile(fp, dtype=dtyp)
    fp.close()

    time = chunk['time'][0]
    u = chunk['u'][0].reshape((im+2,jm+2,km+2),order="F")

    # shift staggerd variables to the cell center
    
    uc=np.zeros([im+2,jm+2,km+2])
    uc[:,1:jm,:]=0.5*(u[:,0:jm-1,:]+u[:,1:jm,:])
    uc[:,0,:]=-uc[:,1,:]
    uc[:,jm+1,:]=-uc[:,jm,:]
    UC=uc[50,:,:]*100

    plt.figure(figsize=(8, 6))
    plt.pcolormesh(YC,ZC,UC, cmap='RdBu_r',shading="gouraud",vmin=-70,vmax=70)
    plt.colorbar(label='velocity [cm/s]')
    plt.title('')
    plt.xlabel('y[m]')
    plt.ylabel('z[m]')
    plt.ylim(-600,0)
    plt.savefig(filename+'/figure_yz.png')
    plt.show()


  
