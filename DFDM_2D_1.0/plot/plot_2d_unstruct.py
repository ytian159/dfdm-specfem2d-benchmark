import pandas as pd
import matplotlib.pyplot as plt 
import sys
import numpy as np

total_elem = 48
time_steps = 10000
gz_max = 0
gz_min = 0
for i in range(0,time_steps,50):
    for elm in range(0,total_elem):
        file_step = str("elem_"+str(elm)+"_"+str(i)+".out")
        step_in = pd.read_csv(file_step, header=None)
        step_in = np.array(step_in)
        zmax = step_in.max()
        zmin = step_in.min()
        if zmax > gz_max:
            gz_max = zmax
        if zmin < gz_min:
            gz_min = zmin
# gz_min = 0


for i in range(0,time_steps,50):
    plt.clf()
    print(f"time_step:{i}")
    for elm in range(0,total_elem):
        file_step = str("elem_"+str(elm)+"_"+str(i)+".out")
        Umid =  np.array(pd.read_csv(file_step, header=None))
        x2d =  np.array(pd.read_csv("grid_x_"+str(elm), header=None))
        z2d =  np.array(pd.read_csv("grid_z_"+str(elm), header=None))

        nx, nz = x2d.shape
        # print(x2d.shape)
        id1 = np.arange(0, nx*nz, nx)
        id2 = np.arange(nx-1, nx*nz, nx)
        id3 = np.arange(0, nx)
        id4 = np.arange(nx*(nz-1), nx*nz)
        x2d_F = x2d.flatten()
        z2d_F = z2d.flatten()

        plt.plot(x2d[0]/1e3, z2d[0]/1e3, '--k', linewidth=1)
        plt.plot(x2d[nx-1]/1e3, z2d[nx-1]/1e3, '--k', linewidth=1)
        plt.plot(x2d[:,0]/1e3, z2d[:,0]/1e3, '--k', linewidth=1)
        plt.plot(x2d[:,nz-1]/1e3, z2d[:,nz-1]/1e3, '--k', linewidth=1)
        cmap = plt.colormaps["jet"]
        
        plt.pcolormesh(x2d/1e3, z2d/1e3, Umid, shading='gouraud', cmap="jet")
        plt.clim([-8e-8, 8e-8])
    plt.pause(0.5)
