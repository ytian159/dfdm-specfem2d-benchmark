import pandas as pd
import matplotlib.pyplot as plt 
import sys
import numpy as np
import glob

total_elem=len(glob.glob("grid_z_*"))

gz_max = 0
gz_min = 0
for i in range(0,1000,50):
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

print(f"gz_max:{gz_max}, gz_min:{gz_min}")





for i in range(0,1000,50):
    # fig=plt.figure()
    # ax=plt.gca()
    print(f"step:{i}")
    for elm in range(0,total_elem):
        fig=plt.figure()
        ax=plt.gca()
        file_step = str("elem_"+str(elm)+"_"+str(i)+".out")
        step_in = pd.read_csv(file_step, header=None)
        x2d = pd.read_csv("grid_x_"+str(elm), header=None).transpose()
        z2d = pd.read_csv("grid_z_"+str(elm), header=None).transpose()
        x,z = x2d, z2d 

        step_in = np.array(step_in)
        cmap = plt.colormaps["jet"]

        ax.pcolormesh(x,z,step_in.transpose(), cmap=cmap, shading='nearest') #vmin=gz_min, vmax=gz_max,
        plt.savefig('out_'+str(i)+'_'+str(elm)+'.png')
    # plt.pause(1)