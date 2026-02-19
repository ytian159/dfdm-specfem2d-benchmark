import pandas as pd
import matplotlib.pyplot as plt 
import sys
import numpy as np
import argparse
import os
def calculate_global_min_max_param(data_dir, param, total_elem):
    gz_min = 0
    gz_max = 0
    for elm in range(0, total_elem):
        file_step = f"{data_dir}/{param}_{elm}"
        step_in = pd.read_csv(file_step, header=None)
        step_in = np.array(step_in)
        zmax = step_in.max()
        zmin = step_in.min()
        # print(zmin,zmax)
        if zmax > gz_max:
            gz_max = zmax
        if zmin < gz_min:
            gz_min = zmin
    return gz_min, gz_max

def calculate_min_max(data_dir, time_step, total_elem):
    gz_min = 0
    gz_max = 0
    i=time_step
    for elm in range(0, total_elem):
        file_step = f"{data_dir}/elem_{elm}_{i}.out"
        step_in = pd.read_csv(file_step, header=None)
        step_in = np.array(step_in)
        zmax = np.amax(step_in)
        zmin = np.amin(step_in)
        if zmax > gz_max:
            gz_max = zmax
        if zmin < gz_min:
            gz_min = zmin
        # print(gz_min, gz_max)
    return gz_min, gz_max

def plot_grid(x2d, z2d, elm, show_element_ids):
    nx, nz = x2d.shape
    plt.plot(x2d[0] / 1e3, z2d[0] / 1e3, '--k', linewidth=1)
    plt.plot(x2d[nx - 1] / 1e3, z2d[nx - 1] / 1e3, '--k', linewidth=1)
    plt.plot(x2d[:, 0] / 1e3, z2d[:, 0] / 1e3, '--k', linewidth=1)
    plt.plot(x2d[:, nz - 1] / 1e3, z2d[:, nz - 1] / 1e3, '--k', linewidth=1)
    if show_element_ids:
        x_center = np.mean(x2d)  # Average of all x-coordinates
        z_center = np.mean(z2d)  # Average of all z-coordinates
        plt.text(x_center / 1e3, z_center / 1e3, f"{elm}", color="red", fontsize=10, ha="center", va="center")
    
    

def plot_umid(x2d, z2d, Umid, vmin, vmax):
    plt.pcolormesh(x2d / 1e3, z2d / 1e3, np.abs(Umid), shading='gouraud', cmap="jet")
    # print(vmax)
    # plt.clim([5e10, 1e11])
    plt.clim([vmin, vmax])

def main():
    parser = argparse.ArgumentParser(description="Plot 2D unstructured data.")
    parser.add_argument("--data_dir", type=str, required=True, help="Directory containing the grid and element files.")
    parser.add_argument("--total_elem", type=int, required=True, help="Total number of elements.")
    args = parser.parse_args()

    data_dir = args.data_dir
    total_elem = args.total_elem

    # if not args.plot_grid_only:
    #     vmin, vmax = calculate_global_min_max(data_dir, time_steps, total_elem)
    # vmax=max(abs(gz_min), abs(gz_max))
    # print(vmax)
    # vmax=0.
    # vmin=0.
    
    os.makedirs(data_dir+"/plots/", exist_ok=True)
    for param in ['mu11']: #, 'mu22', 'rho12', 'rho21'
        plt.clf()
        vmin, vmax = calculate_global_min_max_param(data_dir, param, total_elem)
        print(vmin,vmax)
        
        # vmax=max(abs(gz_min), abs(gz_max))
        # if vmax==0.:
        #     vmax=1.
        # print(vmax)
        for elm in range(0, total_elem):
            x2d = np.array(pd.read_csv(f"{data_dir}/grids/grid_x_{elm}", header=None))
            z2d = np.array(pd.read_csv(f"{data_dir}/grids/grid_z_{elm}", header=None))

            show_element_ids=False
            plot_grid(x2d, z2d, elm, show_element_ids)

            file_step = f"{data_dir}/{param}_{elm}"
            Umid = np.array(pd.read_csv(file_step, header=None))
            # print(np.amax(np.abs(Umid)))
            plot_umid(x2d, z2d, Umid, vmin, vmax)
        # if vmax>0.:
        plt.gca().axis('equal')
        plt.colorbar()
        print(data_dir+'/plots/plot_'+param+'.png')
        plt.savefig(data_dir+'/plots/plot_'+param+'.png', dpi=300)
        # plt.pause(0.5)


if __name__ == "__main__":
    main()
