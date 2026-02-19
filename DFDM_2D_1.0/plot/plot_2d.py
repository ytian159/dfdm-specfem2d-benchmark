import pandas as pd
import matplotlib.pyplot as plt 
import sys
import numpy as np

dt = 0.001
if len(sys.argv) == 1:
    print("Simulation Inputs not provided\n")
elif len(sys.argv) == 3:
    # step_file = sys.argv[1]
    x2d_file = sys.argv[1]
    z2d_file = sys.argv[2]

    x2d = pd.read_csv(x2d_file, header=None)
    z2d = pd.read_csv(z2d_file, header=None)
    # row_ = z2d.iloc[0,:]
    # row_ = [i for i in range(141)]
    # x, y = np.meshgrid(row_,row_)
    # print(x)
    for i in range(500,10000,500):
        step_in = pd.read_csv(str("step"+str(i)+".out"), header=None)
        # step_in = pd.read_csv(str("step1000.out"), header=None)
        step_in = np.array(step_in)
        fig, (ax0) = plt.subplots()
        cmap = plt.colormaps["jet"]
        ax0.pcolormesh(x2d,z2d,step_in, cmap=cmap)
    
        plt.pause(1)
        plt.clf()
  
