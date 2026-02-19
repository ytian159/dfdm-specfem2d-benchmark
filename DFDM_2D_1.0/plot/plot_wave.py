import pandas as pd
import matplotlib.pyplot as plt 
import sys

dt = 0.001
if len(sys.argv) == 1:
    print("Simulation Input not provided\n")
elif len(sys.argv) == 2:
    file_name = sys.argv[1]
    df = pd.read_csv(file_name)

    plt.figure()
    for i in range(1, 40001, 100):
        time = dt*i
        plt.ylim(-2500,2500)
        plt.plot(df.iloc[:,0], df.iloc[:,i])
        plt.title(f"{time:.2f} seconds")
        plt.pause(0.0001)
        plt.clf()

    plt.show()
else:
    file_name = sys.argv[1]
    df = pd.read_csv(file_name)
    i= int(sys.argv[2])
    time = dt*i
    print(f"Saving plot for time step {time:.2f}\n")
    plt.ylim(-2500,2500)
    plt.plot(df.iloc[:,0], df.iloc[:,i])
    plt.title(f"{time:.2f} seconds")
    plt.savefig("wave.jpg")
