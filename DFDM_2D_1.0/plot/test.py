import matplotlib.pyplot as plt
import numpy as np


fig, ax = plt.subplots(2, 3, sharex='col', sharey='row')

# axes are in a two-dimensional array, indexed by [row, col]
for i in range(2):
    for j in range(3):
        ax[i, j].text(0.5, 0.5, str((i, j)),
                      fontsize=18, ha='center')
plt.show()