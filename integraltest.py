import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys
import os.path
from matplotlib.colors import LogNorm
from matplotlib.ticker import FuncFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.gridspec as gridspec

fmt = lambda x, pos: '{:.6f}'.format(x)

# Set font and font size
plt.rcParams.update({
    #'font.family': 'serif',  # Use Times New Roman font
    'font.size': 14,         # Set the font size for all text in the plot
})

# Load data from the file
datafile = "integraltest.out"
if(len(sys.argv)>1 and os.path.isfile(sys.argv[1])): datafile=sys.argv[1]

column_names = ["NMaxIter", "LogTolerance", "Integral", "NFuncCalls"]
df = pd.read_csv(datafile, sep="\t", header=None, names=column_names)

# Shift the values in "NMaxIter" and "LogTolerance" column by 0.5, to have centered bins
df["LogTolerance"] -= 0.5
df["NMaxIter"] -= 0.5

# Create the matrix plot for "Integral"
integral_matrix = df.pivot_table(index="NMaxIter", columns="LogTolerance", values="Integral", aggfunc=np.mean)

# Create the matrix plot for "NFuncCalls" with logarithmic scale
nfunccalls_matrix = df.pivot_table(index="NMaxIter", columns="LogTolerance", values="NFuncCalls", aggfunc=np.mean)

# Create the two-panel plot with overlapping y-axes
fig, axs = plt.subplots(1, 2, figsize=(16, 8), gridspec_kw={'wspace': 0}) #gridspec_kw={'width_ratios': [0.8, 0.05]})

# Plot for "Integral" with color bar on the left side
im1 = axs[0].imshow(integral_matrix, cmap="viridis", origin="lower", aspect="auto", extent=[min(df["LogTolerance"]), max(df["LogTolerance"]), min(df["NMaxIter"]), max(df["NMaxIter"])])
axs[0].set_title(r'$C(Q)-1$ for' "\n" r'$\lambda={0}$, $R={1}$ fm, $\alpha={2}$, $Q={3}$ GeV/c'.format(1.0, 6.0, 1.4, 0.04))
axs[0].set_xlabel(r"-log$_{10}$(tolerance)")
axs[0].set_ylabel("# of max iterations")
cbar0 = plt.colorbar(im1, ax=axs[0], location='left', format=FuncFormatter(fmt))
cbar0.set_label("C(Q)-1", rotation=90, labelpad=10)  # Rotate the color bar label and increase padding

# Plot for "NFuncCalls" with logarithmic scale and color bar on the right side
im2 = axs[1].imshow(nfunccalls_matrix, cmap="viridis", origin="lower", aspect="auto", extent=[min(df["LogTolerance"]), max(df["LogTolerance"]), min(df["NMaxIter"]), max(df["NMaxIter"])], norm=LogNorm())
axs[1].yaxis.tick_right()
axs[1].set_title(r'# of function calls for' "\n" r'$\lambda={0}$, $R={1}$ fm, $\alpha={2}$, $Q={3}$ GeV/c'.format(1.0, 6.0, 1.4, 0.04))
axs[1].set_xlabel(r"-log$_{10}$(tolerance)")
axs[1].set_ylabel("# of max iterations")  # Adjust labelpad to position the title on the right
axs[1].yaxis.set_label_position("right")  # Set the y-axis title to appear on the right side
cbar1 = plt.colorbar(im2, ax=axs[1], pad=0.11)
cbar1.set_label("# of function calls", rotation=270, labelpad=15) 

plt.subplots_adjust(right=0.85)  # You can adjust the value (0.85 in this case) to create enough space between the y-axis label and color bar

plt.tight_layout()

pngfile = "integraltest.png"
if(len(sys.argv)>2): pngfile=sys.argv[2]
plt.savefig(pngfile)
