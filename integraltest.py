import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm

# Set font and font size
plt.rcParams.update({
    'font.family': 'serif',  # Use Times New Roman font
    'font.size': 14,         # Set the font size for all text in the plot
})

# Load data from the file
datafile = "integraltest.out"
column_names = ["NMaxIter", "LogTolerance", "Integral", "NFuncCalls"]
df = pd.read_csv(datafile, sep="\t", header=None, names=column_names)

# Shift only the values in "NMaxIter" column by 0.5
df["LogTolerance"] -= 0.5
df["NMaxIter"] -= 0.5

# Create the matrix plot for "Integral"
integral_matrix = df.pivot_table(index="NMaxIter", columns="LogTolerance", values="Integral", aggfunc=np.mean)

# Create the matrix plot for "NFuncCalls" with logarithmic scale
nfunccalls_matrix = df.pivot_table(index="NMaxIter", columns="LogTolerance", values="NFuncCalls", aggfunc=np.mean)

# Create the two-panel plot with overlapping y-axes
fig, axs = plt.subplots(1, 2, figsize=(16, 8), gridspec_kw={'wspace': 0.05})

# Plot for "Integral" with color bar on the left side
im1 = axs[0].imshow(integral_matrix, cmap="viridis", origin="lower", aspect="auto", extent=[min(df["LogTolerance"]), max(df["LogTolerance"]), min(df["NMaxIter"]), max(df["NMaxIter"])])
axs[0].set_title("Matrix Plot of Integral")
axs[0].set_xlabel("-log(tolerance)")
axs[0].set_ylabel("NMaxIter")
cbar1 = plt.colorbar(im1, ax=axs[0], label="Integral", location='left')
cbar1.set_label("Integral", rotation=270, labelpad=20)  # Rotate the color bar label and increase padding

# Plot for "NFuncCalls" with logarithmic scale and color bar on the right side
im2 = axs[1].imshow(nfunccalls_matrix, cmap="viridis", origin="lower", aspect="auto", extent=[min(df["LogTolerance"]), max(df["LogTolerance"]), min(df["NMaxIter"]), max(df["NMaxIter"])], norm=LogNorm())
axs[1].set_title("Matrix Plot of NFuncCalls")
axs[1].set_xlabel("-log(tolerance)")
axs[1].set_ylabel("")  # Remove y-axis label for the right plot
axs[1].set_yticks([])  # Remove y-axis ticks for the right plot
cbar2 = plt.colorbar(im2, ax=axs[1], label="NFuncCalls (log scale)")
cbar2.set_label("NFuncCalls (log scale)", rotation=270, labelpad=20)  # Rotate the color bar label and increase padding

plt.tight_layout()
plt.subplots_adjust(wspace=0)  # Adjust the space between the plots
plt.show()
