import numpy as np
import matplotlib.pyplot as plt

# Read the data from the file
data = np.loadtxt('coulcorrtest.out')

# Plot the data
plt.plot(data[:,0], data[:,1], label='$C_2^{(0)}(Q)$', linestyle='solid')
plt.plot(data[:,0], data[:,2], label='$C_2^{(full)}(Q)$', linestyle='dashed')
plt.plot(data[:,0], data[:,3], label='$K_{Coul}(Q)$', linestyle='dotted')
plt.title(r'$\lambda={0}$, $R={1}$, $\alpha={2}$'.format(0.8, 5.3, 1.2))
plt.legend()
#plt.show()

# Set the x axis to start from 0
plt.xlim(0, None)

# Add an axis title for the x axis
plt.xlabel('Q [GeV/c]')

# Save the plot as a png file
plt.savefig('coulcorrtest.png')
