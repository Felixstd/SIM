import numpy as np
import matplotlib.pyplot as plt

# Parameters for the new angle of 25 degrees
theta_25 = 25.0  # Angle in degrees
N = 20.0          # Distance between lines

# Convert angle to radians and calculate the slope
theta_25_rad = np.radians(theta_25)
slope_25 = np.tan(theta_25_rad)

# Set the matrix size
M = 100
center = M / 2

# Calculate the intercepts for the two parallel lines, symmetric about (0,0)
intercept1 = -N / 2.0 / np.cos(theta_25_rad)
intercept2 = intercept1 + N / np.cos(theta_25_rad)

# Initialize the matrix to zeros
matrix_25 = np.zeros((M, M))

# Fill the matrix using the centered intercepts for 25 degrees
for i in range(M):
    for j in range(M):
        x = i - center  # Adjust x relative to the center of the matrix
        y = j - center  # Adjust y relative to the center of the matrix
        if slope_25 * x + intercept1 <= y <= slope_25 * x + intercept2:
            matrix_25[i, j] = 1

# Display the resulting matrix for the angle of 25 degrees
plt.imshow(matrix_25, cmap='gray', origin='lower')
plt.title(f'Matrix centered at (0,0) between two parallel lines at {theta_25}°')
plt.colorbar(label='Value')
plt.savefig('test.png')

