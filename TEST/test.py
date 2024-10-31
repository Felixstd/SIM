import numpy as np
import matplotlib.pyplot as plt


theta = np.deg2rad(45)
nx = 502
ny = 502
d = 282
thickness = d
maskC = np.zeros((nx, ny))
# matrix_thick = np.zeros((M, M))

# # Set the top-left corner and bottom-right corner with a "thickness"
# maskC[0:thickness, 0:thickness] = 1  # Top-left corner
# maskC[nx-thickness:nx, ny-thickness:ny] = 1


A = np.zeros((nx, ny))
h = np.zeros((nx, ny))
# Parameters for the new angle of 25 degrees
slope = np.tan(theta)

center_x = nx / 2.0
center_y = ny / 2.0

# ! intercept_1 = - d / (2 * cos(theta))
# ! intercept_2 =  d / (2 * cos(theta))

intercept_2 = 200
intercept_1 = intercept_2 - d*np.sqrt(slope**2 + 1)
print(intercept_1, intercept_2)

y_end_1 = slope * nx + intercept_1
x_end_2 = (ny - intercept_2) / slope 
slope_2 = (y_end_1 - ny)/(nx - x_end_2)
intercept_3 = y_end_1 - slope_2*nx

for i in range(0, nx):
    for j in range(0, ny):

        # x = i - center_x
        # y = j - center_y
        x = i
        y = j

        if (y < (- x * slope + intercept_2)):
        
            maskC[i, j] = 1.0
            
        if y > ( x * slope_2 + intercept_3):
            # print('here')
            maskC[i, j] = 1
            
        # if (y >= (slope * x + intercept_1) and y <=( slope * x + intercept_2)) : 
        #     # maskC[i, j] = 1
            
        if (y >= (slope * x + intercept_1) and y <= ( slope * x + intercept_2)) :
            print('here')
            A[i, j] = 1
            h[i, j] = 1
                    
        
            
        # if (y == slope * x + intercept_2):
        #     maskC[i, j] = 3
            
        # else:
        #     maskC[i, j] = 1
        
        # ! elseif (y >= (slope * x + intercept_1) .and. y <=( slope * x + intercept_2)) then
        # !     maskC(i, j) = 1.0
    
#         endif
#     end do
# end do

# for i in range(thickness):
#     for j in range(i+1):  # j goes up to i to form the right triangle in top-left
#         maskC[i, j] = 1

# # Set the bottom-right triangle (right angle in the corner)
# for i in range(nx-thickness, nx):
#     for j in range(ny-1, ny-(ny-i)-1, -1):  # j starts from the end to form the right triangle
#         maskC[i, j] = 1

A[maskC == 1] = 0


# Display the resulting matrix for the angle of 25 degrees
plt.imshow(A, cmap='viridis', origin='lower')
# plt.title(f'Matrix centered at (0,0) between two parallel lines at {theta_25}°')
# plt.colorbar(label='Value')
plt.savefig('test.png')

