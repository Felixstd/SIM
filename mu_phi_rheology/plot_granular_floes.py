import numpy as np
import matplotlib.pyplot as plt
# Function to check if two circles overlap
def circles_overlap(x1, y1, r1, x2, y2, r2):
    distance = np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)
    return distance < (r1 + r2)

# Function to find a non-overlapping position for a new circle
def find_position(existing_circles, radius, start_x=0, start_y=0):
    x, y = start_x, start_y
    while any(circles_overlap(x, y, radius, xc, yc, rc) for xc, yc, rc in existing_circles):
        x += 0.1  # Increment x slightly to avoid overlap
        if x > 1000:  # Prevent infinite loop in extreme cases
            raise RuntimeError("Could not find a non-overlapping position.")
    return x, y
circle_radii_grid = [
    [0.5, 0.8, 0.6],
    [1.2, 0.7, 0.9],
    [0.6, 0.8, 1.0]
]


# List to store all circles' positions and radii
placed_circles = []

# Start placing circles
x_start, y_start = 0, 0
for row_index, row in enumerate(circle_radii_grid):
    for radius in row:
        x, y = find_position(placed_circles, radius, start_x=x_start, start_y=y_start)
        placed_circles.append((x, y, radius))
        x_start = x + radius + 0.1  # Prepare x_start for the next circle in the row
    y_start += max(row) * 2  # Prepare y_start for the next row

# Reinitialize the plot for tightly packed non-overlapping circles
fig, ax = plt.subplots()

# Plot each circle using the adjusted positions and radii
for x, y, r in placed_circles:
    circle = plt.Circle((x, y), r, edgecolor='blue', facecolor='none', linewidth=2)
    ax.add_patch(circle)

# Adjust the plot limits to ensure all circles are visible
x_positions = [x for x, y, r in placed_circles]
y_positions = [y for x, y, r in placed_circles]
radii = [r for x, y, r in placed_circles]

x_max = max(x_positions) + max(radii) + 1
y_max = max(y_positions) + max(radii) + 1
ax.set_xlim(-1, x_max)
ax.set_ylim(-1, y_max)
ax.set_aspect('equal')  # Equal aspect ratio

# Add labels and grid for clarity
plt.title("Non-Overlapping Packed Circles")
plt.xlabel("X")
plt.ylabel("Y")
plt.grid(True)

plt.savefig('circles.png')