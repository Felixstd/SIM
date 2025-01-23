import numpy as np
import matplotlib.pyplot as plt
import scienceplots

plt.style.use('science')

phi_min = 0.0
phi_max = 0.9

A = np.linspace(0,1+1e-4,1000)

def phiofA(A,phi_min,phi_max,  tau = 0.2):
    
    phi_s = (1-phi_min) * A + phi_min
    
    phi_g = - (1-phi_max) * np.tanh((A-1)/tau) + (phi_max)
     
    phi = phi_s * phi_g
    
    return phi, phi_s, phi_g

def phiofA2(A,phi_min,phi_max):
    
    phi_s = (1-phi_min) * A + phi_min
    
    phi = (phi_max-phi_min) * np.sin(np.pi*A/2) + phi_min
         
    phi_g = phi / phi_s 
    
    return phi, phi_s, phi_g

def phiofA3(A, phi_min, phi_max, phi_0 = 0.01):
    
    phi_s = (1-phi_min) * A + phi_min
    
    phi = phi_min + (phi_max - phi_min)/(phi_0/A + 1)
    
    phi_g = phi / phi_s
    
    return phi, phi_s, phi_g

def phiofA4(A, phi_min, phi_max, b = 2):
    
    phi_s = (1-phi_min) * A + phi_min
    
    phi = (phi_max - phi_min) /(1 - b) * (1-b**A) + phi_min
    
    phi_g = phi  / phi_s
    
    return phi, phi_s, phi_g

def Aofphi2(phi,phi_min,phi_max):
    
    A = 2/np.pi * np.arcsin((phi-phi_min)/(phi_max-phi_min))
        
    return A

def Aofphi3(phi, phi_min, phi_max, phi_0 = 2):
    
    A = phi_0/((phi_max - phi_min)/(phi - phi_min)-1)
    
    return A

def Aofphi4(phi, phi_min, phi_max, b = 2):
    
    A = np.log(1 - (phi*(1-b) - phi_min)/(phi_max - phi_min)) / np.log(b)
    
    return A

def phiofA5(A, phi_min, phi_max, b = 2):
    
    phi_s = (1-phi_min) * A + phi_min
    
    phi = (phi_max - phi_min) * A + phi_min
    
    phi_g = phi  / phi_s
    
    return phi, phi_s, phi_g

def plot_phiofA(phi, phi_s, phi_g, phi_min, phi_max, namefig):
    
    fig = plt.figure(figsize=(4,4))
    plt.plot(A,phi,label='corrected')
    plt.plot(A,phi_s, ls='dotted', label='linear')
    plt.plot(A,phi_g,label='correction')
    plt.plot()
    
    if phi_max != 1.0 :
        plt.axhline(phi_max,color='k',ls='dashed')
    if phi_min != 0.0 :
        plt.axhline(phi_min,color='k',ls='dashed')

    plt.axhline(0,color='k')
    plt.axvline(1,color='k')
    plt.axvline(0,color='k')
    plt.axhline(1,color='k')
    plt.grid()
    plt.legend(loc='lower right', bbox_to_anchor=(1.4, 0.5))
    plt.axis('equal')
    plt.ylim([0-0.1,1+0.1])
    plt.xlim([0-0.1,1+0.1])
    plt.xlabel('$A$')
    plt.ylabel(r'$\phi$')
    
    plt.savefig(namefig)
    

phi_1, phi_s_1, phi_g_1 = phiofA(A,phi_min,phi_max, tau = 0.2)
plot_phiofA(phi_1, phi_s_1, phi_g_1, phi_min, phi_max, 'phiofA_1.png')

phi_2, phi_s_2, phi_g_2 = phiofA2(A,phi_min,phi_max)
plot_phiofA(phi_2, phi_s_2, phi_g_2, phi_min, phi_max, 'phiofA_2.png')

phi_3, phi_s_3, phi_g_3 = phiofA3(A,phi_min,phi_max)
plot_phiofA(phi_3, phi_s_3, phi_g_3, phi_min, phi_max, 'phiofA_3.png')

phi_4, phi_s_4, phi_g_4 = phiofA4(A,phi_min,phi_max)
plot_phiofA(phi_4, phi_s_4, phi_g_4, phi_min, phi_max, 'phiofA_4.png')

phi_5, phi_s_5, phi_g_5 = phiofA5(A,phi_min,phi_max)
plot_phiofA(phi_5, phi_s_5, phi_g_5, phi_min, phi_max, 'phiofA_5.png')

plt.figure(figsize=(4,4))
plt.plot(A, phi_1, label = r'$\phi_1$')
plt.plot(A, phi_2, label = r'$\phi_2$')
plt.plot(A, phi_3, label = r'$\phi_3$')
plt.plot(A, phi_4, label = r'$\phi_4$')
plt.plot(A, phi_5, label = r'$\phi_l$')
if phi_max != 1.0 :
    plt.axhline(phi_max,color='k',ls='dashed')
if phi_min != 0.0 :
    plt.axhline(phi_min,color='k',ls='dashed')
plt.axhline(0,color='k')
plt.axvline(1,color='k')
plt.axvline(0,color='k')
plt.axhline(1,color='k')
plt.grid()
plt.legend(loc='lower right', bbox_to_anchor=(1.22, 0.5))
plt.axis('equal')
plt.ylim([0-0.1,1+0.1])
plt.xlim([0-0.1,1+0.1])
plt.xlabel('$A$')
plt.ylabel(r'$\phi$')
plt.savefig('PhiofA_all.png')



def hexagonal_packing(R, W, H):
    # Circle diameter
    D = 2 * R
    
    # Hexagonal grid spacing
    dx = D
    dy = np.sqrt(3) * R  # Height of equilateral triangle
    
    # Generate grid points
    x_coords = np.arange(R, W, dx)
    y_coords = np.arange(R, H, dy)
    
    circles = []
    for i, y in enumerate(y_coords):
        for x in x_coords:
            if i % 2 == 0:
                # Offset alternate rows for hexagonal packing
                cx = x
            else:
                cx = x + R
            if cx + R <= W and y + R <= H:  # Ensure circles are within bounds
                circles.append((cx, y))
    
    return circles

def plot_circles(circles, R, W, H):
    fig, ax = plt.subplots()
    ax.set_xlim(0, W)
    ax.set_ylim(0, H)
    ax.set_aspect('equal', adjustable='box')
    
    # Plot the domain
    rect = plt.Rectangle((0, 0), W, H, linewidth=1, edgecolor='black', facecolor='none')
    ax.add_patch(rect)
    
    # Plot circles
    for (cx, cy) in circles:
        circle = plt.Circle((cx, cy), R, edgecolor='blue', facecolor='lightblue', alpha=0.6)
        ax.add_patch(circle)
    
    plt.savefig('circles.png')

# Parameters
R = 10  # Radius of circles
W = 5000  # Width of domain
H = 1000  # Height of domain

# Compute circle packing
circles = hexagonal_packing(R, W, H)

# Calculate packing fraction
n = len(circles)
packing_fraction = (n * np.pi * R**2) / (W * H)
print(f"Number of circles: {n}")
print(f"Packing fraction: {packing_fraction:.4f}")

plot_circles(circles, R, W, H)