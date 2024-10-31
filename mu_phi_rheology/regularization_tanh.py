
#%%
import numpy as np
import matplotlib.pyplot as plt
import scienceplots

plt.style.use('science')

#%%
def pressure_hibler(h, A): 
    
    return Pstar*h*np.exp(-C*(1-A))


def pressure_muphi(h, A, max_diff_A = 1e-10):
    
    conc_max_diff = 1-A
    conc_max_diff[conc_max_diff < max_diff_A] = max_diff_A
    
    p = rhoice*h*(shear*d_mean / conc_max_diff)**2
    
    # p[p>Pstar] = Pstar
    
    return p


def pressure_tanh_1(h, A, shear, I_0, x_0, max_diff_A = 1e-20) : 
    
    conc_max_diff = 1-A
    conc_max_diff[conc_max_diff < max_diff_A] = max_diff_A
    
    p_reg = pressure_hibler(h, A) * np.tanh(rhoice*h*(d_mean*shear/(conc_max_diff))**2 * I_0 / (pressure_hibler(h, A)))
    
    return p_reg

def pressure_tanh_2(h, A, A_0) : 
    
    weight2 = np.tanh((1-A)/A_0)
    # print(weight2[10])
    p_reg = pressure_hibler(h, A)*(1-weight2) + pressure_muphi(h, A)*(weight2)
    # p_reg =  pressure_muphi(h, A)*(weight2)
    return p_reg

#%%
rhoice = 900
d_mean = 1e3
I_0    = 1e-1
Pstar  = 27.5e3
C      = 20
h      = 1


# A = np.arange(0, 1, 0.000001)
oneminus_A = np.linspace(0, 0.03, 100000)
A = 1-oneminus_A

A_0 = np.linspace(0.0001, 0.01, 10)
x_0 = np.linspace(0, 1, 10)
# shear = np.linspace(2e-7, 1e-5, len(A))

shear  = 1e-5

for c, x in enumerate(x_0):
    Preg    = pressure_tanh_1(h, A, 0, x, x)
    Preg2    = pressure_tanh_2(h, A, A_0[c])
    # print(Preg2[-1]/Pstar)
    Phibler = pressure_hibler(h, A)
    Pmuphi  = pressure_muphi(h, A)

    plt.figure()
    ax = plt.axes()
    # ax.set_box_aspect(aspect=1)
    # plt.plot(1-A, Preg2/Phibler, label = r'$P_H(1-w) + P_\mu(w)$')
    plt.plot(1-A, Phibler/Phibler,color = 'orange', label = r'$P_H$')
    plt.plot(1-A, Pmuphi/Phibler, color = 'r', label = r'$P_\mu$')
    plt.plot(1-A, Preg/Phibler, color = 'b', label = r'$\tanh(P_\mu * I_0 /P_H)$')
    # plt.plot(1-A, np.tanh((A-x)/10), label = r'$\mu-\Phi$')
    plt.xlabel('1-A')
    plt.ylabel(r'$P/P_H$ (N/m)')
    plt.title('Reg(shear = {:.2}, $I_0$ = {:.2})'.format(shear, x))
    plt.grid()
    plt.ylim(0,1.1)
    plt.legend()
    plt.savefig('preg1_{}.png'.format(c))
    
    max_p = np.array([np.max(Preg), np.max(Preg2), np.max(Phibler)])
    
    
    plt.figure()
    ax = plt.axes()
    # ax.set_box_aspect(aspect=1)
    plt.scatter(np.arange(0, len(max_p)), max_p/Pstar)
    # plt.plot(1-A, np.tanh((A-x)/10), label = r'$\mu-\Phi$')
    plt.xlabel('Pressures')
    plt.ylabel(r'Max P/Pstar')
    # plt.title('Reg({}, {}, {:.2},{:.2f})'.format(I_0, shear, x, A_0[c]))
    # plt.grid()
    # plt.ylim(0,2)
    # plt.legend()
    plt.savefig('maxp_{}.png'.format(c))
    
    
    
    
    # plt.figure()
    # ax = plt.axes()
    # # ax.set_box_aspect(aspect=1)
    # plt.plot(A, Preg2, label = 'tanh')
    # plt.plot(A, Phibler, label = 'Hibler')
    # plt.plot(A, Pmuphi, label = r'$\mu-\Phi$')
    # plt.xlabel('A')
    # plt.ylabel('P (N/m)')
    # plt.title('Reg({}, {}-{}, {:.2})'.format(I_0, shear[0], shear[-1], x))
    # plt.grid()
    # plt.legend()
    # plt.savefig('preg1_{}.png'.format(c))


#%%
def w(x, a, b):
    return (np.tanh(a * x) + np.tanh(b * x)) / 2

def f1_new(x):
    return np.exp(-20 * x)

def f2_modified(x, k=10):
    return (100 / (1 + k * (x + 0.75)))**2

# Update the combined function F(x)
def F_modified(x, a, b, k):
    return w(x, a, b) * f2_modified(x, k) + (1 - w(x, a, b)) * f1_new(x)
def f2_new(x):
    return (1 / (1 - x))**2
x_new = np.linspace(-0.9, 0.9, 500)
a = 1.0
b = 0.5
# Compute F(x) with the modified f2 function
Fx_modified = F_modified(x_new, a, b, k=10)
def w_modified(x, a, b, A=1.5, sigma=0.1):
    base_weight = (np.tanh(a * x) + np.tanh(b * x)) / 2
    bump = A * np.exp(-((x + 0.75) / sigma)**2)
    return base_weight + bump

# Update the combined function F(x) with the modified weight
def F_with_bump(x, a, b, A, sigma):
    return w_modified(x, a, b, A, sigma) * f2_new(x) + (1 - w_modified(x, a, b, A, sigma)) * f1_new(x)

# Compute F(x) with the modified weight
Fx_bump = F_with_bump(x_new, a, b, A=2, sigma=0.05)

# Plot F(x)
plt.figure(figsize=(8, 6))
plt.plot(x_new, Fx_bump, label=r'$F(x)$ with bump near $x = -0.75$')
plt.title(r'Modified Weight to Reach 100 Near $x = -0.75$')
plt.xlabel(r'$x$')
plt.ylabel(r'$F(x)$')
plt.grid(True)
plt.legend()
plt.show()

# Plot F(x)
plt.figure(figsize=(8, 6))
plt.plot(x_new, Fx_modified, label=r'$F(x)$ (Modified)')
plt.title(r'Modified $(\frac{100}{1+k(x+0.75)})^2$ to reach 100 near $x = -0.75$')
plt.xlabel(r'$x$')
plt.ylabel(r'$F(x)$')
plt.grid(True)
plt.legend()
plt.savefig('test.png')















# %%
