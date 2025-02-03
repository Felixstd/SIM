import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import skewnorm
from scipy.integrate import quad
import scienceplots

plt.style.use('science')


def methode_simpson(func, a, b) : 

    middle = (b+a)/2
    h = abs(b-a)/2

    integral = (h/3) * (func(a) + 4*func(middle) + func(b)) #sum of the intervals

    return integral

def recursion_simpson(func, a, b, erreur, total) : 

    """Implementation of the recursion method associated to the simpson method 

    Inputs :
        func : function to integrate
        a : lower bound
        b : upper bound
        erreur : error that we want 
        total : the left interval + the right interval

    Returns: 
        if the error associated with the total is smaller than 15*erreur (15 is from Prof. Mousseau)
            return the integral value

        if not : 
            take the left and right interval and call the recursion function to have a smaller error. 
    """

    middle = (b+a)/2

    intervalle_gauche = methode_simpson(func, a, middle)
    intervalle_droit = methode_simpson(func, middle, b)

    #If yes : end
    if abs(intervalle_gauche + intervalle_droit - total) <= erreur : 


        integral_value = intervalle_gauche + intervalle_droit + (intervalle_gauche + intervalle_droit - total)/15
        return integral_value

    #If no : start again with smaller error and more intervals

    return recursion_simpson(func, a, middle, erreur/2, intervalle_gauche) + recursion_simpson(func, middle, b, erreur/2, intervalle_droit)


def methode_adaptative_simpson(func,a, b, erreur) : 

    """Implementation of the adaptative method 

    Inputs : 
        func : function to integrate
        a : lower bound
        b : upper bound
        erreur : error 

    Returns:
       Integral : value of the integral
    """

    #Calculate one integral to feed the recursion loop
    integral = methode_simpson(func, a, b) 
    
    return recursion_simpson(func, a, b, erreur, integral) #recursion loop


def g_of_h(h, a = 8, loc = 1, scale = 2):
    
    # return skewnorm.pdf(h, a=a, loc=loc, scale=scale)
    A = scale*np.exp(-scale*loc)
    return A*np.exp(-scale*(h-loc))

def b_of_h(h, G_star = 0.15):
    
    G_of_h = quad(g_of_h, 0, h)
    print(G_of_h)
    G_of_h = G_of_h[0]
    if G_of_h <= G_star:
        bh = 2/G_star*(1-G_of_h/G_star)
    
    elif (G_of_h > G_star) and (G_of_h <= 1):
        bh = 0
        
    return bh

def a_of_h(h) : 
    
    return b_of_h(h)*g_of_h(h)

def integrand_a_times_h(h):
    print('a of h', a_of_h(H))
    return h**2*a_of_h(H)

def icestrength(h, cf, cp, k):
    
    pstar = (k*cp + k*cf/(k-1))*methode_adaptative_simpson(integrand_a_times_h, 0, 5000, 1e-5)
    # pstar = (k*cp + k*cf/(k-1))*quad(integrand_a_times_h, 0, 10000)[0]
    return pstar


g = 9.80
rhow = 1026
rhoi = 900
tanphi = 0.8
k = 5
mu = 0.7
H=2

cp = rhoi*(rhow-rhoi)*g/(2*rhow)
print(rhoi*(rhow-rhoi)*g)
cf = mu*(rhow - rhoi)*g*(rhoi*(k-1)/rhow)**2/(2*tanphi)

a = 8
h = np.linspace(0, 10, 10000)

gh = g_of_h(h, a =0.1, loc = 3, scale = 0.6)
# print(gh)

print(quad(g_of_h, 0, 10000), quad(g_of_h, 0, 2), quad(g_of_h, 0, 3))
plt.figure()
plt.plot(h, gh)
plt.grid()
plt.axhline(0.15)
plt.xlabel('Ice Thickness (m)')
plt.ylabel(r'$g(h)$')
plt.savefig("gh.png")


print(icestrength(1, cf, cp, k)/1e3)
# 