# Assignment-5
Numerical Integration, Numerical Differentiation, and the Stefan-Boltzmann Constant 


#Elizabeth Rodriguez
#Computational Physics
#Assignment 5


# 1: Numerical Integration 
# (A): Consider the integral:
# I = int from 0-1 sin^2[sqrt(100x)]dx
# Write a program that uses the adaptive trapezoidal rule method described in
#Section 5.3 of the book particularly Eq. (5.34):
# I(a,b) = h[1/2f(a)+1/2f(b)+ summation f(a+kh)]
#to calculate the value of this integral to an approximate accuracy of e = 10^-6
# (i.e., correct to six digits after the decimal point). Start with one single 
#integration slice and work up from there to two, four, eight, and so forth. Have
#your program print out the number of slices, its estimate of the integral, and its
#estimate of the error on the integral, for each value of the number of slices N,
#until the target accuracy is reached. (Hint: You should find the result is around
# I = 0.45)

from math import cos, sqrt
from operator import truediv # Helps for division
from scipy.misc import derivative

def f(x):
    return truediv(1,2)-truediv(1,2)*cos(2*sqrt(100*x)) # Using trig identities for sin^2(x)

N = 37 # Number of slices
a = 0.0
b = 1.0
h = truediv((b-a),N)
print 'The number of slices is:', N

# Using Trapezoidal Rule

s = 0.5*f(a) + 0.5*f(b)
for k in range(1,N):
    s += f(a+k*h)

print 'The estimate of the integral is:', round(h*s,6)

# Calculating Error
c = derivative(f, 1.0, dx=1e-6,  n=2) # Second Derivative 
error = truediv((b-a)**3, (12*(N)**2))*c
print 'The estimated error is:', round(error,6)

# (B): Now modify your program to evaluate the same integral using Romberg integration.
# Have your program print out a triangular table of values, as on page 161 of the 
#book, of all the Romberg estimates of the inttegral. Calculate the error on your
#estimates using Eq. (5.49):
# c_m h_i^2m = (1/(4^m -1))(R_i - R_i-1)+ O(h^(2m+2))
#and again continue the calculation until you reach an accuracy of e = 10^-6. You
#should find that the Romberg method reaches the required accuracy considerably 
#faster than the trapezoidal rule alone. 

import numpy as np
from math import cos, sqrt
from operator import truediv # Helps for division
 
def trapezcomp(f, a, b, n): # Composite trapezoidal function integration 
    # Initialization
    h = truediv((b - a),n)
    x = a
    n = 2
    
    # Composite rule
    In = f(a)
    for k in range(1, n):
        x  = x + h
        In += 2*f(x)
 
    return (In + f(b))*h*0.5
 
def romberg(f, a, b, p): # Romberg integration 
    I = np.zeros((p, p)) # p is number of rows
    for k in range(0, p):
        # Composite trapezoidal rule for 2^k panels
        I[k, 0] = trapezcomp(f, a, b, 2**k)
 
        # Romberg recursive formula
        for j in range(0, k):
            I[k, j+1] = (4**(j+1) * I[k, j] - I[k-1, j]) / (4**(j+1) - 1)
 
        print(I[k,0:k+1])   # display intermediate results
 
    return I
 
 
if __name__ == '__main__':
    def func(x):
        return truediv(1,2)-truediv(1,2)*cos(2*sqrt(100*x))
 
    p_rows = 4
    I = romberg(func, 0, 1.0, p_rows)
    solution = I[p_rows-1, p_rows-1]
    print(solution) 
    
# For full credit: Email your program from part (b) and printouts of the output of 
#the programs from parts (a) and (b), showing that the Romberg method reaches the 
#required accuracy in fewer steps than the adaptive trapezoidal rule alone.



# 2: Numerical Differentiation
# Create a user-defined function f(x) that returns the value 1 + 1 tanh(2x), then use
#a central difference to calculate the derivative of the function in the range
# -2<= x <= 2. Calculate an analytic formula for the derivative and make a graph
#with your numerical result and the analytic answer on the same plot. It may help 
#to plot the exact answer as lines and the numerical one as dots. (Hint: In Python
#the tanh function is found in the math package, and its called simply tanh).

# (I):

from math import tanh
from operator import truediv # Helps for division

# Actual integration process 

def f(x):
    return 1+1*tanh(2*x)
    
a = -2
b = 2
    
def derivative(x):   
    h = truediv(1,10)
    rise = f(a+h)-f(b+x)
    run = 2*h
    slope = truediv(rise,run)
    return slope
    
print 'The calculated derivative is:', derivative(a-b)

# (II):

from math import tanh
from sympy import *
    
for x in range(a, b+1):
    print f(x)

# Plotting

plot(derivative(x), 'bo', f(x), 'ro' )
title('Derivative of F(x) as a function of x')
xlabel('x')
ylabel('F(x)')
show()


# For full credit: Email your program plus the plotted output.



# 3: The Stefan-Boltzman Constant
# The Planck theory of thermal radiation tells us that in the (angular) frequency
#interval w to w+dw, a blackbody of unit area radiates electromagnetically an amount
#of thermal energy per second equal to I(w)dw, where
# I(w) = (h-bar/(4*pi^2*c^2))*(w^3/(e^(h-bar*w/k*T)-1))
#here h-bar is Planck's constant over 2pi, c is the speed of light, and k is 
#Boltzmann's constant.

# (A): Show that the total energy per unit area radiated by a blackbody is:
# W = ((k^4*T^4)/(4*pi^2*c^2*h-bar^3)) int from 0-infinity x^3/(e^x -1) dx

# In a Picture File

# (B): Write a program to evaluate the integral in this expression. Explain what
#method you used, and how accurate you think your answer is. 

from sympy import *
from math import exp, pi
import sympy.physics.units as units
from operator import truediv # Helps for division

def plancks_law(wavelength,temperature):
    T = temperature
    h_bar = units.hbar
    k = units.boltzmann
    c = units.c
    freq = (2*pi*c)/wavelength
    I = (h_bar/(4*pi**2*c**2))*((freq**3)/(exp(h_bar*freq/(k*T))-1)) # Equation
    return I.evalfreq() # Evaluate 
   
# Integrate 
   
N = 37 # Number of slices
a = 0.0
b = 1.0
h = truediv((b-a),N)

s = 0.5*f(a) + 0.5*f(b)
for k in range(1,N):
    s += f(a+k*h)

print 'The estimate of the integral is:', round(h*s,6)

# Further explanation is located in an LaTex file 

# (C): Even before Planck gave his theory of thermal radiation around the turn of
#the 20th century, it was known that the total energy W given off a blackbody per
#unit area per second followed Stefans law: W= sigma*T^4, where sigma is the Stefan-
#Boltzmann constant. Use your value for the integral above to compute a value for the
#Stefan-Boltzmann constant (in SI units) to three significant figures. Check your
#result against the known value, which you can find in books or online. You should 
#get good agreement.


# For full credit: Email your program plus the comparison against the known value
