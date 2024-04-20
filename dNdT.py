from sympy import symbols, diff, exp, pi, Pow
from sympy.functions.special.bessel import besselk

degeneracy, V, T, pi, statistics, k, mass, muB = symbols('degeneracy V T pi statistics k mass muB')

# Define the function
f = degeneracy*V/(2*Pow(pi, 2)) * Pow(-1*statistics, k+1) * Pow(mass, 2)*T/k * besselk(2, k*mass/T) * exp(k*muB/T)

# Differentiate with respect to T and print the result
dfdT = diff(f, T)
print(dfdT)