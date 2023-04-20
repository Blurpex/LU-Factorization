# Project 2
#
# Group 4:
# Shiza Qureshi
# Haris Ahmad
# Zakaria Coulibaly
# Lilac Dixon
# Quazi Maliha

import math
import matplotlib.pyplot as plt
import numpy as np
        
# x and y represent arrays for the base points
# val is the value being interpolated for
# returns the interpretted value at val
def langrangeInterpolation(x, y, val):
    n = len(x)  # get the number of total base points
    sum = 0     # the total sum of the lagrange polynomial
    for i in range(n):
        product = y[i]      # set the starting value of product to the y value of the current base point.
        for j in range(n):  
            if i != j:      # skip the term for the current base point
                product = product * (val - x[j]) / (x[i] - x[j])   # evaluate the current term's product.
        sum = sum + product     # add the current term to the total sum
    return sum  # return the interpolated value of the inputed value.
    

n = int(input("\nPlease enter the number of base points, n, within [-1,1]:\n"))

chebyshevX = []
chebyshevY = []
equalSpaceX = []
equalSpaceY = []

# Create n Chebyshev and equally spaced base points
for i in range(1, n + 1):
    chebyshevX.append(math.cos(((2 * i - 1) * math.pi) / (2 * n)))
    equalSpaceX.append(-1 + ((2 / (n - 1)) * (i - 1)))
    chebyshevY.append(math.exp(abs(chebyshevX[i - 1])))
    equalSpaceY.append(math.exp(abs(equalSpaceX[i - 1])))

# Create an x-axis from -1 to 1 with 200 steps, so step-size = 0.01.
xAxis = np.linspace(-1, 1, 200)
# Calculate the Chebyshev interpolated values using Lagrange Interpolation
chebyshevInterpolated = [langrangeInterpolation(chebyshevX, chebyshevY, i) for i in xAxis]
# Calculate the Equal Space interpolated values using Lagrange Interpolation
equalSpaceInterpolated = [langrangeInterpolation(equalSpaceX, equalSpaceY, i) for i in xAxis]
# Calculate the actual values of the function e^|x|
actualValues = np.exp(abs(xAxis))

# Add grid lines to the plot
plt.grid(True, which='both')
plt.axvline(x=0, color='k')
plt.axhline(y=0, color='k')

# Plot the actual values in black.
plt.plot(xAxis, actualValues, label = "f(x) = e^|x|", color = "black")
# Plot the Lagrange Interpolated polynomial with Chebyshev base points in purple
plt.plot(xAxis, chebyshevInterpolated, label = "Chebyshev Points", color = "purple")
# Plot the Lagrange Interpolated polynomial with equally spaced base points in red
plt.plot(xAxis, equalSpaceInterpolated, label = "Equally Spaced Points", color = "red")

plt.show()