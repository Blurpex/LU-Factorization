import numpy as np
import matplotlib.pyplot as plt


# Define the function f(x) = e|x|
def f(x):
    return np.exp(np.abs(x))


# Define the evenly spaced interpolation nodes
def evenly_spaced(n):
    return np.linspace(-1, 1, n)


# Define the Chebyshev interpolation nodes
def chebyshev(n):
    return np.cos((np.arange(1, 2 * n + 1, 2) * np.pi) / (2 * n))


# Define the Newton divided difference formula
def newtdd(x, y):
    n = len(x)
    coef = np.zeros((n, n))
    coef[:, 0] = y
    for j in range(1, n):
        for i in range(n - j):
            coef[i][j] = (coef[i + 1][j - 1] - coef[i][j - 1]) / (x[i + j] - x[i])
    return coef[0]


# Define the nested multiplication formula
def nest(coef, x, r):
    n = len(coef)
    ans = coef[n - 1]
    for i in range(n - 2, -1, -1):
        ans = ans * (r - x[i]) + coef[i]
    return ans


# Define the interpolation function
def interpolate(x, y, n, domain):
    x_interp = np.arange(domain[0], domain[1] + 0.01, 0.01)
    y_true = f(x_interp)
    x_nodes = x(n)
    y_nodes = f(x_nodes)
    coef = newtdd(x_nodes, y_nodes)
    y_interp = nest(coef, x_nodes, x_interp)
    error = y_true - y_interp
    return x_interp, y_true, y_interp, error


# Set the parameters
n_values = [10, 20]
domain = [-1, 1]

# Plot the results
plt.figure(figsize=(12, 8))
for i, n in enumerate(n_values):
    plt.subplot(2, 2, i + 1)
    plt.title("Degree {} polynomials".format(n))
    x_interp1, y_true1, y_interp1, error1 = interpolate(evenly_spaced, f, n, domain)
    x_interp2, y_true2, y_interp2, error2 = interpolate(chebyshev, f, n, domain)
    plt.plot(x_interp1, y_true1, label="True function")
    plt.plot(x_interp1, y_interp1, label="Evenly spaced")
    plt.plot(x_interp2, y_interp2, label="Chebyshev")
    plt.legend()

    plt.subplot(2, 2, i + 3)
    plt.title("Empirical errors (n={})".format(n))
    plt.plot(x_interp1, error1, label="Evenly spaced")
    plt.plot(x_interp2, error2, label="Chebyshev")
    plt.legend()

plt.tight_layout()
plt.show()

"""
Presentation
Title: Comparing Evenly Spaced Interpolation with Chebyshev Interpolation

Today, I will be presenting a project that compares the performance of two methods of polynomial interpolation - 
Evenly Spaced Interpolation and Chebyshev Interpolation. The aim of this project is to investigate the effectiveness 
of these two methods in approximating a function using a polynomial of a given degree.

Firstly, let's take a brief look at the theory behind polynomial interpolation. Polynomial interpolation is the process 
of finding a polynomial function that passes through a given set of data points. The degree of the polynomial is equal 
to the number of data points minus one. The two methods of polynomial interpolation that we will be discussing today are 
Evenly Spaced Interpolation and Chebyshev Interpolation.

EVEN SPACED INTERPOLATION

Evenly Spaced Interpolation is a method that uses equally spaced nodes as data points to construct a polynomial function 
that approximates a given function. 

CHEBYSHEV INTERPOLATION

Chebyshev Interpolation, on the other hand, uses Chebyshev nodes as data points to 
construct the polynomial function.

COMPARISON

To compare the performance of these two methods, we used the function f(x) = e|x|, which we approximated using polynomial
functions of degree 10 and 20. We then plotted the true function, the polynomial function obtained using Evenly Spaced 
Interpolation, and the polynomial function obtained using Chebyshev Interpolation for each degree.

We also plotted the empirical errors for each method, which represent the difference between the true function and the 
approximated function for each value of x. These errors were plotted against the x values for each method.

Our results show that Chebyshev Interpolation performed better than Evenly Spaced Interpolation in approximating the 
function, especially for higher degrees of the polynomial. Additionally, we observed that the Runge phenomenon was 
present in the Evenly Spaced Interpolation method, which resulted in higher errors for certain values of x.

In conclusion, our project highlights the importance of choosing an appropriate method for polynomial interpolation 
based on the degree of the polynomial and the distribution of the data points. Chebyshev Interpolation is a more 
effective method than Evenly Spaced Interpolation for approximating functions using high degree polynomials.

Thank you for listening to my presentation.


"""


"""

Answer for the question: Can the Runge phenomenon be observed in this problem?

Yes, the Runge phenomenon can be observed in this problem. 
The Runge phenomenon is a phenomenon in numerical analysis 
where the error in polynomial interpolation can oscillate 
widely near the edges of the interpolation interval, even as 
the number of interpolation points increases. In this problem, 
the function being interpolated is the exponential function of 
the absolute value of x, which has a sharp change in curvature 
near x=0. When interpolating this function using high-degree 
polynomials and evenly spaced interpolation nodes, the error 
can oscillate widely near the edges of the interpolation interval, 
as can be seen in the plots of the empirical errors for n=10 and n=20. 
This oscillation is a manifestation of the Runge phenomenon. However, 
when using Chebyshev interpolation nodes, the error near the edges of 
the interpolation interval is much smaller and does not exhibit the 
oscillation observed with evenly spaced nodes. 
This is because Chebyshev nodes are spaced more densely near the edges of 
the interval, which helps to reduce the effect of the sharp change 
in curvature near x=0.

"""