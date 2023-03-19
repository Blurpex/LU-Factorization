"""
This line creates a new Python class named LUFactorization.
"""


class LUFactorization:
    """
    This is the constructor method of the LUFactorization class.
    It takes a single argument A which is a matrix. The method
    initializes some instance variables including the matrix A,
    the size of the matrix n, and the matrices L and U.
    It then calls the decompose() method to perform the LU decomposition.
    """

    def __init__(self, A):
        self.A = A
        self.n = len(A)
        self.L = [[0.0] * self.n for i in range(self.n)]
        self.U = [[0.0] * self.n for i in range(self.n)]
        self.decompose()

    """
    This method performs the LU decomposition of the matrix A. 
    It first initializes the first row of the matrix U and the first 
    column of the matrix L. It then loops over the remaining rows and 
    columns of U and L, computing the values for each element based on 
    the previous elements. The sum() function is used to compute the dot products.
       
    """

    def decompose(self):
        for j in range(self.n):
            self.U[0][j] = self.A[0][j]
            self.L[j][0] = self.A[j][0] / self.U[0][0]

        for i in range(1, self.n):
            for j in range(i, self.n):
                s1 = sum(self.U[k][j] * self.L[i][k] for k in range(i))
                self.U[i][j] = self.A[i][j] - s1

                s2 = sum(self.U[k][i] * self.L[j][k] for k in range(i))
                self.L[j][i] = (self.A[j][i] - s2) / self.U[i][i]

    """
    This method multiplies the matrices L and U together to produce 
    the original matrix A. It initializes a new matrix LU to hold 
    the result and then loops over each element, computing the dot 
    product of the corresponding rows and columns of L and U.
    """

    def multiply_LU(self):
        LU = [[0.0] * self.n for i in range(self.n)]
        for i in range(self.n):
            for j in range(self.n):
                s = 0.0
                for k in range(self.n):
                    s += self.L[i][k] * self.U[k][j]
                LU[i][j] = s
        return LU

    """
    This method solves Lc = b for c of the lower triangular matrix and then use the 
    result to solve Ux = c for x  and then return the solution
    """

    def solve_system(self, b):
        # solve Lc = b for b
        y = [0.0] * self.n
        y[0] = b[0] / self.L[0][0]
        for i in range(1, self.n):
            s = sum(self.L[i][j] * y[j] for j in range(i))
            y[i] = (b[i] - s) / self.L[i][i]

        # solve Ux = c for x
        x = [0.0] * self.n
        x[self.n - 1] = y[self.n - 1] / self.U[self.n - 1][self.n - 1]
        for i in range(self.n - 2, -1, -1):
            s = sum(self.U[i][j] * x[j] for j in range(i + 1, self.n))
            x[i] = (y[i] - s) / self.U[i][i]

        return x


"""
The time complexity of the LU decomposition algorithm is O(n^3), where n is the 
dimension of the matrix. This is because the decomposition involves three nested 
loops: two for iterating through the rows and columns of the matrix, and one for
iterating through the previously computed values. Since each of these loops 
runs for n iterations, the total time complexity is O(n^3).

The time complexity of the solve function is O(n^2), because it involves two nested
loops, one for solving Ly = b and the other for solving Ux = y. Since each of these
loops runs for n iterations, the total time complexity is O(n^2).

The time complexity of the multiply_LU function is also O(n^3), because it involves
three nested loops, one for iterating through the rows of the L matrix, one for iterating
through the columns of the U matrix, and one for iterating through the values of the 
L and U matrices being multiplied. Since each of these loops runs for n iterations, the 
total time complexity is O(n^3).

"""

# ------------------------ Testing the algorithm with some input-------------
n = int(input("\nPlease input the n value of the nxn matrix A:  "))

A = []

print("\nPlease input matrix A:\n")

for row in range(0, n):
    A.append([])
    for col in range(0, n):
        A[row].append(int(input("a" + str(row + 1) + str(col + 1) + ":  ")))

b = []

print("\nPlease input vector b:\n")
for row in range(0, n):
    b.append(int(input("b" + str(row + 1) + ":  ")))

# :::test cases:::

#A = [[1, 2, -1],
#     [2, 1, -2],
#     [-3, 1, 1]]
#
#b = [3, 3, -6]  # b represents the right hand side matrix

# A = [[1, 1],
#      [3, -4]]
#
# b = [3, 2]

# Create a new LUFactorization object passing as a parameter the matrix A
# To decompose it in Lower(L) and upper(U) Triangular matrix
lufactored = LUFactorization(A)

print("\nL:\n" + str(lufactored.L))
print("\nU:\n" + str(lufactored.U))

# Compute the product of LU to check if LU = A the original matrix
print("\nMultiplying L and U.\nShould equal A:\n" + str(lufactored.multiply_LU()))

# Solve the system and print the result.
print("\nSolution for x:\n" + str(lufactored.solve_system(b)))
