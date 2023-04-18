# Project 2
#
# Group 4:
# Shiza Qureshi
# Haris Ahmad
# Zakaria Coulibaly
# Lilac Dixon
# Quazi Maliha

import math

class interpolate:
    # Constructor
    def __init__(self, n):
        self.n = n
        
        self.chebyshev = []
        self.equalSpace = []
        
        # Create n Chebyshev and equally spaced base points
        for i in range(1, n + 1):
            self.chebyshev.append(math.cos(((2 * i - 1) * math.pi) / 2 * n))
            self.equalSpace.append(-1 + ((2 / n - 1) * (i - 1)))
        
        
    @staticmethod
    def function(x):
        return math.exp(abs(x))
 
            