import math
##### Problem 3d #####
def algo(x, num_terms):
    y = 0
    for i in range(1, num_terms):
        y += (x**i / math.factorial(i))
    
    return y

#x = 9.999999995000000 * 10**-10
x = 10**-9
print(algo(x, 5))
##### END #####
