import numpy as np

def find_points(xeval, xint, idx):
    """
    xeval is the larger array of points
    xint is the smaller array of intervals
    idx is the interval index (counting from 1)
    """
    indices = np.where(np.logical_and(xeval <= xint[idx], xeval >= xint[idx-1]))[0]
    return xeval[indices]

# Example run for #1
xeval = np.linspace(0,10,1001)
xint = np.linspace(0,10,11)
# Should output 101 values between 2 and 3
print(find_points(xeval, xint, 3))


def eval_line(x0, fx0, x1, fx1, xeval):
    m = (fx1 - fx0)/(x1 - x0)
    return m*xeval - m*x0 + fx0

# Example run for #2
x0 = 5; fx0 = 2; x1 = 3; fx1 = 7
# Should output 4.5
print(eval_line(x0, fx0, x1, fx1, 4))