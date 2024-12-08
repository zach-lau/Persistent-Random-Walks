"""
Shortcut algorithm to find the stationary distribution of an NRST chain based
on transition probabilities
"""

import numpy as np
from test import Test

def stationary_probs(a_up, a_down):
    """
    Find stationary probs for NRST given up acceptance rates and down 
    acceptance rates. These correspond to a
    """
    return np.array([[0,0,0],[0,0,0]])

def naive_probs(a_up: np.ndarray, a_down: np.ndarray):
    n = len(a_up) + len(a_down) # size of transition matrix
    P = np.zeros((n,n)) # transition probabilites
    # The forward state for level i is at 2*i
    # the backward state is at 2*i+1
    for i,a in enumerate(a_up):
        if 2*i < n-2: # can actually go up
            P[2*i,2*i+2] = a # accept
        P[2*i, 2*i+1] = 1-a # can always move to down state
    for i,a in enumerate(a_down):
        if i > 0:
            P[2*i+1,2*i-1] = a # if we can go down
        P[2*i+1,2*i] = 1-a # we can always go to our up state
    # Last row sums to all one for stationary probs
    P -= np.eye(n) 
    P[:,-1] = np.ones(n)
    # Right hand side of the equation
    b = np.zeros(n)
    b[-1] = 1
    # print(P)
    x = np.linalg.solve(np.transpose(P),b)
    # print(x)
    # Now turn that back into our desired for
    return x.reshape((2,n//2),order="F")

def test():
    tests = [
        Test(
                "",
                np.allclose,
                naive_probs(np.array([0.5,0.5,0]),np.array([0,0.5,0.5])),
                np.array([[1/6,1/6,1/6],[1/6,1/6,1/6]])
            )
    ]
    for i,t in enumerate(tests):
        print(f"[{i+1}/{len(tests)}]")
        t.check()

if __name__ == "__main__":
    test()