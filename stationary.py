"""
Shortcut algorithm to find the stationary distribution of an NRST chain based
on transition probabilities
"""

import numpy as np
from mytest import Test

# def _stationary_probs_helper(a_up, a_down, ans, i, j, total_mass):
#     """Recursive helper for stationary probs that solves for the chain
#     just in the region [i,j]. Solves for the restriction to this region
    
#     a_up: list of up acceptance rates
#     a_down: list of down acceptance rates
#     ans: 1xn array to store answer
#     i: first level we are concerned with
#     j: one more than last level we are concerned with
#     """
#     # Base case
#     if (i > j): # empty array
#         return
#     if i == j: # singleton
#         ans[i] = total_mass
#     # Recursive stage
#     mid_point = (i+j)//2
#     _stationary_probs_helper(a_up, a_down, ans, i, mid_point)
#     _stationary_probs_helper(a_up, a_down, ans, mid_point+1, j)
#     # on the transition rates up and down


# def probs_stable(a_up, a_down):
#     """ Stable but reasonably fast algorithm to find stsationary probs"""


def stationary_probs(a_up, a_down):
    """
    Find stationary probs for NRST given up acceptance rates and down 
    acceptance rates in linear time
    """
    assert(len(a_up) == len(a_down))
    n = len(a_up)
    # This might have poor numerical stability but it is fast
    ans = np.zeros(n)
    ans[0] = 1
    for i in range(1,n):
        ans[i] = ans[i-1]*a_up[i-1]/a_down[i]
    ans = ans/sum(ans)/2
    return np.array([ans, ans])

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
            ),
        Test(
            "Faster algorithm",
            np.allclose,
            stationary_probs(np.array([0.5,0.5,0]),np.array([0,0.5,0.5])),
            np.array([[1/6,1/6,1/6],[1/6,1/6,1/6]])
        ),
        Test(
            "",
            np.allclose,
            stationary_probs(np.array([0.3,0.2,0]),np.array([0,0.8,0.9])),
            naive_probs([0.3,0.2,0],[0,0.8,0.9])
        )
    ]
    for i,t in enumerate(tests):
        print(f"[{i+1}/{len(tests)}]")
        t.check()

if __name__ == "__main__":
    test()
