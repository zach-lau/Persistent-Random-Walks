"""
Shortcut algorithm to find the stationary distribution of an NRST chain based
on transition probabilities
"""

import numpy as np
from test import Test

def _stationary_probs_helper(a_up, a_down, ans, i, j):
    """Recursive helper for stationary probs that solves for the chain
    just in the region [i,j]. Solves for the restriction to this region
    
    a_up: list of up acceptance rates
    a_down: list of down acceptance rates
    ans: 2xn array to store answer
    i: first level we are concerned with
    j: one more than last level we are concerned with
    """
    # Base case
    if (i > j): # empty array
        return
    if i == j: # singleton
        # Just populate accordingn to balance equations truncate to these two
        # states
        rup = 1-a_up[i]
        rdown = 1-a_down[i]
        ans[0,i] = rdown/(rdown+rup)
        ans[1,i] = rup/(rdown+rup)
        return
    # Recursive stage
    mid_point = (i+j)//2
    _stationary_probs_helper(a_up, a_down, ans, i, mid_point)
    _stationary_probs_helper(a_up, a_down, ans, mid_point+1, j)
    # Would need to work out this recursive step but seems like too much
    # linear algebra to be worth it
    # on the transition rates up and down
    # rup = ans[0,mid_point-1]*a_up[mid_point-1] # rate going up
    # rdown = ans[1,mid_point]*a_down[mid_point]
    # pleft = rdown/(rup+rdown)
    # pright = rup/(rup+rdown)
    # # Scale proportionally
    # ans[:,i:mid_point] *= pleft
    # ans[:,mid_point+1:j] *= pright

def stationary_probs(a_up, a_down):
    """
    Find stationary probs for NRST given up acceptance rates and down 
    acceptance rates in nln n time
    """
    assert(len(a_up) == len(a_down))
    n = len(a_up)
    ans = np.zeros((2,n))
    _stationary_probs_helper(a_up, a_down, ans, 0, n-1)
    return ans

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
        # Test(
        #         "Faster algorithm",
        #         np.allclose,
        #         stationary_probs(np.array([0.5,0.5,0]),np.array([0,0.5,0.5])),
        #         np.array([[1/6,1/6,1/6],[1/6,1/6,1/6]])
        # )
    ]
    for i,t in enumerate(tests):
        print(f"[{i+1}/{len(tests)}]")
        t.check()

if __name__ == "__main__":
    test()