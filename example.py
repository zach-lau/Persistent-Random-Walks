"""
Example calculations for the TE of multimodal distribution provided by
Syed et al.
"""

import numpy as np
from mytest import Test
from matplotlib import pyplot as plt

def find_one(dc, db, a, k, beta):
    """
    Calculate acceptance rates for multimodal example in Syed et al.
    dc : c'-c
    db: beta'-beta. In (0,1)
    a: constant defining distribution. Requires a > 1
    k: integer defining target
    beta: current level. Between 0 and 1
    """
    podd = k/(k+a**beta*(k+1))
    peven = a**beta*(k+1)/(k+a**beta*(k+1))
    # Fine one instance of alpha on a case by case basis
    if db > 0: # going up
        if dc > 0:
            return 1
        elif dc > -db*np.log(a):
            return podd*np.exp(dc)+peven
        else:
            return podd*np.exp(dc)+peven*np.exp(db*np.log(a)+dc)
    else: # going down
        if dc < 0:
            return peven*np.exp(db*np.log(a)+dc) + podd*np.exp(dc)
        elif dc < -db*np.log(a):
            return peven*np.exp(db*np.log(a)+dc) + podd
        else:
            return 1

def find_alphas(dc, a, k, beta):
    """
    dc: vector of differences in level affinity tuning parameters.
    The differences uniquely specify the system so we use these
    k: the size of the system
    a: tuning parameter > 1
    
    return: a matrix of alpha values corresponding to the forward and
    backwards acceptance rates
    """
    n = len(dc)
    alpha = np.zeros((2,n+1)) # 0 to n inclusive

    for i in range(1,n): # interior points
        alpha[0,i] = find_one(dc[i], beta[i+1]-beta[i], a, k, beta[i])
        alpha[1,i] = find_one(-dc[i-1], beta[i-1]-beta[i], a, k, beta[i]) 
    # End points
    alpha[1,0] = 0 # down from 0
    alpha[0,n] = 0 # up from n 
    alpha[0,0] = find_one(dc[0], beta[1]-beta[0], a, k, beta[0]) # up from 0
    alpha[1,n] = find_one(-dc[n-1], beta[n-1]-beta[n], a, k, beta[i])  # down from n
    return alpha
 
def test():
    tests = [
        Test("", np.allclose, find_one(1, 0.1, 2, 10, 0.1), 1)
    ]
    for i,t in enumerate(tests):
        print(f"[{i+1}/{len(tests)}]")
        t.check()

def graph():
    # Test to graph the up and down probabilities as a function of c
    db = 0.1
    beta = 0.5
    a = 2
    k = 10
    dcs = np.linspace(-2*db*np.log(a), 2*db*np.log(a))
    up = [find_one(dc, db, a, k, beta) for dc in dcs]
    down = [find_one(-dc, -db, a, k, beta+db) for dc in dcs]
    # plt.plot(up, down)
    plt.plot(dcs, up)
    plt.plot(dcs, down)
    plt.show()

if __name__ == "__main__":
    test()
    # print(find_one(-1, -0.1, 2, 10, 0.1)) 
    # graph()
    dc = [0,0]
    betas = [0, 0.5, 1]
    a = 2
    k = 10
    print(find_alphas(dc, a, k, betas))
