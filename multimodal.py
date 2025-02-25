"""
Implement the multimiodal exapmle from Syed et al.
"""

import numpy as np
from mytest import Test

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

def find_alphas(dc, a, k, beta, find_one=find_one):
    """
    dc: vector of differences in level affinity tuning parameters.
    The differences uniquely specify the system so we use these
    k: the size of the system
    a: tuning parameter > 1
    find_one: function to find a single upwards or downwards transition rate
    
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
    alpha[1,n] = find_one(-dc[n-1], beta[n-1]-beta[n], a, k, beta[n])  # down from n
    return alpha
 
def unit_test():
    tests = [
        Test("", np.allclose, find_one(1, 0.1, 2, 10, 0.1), 1)
    ]
    for i,t in enumerate(tests):
        print(f"[{i+1}/{len(tests)}]")
        t.check()

if __name__ == "__main__":
    unit_test()
