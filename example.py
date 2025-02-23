"""
Example calculations for the TE of multimodal distribution provided by
Syed et al.
"""

import numpy as np
from mytest import Test
from matplotlib import pyplot as plt
import te
from scipy.optimize import minimize, fsolve
import stationary

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
    alpha[1,n] = find_one(-dc[n-1], beta[n-1]-beta[n], a, k, beta[n])  # down from n
    return alpha
 
def unit_test():
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

def two_layer_test(dc=[0,0], a = 2):
    # print(find_one(-1, -0.1, 2, 10, 0.1)) 
    # graph()
    # dc = [0,0]
    betas = [0, 0.5, 1]
    k = 10
    alpha = find_alphas(dc, a, k, betas)
    # print(alpha)
    tour_eff = te.theory_te(alpha[0,],alpha[1,])
    print(tour_eff)

def n_layer_test(n):
    dc = np.zeros(n)
    betas = np.linspace(start=0, stop=1, num=n+1)
    a = 2
    k = 10
    alpha = find_alphas(dc, a, k, betas)
    # print(alpha)
    tour_eff = te.theory_te(alpha[0,],alpha[1,])
    print(tour_eff)

def optimize_grid():
    betas = [0, 0.5, 1]
    a = 10
    k = 10
    xs = np.linspace(-2*np.log(a), 2*np.log(a))
    ys = xs.copy() # deep copy
    zs = np.zeros((len(xs), len(ys)))
    for i,x in enumerate(xs):
        for j,y in enumerate(ys):
            alpha = find_alphas([x,y], a, k, betas)
            zs[j,i] = te.theory_te(alpha[0,],alpha[1,])
    cs = plt.contour(xs, ys, zs)
    plt.clabel(cs)
    uniform_point = (-0.757, -0.986)
    optimal_point = (0, -1.151)
    plt.plot(*uniform_point, 'ro')
    plt.plot(*optimal_point, 'bo')
    plt.xlabel("c1-c0")
    plt.ylabel("c2-c1")
    plt.legend(["Uniform", "Optimal"])
    plt.show()

def optimize_test():
    dcs = [-1,-2]
    betas = [0, 0.5, 1]
    a = 10
    k = 10
    # def to_opt
    def f(x):
        alpha = find_alphas(x, a, k, betas)
        tour_eff = te.theory_te(alpha[0,],alpha[1,])
        return tour_eff
    # Nelder-mead often flops
    # res = minimize(lambda x : -f(x), dcs, method="Nelder-Mead")
    # Default is BFGS : performs reasonably well in this case
    res = minimize(lambda x : -f(x), dcs)
    print(res)
    print(f(res.x))

def find_uniform_affinities():
    print("==== Beginning calculation for uniform affinity tuning params === ")
    a = 10
    k = 10
    betas = [0, 0.5, 1]
    # We could solve a systme of equations, but its easier to just match
    # upwards and dowanrds rates
    n = len(betas)-1 # number of things to solve for
    cs = [0]*n
    for i in range(n):
        db = betas[i+1]-betas[i]
        def f(dc):
            up = find_one(dc, db, a, k, betas[i])
            down = find_one(-dc, -db, a, k, betas[i+1])
            return up-down
        cs[i] = fsolve(f, -1)[0]
        # print(f(cs[i]))
        # print(find_one(cs[i], db, a, k, betas[i]))
        # print(find_one(-cs[i], -db, a, k, betas[i+1]))
    print(f"The uniform level affinites are given by: {cs}")
    alphas = find_alphas(cs, a, k, betas)
    print(f"We get acceptance rates: {alphas}")
    sprobs = stationary.stationary_probs(alphas[0,], alphas[1,]).sum(axis=0)
    print(f"The stationary probabilties are: {sprobs}")
    myte = te.theory_te(alphas[0,], alphas[1,])
    print(f"The te is: {myte}")

def explore_optimized_region():
    dcs = [0, -1.151]
    a = 10
    k = 10
    betas = [0, 0.5, 1]
    alphas = find_alphas(dcs, a, k, betas)
    print(f"The acceptance rates are: {alphas}")
    sprobs = stationary.stationary_probs(alphas[0,], alphas[1,]).sum(axis=0)
    print(f"The level affinities are: {sprobs}") # these are non statinoary!

if __name__ == "__main__":
    unit_test()
    # graph()
    # two_layer_test([0,0]) # 0.744
    # two_layer_test([-0.303, -0.055]) # 0.712
    # two_layer_test([0,0], a = 10) # 0.548
    # two_layer_test([4,4], a = 10) # 0.5 (i.e. p up = 1 p down = 0)
    # two_layer_test([0, -1], a = 10) # 0.596
    # n_layer_test(10) # 0.744
    optimize_test() # gives best result at 0, -1.151 with value 0.619
    explore_optimized_region()
    find_uniform_affinities() # -0.757, -0.986
    optimize_grid()