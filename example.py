"""
Example calculations for the TE of multimodal distribution provided by
Syed et al.
"""

import numpy as np
from matplotlib import pyplot as plt
import te
from scipy.optimize import minimize, fsolve
import stationary

from multimodal import find_one, find_alphas

# Simple tests to check everything works
def graph_one_transition():
    # Graph the up and down acceptance probabilities as a function of c
    # for a test case
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
    """
    Find te for 2 layers
    """
    # print(find_one(-1, -0.1, 2, 10, 0.1)) 
    # graph()
    # dc = [0,0]
    betas = [0, 0.5, 1]
    k = 10
    alpha = find_alphas(dc, a, k, betas)
    # print(alpha)
    tour_eff = te.theory_te(alpha[0,],alpha[1,])
    return tour_eff

def n_layer_test(n):
    """
    Find te for n layers evenly spaced
    """
    dc = np.zeros(n)
    betas = np.linspace(start=0, stop=1, num=n+1)
    a = 2
    k = 10
    alpha = find_alphas(dc, a, k, betas)
    # print(alpha)
    tour_eff = te.theory_te(alpha[0,],alpha[1,])
    print(tour_eff)

# Optimization tests

def contour_plot(uniform_point, optimal_point, a = 10, k = 10, f=2):
    """
    Show a contour plot for the two layer multimodal example
    """
    betas = [0, 0.5, 1]
    a = 10
    k = 10
    xs = np.linspace(-2*np.log(a), 2*np.log(a))
    ys = xs.copy() # deep copy
    zs = np.zeros((len(xs), len(ys)))
    for i,x in enumerate(xs):
        for j,y in enumerate(ys):
            alpha = find_alphas([x,y], a, k, betas)
            zs[j,i] = te.theory_te(alpha[0,],alpha[1,], f)
    cs = plt.contour(xs, ys, zs)
    plt.clabel(cs)
    # uniform_point = (-0.757, -0.986)
    # optimal_point = (0, -1.151)
    plt.plot(*uniform_point, 'ro')
    plt.plot(*optimal_point, 'bo')
    plt.xlabel("c1-c0")
    plt.ylabel("c2-c1")
    plt.legend(["Uniform", "Optimal"])
    plt.show()

def optimize_two_layer(a=10, k=10, betas=[0, 0.5, 1], f=2):
    """
    Find the optimal parametrs in the two layer case
    """
    print("=== Optimizing in the two layer case ===")
    dcs = [-1,-2] # starting point
    def get_te(x):
        alpha = find_alphas(x, a, k, betas)
        tour_eff = te.theory_te(alpha[0,],alpha[1,], f)
        return tour_eff
    # Nelder-mead often flops
    # res = minimize(lambda x : -f(x), dcs, method="Nelder-Mead")
    # Default is BFGS : performs reasonably well in this case
    res = minimize(lambda x : -get_te(x), dcs)
    print(res)
    return res.x

def find_uniform_affinities(a=10, k=10, betas=[0, 0.5, 1]):
    """
    Find acceptance rates for uniform level affinities. This should be the
    same for NRST and ST
    """
    print("==== Beginning calculation for uniform affinity tuning params === ")
    # We could solve a systme of equations, but its easier to just match
    # upwards and dowanrds rates if we know them
    n = len(betas)-1 # number of things to solve for
    cs = [0]*n
    for i in range(n):
        db = betas[i+1]-betas[i]
        def f(dc):
            up = find_one(dc, db, a, k, betas[i])
            down = find_one(-dc, -db, a, k, betas[i+1])
            return up-down
        cs[i] = fsolve(f, -1)[0]
    return cs

def explore_optimized_region(dcs, a=10, k=10, betas = [0, 0.5, 1], f=2):
    """
    Exploring optimized region
    """
    print(f"=== Exploring the region near {dcs} ===")
    alphas = find_alphas(dcs, a, k, betas)
    print(f"The acceptance rates are: {alphas}")
    sprobs = stationary.stationary_probs(alphas[0,], alphas[1,]).sum(axis=0)
    print(f"The level affinities are: {sprobs}") # these are non statinoary!
    myte = te.theory_te(alphas[0,], alphas[1,],f)
    print(f"The te is: {myte}")

if __name__ == "__main__":
    # graph_one_transition()

    # print(two_layer_test([0,0])) # 0.744
    # print(two_layer_test([-0.303, -0.055])) # 0.712
    # print(two_layer_test([0,0], a = 10)) # 0.548
    # print(two_layer_test([4,4], a = 10)) # 0.5 (i.e. p up = 1 p down = 0)
    # print(two_layer_test([0, -1], a = 10)) # 0.596

    # n_layer_test(10) # 0.744
    
    # First for NRST
    optimal_nrst = optimize_two_layer() # gives best result at 0, -1.151 with value 0.619
    explore_optimized_region(optimal_nrst)
    uniform_dcs = find_uniform_affinities() # -0.757, -0.986
    explore_optimized_region(uniform_dcs)
    contour_plot(uniform_dcs, optimal_nrst)

    # Then for ST
    optimal_st = optimize_two_layer(f=1) # gives best result at 8.88, 4.3?
    explore_optimized_region(optimal_st, f=1)
    explore_optimized_region(uniform_dcs, f=1)
    contour_plot(uniform_dcs, optimal_st, f=1)