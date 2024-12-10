import random 
import matplotlib.pyplot as plt
import numpy as np
from stationary import stationary_probs
from te import *
    
def nrst_random_walk(ap1, am1, b = 100):
    assert(len(ap1) == len(am1))
    assert(ap1[-1] == 0)
    assert(am1[0] == 0)

    # Initial values
    l = 0
    epsilon = 1
    l_values = []
    e_values = []
    for _ in range(b):
        if (epsilon == 1):
            # Moving up
            if (random.random() < ap1[l]):
                # Success
                l = l + 1
            else:
                epsilon = -1
        else:
            # Going down
            if (random.random() < am1[l]):
                # Success
                l = l-1
            else:
                epsilon = 1
        l_values.append(l)
        e_values.append(epsilon)
    return (l_values, e_values)


def st_random_walk(ap1, am1, b = 100):
    assert(len(ap1) == len(am1))
    assert(ap1[-1] == 0)
    assert(am1[0] == 0)
    levels = len(ap1)

    # Initial values
    l = 0
    epsilon = 1
    l_values = []
    e_values = []
    for _ in range(b):
        if random.random() < 1/2:
            # Go up
            if random.random() < ap1[l]:
                l += 1
                epsilon = 1
        else:
            # Go down
            if random.random() < am1[l]:
                l -= 1
                epsilon = -1
        l_values.append(l)
        e_values.append(epsilon)
    return (l_values, e_values)

def main():
    levels = 3
    np.random.seed(547)
    ap1 = np.random.rand(levels)
    ap1[-1] = 0
    am1 = np.random.rand(levels)
    am1[0] = 0
    levels = len(am1)
    b = 100000

    # NRST
    l_values, e_values = nrst_random_walk(ap1, am1, b)
    # plt.figure()
    # plt.plot(l_values)
    # plt.hist(l_values)
    affinities = stationary_probs(ap1, am1).sum(axis=0)
    # plt.plot(range(levels), affinities*b)
    print("NRST")
    print(f"TE is {te(l_values, e_values, levels-1)}")
    print(f"Theoretical TE is {theory_te(ap1, am1)}")
    # print(f"Emprical level affinities are {affinities}")

    # ST
    l_values, e_values = st_random_walk(ap1, am1, b)
    # plt.figure()
    # plt.plot(l_values)
    # plt.hist(l_values)
    # affinities = stationary_probs(ap1, am1).sum(axis=0)
    # plt.plot(range(levels), affinities*b)
    print(f"TE is {te(l_values, e_values, levels-1)}")
    # print(f"Emprical level affinities are {affinities}")
    # plt.show()

if __name__ == "__main__":
    main()
