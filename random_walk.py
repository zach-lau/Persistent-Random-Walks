import random 
import matplotlib.pyplot as plt
import numpy as np
from stationary import naive_probs
from te import te
    
def basic_random_walk(ap1, am1):
    assert(len(ap1) == len(am1))
    assert(ap1[-1] == 0)
    assert(am1[0] == 0)
    levels = len(ap1)

    b = 10000
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
    plt.figure()
    plt.plot(l_values)
    plt.figure()
    plt.hist(l_values)
    affinities = naive_probs(ap1, am1).sum(axis=0)
    plt.plot(range(levels), affinities*b)
    print(f"TE is {te(l_values, e_values, levels-1)}")
    print(f"Emprical level affinities are {affinities}")
    plt.show()

def main():
    levels = 5
    np.random.seed(547)
    ap1 = np.random.rand(levels)
    ap1[-1] = 0
    am1 = np.random.rand(levels)
    am1[0] = 0
    basic_random_walk(ap1, am1)

if __name__ == "__main__":
    main()
