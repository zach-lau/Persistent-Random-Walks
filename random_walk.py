import random 
import matplotlib.pyplot as plt
import numpy as np

def te(level, momentum, target):
    if (target == 0):
        raise ValueError("Cannot have target level 0")
    tl = [] # tour lengths
    visits = [] # number of vists

    this_visits = 0 # Visits this tour
    start = 0 # start of the tour
    for i,state in enumerate(zip(level,momentum)):
        l,e = state
        if l == target:
            this_visits = this_visits + 1
        if l == 0 and e == -1: # end of tour
            tl.append(i-start+1) # inclusive
            visits.append(this_visits)
            start = i+1 # start on next one
            this_visits = 0
    k = len(tl)
    if (k == 0):
        return None # NA basically
    return sum(visits)**2/k/sum([x**2 for x in visits])
    
def main():
    levels = 10
    # np.random.seed(547)
    # ap1 = np.random.rand(levels)
    # ap1[-1] = 0
    # am1 = np.random.rand(levels)
    # am1[0] = 0
    ap1 = [1]*levels
    ap1[-1] = 0
    am1 = [1]*levels
    am1[0] = 0

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
    print(f"TE is {te(l_values, e_values, levels-1)}")
    plt.show()

if __name__ == "__main__":
    main()
