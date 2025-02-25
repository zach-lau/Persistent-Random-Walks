import numpy as np
from mytest import Test

def te(level, momentum, target):
    if (target == 0):
        raise ValueError("Cannot have target level 0")
    tl = [] # tour lengths
    visits = [] # number of vists

    this_visits = 0 # Visits this tour
    start = -1 # start of the tour according to Biron-Lattes et al
    for i,state in enumerate(zip(level,momentum)):
        l,e = state
        if l == target:
            this_visits = this_visits + 1
        if l == 0 and e == -1: # atom
            tl.append(i-start) # inclusive
            visits.append(this_visits)
            start = i # start on next one
            this_visits = 0
    k = len(tl)
    if (k == 0): # Didn't finish a single tour 
        return None # NA basically 
    if (sum(visits) == 0): # no visits at all
        return 0
    return sum(visits)**2/k/sum([x**2 for x in visits])

def end_probs(a_up, a_down):
    """
    Chance of transitioning to the target layer before transitioning back to
    the atom from the given state
    """
    n = len(a_up)
    assert(len(a_down) == n)
    P = np.zeros((2*n,2*n))
    for i, (up, down) in enumerate(zip(a_up, a_down)):
        # 2*i states are up
        # 2*i+1 states are down
        if i < n-1: # not last layer
            P[2*i,2*i+2] = up
        P[2*i,2*i+1] = 1-up
        if i > 0: # not first layer
            P[2*i+1,2*i-1] = down
        P[2*i+1,2*i] = 1-down
    # print(P)
    # Make our atom and our end absorbing states
    # print(P)
    b = P[:,-2].copy()
    P[:,1] = 0 # fail if we transition to start again
    P[:,-2] = 0 # dont double count tranisitions into target
    return np.linalg.solve(np.eye(2*n)-P,b)

def top_bottom_probs(a_up, a_down):
    x = end_probs(a_up, a_down)
    pup = x[1]
    pdown = 1-x[-1]
    return (pup, pdown)

def fast_tb_probs(a_up, a_down, f=2):
    """
    a_up: upwards transition acceptance rates
    a_down: downwarsd transition acceptance rates
    f: nrst or st. 1 for st 2 for nrst
    """
    # Coefficients for nrst vs st
    if f==1:
        k0 = 2 # constant
        k1 = 2
    elif f==2:
        k0 = 0
        k1 = 1
    # TODO deal with invalid f values

    # Calculate the probability of transitions up and down 
    n = len(a_up)
    assert(n > 1)
    assert(len(a_down) == n)
    pdinv = 1/a_down[1]
    # Start at the bottom and work up
    for i in range(2,n):
        # pdinv = a_up[i-1]/a_down[i]*((1-a_up[i-1])/a_up[i-1]+pdinv)
        pdinv = a_up[i-1]/a_down[i]*(pdinv+k0) + (1-a_up[i-1])/a_down[i]*k1
    # Start at the to and work down
    puinv = 1/a_up[n-2]
    for j in range(n-3,-1,-1): # work down to 0
        # puinv = a_down[j+1]/a_up[j]*((1-a_down[j+1])/a_down[j+1]+puinv)    
        puinv = a_down[j+1]/a_up[j]*(puinv+k0) + (1-a_down[j+1])/a_up[j]*k1

    return (1/puinv, 1/pdinv)

def theory_te(a_up, a_down, f=2):
    """
    Find the theoretical te for given acceptance rates

    f = 2 for nrst and 1 for st
    """
    # p = end_probs(a_up, a_down) # probability of reaching the end first
    # pup = p[1] # chance of making it to last level before coming back
    # pdown = 1-p[-2] # chance of making it down before coming back
    # ev = f*pup/pdown # bernouilli x geometric
    # ev2 = (f**2)*pup*(1/pdown**2 + (1-pdown)/pdown**2)
    # return ev**2/ev2
    pup, pdown = fast_tb_probs(a_up, a_down, f)
    return pup/(2-pdown) # s    

def test():
    
    # print(end_probs([0.3,0],[0,0.5])) # easy to solve analytically for 2 states
    # print(end_probs([1,1,0],[0,1,1])) # also easy for a deterministic tour
    # print(end_probs([0.5,0.3,0.2,0],[0,0.7,0.8,0.9])) # harder
    # print(end_probs([0.5,0.5,0],[0,0.5,0.5]))
    # print(theory_te([0.5,0.5,0],[0,0.5,0.5])) # we know from the symmetric case we have 0.2
    # print(theory_te([1,1,0],[0,1,1])) # should be 1 for a deterministic walk
    # print(theory_te([1,1,0],[0,1,1],1)) # pup = pdown = 1/4, should be 1/5

    tests = [
        # 2x2 case is easy
        Test("2x2", np.allclose, end_probs([0.3,0],[0,0.5]), [0.3, 0.3, 0.5, 0.5]),
        # Deterministic tour
        Test("", np.allclose, end_probs([1,1,0],[0,1,1]), [1,1,1,0,0,0]),
        Test("", np.allclose, theory_te([1,1,0],[0,1,1]), 1),
        # We also know the answer for the symmetric case
        Test("Symmetric", np.allclose,theory_te([0.5,0.5,0],[0,0.5,0.5]),.2),
        # Isolate top and bottom probs
        Test("Top bottom", np.allclose, top_bottom_probs([1,1,0],[0,1,1]), (1,1)),
        # Fast version
        Test("Deterministic",
             np.allclose, fast_tb_probs([1,1,0],[0,1,1]), (1,1)),
        Test("Matches with full solution", np.allclose,
             fast_tb_probs([0.5,0.3,0],[0,0.7,0.8]),
             top_bottom_probs([0.5,0.3,0],[0,0.7,0.8])),
        Test("Big test", np.allclose,
             fast_tb_probs([0.5,0.5,0.5,0.5,0.3,0],[0,0.7,0.7,0.7,0.7,0.8]),
             top_bottom_probs([0.5,0.5,0.5,0.5,0.3,0],[0,0.7,0.7,0.7,0.7,0.8])),
        Test("ST", np.allclose,theory_te([1,1,0],[0,1,1],1),0.2)
    ]
    for i,t in enumerate(tests):
        print(f"[{i+1}/{len(tests)}] ({t.get_description()})")
        t.check()

if __name__ == "__main__":
    test()
