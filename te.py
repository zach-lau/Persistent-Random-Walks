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
    if (visits == 0):
        return 0
    return sum(visits)**2/k/sum([x**2 for x in visits])
