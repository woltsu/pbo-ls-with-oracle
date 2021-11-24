import os
from numpy.random import choice


def weighted_shuffle(arr, weights):
    return list(choice(a=arr, p=weights, replace=False, size=len(arr)))


def cost(T, model, C_map):
    r = 0
    for l in T:
        if model[l - 1] > 0:
            r += C_map[l]
    return r


def is_good(c, l):
    if c < 0 and l > 0:
        return True
    elif c < 0 and l < 0:
        return False
    elif c > 0 and l < 0:
        return True
    elif c > 0 and l > 0:
        return False


def is_sat(result):
    return result[2] == 1


# Put the literals with the heaviest weights to the front
# TODO: Maybe test other way around
def calculate_weights(T, C_map):
    weights = [0] * len(T)
    for (i, l) in enumerate(T):
        weights[i] = abs(C_map[l])

    # Normalize weights
    total_weights = 0
    for w in weights:
        total_weights += w

    for (i, w) in enumerate(weights):
        weights[i] = w / total_weights

    return weights


def split_to_good_and_bad(T, C_map, model):
    good = []
    bad = []
    for l in T:
        l_idx = l - 1
        if is_good(C_map[l], model[l_idx]):
            good.append(l)
        else:
            bad.append(l)
    return [good, bad]


def load_input(in_file):
    if not os.path.exists(in_file):
        print("cannot open path \""+in_file+"\"")
        exit(1)

    print("c solving file: {}".format(in_file))

    lines = []
    constraints = []
    objectiveline = None
    nextVar = 0
    origMaxVar = 0
    with open(in_file, 'r') as file:
        prev = ""
        for l in file:
            if l.startswith("*"):
                # comments/header doesn't end with ;
                lines.append(l)
                continue
            for kw in l.split():
                if kw[0] == 'x':
                    nextVar = max(nextVar, int(kw[1:]))
            if len(prev) > 0:
                # append to previous line
                l = prev + l
            if l.strip()[-1] != ";":
                # line continues
                prev = l.strip() + " "
                continue
            if l.startswith("min") or l.startswith("max"):
                objectiveline = l
            else:
                lines.append(l)
                sides = l[:-1].split("=")
                equality = sides[0][-1] != ">"
                rhs = int(sides[1][:-1])
                wrds = sides[0].split()
                if wrds[-1] == '>':
                    wrds = wrds[:-1]
                ind = list(map(lambda x: int(x[1:]), wrds[1::2]))
                coefs = list(map(int, wrds[::2]))
                constraints.append((dict(zip(ind, coefs)), rhs))
                if equality:
                    constraints.append(
                        (dict(zip(ind, [-v for v in coefs])), -rhs))
            prev = ""
    nextVar += 1

    # Parse objective line
    objectiveline = objectiveline[4:-2].strip().split()
    origMaxVar = nextVar
    objvars = [int(t[1:]) for t in objectiveline[1::2]]
    objcoefs = list(map(int, objectiveline[0::2]))

    return origMaxVar, objvars, objcoefs, constraints
