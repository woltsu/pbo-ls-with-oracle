import os
from numpy.random import choice
from time import time
import argparse


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
    weights_per_l = {}
    weights = [0] * len(T)
    for (i, l) in enumerate(T):
        weights[i] = abs(C_map[l])
        weights_per_l[l] = abs(C_map[l])

    # Normalize weights
    total_weights = 0
    for w in weights:
        total_weights += w

    for key in weights_per_l.keys():
        weights_per_l[key] = weights_per_l[key] / total_weights

    return weights_per_l


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
                # Only assumes normalized instances?
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


def load_args():
    parser = argparse.ArgumentParser(
        description='PBO-#oracle solver', formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('instance', help='instance file name')
    parser.add_argument('--timeout', type=int, default=60,
                        help='Timeout after given seconds')
    parser.add_argument('--verbose', action='store_true',
                        help='Print debug info')
    parser.add_argument('--baseline', action='store_true',
                        help='Run without any semantics')
    parser.add_argument('--print-solution', action='store_true',
                        help='Print solution if any')
    args = parser.parse_args()

    seed = int(time())

    print("c PBO-#with-oracle")
    print(f"c using seed {seed}")
    for k, v in args.__dict__.items():
        if k in ["instance"]:
            continue
        print("c {}: {}".format(k, v))
    print("c ---------------------------------")

    in_file = args.instance
    verbose = args.verbose
    baseline = args.baseline
    timeout = args.timeout
    print_solution = args.print_solution

    return in_file, verbose, baseline, timeout, print_solution, seed


def init_timer():
    start = time()

    def time_print(*out):
        print("t {:.3f}s - {}".format(time() -
                                      start, " ".join(list(map(str, out)))))

    return time_print
