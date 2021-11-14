#!/usr/bin/env python3
import os
import signal
from time import time
import roundingsat
import argparse
import random
import copy
from numpy.random import choice

start = 0
total_hs_time = 0
total_core_time = 0
total_minimize_time = 0
total_abs_time = 0

coeff_sum = 0

UB = float('inf')
best_model = []
abs_cores = []
countVarMap = {}


def time_print(*out):
    print("c {:8.3f}s - {}".format(time() -
                                   start, " ".join(list(map(str, out)))))


class TimeoutHandler:
    def __init__(self, rsat):
        self.rsat = rsat

    def __call__(self, signo, frame):
        if len(best_model) > 0:
            print("s SATISFIABLE")
            print("v ", end='')
            for mm in best_model:
                print("{} ".format(mm), end='')
        else:
            print("s UNKNOWN")


def sign(x): return -1 if x < 0 else 1


def weighted_shuffle(arr, weights):
    return list(choice(a=arr, p=weights, replace=False, size=len(arr)))


def cost(T, model, C_map):
    r = 0
    for l in T:
        if model[l - 1] > 0:
            r += C_map[l]
    return r


# TODO: Handle max instances other way around?
def is_good(c, x):
    if c < 0 and x > 0:
        return False
    elif c < 0 and x < 0:
        return True
    elif c > 0 and x < 0:
        return True
    elif c > 0 and x > 0:
        return False


def is_sat(result):
    return result[2] == 1


# We want to have "good" literals in front and "bad" literals at back
def calculate_weights(T, C_map):
    weights = [0] * len(T)
    for (i, l) in enumerate(T):
        weights[i] = 1 / C_map[l]

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


def solve(model, T, C_map, best):
    best_copy = copy.copy(best)
    final_result = None
    assumptions = []

    for i in range(len(T)):
        l = T[i]
        l_idx = l - 1

        if is_good(C_map[l], best_copy[l_idx]):
            assumptions.append(best_copy[l_idx])
        else:
            tmp_assumptions = copy.copy(assumptions)
            tmp_assumptions.append(-best_copy[l_idx])
            model.solve(tmp_assumptions, 0)
            final_result = model.getResult()

            if not is_sat(final_result):
                assumptions.append(best_copy[l_idx])
            else:
                best_copy = final_result[1]
                (good, bad) = split_to_good_and_bad(T[i+1:], C_map, best_copy)
                T = T[:i+1] + good + bad
                assumptions = tmp_assumptions

    if final_result is not None and not is_sat(final_result):
        return best

    return best_copy


# TODO: Change to time limit
K = 2500


def solve_inc(model, T, C_map, debug=True):
    model.solve([], 0)
    result = model.getResult()
    weights = calculate_weights(T, C_map)
    T_copy = copy.copy(T)

    if not is_sat(result):
        return "UNSAT"

    best = copy.copy(result[1])

    for _ in range(K):
        shuffled_T = weighted_shuffle(T_copy, weights)
        tmp_best = solve(model, shuffled_T, C_map, best)
        if debug:
            print("cost(T, tmp_best)", cost(T, tmp_best, C_map),
                  ", cost(T, best)", cost(T, best, C_map))
        if cost(T, tmp_best, C_map) < cost(T, best, C_map):
            best = copy.copy(tmp_best)

    return best


# TODO: Max to min -> flip coefficents (* -1) and then flipbits. Deflip result.
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='PBO-#oracle solver', formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('instance', help='instance file name')
    args = parser.parse_args()

    start = time()

    print("c PBO-#ihs")
    in_file = args.instance

    if not os.path.exists(in_file):
        print("cannot open path \""+in_file+"\"")
        exit(1)

    print("c solving file: {}".format(in_file))
    for k, v in args.__dict__.items():
        if k in ["instance"]:
            continue
        print("c {}: {}".format(k, v))
    print("c ---------------------------------\n")

    #
    # Read opb file and strip out objective function
    #
    rsat = roundingsat.Roundingsat()

    signal.signal(signal.SIGINT, TimeoutHandler(rsat))
    signal.signal(signal.SIGTERM, TimeoutHandler(rsat))
    signal.signal(signal.SIGXCPU, TimeoutHandler(rsat))

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
    origMaxVar = nextVar
    transformedObjectiveline = objectiveline[:-1].replace('-', '+')

    #
    # Parse objective line
    #
    objectiveline = objectiveline[4:-2].strip().split()

    objvars = [int(t[1:]) for t in objectiveline[1::2]]
    objcoefs = list(map(int, objectiveline[0::2]))
    negative_coeffs = []
    smallest_coeff = min(objcoefs)
    # TODO: Does this turn max instances to min instances?
    # If --fliplits is on, each literal that was flipped in objective function is replaced: l -> (1 - l)
    """ if args.fliplits:
        for i in range(len(objvars)):
            if objcoefs[i] < 0:
                negative_coeffs.append(objvars[i])
                objcoefs[i] = abs(objcoefs[i])
                coeff_sum += objcoefs[i]
        for idx,constraint in enumerate(constraints):
            offset = constraint[1]
            newconstraint = {}
            for var,coef in constraint[0].items():
                if var in negative_coeffs:
                    offset += -coef
                    newconstraint[var] = -coef
                else:
                    newconstraint[var] = coef
            constraints[idx] = (newconstraint, offset)
        print("c PBO instance transformed, objective offset {}".format(coeff_sum)) """
    obj_var_coef = {var: coef for (var, coef) in zip(objvars, objcoefs)}

    rsat.init(origMaxVar)
    set_objvars = set(objvars)
    for idx, c in enumerate(constraints):
        rsat.addConstraint(c[0], c[1])
        constraints[idx] = (c[0], c[1])

    # rsat.print()
    C_map = {}
    for i in range(len(objvars)):
        C_map[objvars[i]] = objcoefs[i]

    result = solve_inc(rsat, objvars, C_map, debug=True)
    print("RESULT:", cost(objvars, result, C_map))
    # print(result)
    # TODO: Print cost / calculate cost
    # TODO: Exclude non variables from result

    # rsat.print()
    # rsat.solve([], 0)
    # result = rsat.getResult()
    # print(result)
