#!/usr/bin/env python3
import os
import signal
from time import time
import roundingsat
import argparse
import random
import copy

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


def cost(T, model):
    r = 0
    for l in T:
        if model[l] > 0:
            r += 1
    return r


def solve(model, T, best):
    final_result = None
    assumptions = []
    for i in range(len(T)):
        # TODO: If max, then > 0
        if best[T[i] - 1] < 0:
            assumptions.append(T[i])
        else:
            tmp_assumptions = copy.copy(assumptions)
            tmp_assumptions.append(-T[i])
            model.solve(tmp_assumptions, 0)
            final_result = model.getResult()
            if final_result[2] == 2:
                assumptions.append(T[i])
            else:
                assumptions = tmp_assumptions
    return final_result[1]


K = 1000


def solve_inc(model, T):
    model.solve([], 0)
    result = model.getResult()

    if result[2] == 2:
        return "UNSAT"

    best = copy.copy(result[1])

    for _ in range(K):
        tmp_best = solve(model, T, best)
        # print("cost(T, tmp_best)", cost(T, tmp_best),
        #      ", cost(T, best)", cost(T, best))
        if cost(T, tmp_best) < cost(T, best):
            best = copy.copy(tmp_best)
        random.shuffle(T)
        # model.clearLearnedConstraints()  # Is this wanted?

    return best


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='PBO-#ihs solver', formatter_class=argparse.RawTextHelpFormatter)
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

    # LOOPPAA NÄIDEN YLI
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

    rsat.print()
    result = solve_inc(rsat, objvars)

    # rsat.print()
    #rsat.solve([], 0)
    #result = rsat.getResult()
    # print(result)
