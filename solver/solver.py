import sys
import signal
import roundingsat
import argparse
import copy
import random
import numpy as np
from solver_util import weighted_shuffle, cost, is_good, is_sat, calculate_weights, split_to_good_and_bad, load_input


# Global variables
best_model = []
T = []
C_map = {}
betterments = 0
verbose = False
baseline = False


# Outputs result after given timeout
class TimeoutHandler:
    SIGINT = False

    def log_exit(self):
        print("cost:", cost(T, best_model, C_map))
        print("betterments:", betterments)
        sys.exit()

    def __call__(self, signo, frame):
        self.SIGINT = True
        self.log_exit()


def log(*args):
    global verbose

    if verbose:
        print(*args)


def solve(model, T, C_map, best):
    global baseline

    best_copy = copy.copy(best)
    final_result = None
    assumptions = []

    for i in range(len(T)):
        l = T[i]
        l_idx = l - 1

        # i.e. do we want to try to flip the literal or not?
        if is_good(C_map[l], best_copy[l_idx]):
            assumptions.append(best_copy[l_idx])
        else:
            # Solve with new assumption
            tmp_assumptions = copy.copy(assumptions)
            tmp_assumptions.append(-best_copy[l_idx])
            model.solve(tmp_assumptions, 0)
            tmp_result = model.getResult()

            if is_sat(tmp_result):
                best_copy = tmp_result[1]

                if not baseline:
                    # Move good literals to front for skipping
                    # TODO: Maybe test other way around (bad on front)
                    (good, bad) = split_to_good_and_bad(
                        T[i+1:], C_map, best_copy)
                    T = T[:i+1] + good + bad

                assumptions = tmp_assumptions
                final_result = copy.copy(tmp_result)
            else:
                assumptions.append(best_copy[l_idx])

    if final_result is None or not is_sat(final_result):
        return best

    return best_copy


def solve_inc(model, T, C_map):
    global best_model
    global betterments
    global verbose
    global baseline

    model.solve([], 0)
    result = model.getResult()
    weights = calculate_weights(T, C_map)

    if not is_sat(result):
        print("UNSAT")
        exit(1)

    best_model = copy.copy(result[1])

    while True and not TimeoutHandler.SIGINT:
        r = np.random.random(1)[0]
        if r < 0.05 or baseline:
            shuffled_T = random.sample(T, len(T))
        else:
            shuffled_T = weighted_shuffle(T, weights)

        tmp_best = solve(model, shuffled_T, C_map, best_model)

        if cost(T, tmp_best, C_map) < cost(T, best_model, C_map):
            log("betterment:", cost(T, best_model, C_map))
            betterments += 1
            best_model = tmp_best

    return best_model


if __name__ == "__main__":
    # Seed the random generator
    # TODO: In the future, test with different seed values and see whether it affects the results
    random.seed(42)
    np.random.seed(42)

    # Load arguments
    parser = argparse.ArgumentParser(
        description='PBO-#oracle solver', formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('instance', help='instance file name')
    parser.add_argument('--timeout', type=int, default=60,
                        help='Timeout after given seconds')
    parser.add_argument('--verbose', action='store_true',
                        help='Print debug info')
    parser.add_argument('--baseline', action='store_true',
                        help='Run without any semantics')
    args = parser.parse_args()

    print("c PBO-#ihs")
    for k, v in args.__dict__.items():
        if k in ["instance"]:
            continue
        print("c {}: {}".format(k, v))
    print("c ---------------------------------\n")

    in_file = args.instance
    verbose = args.verbose
    baseline = args.baseline
    timeout = args.timeout

    # Init timeout hanlder
    signal.signal(signal.SIGINT, TimeoutHandler())
    signal.signal(signal.SIGTERM, TimeoutHandler())
    signal.signal(signal.SIGXCPU, TimeoutHandler())
    signal.signal(signal.SIGALRM, TimeoutHandler())
    signal.alarm(timeout)

    # Init solver
    origMaxVar, objvars, objcoefs, constraints = load_input(in_file)
    rsat = roundingsat.Roundingsat()
    rsat.init(origMaxVar)

    for idx, c in enumerate(constraints):
        rsat.addConstraint(c[0], c[1])
        constraints[idx] = (c[0], c[1])

    # Init data structures
    T = objvars
    for i in range(len(objvars)):
        C_map[objvars[i]] = objcoefs[i]

    # Start solver
    solve_inc(rsat, objvars, C_map)
