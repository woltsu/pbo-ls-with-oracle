import sys
import signal
import roundingsat
import copy
import random
import numpy as np
from math import inf
from solver_util import weighted_shuffle, cost, is_good, is_sat, calculate_weights, split_to_good_and_bad, load_input, load_args, init_timer
from sanity_check import check

# Global variables
SEED = 42
best_model = []
T = []
C_map = {}
best_cost = inf
betterments = 0
verbose = False
baseline = False
print_solution = False
constraints = None
timer = None


# Outputs result after given timeout
class TimeoutHandler:
    SIGINT = False

    def handle_exit(self, signo, frame):
        self.SIGINT = True
        if (check(model=best_model, constraints=constraints)):
            timer(cost(T, best_model, C_map))
            print("r SAT")
            if print_solution:
                print(f"o {' '.join(map(str, best_model))}")
        else:
            print("r UNSAT")

        sys.exit()

    def init(self, timeout):
        signal.signal(signal.SIGINT, self.handle_exit)
        signal.signal(signal.SIGTERM, self.handle_exit)
        signal.signal(signal.SIGXCPU, self.handle_exit)
        signal.signal(signal.SIGALRM, self.handle_exit)
        signal.alarm(timeout)


def shuffle_literals(T, model, weights_per_l):
    global baseline
    r = np.random.random(1)[0]

    (good, bad) = split_to_good_and_bad(
        T, C_map, model)

    if baseline:
        return random.sample(T, len(T))

    if r < 0.2:
        return good + bad
    elif r < 0.4:
        return bad + good
    elif r < 0.6:
        return random.sample(T, len(T))
    elif r < 0.8:
        weights = np.array([weights_per_l[l] for l in T])
        return weighted_shuffle(T, weights / weights.sum())
    else:
        return T


def solve(model, T, C_map, best, weights_per_l):
    global best_cost

    best_copy = copy.copy(best)
    final_result = None
    assumptions = []
    T_len = len(T)

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

                if i != T_len - 1:
                    T = T[:i+1] + \
                        shuffle_literals(T[i+1:], best_copy, weights_per_l)

                assumptions = tmp_assumptions
                final_result = copy.copy(tmp_result)
            else:
                assumptions.append(best_copy[l_idx])

    if final_result is None or not is_sat(final_result):
        return best

    return best_copy


def solve_inc(model, T, C_map):
    global best_model
    global best_cost
    global betterments
    global verbose

    model.solve([], 0)
    result = model.getResult()
    weights_per_l = calculate_weights(T, C_map)

    if not is_sat(result):
        print("r UNSAT")
        exit(0)

    best_model = copy.copy(result[1])
    best_cost = cost(T, best_model, C_map)

    if verbose:
        timer(best_cost)

    while True and not TimeoutHandler.SIGINT:
        shuffled_T = shuffle_literals(T, best_model, weights_per_l)
        tmp_best = solve(model, shuffled_T, C_map, best_model, weights_per_l)
        tmp_cost = cost(T, tmp_best, C_map)

        if tmp_cost < best_cost:
            if verbose:
                timer(tmp_cost)
            betterments += 1
            best_model = tmp_best
            best_cost = tmp_cost

    return best_model


if __name__ == "__main__":
    # Seed the random generator
    # TODO: In the future, test with different seed values and see whether it affects the results
    random.seed(SEED)
    np.random.seed(SEED)

    # Load arguments
    in_file, verbose, baseline, timeout, print_solution = load_args()

    # Init timeout hanlder
    timeout_handler = TimeoutHandler()
    timeout_handler.init(timeout)

    # Init timer
    timer = init_timer()

    # Init solver
    origMaxVar, objvars, objcoefs, constraints = load_input(in_file)
    rsat = roundingsat.Roundingsat()
    rsat.init(origMaxVar)

    for idx, c in enumerate(constraints):
        rsat.addConstraint(c[0], c[1])

    # Init data structures
    T = objvars
    for i in range(len(objvars)):
        C_map[objvars[i]] = objcoefs[i]

    # Start solver
    solve_inc(rsat, objvars, C_map)
