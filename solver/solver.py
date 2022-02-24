import sys
import signal
import roundingsat
import copy
import random
import numpy as np
from time import time_ns
from math import inf
from solver_util import weighted_shuffle, cost, is_good, is_sat, calculate_weights, split_to_good_and_bad, load_input, load_args, init_timer
from sanity_check import check

# Global variables
best_model = []
T = []
C_map = {}
best_cost = inf
betterments = 1
verbose = False
baseline = False
print_solution = False
constraints = None
timer = None
current_solver_calls = 1
total_solver_calls = 1
total_time = 0


# Outputs result after given timeout
class TimeoutHandler:
    SIGINT = False

    def handle_exit(self, signo, frame):
        self.SIGINT = True
        if (check(model=best_model, constraints=constraints)):
            timer(cost(T, best_model, C_map),
                  current_solver_calls, total_solver_calls, betterments)
            print(f"T {total_time / total_solver_calls}")
            print("r SAT")
            if print_solution:
                print(
                    f"o {' '.join(map(str, best_model))}")
        else:
            print("r UNSAT")

        sys.exit()

    def init(self, timeout):
        signal.signal(signal.SIGINT, self.handle_exit)
        signal.signal(signal.SIGTERM, self.handle_exit)
        signal.signal(signal.SIGXCPU, self.handle_exit)
        signal.signal(signal.SIGALRM, self.handle_exit)
        signal.alarm(timeout)


def shuffle_literals(T, model, weights):
    r = np.random.random(1)[0]
    if r < 0.05:
        return random.sample(T, len(T))
    else:
        return weighted_shuffle(T, weights)


def to_good_and_bad(T, C_map, model):
    (good, bad) = split_to_good_and_bad(
        T, C_map, model)
    return good + bad


def model_solve(model, *args):
    global total_time
    start = time_ns() // 1000000
    model.solve(*args)
    total_time += (time_ns() // 1000000) - start


def solve(model, T, C_map, best):
    global best_cost
    global current_solver_calls
    global current_solver_calls
    global total_solver_calls

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
            model_solve(model, tmp_assumptions, 0)
            current_solver_calls += 1
            total_solver_calls += 1

            tmp_result = model.getResult()

            if is_sat(tmp_result):
                best_copy = tmp_result[1]

                if i != T_len - 1:
                    T = T[:i+1] + \
                        to_good_and_bad(T[i+1:], C_map, best_copy)

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
    global current_solver_calls
    global total_solver_calls

    model_solve(model, [], 0)
    result = model.getResult()
    weights_per_l = calculate_weights(T, C_map)
    weights = np.array([weights_per_l[l] for l in T])
    weights = weights / weights.sum()

    if not is_sat(result):
        print("r UNSAT")
        exit(0)

    best_model = copy.copy(result[1])
    best_cost = cost(T, best_model, C_map)

    if verbose:
        timer(best_cost, current_solver_calls, total_solver_calls, betterments)

    while True and not TimeoutHandler.SIGINT:
        shuffled_T = shuffle_literals(T, best_model, weights)

        tmp_best = solve(model, shuffled_T, C_map, best_model)
        tmp_cost = cost(T, tmp_best, C_map)

        if tmp_cost < best_cost:
            betterments += 1
            if verbose:
                timer(tmp_cost, current_solver_calls,
                      total_solver_calls, betterments)
            current_solver_calls = 0
            best_model = tmp_best
            best_cost = tmp_cost

    return best_model


if __name__ == "__main__":
    # Load arguments
    in_file, verbose, baseline, timeout, print_solution, seed = load_args()

    # Seed the random generator
    random.seed(seed)
    np.random.seed(seed)

    # Init timeout hanlder
    timeout_handler = TimeoutHandler()
    timeout_handler.init(timeout)

    # Init timer
    timer = init_timer()

    # Init solver
    origMaxVar, objvars, objcoefs, constraints = load_input(in_file)
    rsat = roundingsat.Roundingsat()
    rsat.init(origMaxVar)

    M = 0
    for c in objcoefs:
        if c < 0:
            M += c

    print(f"M {M}")
    print(f"S {origMaxVar - 1} {len(constraints)}")

    for idx, c in enumerate(constraints):
        rsat.addConstraint(c[0], c[1])

    # Init data structures
    T = objvars
    for i in range(len(objvars)):
        C_map[objvars[i]] = objcoefs[i]

    # Start solver
    solve_inc(rsat, objvars, C_map)
