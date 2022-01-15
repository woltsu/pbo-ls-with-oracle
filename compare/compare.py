import glob
import argparse
import numpy as np
import matplotlib.pyplot as plt


def load_args():
    parser = argparse.ArgumentParser(
        description='compare-ls-to-oracle', formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--ls',
                        help='Path to ls results', required=True)
    parser.add_argument('--oracle',
                        help='Path to oracle results', required=True)
    parser.add_argument(
        "--print", help="Print oracle instances", action='store_true')
    args = parser.parse_args()

    path_to_ls = args.ls
    path_to_oracle = args.oracle
    print_oracle_instances = args.print

    return path_to_ls, path_to_oracle, print_oracle_instances


def find_log_files(path):
    files = []

    for log_file in glob.iglob(f"{path}/**/*.log", recursive=True):
        files.append(log_file)

    return files


def parse_instance_name(path):
    return path.split(".")[0].split("/")[-1]


def parse_ls_results(path_to_ls):
    results = {}
    ls_log_files = find_log_files(path_to_ls)

    for log_file in ls_log_files:
        instance_name = parse_instance_name(log_file)
        results[instance_name] = {}

        with open(log_file) as log:
            data = log.readlines()
            score = data[0].split("\t")[0]

            try:
                results[instance_name]["score"] = int(score)
            except:
                # No feasible solution
                results[instance_name]["score"] = float('inf')

    return results


def parse_oracle_results(path_to_oracle):
    results = {}
    oracle_log_files = find_log_files(path_to_oracle)

    for log_file in oracle_log_files:
        instance_name = parse_instance_name(log_file)
        results[instance_name] = {}

        with open(log_file) as log:
            data = log.readlines()
            for l in reversed(data):
                if l[0] == "t":
                    score = int(l.split(" - ")[-1])
                    results[instance_name]["score"] = int(score)
                    break

            # No feasible solution
            if not "score" in results[instance_name]:
                results[instance_name]["score"] = float('inf')

    return results


def compare_ls_and_oracle_results(ls_results, oracle_results):
    ls_points = 0
    oracle_points = 0
    equal_points = 0

    ls_scores = []
    oracle_scores = []
    oracle_instances = []

    for instance, ls_result in ls_results.items():
        oracle_result = oracle_results[instance]

        ls_score = ls_result["score"]
        oracle_score = oracle_result["score"]

        ls_scores.append(ls_score)
        oracle_scores.append(oracle_score)

        if ls_score < oracle_score:
            ls_points += 1
        elif oracle_score < ls_score:
            oracle_points += 1
            oracle_instances.append(instance)
        else:
            equal_points += 1

    plt.scatter(ls_scores, oracle_scores, alpha=0.5)
    # plt.show()

    return ls_points, oracle_points, equal_points, oracle_instances


if __name__ == "__main__":
    path_to_ls, path_to_oracle, print_oracle_instances = load_args()

    ls_results = parse_ls_results(path_to_ls)
    oracle_results = parse_oracle_results(path_to_oracle)

    ls_points, oracle_points, equal_points, oracle_instances = compare_ls_and_oracle_results(
        ls_results, oracle_results)

    print(
        f"ls points: {ls_points}, oracle points: {oracle_points}, equal points: {equal_points}")

    if print_oracle_instances:
        for instance in oracle_instances:
            print(instance)