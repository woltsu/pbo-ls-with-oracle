import glob
import argparse
import matplotlib.pyplot as plt
import matplotlib.patches as mpl_patches
import numpy as np
import pandas as pd
import os


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
    domain = path.split(".")[0].split("/")[-2]
    name = path.split(".")[0].split("/")[-1]
    return domain, name


def parse_ls_results(path_to_ls):
    results = {}
    ls_log_files = find_log_files(path_to_ls)

    for log_file in ls_log_files:
        instance_domain, instance_name = parse_instance_name(log_file)
        results[instance_name] = {}
        results[instance_name]["domain"] = instance_domain
        results[instance_name]["solved"] = True

        with open(log_file) as log:
            try:
                data = log.readlines()
                score = data[0].split("\t")[0]
                results[instance_name]["score"] = int(score)
            except:
                # No feasible solution
                results[instance_name]["score"] = float('inf')
                results[instance_name]["solved"] = False

    return results


def parse_oracle_results(path_to_oracle):
    results = {}
    oracle_log_files = find_log_files(path_to_oracle)

    for log_file in oracle_log_files:
        instance_domain, instance_name = parse_instance_name(log_file)
        results[instance_name] = {}
        results[instance_name]["domain"] = instance_domain
        results[instance_name]["solved"] = True
        results[instance_name]["scores"] = []

        with open(log_file) as log:
            data = log.readlines()
            if "r SAT\n" in data:
                for l in reversed(data):
                    if l[0] == "t":
                        splitted = l.split(" - ")[-1].split(" ")
                        score = int(splitted[0])
                        n_betterments = int(splitted[-1])
                        n_calls = int(splitted[-2])
                        results[instance_name]["score"] = score
                        results[instance_name]["n_betterments"] = n_betterments
                        results[instance_name]["n_calls"] = n_calls
                        break

                for l in data:
                    if l[0] == "M":
                        results[instance_name]["M"] = int(l.split(" ")[-1])
                    elif l[0] == "T":
                        results[instance_name]["T"] = float(l.split(" ")[-1])
                    elif l[0] == "S":
                        a = int(l.split(" ")[-1])
                        b = int(l.split(" ")[-2])
                        results[instance_name]["n_vars"] = b
                        results[instance_name]["n_cons"] = a
                    elif l[0] == "t":
                        splitted = l.split(" - ")
                        time = int(
                            float(splitted[0].split(" ")[1].replace("s", "")))
                        score = int(splitted[1].split(" ")[0])
                        results[instance_name]["scores"].append((time, score))

            # No feasible solution
            if not "score" in results[instance_name]:
                results[instance_name]["score"] = float('inf')
                results[instance_name]["solved"] = False

    return results


def compare_ls_and_oracle_results(ls_results, oracle_results):
    ls_points = 0
    oracle_points = 0
    equal_points = 0
    missing = 0

    ls_scores = []
    oracle_scores = []
    oracle_instances = []
    scatter_x = []
    scatter_y = []

    domains = []
    # TODO: #solver calls made
    # TODO: %solver calls with betterments
    results = {}

    for instance, ls_result in ls_results.items():
        if instance in oracle_results:
            oracle_result = oracle_results[instance]
        else:
            missing += 1
            continue

        # TODO: Handle?
        if not oracle_result["solved"] and not ls_result["solved"]:
            continue

        domain = oracle_result["domain"]
        if not domain in domains:
            domains.append(domain)
            results[domain] = {"wins": 0, "equal": 0,
                               "instances": 0, "solved": 0, "score": 0, "ls_score": 0, "solver_time": 0, "solver_times": [], "solver_calls": 0, "solver_betterments": 0,
                               "n_vars": 0, "n_cons": 0, "solved_instances": 0}

        results[domain]["instances"] += 1
        if not oracle_result["solved"]:
            ls_points += 1
            results[domain]["score"] += 0
            continue

        results[domain]["n_vars"] += oracle_result["n_vars"]
        results[domain]["n_cons"] += oracle_result["n_cons"]

        results[domain]["solved"] += 1
        T = oracle_result["T"]
        results[domain]["solver_time"] += T
        results[domain]["solver_times"].append(T)

        results[domain]["solver_calls"] += oracle_result["n_calls"]
        results[domain]["solver_betterments"] += oracle_result["n_betterments"]

        M = oracle_result["M"]
        ls_score = ls_result["score"] + M
        oracle_score = oracle_result["score"]

        ls_scores.append(ls_score)
        oracle_scores.append(oracle_score)

        if ls_score < oracle_score:
            ls_points += 1
        elif oracle_score < ls_score:
            oracle_points += 1
            oracle_instances.append(instance)
            results[domain]["wins"] += 1
        else:
            results[domain]["equal"] += 1
            equal_points += 1

        s_min = min(ls_score, oracle_score)
        score = (abs(s_min - M) + 1) / (abs(oracle_score - M) + 1)

        if not ls_result["solved"]:
            score_ls = 0
        else:
            score_ls = (abs(s_min - M) + 1) / (abs(ls_score - M) + 1)

        results[domain]["score"] += score
        results[domain]["ls_score"] += score_ls

    df = pd.DataFrame(
        columns=["domain", "#instances", "#wins", "#equal", "#calls", "#betterments", "score", "score (ls)", "solve time (ms)", "median solve time (ms)"], index=domains)

    meta = pd.DataFrame(
        columns=["domain", "average number of variables",
                 "average number of constraints"], index=domains
    )

    for domain in results.keys():
        n_solved_instances = results[domain]["solved"] or 1
        n_instances = results[domain]["instances"]
        results[domain]["n_vars"] /= n_solved_instances
        results[domain]["n_cons"] /= n_solved_instances

        meta.loc[domain, "domain"] = domain
        meta.loc[domain, "average number of variables"] = results[domain]["n_vars"]
        meta.loc[domain, "average number of constraints"] = results[domain]["n_cons"]

        df.loc[domain, "domain"] = domain
        df.loc[domain, "#wins"] = results[domain]["wins"]
        df.loc[domain, "#equal"] = results[domain]["equal"]
        df.loc[domain, "#instances"] = results[domain]["instances"]
        df.loc[domain, "score"] = results[domain]["score"] / n_instances
        df.loc[domain,
               "score (ls)"] = results[domain]["ls_score"] / n_instances
        df.loc[domain, "#calls"] = results[domain]["solver_calls"] / \
            n_solved_instances
        df.loc[domain, "#betterments"] = results[domain]["solver_betterments"] / \
            n_solved_instances
        df.loc[domain,
               "solve time (ms)"] = results[domain]["solver_time"] / n_solved_instances
        df.loc[domain, "median solve time (ms)"] = np.median(
            np.array(results[domain]["solver_times"]))

        scatter_x.append(df.loc[domain, "score"])
        scatter_y.append(df.loc[domain, "score (ls)"])

    df["score"] = df["score"].astype("float").round(2)
    df["score (ls)"] = df["score (ls)"].astype("float").round(2)
    df["#calls"] = df["#calls"].astype("int64")
    df["#betterments"] = df["#betterments"].astype("float").round(2)
    df["solve time (ms)"] = df["solve time (ms)"].astype("float").round(2)
    df["median solve time (ms)"] = df["median solve time (ms)"].astype(
        "float").round(2)

    highlighted = df[abs(df["score"] - df["score (ls)"]) > 0.1]
    highlighted = highlighted.style.highlight_max(
        axis=1, props='bfseries: ;', subset=highlighted.columns[6:8])
    highlighted = highlighted.format(precision=2)
    highlighted = highlighted.hide_index()
    highlighted = highlighted.format(
        "\\textbf{{{}}}", escape="latex", subset=highlighted.columns[0:1])
    print(highlighted.to_latex())

    df = df.style.highlight_max(
        axis=1, props='bfseries: ;', subset=highlighted.columns[6:8])
    df = df.format(precision=2)
    df = df.hide_index()
    df = df.format("\\textbf{{{}}}", escape="latex",
                   subset=highlighted.columns[0:1])
    print(df.to_latex())

    meta = meta.style.format(precision=0)
    meta = meta.hide_index()
    meta = meta.format("\\textbf{{{}}}", escape="latex",
                       subset=highlighted.columns[0:1])
    print(meta.to_latex())

    graph_data = {
        "scatter_x": scatter_x,
        "scatter_y": scatter_y
    }

    return ls_points, oracle_points, equal_points, oracle_instances, missing, graph_data


def scatter(graph_data):
    plt.scatter(graph_data["scatter_x"], graph_data["scatter_y"], )
    plt.xlabel("LS-ORACLE-PBO score")
    plt.ylabel("LS-PBO score")
    plt.show()


def lines(oracle_results):
    i = 0
    for instance, data in oracle_results.items():
        i += 1
        scores = data["scores"]
        domain = data["domain"]

        results_dir = f"plots/{domain}"
        if not os.path.isdir(results_dir):
            os.makedirs(results_dir)

        x_vals = []
        y_vals = []
        for (x, y) in scores:
            x_vals.append(x)
            y_vals.append(y)
        plt.plot(x_vals, y_vals)
        plt.xlabel("Time (s)")
        plt.ylabel("Solution cost")
        plt.title(instance)

        if "n_vars" in data and "n_cons" in data:
            # create a list with two empty handles (or more if needed)
            handles = [mpl_patches.Rectangle((0, 0), 1, 1, fc="white", ec="white",
                                            lw=0, alpha=0)] * 2

            # create the corresponding number of labels (= the text you want to display)
            labels = []
            labels.append(f"No. variables = {data['n_vars']}")
            labels.append(f"No. constraints = {data['n_cons']}")
            plt.legend(handles, labels, loc='best', fontsize='small',
                    fancybox=True, framealpha=0.7,
                    handlelength=0, handletextpad=0)

        plt.savefig(f"{results_dir}/{instance}")
        plt.close()
    exit(0)


if __name__ == "__main__":
    path_to_ls, path_to_oracle, print_oracle_instances = load_args()

    oracle_results = parse_oracle_results(path_to_oracle)
    ls_results = parse_ls_results(path_to_ls)

    lines(oracle_results)

    ls_points, oracle_points, equal_points, oracle_instances, missing, graph_data = compare_ls_and_oracle_results(
        ls_results, oracle_results)

    print(f"missing: {missing}")

    print(
        f"ls points: {ls_points}, oracle points: {oracle_points}, equal points: {equal_points}")

    if print_oracle_instances:
        for instance in oracle_instances:
            print(instance)

    # scatter(graph_data)
