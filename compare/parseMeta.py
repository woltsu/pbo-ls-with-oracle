

import glob
import os
import pandas as pd


def load_input(in_file):
    if not os.path.exists(in_file):
        print("cannot open path \""+in_file+"\"")
        exit(1)

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


def parse_instance_name(path):
    domain = path.split(".")[0].split("/")[-2]
    name = path.split(".")[0].split("/")[-1]
    return domain, name


def find_files(path):
    files = []

    for log_file in glob.iglob(f"{path}/**/*.opb", recursive=True):
        files.append(log_file)

    return files


if __name__ == "__main__":
    INSTANCES_PATH = "instances/"
    domains = []
    results = {}

    for f in find_files(INSTANCES_PATH):
        domain, name = parse_instance_name(f)
        if not domain in domains:
            domains.append(domain)

        if not domain in results:
            results[domain] = {
                "n_vars": 0,
                "n_cons": 0,
                "n_instances": 0
            }

        origMaxVar, objvars, objcoefs, constraints = load_input(f)
        n_vars = origMaxVar - 1
        n_cons = len(constraints)

        results[domain]["n_instances"] += 1
        results[domain]["n_vars"] += n_vars
        results[domain]["n_cons"] += n_cons

    meta = pd.DataFrame(
        columns=["domain", "average number of variables",
                 "average number of constraints"], index=domains
    )

    for domain in results.keys():
        n_instances = results[domain]["n_instances"]
        results[domain]["n_vars"] /= n_instances
        results[domain]["n_cons"] /= n_instances
        meta.loc[domain, "domain"] = domain
        meta.loc[domain, "average number of variables"] = results[domain]["n_vars"]
        meta.loc[domain, "average number of constraints"] = results[domain]["n_cons"]

    meta = meta.style.format(precision=0)
    meta = meta.hide_index()
    meta = meta.format("\\textbf{{{}}}", escape="latex",
                       subset=meta.columns[0:1])
    print(meta.to_latex())
