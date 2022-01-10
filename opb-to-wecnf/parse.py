import os
import argparse


def load_input(in_file):
    if not os.path.exists(in_file):
        print("c cannot open path \""+in_file+"\"")
        exit(1)

    print("c parsing file: {}".format(in_file))

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

    return origMaxVar - 1, objvars, objcoefs, constraints


def load_args():
    parser = argparse.ArgumentParser(
        description='opb-to-wecnf', formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('instance', help='instance file name')
    args = parser.parse_args()

    print("c opb-to-wecnf")
    for k, v in args.__dict__.items():
        if k in ["instance"]:
            continue
        print("c {}: {}".format(k, v))
    print("c ---------------------------------")

    in_file = args.instance

    return in_file


def parse_instance(nbvar, objvars, objcoefs, constraints):
    top = 0
    soft = []

    for i, c in enumerate(objcoefs):
        var = objvars[i]
        weight = abs(c)
        top += weight

        degree = 1
        if c > 0:
            var = f"-{var}"
            degree = 0

        soft.append(f"{weight} {degree} 1 {var} 0")

    print("c ----------- PARAMETERS -----------")
    print(f"p wcnf {nbvar} {len(soft) + len(constraints)} {top}")

    print("c ----------- HARD CONSTRAINTS -----------")
    for coefs, degree in constraints:
        """ v = " ".join(f"{abs(c)} {var if c > 0 else -var}" for var,
                     c in coefs.items()) """
        v = " ".join(f"{c} {var}" for var,
                    c in coefs.items())
        print(f"{top} {degree} {v} 0")

    print("c ----------- SOFT CONSTRAINTS -----------")
    for c in soft:
        print(c)


if __name__ == "__main__":
    in_file = load_args()
    nbvar, objvars, objcoefs, constraints = load_input(in_file)
    parse_instance(nbvar, objvars, objcoefs, constraints)
