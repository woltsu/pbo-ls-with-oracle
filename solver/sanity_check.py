def check(model, constraints):
    is_sat = True

    for (coefs, rhs) in constraints:
        c_sum = 0
        for var in model:
            if var < 0:
                continue

            l = abs(var)
            if l in coefs:
                c_sum += coefs[l]

        if c_sum < rhs:
            is_sat = False

    return is_sat
