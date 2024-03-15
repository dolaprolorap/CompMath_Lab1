from sympy import symbols, diff, maximum, minimum, cos, Interval

# Input values
x, eps = symbols("x eps")
func = 2*x - cos(x)
#func = x * 0 + 1
borders = [0.1, 0.6]
arr_x = [0.44, 0.13, 0.58, 0.37]

# 1 point
h = (borders[1] - borders[0]) / 10
net = [borders[0] + h * i for i in range(11)]
func_in_net = [[val, func.subs(x, val)] for val in net]

print("Function " + str(func) + " in net:")
for val in func_in_net:
    print("{0:0.2f}: {1}".format(val[0], val[1]))

def exec(f, x_):
    def L_1(f, x_1, x_2):
        return f.subs(x, x_1) * (x - x_2) / (x_1 - x_2) + f.subs(x, x_2) * (x - x_1) / (x_2 - x_1)

    def L_2(f, x_1, x_2, x_3):
        return f.subs(x, x_1) * (x - x_2) * (x - x_3) / ((x_1 - x_2) * (x_1 - x_3)) + \
            f.subs(x, x_2) * (x - x_1) * (x - x_3) / ((x_2 - x_1) * (x_2 - x_3)) + \
            f.subs(x, x_3) * (x - x_1) * (x - x_2) / ((x_3 - x_1) * (x_3 - x_2))

    def find_interval(net, x):
        cnt = len(net)
        for i in range(cnt - 1):
            if (net[i] <= x) and (net[i + 1] >= x):
                return (net[i], net[i + 1])
            
    def find_double_interval(net, x):
        prev_dist = -1
        cnt = len(net)
        for i in range(cnt - 2):
            if (net[i] <= x) and (net[i + 2] >= x):
                if(prev_dist < 0):
                    prev_dist = abs(net[i] - x) + abs(net[i + 2] - x)
                else:
                    cur_dist = abs(net[i] - x) + abs(net[i + 2] - x)
                    if(cur_dist < prev_dist):
                        return (net[i], net[i + 1], net[i + 2])
                    else:
                        return (net[i - 1], net[i], net[i + 1])
        return (net[cnt - 3], net[cnt - 2], net[cnt - 1])

    def divided_diff(f, *args):
        s = 0
        for i in args:
            for j in args:
                mul = 1
                if args.index(i) != args.index(j):
                    mul *= (j - i)
                s += f.subs(x, i) / mul
        return s

    def N_1(f, x_1, x_2):
        return f.subs(x, x_1) + divided_diff(f, x_1, x_2) * (x - x_1)

    def N_2(f, x_1, x_2, x_3):
        return N_1(f, x_1, x_2) + divided_diff(f, x_1, x_2, x_3) * (x - x_1) * (x - x_2)
    
    x_1, x_2 = find_interval(net, x_)

    print("\nFunction:")
    print("{0}: {1}".format(x_, f.subs(x, x_)))

    print("\nFirst order Lagrange: ")
    print("{0}: {1}".format(x_, L_1(f, x_1, x_2).subs(x, x_)))

    sec_der = diff(f, x, 2).subs(x, eps)
    max_sec_der = maximum(sec_der, eps, Interval(x_1, x_2))
    min_sec_der = minimum(sec_der, eps, Interval(x_1, x_2))

    print("\nSecond derivative:")
    print("Maximum: {0}".format(max_sec_der))
    print("Minimum: {0}".format(min_sec_der))

    w_2 = (x - x_1) * (x - x_2)
    max_w_2 = maximum(w_2, x, Interval(x_1, x_2))
    min_w_2 = minimum(w_2, x, Interval(x_1, x_2))
    max_R1 = max_sec_der * max_w_2 / 2
    min_R1 = min_sec_der * min_w_2 / 2
    print("\nR_1(x):")
    print("Maximum: {0}".format(max_R1))
    print("Minimum: {0}".format(min_R1))

    print("\nR_1(x*):")
    print("{0}: {1}".format(x_, f.subs(x, x_) - L_1(f, x_1, x_2).subs(x, x_)))

    x_1, x_2, x_3 = find_double_interval(net, x_)

    print("\nSecond order Lagrange:")
    print("{0}: {1}".format(x_, L_2(f, x_1, x_2, x_3).subs(x, x_)))

    third_der = diff(func, x, 3).subs(x, eps)
    max_third_der = maximum(third_der, eps, Interval(x_1, x_3))
    min_third_der = minimum(third_der, eps, Interval(x_1, x_3))

    print("\nThird derivative:")
    print("Maximum: {0}".format(max_third_der))
    print("Minimum: {0}".format(min_third_der))

    w_3 = (x - x_1) * (x - x_2) * (x - x_3)
    max_w_3 = maximum(w_3, x, Interval(x_1, x_3))
    min_w_3 = minimum(w_3, x, Interval(x_1, x_3))
    max_R2 = max_third_der * max_w_3 / 6
    min_R2 = min_third_der * min_w_3 / 6
    print("\nR_2(x):")
    print("Maximum: {0}".format(max_R2))
    print("Minimum: {0}".format(min_R2))

    print("\nR_2(x*):")
    print("{0}: {1}".format(x_, f.subs(x, x_) - L_2(f, x_1, x_2, x_3).subs(x, x_)))

    print("\nDivided diffs:")
    print("{0}: {1}".format("x_1 x_2", divided_diff(f, x_1, x_2)))
    print("{0}: {1}".format("x_1 x_3", divided_diff(f, x_1, x_3)))
    print("{0}: {1}".format("x_2 x_3", divided_diff(f, x_2, x_3)))
    print("{0}: {1}".format("x_1 x_2 x_3", divided_diff(f, x_1, x_2, x_3)))

    print("\nNewton:")
    print("N_1(x) | {0}: {1}".format(x_, N_1(f, x_1, x_2).subs(x, x_)))
    print("N_2(x) | {0}: {1}".format(x_, N_2(f, x_1, x_2, x_3).subs(x, x_)))

exec(func, arr_x[0])