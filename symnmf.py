import sys
import pandas as pd
import numpy as np
np.random.seed(1234)

def print_error():
    print("An Error Has Occured.")

def initialize_H(W, n, k):
    m = np.mean(W)
    H = np.random.uniform(0, 2 * ((m/k)**0.5), size=(n, k))
    return H

def update_H(W, H):
    numerator = W@H
    denominator = H@H.T@H
    H = H*(0.5 + 0.5*numerator/denominator)  # notice that there isn't a check for dividing by 0 
    return H

def symnmf(H, W):
    epsilon = 0.0001
    prevH = H
    H = update_H(W, H)
    diff_sum = np.sum((H - prevH)**2)
    iteration  = 0
    while (diff_sum > epsilon and iter < 300):
        iteration  += 1
        prevH = H
        H = update_H(W, H)
        diff_sum = np.sum((H - prevH)**2)
    return H

def load_datapoints(filename):
    datapoints = pd.read_csv(filename, sep = ",", header=None)
    lst = datapoints.values.tolist()
    return lst

def print_mat(mat):
    for row in mat:
        print(",".join(row))

def main():
    try:
        k = int(sys.argv[1])
        goal, filename = sys.argv[2:]

        datapoints = load_datapoints(filename)
        N = len(datapoints)
        d = len(datapoints[0])

        if not 0 < k < N:
            raise Exception()

        if goal == "sym":
            # result = sym(datapoints, N, d)
            pass
        elif goal == "ddg":
            # result = ddg(datapoints, N, d)
            pass
        elif goal == "norm":
            # result = norm(datapoints, N, d)
            pass
        elif goal == "symnmf":
            # W = norm(datapoints, N, d)
            # H0 = initialize_H(W, N, k)
            # result = symnmf(H0, W, N, k)
            pass
        else: # invalid goal
            raise Exception()
        # print_mat(result)
    
    except Exception:
        print_error()
        return 0

if __name__ == "__main__":
    main()