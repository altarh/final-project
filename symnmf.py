import sys
import pandas as pd
import numpy as np
from numpy import inf
 
import symnmfmodule

np.random.seed(1234)

def print_error():
    '''Print an error message in case of an exception.'''
    print("An Error Has Occured.")

def initialize_H(W, N, k):
    '''Create initial decomposition matrix H.'''
    m = np.mean(W)
    H = np.random.uniform(0, 2 * ((m/k)**0.5), size=(N, k))
    return H

def update_H(W, H):
    '''Perform a single iteration of updating decomposition matrix H.'''
    numerator = W@H
    denominator = H@H.T@H
    H = H*(0.5 + 0.5*numerator/denominator)  # notice that there isn't a check for dividing by 0 
    return H

def symnmf(H, W):
    '''
    Given initial decomposition matrix H and normalized similarity matrix W,
    Update H until max iteration or convergence.
    '''
    epsilon = 1e-4
    iteration = 0
    diff_sum = inf # initial squared norm to allow first iteration
    while (diff_sum > epsilon and iteration < 300):
        prevH = H
        H = update_H(W, H)
        diff_sum = np.sum((H - prevH)**2)
        iteration += 1
    return H

def load_datapoints(filename):
    '''Load datapoints from given text file into 2D ndarray.'''
    datapoints = pd.read_csv(filename, sep=",", header=None)
    return datapoints.values.tolist()

def print_mat(mat):
    '''Print given 2D array-like, comma-delimited.'''
    for row in mat:
        print(",".join(["%.4f" % val for val in row]))

def main():
    k = int(sys.argv[1])
    goal, filename = sys.argv[2:]

    datapoints = load_datapoints(filename)
    N = len(datapoints)

    if not 0 < k < N: # invalid k
        raise Exception()

    if goal == "sym":
        result = symnmfmodule.sym(datapoints)
    elif goal == "ddg":
        result = symnmfmodule.ddg(datapoints)
    elif goal == "norm":
        result = symnmfmodule.norm(datapoints)
    elif goal == "symnmf":
        W = symnmfmodule.norm(datapoints)
        H0 = initialize_H(W, N, k)
        result = symnmf(H0, W)
    else: # invalid goal
        raise Exception()
    
    print_mat(result)
    
if __name__ == "__main__":
    try:
        main()
    except Exception:
        print_error()