import numpy as np
np.random.seed(1234)

def initialize_H(W, n, k):
    m = np.mean(W)
    H = np.random.uniform(0, 2 * ((m/k)**0.5), size=(n, k))
    return H

def update_H(W, H, n, k):
    halfs_mat = np.full((n, k), 0.5)
    numerator = W@H
    denominator = H@H.T@H
    H = H*(0.5 + 0.5*numerator/denominator)  # notice that there isn't a check for dividing by 0 
    return H

def symnmf(H, W, n, k):
    epsilon = 0.0001
    prevH = H
    H = update_H(W, H, n, k)
    diff_sum = np.sum((H - prevH)**2)
    iteration  = 0
    while (diff_sum > epsilon and iter < 300):
        iteration  += 1
        prevH = H
        H = update_H(W, H, n, k)
        diff_sum = np.sum((H - prevH)**2)
    return H