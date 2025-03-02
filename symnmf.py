import numpy as np
np.random.seed(1234)

def initialize_H(W, n, k):
    m = np.mean(W)
    H = np.random.uniform(0, 2 * ((m/k)**0.5), size=(n, k))

def update_H(W, H, n, k):
    halfs_mat = np.full((n, k), 0.5)
    numerator = W@H
    denominator = H@H.T@H
    H = H*(0.5 + 0.5*numerator/denominator)  # notice that there isn't a check for dividing by 0 