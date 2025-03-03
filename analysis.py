from numpy import inf
from sklearn import metrics 
import sys
import pandas as pd
import HW1
import symnmfmodule
import symnmf
import numpy as np

def main():
    k = int(sys.argv[1])
    file = sys.argv[2]
    EPSILON = 1e-4
    DEFAULT_ITER = 300

    datapoints = symnmf.load_datapoints(file)
    N = len(datapoints)
    d = len(datapoints[0])

    centroids = HW1.run_kmeans(k, datapoints, HW1.init_centroids(k, datapoints), d, DEFAULT_ITER)
    kmeans_labeling = [0 for i in range(N)]
    for datapoint_idx in range(N):  # computing each datapoint's label
        min_dist = inf
        for centroid_idx in range(len(centroids)):
            curr_dist = HW1.euclidean_distance(datapoints[datapoint_idx], centroids[centroid_idx])
            if curr_dist < min_dist:
                min_dist = curr_dist
                min_idx = centroid_idx
        kmeans_labeling[datapoint_idx] = min_idx
    
    kmeans_silhouette_score = metrics.silhouette_score(datapoints, kmeans_labeling, metric='euclidean')

    W = symnmfmodule.norm(datapoints)
    H = symnmf.initialize_H(W, N, k)
    H = symnmf.symnmf(H, W)

    symnmf_clustering = np.argmax(H, axis=1)  # computing each datapoint's label
    symnmf_silhouette_score = metrics.silhouette_score(datapoints, symnmf_clustering, metric='euclidean')

    print(f"nmf: {symnmf_silhouette_score}")
    print(f"kmeans: {kmeans_silhouette_score}")

if __name__ == "__main__":
    try:
        main()
    except Exception:
        symnmf.print_error()