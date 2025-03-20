from numpy import inf
from sklearn import metrics 
import sys
import symnmfmodule
import symnmf
import numpy as np

# ------------------------ relevant funcs from HW1 ------------------------

EPSILON = 1e-4

def euclidean_distance(point1, point2):
    return (sum([(point1[i] - point2[i])**2 for i in range(len(point1))]))**0.5

def init_centroids(K, datapoints):
    return datapoints[:K]

def run_kmeans(K, datapoints, centroids, d, iter):
    for i in range(iter):
        old_centroids = centroids
        centroids_sums = [[0 for i in range(d)] for j in range(K)]
        centroids_counters = [0 for i in range(K)]
        
        for point in datapoints:
            min_dist_centroid = 0
            min_dist = euclidean_distance(point, centroids[0])

            for ind, centroid in enumerate(centroids):
                dist = euclidean_distance(point, centroid)

                # if the point is closest to this centroid so far:
                if dist < min_dist:
                    # update the info of the closest centroid
                    min_dist = dist
                    min_dist_centroid = ind
            
            # add the point to the closest centroid's sum and up its counter by 1
            for coord in range(d):
                centroids_sums[min_dist_centroid][coord] += point[coord]
            centroids_counters[min_dist_centroid] += 1
        
        # set the new centroids to the mean of all of the points closest to them
        centroids = [[centroids_sums[j][i] / centroids_counters[j] for i in range(d)] for j in range(K)]

        # if the centroids barely moved since the last iteration, stop
        if max([euclidean_distance(centroids[i], old_centroids[i]) for i in range(K)]) < EPSILON:
            break
    
    return centroids

# ------------------------ comparing symnmf & kmeans from HW1 ------------------------

def main():
    k = int(sys.argv[1])
    file = sys.argv[2]
    DEFAULT_ITER = 300

    datapoints = symnmf.load_datapoints(file)
    N = len(datapoints)
    d = len(datapoints[0])

    # running kmeans from HW1
    centroids = run_kmeans(k, datapoints, init_centroids(k, datapoints), d, DEFAULT_ITER)
    kmeans_labeling = [0 for i in range(N)]
    for datapoint_idx in range(N):  # computing each datapoint's label
        min_dist = inf
        for centroid_idx in range(len(centroids)):
            curr_dist = euclidean_distance(datapoints[datapoint_idx], centroids[centroid_idx])
            if curr_dist < min_dist:
                min_dist = curr_dist
                min_idx = centroid_idx
        kmeans_labeling[datapoint_idx] = min_idx # label == index of closest centroid
    
    kmeans_silhouette_score = metrics.silhouette_score(datapoints, kmeans_labeling, metric='euclidean')

    # running symnmf
    W = symnmfmodule.norm(datapoints)
    H = symnmf.initialize_H(W, N, k)
    H = symnmf.symnmf(H, W)

    symnmf_clustering = np.argmax(H, axis=1)  # computing each datapoint's label
    symnmf_silhouette_score = metrics.silhouette_score(datapoints, symnmf_clustering, metric='euclidean')

    # printing comparison
    print(f"nmf: {symnmf_silhouette_score}")
    print(f"kmeans: {kmeans_silhouette_score}")

if __name__ == "__main__":
    try:
        main()
    except Exception:
        symnmf.print_error()