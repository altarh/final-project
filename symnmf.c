#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

const int SUCCESS = 0;
const int ERROR = 1;

#define GOTO_CLEANUP_IF_NULL(x) { if ((x) == NULL) { goto cleanup; } }
#define GOTO_CLEANUP_IF_ERROR(x) { if ((x) == ERROR) { goto cleanup; } }

const char *GENERIC_ERROR_MSG = "An Error Has Occurred\n";

typedef struct coord {
    double coord;
    struct coord *next;
} coord;

typedef struct centroid {
    coord *centroid_coords;
    coord *sum;
    int count;
} centroid;

typedef struct datapoint {
    coord *coords;
    struct datapoint *next;
} datapoint;

/* internal funcs for read_args - initializes memory for datapoints and coordinates */
int init_datapoint(datapoint **dp, coord *first_coord) {
    datapoint *new_datapoint = malloc(sizeof(datapoint));
    if (new_datapoint == NULL) {
        return ERROR;
    }
    new_datapoint->coords = first_coord;
    new_datapoint->next = NULL;

    *dp = new_datapoint;
    return SUCCESS;    
}

int init_coord(coord **c, double n) {
    coord *new_coord = malloc(sizeof(coord));
    if (new_coord == NULL) {
        return ERROR;
    }
    new_coord->coord = n;
    new_coord->next = NULL;

    *c = new_coord;
    return SUCCESS;
}

/**
 * Reads a file and parses it into a linked list of datapoints.
 *
 * @param filename The name of the file to read.
 * @param d The calculated number of coordinates in each datapoint.
 * @param N The calculated number of datapoints in the linked list.
 * @param datapoints The resulting linked list of datapoints.
 *
 * @return SUCCESS if the file was parsed successfully, ERROR otherwise.
 */
int parse_file(const char *filename, int *d, int *N, datapoint **datapoints) {
    int return_code = ERROR;
    double n; /* for the double values */
    char delim; /* commas, \n, ... */
    datapoint **curr_datapoint = datapoints;
    coord *first_coord = NULL;
    coord **curr_coord = NULL;

    FILE *file = fopen(filename, "r");
    GOTO_CLEANUP_IF_NULL(file);

    /* go over first line to get d */
    first_coord = NULL;
    curr_coord = &first_coord;
    do {
        fscanf(file, "%lf%c", &n, &delim);
        GOTO_CLEANUP_IF_ERROR(init_coord(curr_coord, n));
        curr_coord = &(*curr_coord)->next;
        (*d)++;
    } while (delim != '\n');

    /* initialize the first datapoint */
    GOTO_CLEANUP_IF_ERROR(init_datapoint(curr_datapoint, first_coord));
    curr_datapoint = &(*curr_datapoint)->next;
    (*N)++;

    /* go over the rest of the lines to get datapoints and N */
    first_coord = NULL;
    curr_coord = &first_coord;
    while (fscanf(file, "%lf%c", &n, &delim) == 2) {
        GOTO_CLEANUP_IF_ERROR(init_coord(curr_coord, n));
        curr_coord = &(*curr_coord)->next;

        if (delim == '\n') { /* if at the end of the line */
            GOTO_CLEANUP_IF_ERROR(init_datapoint(curr_datapoint, first_coord));
            curr_datapoint = &(*curr_datapoint)->next;
            first_coord = NULL;
            curr_coord = &first_coord;
            (*N)++;
        }
    }

    return_code = SUCCESS;

cleanup:
    if (file != NULL) {
        fclose(file);
    }
    return return_code;
}

void free_coords(coord *c) {
    coord *curr_coord = c;
    while (curr_coord != NULL) {
        coord *next_coord = curr_coord->next;
        free(curr_coord);
        curr_coord = next_coord;
    }
}

void free_datapoints_structs(datapoint *datapoints) {
    datapoint *curr_datapoint;
    datapoint *next_datapoint;

    curr_datapoint = datapoints;
    while (curr_datapoint != NULL) {
        next_datapoint = curr_datapoint->next;
        free_coords(curr_datapoint->coords);
        free(curr_datapoint);
        curr_datapoint = next_datapoint;
    }
}

void free_2D_array(double **array) {
    if (array != NULL) {
        free(array[0]);
        free(array);
        array = NULL;
    }
}

/**
 * Converts a linked list of datapoints into a 2D array.
 *
 * @param filename The name of the file to read.
 * @param d The calculated number of coordinates in each datapoint.
 * @param N The calculated number of datapoints in the linked list.
 * @param result The resulting 2D array of datapoint coordinates.
 *
 * @return SUCCESS if the conversion was successful, ERROR otherwise.
 */
int linked_list_to_2D_array(datapoint* point, int N, int d, double ***result){
    int i;
    int j;
    int return_code = ERROR;
    double *p;
    double **a;
    coord* coord1;

    p = calloc(d * N, sizeof(double));
    GOTO_CLEANUP_IF_NULL(p);
    a = calloc(N, sizeof(double *));
    GOTO_CLEANUP_IF_NULL(a);

    for (i = 0; i < N; i++) {
        coord1 = point->coords;
        a[i] = p + (i * d);
        for (j = 0; j < d; j++) {
            a[i][j] = coord1->coord;
            coord1 = coord1->next;
        }
        point = point->next;
    }

    return_code = SUCCESS;
    *result = a;

cleanup:
    if (return_code == ERROR) {
        /* try to free memory */
        free_2D_array(a);
    }
    return return_code;
}

/**
 * Reads a file and parses its contents into a 2D array of datapoint coordinates.
 *
 * @param filename The name of the file to read.
 * @param d A pointer to the number of coordinates in each datapoint.
 * @param N A pointer to the number of datapoints in the file.
 * @param datapoints_array A pointer to the resulting 2D array of datapoint coordinates.
 *
 * @return SUCCESS if the parsing was successful, ERROR otherwise.
 */
int parse(const char *filename, int *d, int *N, double ***datapoints_array) {
    int return_code = ERROR;
    datapoint *datapoints = NULL;

    GOTO_CLEANUP_IF_ERROR(parse_file(filename, d, N, &datapoints));
    GOTO_CLEANUP_IF_ERROR(linked_list_to_2D_array(datapoints, *N, *d, datapoints_array));
    return_code = SUCCESS;

cleanup:
    free_datapoints_structs(datapoints);

    return return_code;
}

/* Calculates the Euclidean distance between two points in d-dimensional space.  */
double calc_euclidean_distance(double *coord1, double *coord2, int d){
    double sum = 0;
    int i;

    for (i = 0; i < d; i++){
        sum += pow((coord1[i] - coord2[i]), 2);
    }
    return sqrt(sum);
}

/**
 * Calculates the similarity matrix A.
 *
 * @param datapoint_coords A 2D array of the datapoint coordinates.
 * @param N The number of datapoints in the array.
 * @param d The number of coordinates in each datapoint.
 * @param result The resulting similarity matrix.
 *
 * @return SUCCESS if the calculation was successful, ERROR otherwise.
 */
int sym(double** datapoint_coords, int N, int d, double ***result) {
    int i;
    int j;
    int return_code = ERROR;
    double* arr = malloc(N*N*sizeof(double));
    double** A = malloc(N*sizeof(double*));

    GOTO_CLEANUP_IF_NULL(arr);
    GOTO_CLEANUP_IF_NULL(A);
    
    for (i=0; i<N; i++) {
        A[i] = arr + (i*N);
        for (j=0; j<N; j++) {
            if (i == j) {
                A[i][j] = 0;
            }
            else {
                double dist = calc_euclidean_distance(datapoint_coords[i], datapoint_coords[j], d);
                A[i][j] = exp(-dist/2);
            }
        }
    }

    return_code = SUCCESS;
    *result = A;

cleanup:
    if (return_code == ERROR) {
        /* try to free memory */
        free_2D_array(A);
    }
    return return_code;
}

/**
 * Internal function which calculates the diagonal degree matrix D using the similarity matrix A.
 *
 * @param A The similarity matrix - a 2D array of size N x N.
 * @param N The number of rows and columns in A.
 * @param result The resulting diagonal degree matrix D.
 *
 * @return SUCCESS if the calculation was successful, ERROR otherwise.
 */
int _ddg(double** A, int N, double ***result) {
    int i;
    int j;
    int return_code = ERROR;
    double d;
    double* arr = calloc(N*N, sizeof(double));
    double** D = malloc(N*sizeof(double*));

    GOTO_CLEANUP_IF_NULL(arr);
    GOTO_CLEANUP_IF_NULL(D);

    for (i=0; i<N; i++) {
        D[i] = arr + i*N;
        d = 0;
        for (j=0; j<N; j++) {
            d += A[i][j];
        }
        D[i][i] = d;
    }

    return_code = SUCCESS;
    *result = D;

cleanup:
    if (return_code == ERROR) {
        /* try to free memory */
        free_2D_array(D);
    }
    return return_code;
}

/**
 * Calculates the diagonal degree matrix D from datapoint coordinates.
 *
 * @param datapoint_coords A 2D array of the datapoint coordinates.
 * @param N The number of datapoints in the array.
 * @param d The number of coordinates in each datapoint.
 * @param result The resulting diagonal degree matrix D.
 *
 * @return SUCCESS if the calculation was successful, ERROR otherwise.
 */
int ddg(double **datapoint_coords, int N, int d, double ***result) {
    int return_code = ERROR;
    double **A = NULL;
    double **D = NULL;

    GOTO_CLEANUP_IF_ERROR(sym(datapoint_coords, N, d, &A));
    GOTO_CLEANUP_IF_ERROR(_ddg(A, N, &D));

    return_code = SUCCESS;
    *result = D;

cleanup:
    free_2D_array(A);
    if (return_code == ERROR) {
        /* try to free memory */
        free_2D_array(D);
    }
    return return_code;
}

/**
 * Calculates D^(-0.5) for the normalized similarity matrix (in-place).
 *
 * @param D The diagonal degree matrix.
 * @param N The number of rows and columns in D.
 */
void _D_pow(double** D, int N) {
    int i;
    
    for (i=0; i<N; i++) {
        D[i][i] = pow(D[i][i], -0.5);
    }
}

/**
 * Computes the matrix multiplication of two N x N matrices A and B.
 *
 * @param A The first matrix in the multiplication.
 * @param B The second matrix in the multiplication.
 * @param N The number of rows and columns in matrices A and B.
 * @param result A pointer to the resulting matrix from the multiplication.
 *
 * @return SUCCESS if the multiplication was successful, ERROR otherwise.
 */
int _mat_dot(double** A, double** B, int N, double ***result) {
    int i;
    int j;
    int k;
    int return_code = ERROR;
    double* arr = calloc(N*N, sizeof(double));
    double** mat = malloc(N*sizeof(double*));

    GOTO_CLEANUP_IF_NULL(arr);
    GOTO_CLEANUP_IF_NULL(mat);

    for (i=0; i<N; i++) {
        mat[i] = arr + (i*N);
        for (j=0; j<N; j++) {
            for (k=0; k<N; k++) {
                mat[i][j] += A[i][k] * B[k][j];
            }
        }
    }

    return_code = SUCCESS;
    *result = mat;

cleanup:
    if (return_code == ERROR) {
        /* try to free memory */
        free_2D_array(mat);
    }
    return return_code;
}

/**
 * Calculates the normalized similarity matrix W.
 *
 * @param datapoint_coords A 2D array of the datapoint coordinates.
 * @param N The number of datapoints in the array.
 * @param d The number of coordinates in each datapoint.
 * @param result The resulting normalized similarity matrix W.
 *
 * @return SUCCESS if the calculation was successful, ERROR otherwise.
 */
int norm(double** datapoint_coords, int N, int d, double ***result) {
    int return_code = ERROR;
    double** A = NULL;
    double** D = NULL;
    double** W_ = NULL; /* will be: D^(-0.5)*A */
    double** W = NULL; /* will be: W_*D^(-0.5) = D^(-0.5)*A*D^(-0.5) */

    GOTO_CLEANUP_IF_ERROR(sym(datapoint_coords, N, d, &A));
    GOTO_CLEANUP_IF_ERROR(_ddg(A, N, &D));

    _D_pow(D, N); /* D = D^(-0.5) */
    GOTO_CLEANUP_IF_ERROR(_mat_dot(D, A, N, &W_));
    GOTO_CLEANUP_IF_ERROR(_mat_dot(W_, D, N, &W));

    return_code = SUCCESS;
    *result = W;

cleanup:
    free_2D_array(A);
    free_2D_array(D);
    free_2D_array(W_);

    if (return_code == ERROR) {
        /* try to free memory */
        free_2D_array(W);
    }
    return return_code;
}

/* Prints a matrix to the console. */
void print_mat(double** mat, int N) {
    int i;
    int j;

    for (i=0; i<N; i++) {
        for (j=0; j<N; j++) {
            printf("%.4f", mat[i][j]);
            if (j == N-1) {
                printf("\n");
            }
            else {
                printf(",");
            }
        }
    }
}

int main(int argc, char *argv[]) {
    int return_code = ERROR;
    int d = 0;
    int N = 0;
    char *goal;
    char *filename;
    double **result = NULL;
    double **datapoints = NULL;

    if (argc != 3) {
        goto cleanup;
    }

    goal = argv[1];
    filename = argv[2];

    GOTO_CLEANUP_IF_ERROR(parse(filename, &d, &N, &datapoints));

    if (strcmp(goal, "sym") == 0) {
        GOTO_CLEANUP_IF_ERROR(sym(datapoints, N, d, &result));
    }
    else if (strcmp(goal, "ddg") == 0) {
        GOTO_CLEANUP_IF_ERROR(ddg(datapoints, N, d, &result));
    }
    else if (strcmp(goal, "norm") == 0) {
        GOTO_CLEANUP_IF_ERROR(norm(datapoints, N, d, &result));
    }
    else { /* invalid goal */
        goto cleanup;
    }

    return_code = SUCCESS;
    print_mat(result, N);

cleanup:
    free_2D_array(datapoints);
    free_2D_array(result);
    if (return_code == ERROR) {
        printf("%s", GENERIC_ERROR_MSG);
    }
    return return_code;
}
