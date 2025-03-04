#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "symnmf.h"

#define MAX_ITERATIONS (300)

const int SUCCESS = 0;
const int ERROR = 1;
const double BETA = 0.5;

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

/* initializes 2D matrix of num_rows rows and num_cols columns with zeroes */
int init_matrix(matrix *m, int num_rows, int num_cols) {
    int i;
    int return_code = ERROR;
    double *arr = calloc(num_rows * num_cols, sizeof(double));
    matrix mat = calloc(num_rows, sizeof(double *));

    GOTO_CLEANUP_IF_NULL(arr);
    GOTO_CLEANUP_IF_NULL(mat);

    for (i = 0; i < num_rows; i++) {
        mat[i] = arr + (i * num_cols);
    }

    return_code = SUCCESS;
    *m = mat;

cleanup:
    if (return_code == ERROR) {
        /* try to free memory */
        free_2D_matrix(&mat);
    }
    return return_code;
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
    double n;   /* for the double values */
    char delim; /* commas, \n, ... */
    datapoint **curr_datapoint = datapoints;
    coord *first_coord = NULL;
    coord **curr_coord = &first_coord;

    FILE *file = fopen(filename, "r");
    GOTO_CLEANUP_IF_NULL(file);

    do { /* reading the first line of the file to extract d */
        fscanf(file, "%lf%c", &n, &delim);
        GOTO_CLEANUP_IF_ERROR(init_coord(curr_coord, n));
        curr_coord = &(*curr_coord)->next;
        (*d)++; /* counting number of coordinates */
    } while (delim != '\n');

    GOTO_CLEANUP_IF_ERROR(init_datapoint(curr_datapoint, first_coord));
    curr_datapoint = &(*curr_datapoint)->next; /* add first point to datapoint linked list */
    (*N)++;

    curr_coord = &first_coord;                                                  /* reset coordinates linked list */
    while (fscanf(file, "%lf%c", &n, &delim) == 2) {                            /* reading the rest of the file */
        GOTO_CLEANUP_IF_ERROR(init_coord(curr_coord, n));                       /* init coordinate */
        curr_coord = &(*curr_coord)->next;                                      /* add to coordinates linked list */
        if (delim == '\n') {                                                    /* if at the end of the line */
            GOTO_CLEANUP_IF_ERROR(init_datapoint(curr_datapoint, first_coord)); /* init datapoint */
            curr_datapoint = &(*curr_datapoint)->next;                          /* add it to the linked list */
            curr_coord = &first_coord; /* reset coordinates linked list for next datapoint */
            (*N)++;                    /* counting number of datapoints */
        }
    }

    return_code = SUCCESS;

cleanup:
    if (file != NULL)
        fclose(file);
    return return_code;
}

/* funcs for freeing allocated memory */
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

void free_2D_matrix(matrix *m) {
    if (*m != NULL) {
        free(*m[0]);
        free(*m);
        *m = NULL;
    }
}

void free_1D_array(double **array) {
    if (*array != NULL) {
        free(*array);
        *array = NULL;
    }
}

/**
 * Converts a linked list of datapoints into a 2D matrix.
 *
 * @param filename The name of the file to read.
 * @param d The calculated number of coordinates in each datapoint.
 * @param N The calculated number of datapoints in the linked list.
 * @param result The resulting 2D matrix of datapoint coordinates.
 *
 * @return SUCCESS if the conversion was successful, ERROR otherwise.
 */
int linked_list_to_2D_matrix(datapoint *point, int N, int d, matrix *result) {
    int i;
    int j;
    int return_code = ERROR;
    matrix m = NULL;
    coord *coord1 = NULL;

    GOTO_CLEANUP_IF_ERROR(init_matrix(&m, N, d));

    for (i = 0; i < N; i++) {
        coord1 = point->coords;
        for (j = 0; j < d; j++) {
            m[i][j] = coord1->coord;
            coord1 = coord1->next;
        }
        point = point->next;
    }

    return_code = SUCCESS;
    *result = m;

cleanup:
    if (return_code == ERROR) {
        /* try to free memory */
        free_2D_matrix(&m);
    }
    return return_code;
}

/**
 * Reads a file and parses its contents into a 2D matrix of datapoint coordinates.
 *
 * @param filename The name of the file to read.
 * @param d A pointer to the number of coordinates in each datapoint.
 * @param N A pointer to the number of datapoints in the file.
 * @param result A pointer to the resulting 2D matrix of datapoint coordinates.
 *
 * @return SUCCESS if the parsing was successful, ERROR otherwise.
 */
int parse(const char *filename, int *d, int *N, matrix *result) {
    int return_code = ERROR;
    datapoint *datapoints = NULL;

    GOTO_CLEANUP_IF_ERROR(parse_file(filename, d, N, &datapoints));
    GOTO_CLEANUP_IF_ERROR(linked_list_to_2D_matrix(datapoints, *N, *d, result));
    return_code = SUCCESS;

cleanup:
    free_datapoints_structs(datapoints);

    return return_code;
}

/* Calculates the squared Euclidean distance between two points in d-dimensional space.  */
double calc_squared_euclidean_distance(double *coord1, double *coord2, int d) {
    double sum = 0;
    int i;

    for (i = 0; i < d; i++) {
        sum += pow((coord1[i] - coord2[i]), 2);
    }
    return sum;
}

double calc_squared_frobenius_norm(matrix A, int num_rows, int num_cols) {
    int i;
    int j;
    double sum = 0;

    for (i = 0; i < num_rows; i++) {
        for (j = 0; j < num_cols; j++) {
            sum += pow(A[i][j], 2);
        }
    }
    return sum;
}

/**
 * Calculates the similarity matrix A.
 *
 * @param datapoint_coords A 2D matrix of the datapoint coordinates.
 * @param N The number of datapoints.
 * @param d The number of coordinates in each datapoint.
 * @param result The resulting similarity matrix.
 *
 * @return SUCCESS if the calculation was successful, ERROR otherwise.
 */
int sym_C(matrix datapoint_coords, int N, int d, matrix *result) {
    int i;
    int j;
    int return_code = ERROR;
    matrix A = NULL;

    GOTO_CLEANUP_IF_ERROR(init_matrix(&A, N, N));

    for (i = 0; i < N; i++) {
        for (j = i + 1; j < N; j++) {  /* diagonal initialized already to 0 */
            double dist = calc_squared_euclidean_distance(datapoint_coords[i], datapoint_coords[j], d);
            A[i][j] = exp(-dist / 2);
            A[j][i] = A[i][j];  /* symmetry between A[j][i] and A[i][j] */
        }
    }

    return_code = SUCCESS;
    *result = A;

cleanup:
    if (return_code == ERROR) {
        /* try to free memory */
        free_2D_matrix(&A);
    }
    return return_code;
}

/**
 * Internal function which calculates the degree values using the similarity matrix A.
 *
 * @param A The similarity matrix - a 2D matrix of size N x N.
 * @param N The number of rows and columns in A.
 * @param result The resulting array of degrees.
 *
 * @return SUCCESS if the calculation was successful, ERROR otherwise.
 */
int _calc_diag(matrix A, int N, matrix result) {
    int return_code = ERROR;
    int i;
    int j;
    double d;
    double *diag = calloc(N, sizeof(double));

    GOTO_CLEANUP_IF_NULL(diag);

    for (i = 0; i < N; i++) {
        d = 0;
        for (j = 0; j < N; j++) {
            d += A[i][j];
        }
        diag[i] = d;
    }

    return_code = SUCCESS;
    *result = diag;

cleanup:
    if (return_code == ERROR) {
        /* try to free memory */
        free_1D_array(&diag);
    }
    return return_code;
}

/**
 * Internal function which calculates the diagonal degree matrix D using the similarity matrix A.
 *
 * @param A The similarity matrix - a 2D matrix of size N x N.
 * @param N The number of rows and columns in A.
 * @param result The resulting diagonal degree matrix D.
 *
 * @return SUCCESS if the calculation was successful, ERROR otherwise.
 */
int _ddg(matrix A, int N, matrix *result) {
    int i;
    int return_code = ERROR;
    double *diag = NULL;
    matrix D = NULL;

    GOTO_CLEANUP_IF_ERROR(init_matrix(&D, N, N));  /* initialized with zeroes */
    GOTO_CLEANUP_IF_ERROR(_calc_diag(A, N, &diag));

    for (i = 0; i < N; i++) {
        D[i][i] = diag[i];
    }

    return_code = SUCCESS;
    *result = D;

cleanup:
    free_1D_array(&diag);
    if (return_code == ERROR) {
        /* try to free memory */
        free_2D_matrix(&D);
    }
    return return_code;
}

/**
 * Calculates the diagonal degree matrix D from datapoint coordinates.
 *
 * @param datapoint_coords A 2D matrix of the datapoint coordinates.
 * @param N The number of datapoints.
 * @param d The number of coordinates in each datapoint.
 * @param result The resulting diagonal degree matrix D.
 *
 * @return SUCCESS if the calculation was successful, ERROR otherwise.
 */
int ddg_C(matrix datapoint_coords, int N, int d, matrix *result) {
    int return_code = ERROR;
    matrix A = NULL;
    matrix D = NULL;

    GOTO_CLEANUP_IF_ERROR(sym_C(datapoint_coords, N, d, &A));
    GOTO_CLEANUP_IF_ERROR(_ddg(A, N, &D));

    return_code = SUCCESS;
    *result = D;

cleanup:
    free_2D_matrix(&A);
    if (return_code == ERROR) {
        /* try to free memory */
        free_2D_matrix(&D);
    }
    return return_code;
}

/**
 * Calculates the values on the diagonal of D^(-0.5) for the normalized similarity matrix (in-place).
 *
 * @param diag The degrees array.
 * @param N The length of diag.
 */
void _calc_diag_pow(double *diag, int N) {
    int i;
    for (i = 0; i < N; i++) {
        diag[i] = pow(diag[i], -0.5);
    }
}

/**
 * Calculates the normalized similarity matrix W.
 *
 * @param datapoint_coords A 2D matrix of the datapoint coordinates.
 * @param N The number of datapoints.
 * @param d The number of coordinates in each datapoint.
 * @param result The resulting normalized similarity matrix W.
 *
 * @return SUCCESS if the calculation was successful, ERROR otherwise.
 */
int norm_C(matrix datapoint_coords, int N, int d, matrix *result) {
    int return_code = ERROR;
    int i;
    int j;
    matrix A = NULL;
    double *diag = NULL;
    matrix W = NULL;

    GOTO_CLEANUP_IF_ERROR(init_matrix(&W, N, N));
    GOTO_CLEANUP_IF_ERROR(sym_C(datapoint_coords, N, d, &A));
    GOTO_CLEANUP_IF_ERROR(_calc_diag(A, N, &diag)); /* calculating with diag array instead of matrix, for efficiency */
    _calc_diag_pow(diag, N);

    for (i = 0; i < N; i++) {
        /* not setting the diagonal of W as it is always 0 */
        for (j = i; j < N; j++) {
            W[i][j] = A[i][j] * diag[i] * diag[j]; /* multiplication of the matrixes */
            W[j][i] = W[i][j];                     /* symmetry between W[j][i] and W[i][j] */
        }
    }

    return_code = SUCCESS;
    *result = W;

cleanup:
    free_2D_matrix(&A);
    free_1D_array(&diag);
    if (return_code == ERROR) {
        /* try to free memory */
        free_2D_matrix(&W);
    }
    return return_code;
}

/**
 * Computes the matrix multiplication of two matrices A and B.
 * Expecting amount of rows in B to be equal to amount of columns in A.
 *
 * @param A The first matrix in the multiplication.
 * @param B The second matrix in the multiplication.
 * @param N The number of rows and columns in matrices A and B.
 * @param result A pointer to the resulting matrix from the multiplication.
 *
 * @return SUCCESS if the multiplication was successful, ERROR otherwise.
 */
int mat_dot(matrix A, matrix B, int A_nrows, int A_ncols, int B_ncols, matrix *result) {
    int i;
    int j;
    int k;
    int return_code = ERROR;
    int new_nrows = A_nrows;
    int new_ncols = B_ncols;
    matrix mat = NULL;

    GOTO_CLEANUP_IF_ERROR(init_matrix(&mat, new_nrows, new_ncols));

    for (i = 0; i < new_nrows; i++) {
        for (j = 0; j < new_ncols; j++) {
            for (k = 0; k < A_ncols; k++) {
                mat[i][j] += A[i][k] * B[k][j];
            }
        }
    }

    return_code = SUCCESS;
    *result = mat;

cleanup:
    if (return_code == ERROR) {
        /* try to free memory */
        free_2D_matrix(&mat);
    }
    return return_code;
}

/**
 * Computes and returns the transpose of a N x N matrix.
 *
 * @param mat The matrix to transpose.
 * @param num_rows The number of rows in the matrix
 * @param num_cols The number of columns in the matrix
 * @param result A pointer to the resulting matrix.
 *
 * @return SUCCESS if the multiplication was successful, ERROR otherwise.
 */
int transpose(matrix mat, int num_rows, int num_cols, matrix *result) {
    int i;
    int j;
    int return_code = ERROR;
    matrix mat_transpose = NULL;

    GOTO_CLEANUP_IF_ERROR(init_matrix(&mat_transpose, num_cols, num_rows));

    for (i = 0; i < num_cols; i++) {
        for (j = 0; j < num_rows; j++) {
            mat_transpose[i][j] = mat[j][i];
        }
    }

    return_code = SUCCESS;
    *result = mat_transpose;

cleanup:
    if (return_code == ERROR) {
        /* try to free memory */
        free_2D_matrix(&mat);
    }
    return return_code;
}

int update_H(matrix H, matrix W, int N, int k, matrix *result) {
    int i;
    int j;
    int return_code = ERROR;
    matrix temp = NULL;  /* for intermediate matrix multiplication */
    matrix numerator = NULL;
    matrix denumerator = NULL;
    matrix H_transpose = NULL;
    matrix new_H = NULL;

    GOTO_CLEANUP_IF_ERROR(init_matrix(&new_H, N, k));
    GOTO_CLEANUP_IF_ERROR(mat_dot(W, H, N, N, k, &numerator));  /* TODO H is size n X k */
    GOTO_CLEANUP_IF_ERROR(transpose(H, N, k, &H_transpose));
    GOTO_CLEANUP_IF_ERROR(mat_dot(H, H_transpose, N, k, N, &temp));  /* temp is size n X n */
    GOTO_CLEANUP_IF_ERROR(mat_dot(temp, H, N, N, k, &denumerator));

    for (i = 0; i < N; i++) {
        for (j = 0; j < k; j++) {
            new_H[i][j] = H[i][j] * (1 - BETA + (BETA * (numerator[i][j] / denumerator[i][j])));
        }
        
    }

    return_code = SUCCESS;
    *result = new_H;

cleanup:
    free_2D_matrix(&temp);
    free_2D_matrix(&numerator);
    free_2D_matrix(&denumerator);
    free_2D_matrix(&H_transpose);
    if (return_code == ERROR) {
        /* try to free memory */
        free_2D_matrix(&new_H);
    }
    return return_code;
}

int calculate_distance(matrix prev_H, matrix new_H, int N, int k, double *result) {
    int i;
    int j;
    double frobenius_norm;
    int return_code = ERROR;
    matrix diff = NULL;

    GOTO_CLEANUP_IF_ERROR(init_matrix(&diff, N, k));

    for (i = 0; i < N; i++) {
        for (j = 0; j < k; j++) {
            diff[i][j] = new_H[i][j] - prev_H[i][j];
        }
    }

    frobenius_norm = calc_squared_frobenius_norm(diff, N, k);

    return_code = SUCCESS;
    *result = frobenius_norm;

cleanup:
    free_2D_matrix(&diff);
    return return_code;
}

int symnmf_C(matrix H, matrix W, int N, int k, matrix *result) {
    int i;
    int return_code = ERROR;
    double distance = 0;
    const double EPSILON = exp(-4);
    matrix new_H = NULL;

    for (i = 0; i < MAX_ITERATIONS; i++) {
        GOTO_CLEANUP_IF_ERROR(update_H(H, W, N, k, &new_H));
        /* check for convergence */
        GOTO_CLEANUP_IF_ERROR(calculate_distance(H, new_H, N, k, &distance));
        free_2D_matrix(&H);
        H = new_H;

        if (distance < EPSILON) {
            break;
        }
    }

    return_code = SUCCESS;
    *result = H;

cleanup:
    return return_code;
}

/* Prints a matrix to the console. */
void print_mat(matrix mat, int N) {
    int i;
    int j;

    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            printf("%.4f", mat[i][j]);
            if (j == N - 1) {
                printf("\n");
            } else {
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
    matrix result = NULL;
    matrix datapoints = NULL;

    if (argc != 3) {
        goto cleanup;
    }

    goal = argv[1];
    filename = argv[2];

    GOTO_CLEANUP_IF_ERROR(parse(filename, &d, &N, &datapoints));

    if (strcmp(goal, "sym") == 0) {
        GOTO_CLEANUP_IF_ERROR(sym_C(datapoints, N, d, &result));
    }
    else if (strcmp(goal, "ddg") == 0) {
        GOTO_CLEANUP_IF_ERROR(ddg_C(datapoints, N, d, &result));
    }
    else if (strcmp(goal, "norm") == 0) {
        GOTO_CLEANUP_IF_ERROR(norm_C(datapoints, N, d, &result));
    }
    else { /* invalid goal */
        goto cleanup;
    }

    return_code = SUCCESS;
    print_mat(result, N);

cleanup:
    free_2D_matrix(&datapoints);
    free_2D_matrix(&result);
    if (return_code == ERROR) {
        printf("%s", GENERIC_ERROR_MSG);
    }
    return return_code;
}
