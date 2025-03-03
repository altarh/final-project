#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

const int SUCCESS = 0;
const int ERROR = 1;

const int DEFAULT_ITER = 300;

#define GOTO_CLEANUP_IF_NULL(x) { if (x == NULL) { goto cleanup; } }
#define GOTO_CLEANUP_IF_NEGATIVE(x) { if (x < 0) { goto cleanup; } }
#define GOTO_CLEANUP_IF_PYERROR_OCCURED() { if (NULL != PyErr_Occurred()) { goto cleanup; } }

const char *GENERIC_ERROR_MSG = "An Error Has Occurred\n";

typedef struct centroid {
    struct coord *centroid_coords;
    struct coord *sum;
    int count;
} centroid;

typedef struct coord
{
    double coord;
    struct coord *next;
} coord;

typedef struct datapoint
{
    struct coord *coords;
    struct datapoint *next;
} datapoint;

/* internal funcs for read_args */
int init_datapoint(struct datapoint **datapoint, struct coord *first_coord) {
    struct datapoint *new_datapoint = malloc(sizeof(struct datapoint));
    if (new_datapoint == NULL) {
        return ERROR;
    }
    new_datapoint->coords = first_coord;
    new_datapoint->next = NULL;

    *datapoint = new_datapoint;
    return SUCCESS;    
}

int init_coord(struct coord **coord, double n) {
    struct coord *new_coord = malloc(sizeof(struct coord));
    if (new_coord == NULL) {
        return ERROR;
    }
    new_coord->coord = n;
    new_coord->next = NULL;

    *coord = new_coord;
    return SUCCESS;
}

int parse_file(const char *filename, int *d, int *N, struct datapoint **datapoints) {
    double n; /* for the double values */
    char delim; /* commas, \n, ... */
    struct datapoint **curr_datapoint = datapoints;
    struct coord *first_coord = NULL;
    struct coord **curr_coord = NULL;

    FILE *file = fopen(filename, "r");
    GOTO_CLEANUP_IF_NULL(file);

    /* go over first line to get d */
    first_coord = NULL;
    curr_coord = &first_coord;
    do {
        fscanf(file, "%lf%c", &n, &delim);
        if (SUCCESS != init_coord(curr_coord, n)) {
            return ERROR;
        }
        curr_coord = &(*curr_coord)->next;
        (*d)++;
    } while (delim != '\n');

    /* initialize the first datapoint */
    if (SUCCESS != init_datapoint(curr_datapoint, first_coord)) {
        return ERROR;
    }
    curr_datapoint = &(*curr_datapoint)->next;
    (*N)++;

    /* go over the rest of the lines to get datapoints and N */
    first_coord = NULL;
    curr_coord = &first_coord;
    while (fscanf(file, "%lf%c", &n, &delim) == 2) {
        if (SUCCESS != init_coord(curr_coord, n)) {
            return ERROR;
        }
        curr_coord = &(*curr_coord)->next;

        if (delim == '\n') { /* if at the end of the line */
            if (SUCCESS != init_datapoint(curr_datapoint, first_coord)) {
                return ERROR;
            }
            curr_datapoint = &(*curr_datapoint)->next;
            first_coord = NULL;
            curr_coord = &first_coord;
            (*N)++;
        }
    }

cleanup:
    if (file != NULL) {
        fclose(file);
    }
    return SUCCESS;
}

void free_coords(struct coord *coord) {
    struct coord *curr_coord = coord;
    while (curr_coord != NULL) {
        struct coord *next_coord = curr_coord->next;
        free(curr_coord);
        curr_coord = next_coord;
    }
}

void free_datapoints_structs(struct datapoint *datapoints) {
    struct datapoint *curr_datapoint;
    struct datapoint *next_datapoint;

    curr_datapoint = datapoints;
    while (curr_datapoint != NULL) {
        next_datapoint = curr_datapoint->next;
        free_coords(curr_datapoint->coords);
        free(curr_datapoint);
        curr_datapoint = next_datapoint;
    }
}

void free_2D_array(double **array) {
    free(array[0]);
    free(array);
}

double** linked_list_to_2D_array(datapoint* point, int N, int d){
    int i, j;
    double *p;
    double **a;
    coord* coord1;

    p = calloc(d * N, sizeof(double));
    a = calloc(N, sizeof(double *));

    for (i = 0; i < N; i++) {
        coord1 = point->coords;
        a[i] = p + (i * d);
        for (j = 0; j < d; j++) {
            a[i][j] = coord1->coord;
            coord1 = coord1->next;
        }
        point = point->next;
    }

    return a;
}

double** parse(const char *filename, int *d, int *N) {
    struct datapoint *datapoints = NULL;
    double **datapoints_array = NULL;

    parse_file(filename, d, N, &datapoints);
    datapoints_array = linked_list_to_2D_array(datapoints, *N, *d);
    free_datapoints_structs(datapoints);

    return datapoints_array;
}

double calc_euclidean_distance(double *coord1, double *coord2, int d){
    double sum = 0;
    int i;

    for (i = 0; i < d; i++){
        sum += pow((coord1[i] - coord2[i]), 2);
    }
    return sqrt(sum);
}

double** sym(double** datapoint_coords, int N, int d) {
    int i, j; /* TODO: ask if initing more than 1 per line is allowed */
    double* arr = malloc(N*N*sizeof(double));
    double** mat = malloc(N*sizeof(double*));
    
    for (i=0; i<N; i++) {
        mat[i] = arr + (i*N);
        for (j=0; j<N; j++) {
            if (i == j) {
                mat[i][j] = 0;
            }
            else {
                double dist = calc_euclidean_distance(datapoint_coords[i], datapoint_coords[j], d);
                mat[i][j] = exp(-dist/2);
            }
        }
    }

    return mat;
}

double** _ddg(double** mat, int N) {
    int i,j; /* TODO: ask if initing more than 1 per line is allowed */
    double d = 0;

    for (i=0; i<N; i++) {
        for (j=0; j<N; j++) {
            d += mat[i][j];
            if (i != j) {
                mat[i][j] = 0;
            }
        }
        mat[i][i] = d;
    }

    return mat;
}

double** ddg(double** datapoint_coords, int N, int d) {
    double** A = sym(datapoint_coords, N, d);
    double** D = _ddg(A, N);

    /* TODO: free A */

    return D;
}

double** _mat_pow(double** mat, int N) {
    int i,j; /* TODO: ask if initing more than 1 per line is allowed */
    
    for (i=0; i<N; i++) {
        for (j=0; j<N; j++) {
            if (mat[i][j] != 0) {
                mat[i][j] = pow(mat[i][j], -0.5);
            }
        }
    }

    return mat;
}

double** _mat_dot(double** A, double** B, int N) {
    int i,j,k; /* TODO: ask if initing more than 1 per line is allowed */
    double* arr = calloc(N*N, sizeof(double));
    double** mat = malloc(N*sizeof(double*));

    for (i=0; i<N; i++) {
        mat[i] = arr + (i*N);
        for (j=0; j<N; j++) {
            for (k=0; k<N; k++) {
                mat[i][j] += A[i][k] * B[k][j];
            }
        }
    }

    return mat;
}

double** norm(double** datapoint_coords, int N, int d) {
    double** A = sym(datapoint_coords, N, d);
    double** D = _ddg(A, N);
    double** W;
    double ** temp;
    D = _mat_pow(D, N);
    temp = _mat_dot(D, A, N);
    W = _mat_dot(temp, D, N);

    /* A become D (in place) so no need to free A. */
    free_2D_array(D);
    free_2D_array(temp);
    return W;
}

int main(int argc, char *argv[]) {
    int d = 0;
    int N = 0;
    int i,j; 
    char *goal;
    char *filename;
    double **result = NULL;
    double **datapoints = NULL;

    if (argc != 3) {
        printf("An Error Has Occurred\n");
        return ERROR;
    }

    goal = argv[1];
    filename = argv[2];

    printf("Running %s with input %s\n", goal, filename);

    datapoints = parse(filename, &d, &N);

    if (strcmp(goal, "sym") == 0) {
        result = sym(datapoints, N, d);
    }
    else if (strcmp(goal, "ddg") == 0) {
        result = ddg(datapoints, N, d);
    }
    else if (strcmp(goal, "norm") == 0) {
        result = norm(datapoints, N, d);
    }

    for (i=0; i<N; i++) {
        for (j=0; j<N; j++) {
            printf("%.4f", result[i][j]);
            if (j == N-1) {
                printf("\n");
            }
            else {
                printf(", ");
            }
        }
    }

    free_2D_array(result);
    free_2D_array(datapoints);
    return SUCCESS;
}
