#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

const int SUCCESS = 0;
const int ERROR = 1;

#define GOTO_CLEANUP_IF_NULL(x) { if (x == NULL) { goto cleanup; } }
#define GOTO_CLEANUP_IF_ERROR(x) { if (x == ERROR) { goto cleanup; } }
#define GOTO_CLEANUP_IF_NEGATIVE(x) { if (x < 0) { goto cleanup; } }

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
    int return_code = ERROR;
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
    if (array != NULL) {
        free(array[0]);
        free(array);
        array = NULL;
    }
}

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

int parse(const char *filename, int *d, int *N, double ***datapoints_array) {
    int return_code = ERROR;
    struct datapoint *datapoints = NULL;

    GOTO_CLEANUP_IF_ERROR(parse_file(filename, d, N, &datapoints));
    GOTO_CLEANUP_IF_ERROR(linked_list_to_2D_array(datapoints, *N, *d, datapoints_array));
    return_code = SUCCESS;

cleanup:
    free_datapoints_structs(datapoints);

    return return_code;
}

double calc_euclidean_distance(double *coord1, double *coord2, int d){
    double sum = 0;
    int i;

    for (i = 0; i < d; i++){
        sum += pow((coord1[i] - coord2[i]), 2);
    }
    return sqrt(sum);
}

int sym(double** datapoint_coords, int N, int d, double ***result) {
    int i;
    int j;
    int return_code = ERROR;
    double* arr = malloc(N*N*sizeof(double));
    double** mat = malloc(N*sizeof(double*));

    GOTO_CLEANUP_IF_NULL(arr);
    GOTO_CLEANUP_IF_NULL(mat);
    
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

    return_code = SUCCESS;
    *result = mat;

cleanup:
    if (return_code == ERROR) {
        /* try to free memory */
        free_2D_array(mat);
    }
    return return_code;
}

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

int norm(double** datapoint_coords, int N, int d, double ***result) {
    int return_code = ERROR;
    double** A = NULL;
    double** D = NULL;
    double** W_ = NULL; /* D^(-0.5)*A */
    double** W = NULL; /* D^(-0.5)*A*D^(-0.5) */

    GOTO_CLEANUP_IF_ERROR(sym(datapoint_coords, N, d, &A));
    GOTO_CLEANUP_IF_ERROR(_ddg(A, N, &D));

    D = _mat_pow(D, N);
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

void print_mat(double** mat, int N) {
    int i,j;

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
