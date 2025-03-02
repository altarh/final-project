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
}centroid;

typedef struct coord
{
    double coord;
    struct coord *next;
}coord;

typedef struct datapoint
{
    struct coord *coords;
    struct datapoint *next;
}datapoint;

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

double calc_euclidean_distance(double *coord1, double *coord2, int d){
    double sum = 0;
    int i;

    for (i = 0; i < d; i++){
        sum += pow((coord1[i] - coord2[i]), 2);
    }
    return sqrt(sum);
}

double** linked_list_to_2D_array(datapoint* point, int N, int d){
    int i, j;
    double *p;
    double **a;
    coord* coord1 = point->coords;

    p = calloc(d * N, sizeof(double));
    a = calloc(d, sizeof(double *));

    for (i = 0; i < N; i++) {
        a[i] = p + (i * d);
        for (j = 0; j < d; j++) {
            a[i][j] = coord1->coord;
            coord1 = coord1->next;
        }
        point = point->next;
        coord1 = point->coords;
    }

    return a;
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
};

double** _ddg(double** mat, int N) {
    int i,j; /* TODO: ask if initing more than 1 per line is allowed */
    int d = 0;

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

    return D;
};

double** norm(double** datapoint_coords, int N, int d) {
    double** A = sym(datapoint_coords, N, d);
    double** D = _ddg(A, N);
    
};


int main(int argc, char *argv[]) {
    int d = 0;
    int N = 0;
    char *goal;
    char *filename;
    struct datapoint *datapoints = NULL;

    if (argc != 3) {
        printf("An Error Has Occurred\n");
        return ERROR;
    }

    goal = argv[1];
    filename = argv[2];

    printf("Running %s with input %s\n", goal, filename);
    parse_file(filename, &d, &N, &datapoints);

    return SUCCESS;
}

/* TODO: for debugging, remove */
void print_datapoints(struct datapoint *datapoints) {
    while (datapoints) {
        struct coord *current_coord = datapoints->coords;
        while (current_coord) {
            printf("  %.4f", current_coord->coord);
            current_coord = current_coord->next;
        }
        printf("\n");
        datapoints = datapoints->next;
    }
}
