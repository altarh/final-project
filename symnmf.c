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

struct centroid {
    struct coord *centroid_coords;
    struct coord *sum;
    int count;
};

struct coord
{
    double coord;
    struct coord *next;
};

struct datapoint
{
    struct coord *coords;
    struct datapoint *next;
};

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

    /* go over first line to get d */
    first_coord = NULL;
    curr_coord = &first_coord;
    do {
        scanf("%lf%c", &n, &delim);
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
    while (scanf("%lf%c", &n, &delim) == 2) {
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

    return SUCCESS;
}

int parse_integer(char *src, int *dest) {
    int n = 0;
    int chars_read = 0;
    int amount_parsed = 0;

    amount_parsed = sscanf(src, "%d%n", &n, &chars_read);
    if (amount_parsed != 1 || src[chars_read] != '\0') {
        /* failed to parse integer, or trailing characters found */
        if (src[chars_read] == '.') {
            char *trail_string = src + (chars_read + 1);
            int trail = 0;
            int trail_chars_read = 0;
            amount_parsed = sscanf(trail_string, "%d%n", &trail, &trail_chars_read);
            if (!trail && amount_parsed <= 1 && !trail_string[trail_chars_read]) {
                *dest = n;
                return SUCCESS;
            }
        }
        return ERROR;
    }

    *dest = n;
    return SUCCESS;
}

double calc_euclidean_distance(struct coord *coord1, struct coord *coord2, int d){
    double sum = 0;
    int i;

    for (i = 0; i < d; i++){
        sum += pow((coord1->coord - coord2->coord), 2);
        coord1 = coord1->next;
        coord2 = coord2->next;
    }
    return sqrt(sum);
}

double** crate_similarity_matrix(){};

double** crate_diagonal_matrix(){};

double** crate_normalized_similarity_matrix(){};

double** crate_optimized_H(){};



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
