#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

const int SUCCESS = 0;
const int ERROR = 1;

const double EPSILON = 0.001;
const int DEFAULT_ITER = 200;

const char *GENERIC_ERROR_MSG = "An Error Has Occurred\n";
const char *INVALID_K_ERROR_MSG = "Invalid number of clusters!\n";
const char *INVALID_ITER_ERROR_MSG = "Invalid maximum iteration!\n";

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

int parse_file(int *d, int *N, struct datapoint **datapoints) {
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

/* parse the args (K, iter and the datapoints) */
int read_args(int argc, char *argv[], int *K, int *iter, int *d, int *N, struct datapoint **datapoints) {
    /* Read arguments - argc should be 2 if there is not iter arg, 3 if there is */
    if (argc < 2 || argc > 3) {
        printf("%s", GENERIC_ERROR_MSG);
        return ERROR;
    }

    /* read and validate iter, if given */
    if (argc == 3) {
        if (SUCCESS != parse_integer(argv[2], iter) || *iter <= 1 || *iter >= 1000) {
            printf("%s", INVALID_ITER_ERROR_MSG);
            return ERROR;
        }
    }

    /* read and validate K from below */
    if (SUCCESS != parse_integer(argv[1], K) || *K <= 1) {
        printf("%s", INVALID_K_ERROR_MSG);
        return ERROR;
    }

    /* parses datapoints from file and obtains d and N */
    if (SUCCESS != parse_file(d, N, datapoints)) {
        printf("%s", GENERIC_ERROR_MSG);
        return ERROR;
    }

    /* validate K from above */
    if (*K >= *N) {
        printf("%s", INVALID_K_ERROR_MSG);
        return ERROR;
    }

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

double** linked_list_to_array(datapoint* dpoint){
    coord* crd = dpoint->coords
    
};

double** crate_similarity_matrix(){};

double** crate_diagonal_matrix(){};

double** crate_normalized_similarity_matrix(){};

double** crate_optimized_H(){};



int main(int argc, char *argv[]) {
    char *goal;
    char *filename;

    if (argc != 3) {
        printf("An Error Has Occurred\n");
        return ERROR;
    }

    goal = argv[1];
    filename = argv[2];

    printf("Running %s with input %s\n", goal, filename);

    return SUCCESS;
}
