extern const int SUCCESS;
extern const int ERROR;

typedef double ** matrix;

int sym_C(matrix datapoint_coords, int N, int d, matrix *result);
int ddg_C(matrix datapoint_coords, int N, int d, matrix *result);
int norm_C(matrix datapoint_coords, int N, int d, matrix *result);
int symnmf_C(matrix H, matrix W, int N, int k, matrix *result);

/* helper functions */
#define GOTO_CLEANUP_IF_NULL(x) { if ((x) == NULL) { goto cleanup; } }
#define GOTO_CLEANUP_IF_ERROR(x) { if ((x) == ERROR) { goto cleanup; } }
#define GOTO_CLEANUP_IF_NEGATIVE(x) { if ((x) < 0) { goto cleanup; } }
#define GOTO_CLEANUP_IF_PYERROR_OCCURED() { if (NULL != PyErr_Occurred()) { goto cleanup; } }

void free_2D_matrix(matrix array);
