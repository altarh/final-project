extern const int SUCCESS;
extern const int ERROR;

int sym_C(double** datapoint_coords, int N, int d, double ***result);
int ddg_C(double **datapoint_coords, int N, int d, double ***result);
int norm_C(double** datapoint_coords, int N, int d, double ***result);

/* helper functions */
#define GOTO_CLEANUP_IF_NULL(x) { if ((x) == NULL) { goto cleanup; } }
#define GOTO_CLEANUP_IF_ERROR(x) { if ((x) == ERROR) { goto cleanup; } }
#define GOTO_CLEANUP_IF_NEGATIVE(x) { if ((x) < 0) { goto cleanup; } }
#define GOTO_CLEANUP_IF_PYERROR_OCCURED() { if (NULL != PyErr_Occurred()) { goto cleanup; } }

void free_2D_array(double **array);
