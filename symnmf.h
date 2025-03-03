extern const int SUCCESS;
extern const int ERROR;

int sym_C(double** datapoint_coords, int N, int d, double ***result);
int ddg_C(double **datapoint_coords, int N, int d, double ***result);
int norm_C(double** datapoint_coords, int N, int d, double ***result);

/* helper function */
void free_2D_array(double **array);
