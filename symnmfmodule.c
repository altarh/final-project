#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include "symnmf.h"

typedef int handler(double **, int, int, double ***);

static PyObject* sym(PyObject *self, PyObject *args);
static PyObject* ddg(PyObject *self, PyObject *args);
static PyObject* norm(PyObject *self, PyObject *args);

static PyMethodDef symnmfMethods[] = {
    {"sym", (PyCFunction)sym, METH_VARARGS, PyDoc_STR("sym")},
    {"ddg", (PyCFunction)ddg, METH_VARARGS, PyDoc_STR("ddg")},
    {"norm", (PyCFunction)norm, METH_VARARGS, PyDoc_STR("norm")},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef symnmfmodule = {
    PyModuleDef_HEAD_INIT,
    "symnmfmodule",
    NULL,
    -1,
    symnmfMethods
};

PyMODINIT_FUNC PyInit_symnmfmodule(void) {
    PyObject *m;
    m = PyModule_Create(&symnmfmodule);
    if (!m) {
        return NULL;
    }
    return m;
}

int get_datapoints(int N, int d, PyObject *datapoints, double ***datapoints_array_out) {
    int return_code = ERROR;
    int i;
    int j;
    double *p = NULL;
    double **a = NULL;

    p = calloc(d * N, sizeof(double));
    GOTO_CLEANUP_IF_NULL(p);
    a = calloc(N, sizeof(double *));
    GOTO_CLEANUP_IF_NULL(a);

    for (i = 0; i < N; i++) {
        a[i] = p + (i * d);
        PyObject *curr_item = PyList_GetItem(datapoints, i);
        GOTO_CLEANUP_IF_NULL(curr_item);
        for (j = 0; j < d; j++) {
            PyObject *curr_coor = PyList_GetItem(curr_item, j);
            GOTO_CLEANUP_IF_NULL(curr_coor);
            a[i][j] = PyFloat_AsDouble(curr_coor);
            GOTO_CLEANUP_IF_PYERROR_OCCURED();
        }
    }

    return_code = SUCCESS;
    *datapoints_array_out = a;

cleanup:
    if (return_code == ERROR) {
        /* try to free memory */
        free_2D_matrix(a);
    }
    return return_code;
}

int parse_datapoints(PyObject *self, PyObject *args, int *d, int *N, double ***datapoints_array_out) {
    int return_code = ERROR;
    PyObject *curr_item = NULL;
    PyObject *datapoints = NULL;

    if (!PyArg_ParseTuple(args, "O", &datapoints)) {
        goto cleanup;
    }

    *N = PyObject_Length(datapoints);
    GOTO_CLEANUP_IF_NEGATIVE(*N);

    curr_item = PyList_GetItem(datapoints, 0);
    GOTO_CLEANUP_IF_NULL(curr_item);
    *d = PyObject_Length(curr_item);
    GOTO_CLEANUP_IF_NEGATIVE(*d);

    GOTO_CLEANUP_IF_ERROR(get_datapoints(*N, *d, datapoints, datapoints_array_out));

    return_code = SUCCESS;

cleanup:
    return return_code;
}

int build_output_matrix(int N, double **matrix, PyObject **result) {
    int return_code = ERROR;
    PyObject *output = NULL;
    PyObject *curr_coor = NULL;
    int i;
    int j;

    output = PyList_New(N);
    GOTO_CLEANUP_IF_NULL(output);

    for (i = 0; i < N; i++) {
        PyObject *row = PyList_New(N);
        GOTO_CLEANUP_IF_NULL(row);
        for (j = 0; j < N; j++) {
            curr_coor = PyFloat_FromDouble(matrix[i][j]);
            GOTO_CLEANUP_IF_NULL(curr_coor);
            GOTO_CLEANUP_IF_NEGATIVE(PyList_SetItem(row, j, curr_coor));
        }
        GOTO_CLEANUP_IF_NEGATIVE(PyList_SetItem(output, i, row));
    }

    return_code = SUCCESS;
    *result = output;

cleanup:
    return return_code;
}

static PyObject *call_handler(PyObject *self, PyObject *args, handler handler_func) {
    int N;
    int d;
    double **datapoints = NULL;
    double **result = NULL;
    PyObject *matrix_result = NULL;  /* output value */

    GOTO_CLEANUP_IF_ERROR(parse_datapoints(self, args, &d, &N, &datapoints));
    GOTO_CLEANUP_IF_ERROR(handler_func(datapoints, N, d, &result));
    GOTO_CLEANUP_IF_ERROR(build_output_matrix(N, result, &matrix_result));

cleanup:
    free_2D_matrix(datapoints);
    free_2D_matrix(result);
    return matrix_result;
}

static PyObject *sym(PyObject *self, PyObject *args) {
    return call_handler(self, args, sym_C);
}

static PyObject *ddg(PyObject *self, PyObject *args) {
    return call_handler(self, args, ddg_C);
}

static PyObject *norm(PyObject *self, PyObject *args) {
    return call_handler(self, args, norm_C);
}
