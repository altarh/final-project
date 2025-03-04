#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include "symnmf.h"

typedef int handler(matrix , int, int, matrix *);

static PyObject* symnmf(PyObject *self, PyObject *args);
static PyObject* sym(PyObject *self, PyObject *args);
static PyObject* ddg(PyObject *self, PyObject *args);
static PyObject* norm(PyObject *self, PyObject *args);

static PyMethodDef symnmfMethods[] = {
    {"symnmf", (PyCFunction)symnmf, METH_VARARGS, PyDoc_STR("symnmf")},
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

int get_matrix_content(int num_rows, int num_cols, PyObject *obj, matrix *result) {
    int return_code = ERROR;
    int i;
    int j;
    matrix mat = NULL;

    GOTO_CLEANUP_IF_ERROR(init_matrix(&mat, num_rows, num_cols));

    for (i = 0; i < num_rows; i++) {
        PyObject *curr_item = PyList_GetItem(obj, i);
        GOTO_CLEANUP_IF_NULL(curr_item);
        for (j = 0; j < num_cols; j++) {
            PyObject *curr_coor = PyList_GetItem(curr_item, j);
            GOTO_CLEANUP_IF_NULL(curr_coor);
            mat[i][j] = PyFloat_AsDouble(curr_coor);
            GOTO_CLEANUP_IF_PYERROR_OCCURED();
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

int parse_2D_matrix(PyObject *obj, int *num_rows, int *num_cols, matrix *result) {
    int return_code = ERROR;
    PyObject *row = NULL;

    *num_rows = PyObject_Length(obj);
    GOTO_CLEANUP_IF_NEGATIVE(*num_rows);

    row = PyList_GetItem(obj, 0);
    GOTO_CLEANUP_IF_NULL(row);
    *num_cols = PyObject_Length(row);
    GOTO_CLEANUP_IF_NEGATIVE(*num_cols);

    GOTO_CLEANUP_IF_ERROR(get_matrix_content(*num_rows, *num_cols, obj, result));

    return_code = SUCCESS;

cleanup:
    return return_code;
}

int build_output_matrix(int num_rows, int num_cols, matrix m, PyObject **result) {
    int return_code = ERROR;
    PyObject *output = NULL;
    PyObject *curr_coor = NULL;
    int i;
    int j;

    output = PyList_New(num_rows);
    GOTO_CLEANUP_IF_NULL(output);

    for (i = 0; i < num_rows; i++) {
        PyObject *row = PyList_New(num_cols);
        GOTO_CLEANUP_IF_NULL(row);
        for (j = 0; j < num_cols; j++) {
            curr_coor = PyFloat_FromDouble(m[i][j]);
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
    matrix datapoints = NULL;
    matrix result = NULL;
    PyObject *obj = NULL;
    PyObject *matrix_result = NULL;  /* output value */

    if (!PyArg_ParseTuple(args, "O", &obj)) {
        goto cleanup;
    }

    GOTO_CLEANUP_IF_ERROR(parse_2D_matrix(obj, &N, &d, &datapoints));
    GOTO_CLEANUP_IF_ERROR(handler_func(datapoints, N, d, &result));
    GOTO_CLEANUP_IF_ERROR(build_output_matrix(N, N, result, &matrix_result));

cleanup:
    free_2D_matrix(&datapoints);
    free_2D_matrix(&result);
    return matrix_result;
}

static PyObject *symnmf(PyObject *self, PyObject *args) {
    int N;
    int k;
    matrix H = NULL;
    matrix W = NULL;
    matrix result = NULL;
    PyObject *H_obj = NULL;
    PyObject *W_obj = NULL;
    PyObject *matrix_result = NULL;  /* output value */
    
    if (!PyArg_ParseTuple(args, "OO", &H_obj, &W_obj)) {
        goto cleanup;
    }
    GOTO_CLEANUP_IF_ERROR(parse_2D_matrix(H_obj, &N, &k, &H));
    GOTO_CLEANUP_IF_ERROR(parse_2D_matrix(W_obj, &N, &N, &W));
    GOTO_CLEANUP_IF_ERROR(symnmf_C(H, W, N, k, &result));
    GOTO_CLEANUP_IF_ERROR(build_output_matrix(N, k, result, &matrix_result));

cleanup:
    free_2D_matrix(&H);
    free_2D_matrix(&W);
    free_2D_matrix(&result);
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
