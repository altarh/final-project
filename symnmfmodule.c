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

/**
 * Get content from 2D python list into a 2D matrix, given number of rows and columns.
 */
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

/**
 * Reads a 2D python list and parses its contents into a 2D matrix.
 *
 * @param obj The python object to obtain data from.
 * @param num_rows A pointer to the number of rows in the matrix.
 * @param num_cols A pointer to the number of columns in the matrix.
 * @param result A pointer to the resulting 2D matrix.
 *
 * @return SUCCESS if the parsing was successful, ERROR otherwise.
 */
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

/**
 * Builds a python list of lists from a 2D matrix.
 *
 * @param num_rows The number of rows in the matrix.
 * @param num_cols The number of columns in the matrix.
 * @param m The 2D matrix to build the python list from.
 * @param result A pointer to the resulting python list.
 *
 * @return SUCCESS if the construction was successful, ERROR otherwise.
 */
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

/**
 * Generic helper function which calls a handler function with datapoints and its dimensions after parsing a 2D matrix object, 
 * and builds a python matrix from the result.
 */
static PyObject *call_handler(PyObject *args, handler handler_func) {
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

/**
 * Wrapper function for python API which calls symnmf_C and builds a python matrix from the result.
 */
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
    /* H is freed in symnmf_C */
    free_2D_matrix(&W);
    free_2D_matrix(&result);
    return matrix_result;
}

/**
 * Wrapper function for python API which calls sym_C.
 */
static PyObject *sym(PyObject *self, PyObject *args) {
    return call_handler(args, sym_C);
}

/**
 * Wrapper function for python API which calls ddg_C.
 */
static PyObject *ddg(PyObject *self, PyObject *args) {
    return call_handler(args, ddg_C);
}

/**
 * Wrapper function for python API which calls norm_C.
 */
static PyObject *norm(PyObject *self, PyObject *args) {
    return call_handler(args, norm_C);
}
