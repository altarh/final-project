#define PY_SSIZE_T_CLEAN
#include <Python.h>

#define GOTO_CLEANUP_IF_PYERROR_OCCURED() { if (NULL != PyErr_Occurred()) { goto cleanup; } }

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
