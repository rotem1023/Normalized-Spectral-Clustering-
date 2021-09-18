#define PY_SSIZE_T_CLEAN  /* For all # variants of unit formats (s#, y#, etc.) use Py_ssize_t rather than int. */
#include <Python.h>       /* MUST include <Python.h>, this implies inclusion of the following standard headers:
                             <stdio.h>, <string.h>, <errno.h>, <limits.h>, <assert.h> and <stdlib.h> (if available). */
#include <math.h>         /* include <Python.h> has to be before any standard headers are included */
#include "spkmeans.h"


double ** PyObjToMat (PyObject* mat, int dimRow, int dimCol){
    double **output;
    int i,j;
    output=(double**)malloc(sizeof(double*)*dimRow);
    assert(output!=NULL);
    for (i=0;i<dimRow;i++){
        output[i]=(double*)malloc(sizeof(double)*dimCol);
        assert(output[i]!=NULL);
    }
    for ( i = 0; i < dimRow; i++)
        {
            for (j=0;j<dimCol;j++){
                output[i][j]=PyFloat_AsDouble(PyList_GetItem(mat,i*dimCol+j));
            }
        }
    return output;
}
PyObject* MatToPyObj (double** mat, int dimRow, int dimCol){
    PyObject *line;
    PyObject *ret;
    int i,j;
    ret=PyList_New(dimRow);
    for (i = 0; i < dimRow; i++){
        line=PyList_New(dimCol);
        for (j=0;j<dimCol;j++){
            PyList_SetItem(line,j,PyFloat_FromDouble(mat[i][j]));  
        }
        PyList_SetItem(ret,i,line);
    }
    return ret;
}

static PyObject* geo_capi(PyObject *self, PyObject *args)
{
    int dimRow;
    int dimCol;
    PyObject* a;
    PyObject* retPy;
    double ** aa;
    int goal;
    int k;
    int newK;
    double *** ret;
    /* This parses the Python arguments into a double (d)  variable named z and int (i) variable named n*/
    if(!PyArg_ParseTuple(args, "Oiiii", &a,&dimRow,&dimCol,&goal,&k)) {
        return NULL; /* In the CPython API, a NULL value is never valid for a
                        PyObject* so it is used to signal that an error has occurred. */
    }

/* This builds the answer ("d" = Convert a C double to a Python floating point number) back into a python object */
    aa=PyObjToMat(a,dimRow,dimCol);
    if (goal==1){
        ret=geo_c(aa,dimRow,dimCol,goal,k); /*ret[0]=T,ret[1][0][0]=newK */
        newK=ret[1][0][0];
        retPy=PyList_New(2);
        PyList_SetItem(retPy,0,MatToPyObj(ret[0],dimRow,newK));
        PyList_SetItem(retPy,1,PyFloat_FromDouble(newK));
        return Py_BuildValue("O", retPy); /*  Py_BuildValue(...) returns a PyObject*  */
    }
    geo_c(aa,dimRow,dimCol,goal,k);
    return Py_BuildValue("O",PyList_New(1));
}


static PyObject* fit_capi(PyObject *self, PyObject *args)
{
    int K;
    int max_iter;
    int lines;
    int dim;
    int dNew;
    int it,it2;
    int lenInitialCentroids;
    int lenFinalData;
    PyObject *initialCentroids;
    PyObject *finalData;
    PyObject *rawData;
    double* arrInitialCentroids;
    double* arrFinalData;
    double** rawDataMat;
    if(!PyArg_ParseTuple(args, "iiOOOiii", &K, &max_iter,&initialCentroids,&finalData,&rawData,&lines,&dim,&dNew)) {
        return NULL; 
    }
    lenInitialCentroids=PyObject_Length(initialCentroids);
    lenFinalData=PyObject_Length(finalData);
    arrInitialCentroids=(double *)malloc(sizeof(double)*(lenInitialCentroids+1));
    arrFinalData=(double *)malloc(sizeof(double)*(lenFinalData+1));
    for ( it = 0; it < lenInitialCentroids; it++)
    {
            arrInitialCentroids[it]=PyFloat_AsDouble(PyList_GetItem(initialCentroids,it));;

    }
    for ( it2 = 0; it2 < lenFinalData; it2++)
    {
        arrFinalData[it2]=PyFloat_AsDouble(PyList_GetItem(finalData,it2));;

    }
    rawDataMat=PyObjToMat(rawData,lines,dNew);
    return Py_BuildValue("O", MatToPyObj(mainFuncV2(K, max_iter, arrInitialCentroids,arrFinalData, rawDataMat, lines,dim,dNew),K,dNew)); 
}


/*
 * This array tells Python what methods this module has.
 * We will use it in the next structure
 */
static PyMethodDef capiMethods[] = {
    {"geo",                   /* the Python method name that will be used */
      (PyCFunction) geo_capi, /* the C-function that implements the Python function and returns static PyObject*  */
      METH_VARARGS,           /* flags indicating parametersaccepted for this function */
      PyDoc_STR("our main func")}, /*  The docstring for the function */
    {"fit",                   
      (PyCFunction) fit_capi, 
      METH_VARARGS,           
      PyDoc_STR("fit from kmeans++")}, /*  The docstring for the function */
    {NULL, NULL, 0, NULL}     /* The last entry must be all NULL as shown to act as a
                                 sentinel. Python looks for this entry to know that all
                                 of the functions for the module have been defined. */
};


/* This initiates the module using the above definitions. */
static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "capi_demo1", /* name of module */
    NULL, /* module documentation, may be NULL */
    -1,  /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
    capiMethods /* the PyMethodDef array from before containing the methods of the extension */
};


/*
 * The PyModuleDef structure, in turn, must be passed to the interpreter in the module’s initialization function.
 * The initialization function must be named PyInit_name(), where name is the name of the module and should match
 * what we wrote in struct PyModuleDef.
 * This should be the only non-static item defined in the module file
 */
PyMODINIT_FUNC
PyInit_capi_demo1(void)
{
    PyObject *m;
    m = PyModule_Create(&moduledef);
    if (!m) {
        return NULL;
    }
    return m;
}


