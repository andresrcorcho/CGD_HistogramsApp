import numpy as np
cimport numpy as cnp
import cython
from libc.math cimport sqrt, exp

ctypedef double dtype_t
#ctypedef int dtype_I

@cython.boundscheck(False)
@cython.wraparound(False)

cpdef solve_gaussian(cnp.ndarray[dtype_t, ndim=1] x,double sigma):
    cdef cnp.ndarray[dtype_t, ndim=1] gaussian
    
    cdef double pi=np.pi
    cdef double term1= sigma*sqrt(2.*pi)
    cdef double term2= 2. * sigma * sigma
    
    gaussian=divideArray((expArray(divideArray(-multArrays(x,x),term2))),(term1))
    return gaussian
    
#Substract from one array a double value and returns an array
cpdef substractArray(cnp.ndarray[dtype_t, ndim=1] arrayD,double dat):
    cdef cnp.ndarray[dtype_t, ndim=1] subs
    
    cdef int i=0
    cdef int sizeA = arrayD.shape[0]
    
    subs=np.empty((sizeA), dtype=np.double)
    
    for i in range(sizeA):
        subs[i]=arrayD[i]-dat
    return subs
    
    
#sum component by component two arrays of the same size
cpdef sumArrays (cnp.ndarray[dtype_t, ndim=1] arrayD1,cnp.ndarray[dtype_t, ndim=1] arrayD2):
    cdef cnp.ndarray[dtype_t, ndim=1] sumA
    cdef int i=0
    cdef int sizeA = arrayD1.shape[0]
    
    sumA=np.empty((sizeA), dtype=np.double)
    
    for i in range (sizeA):
        sumA[i]=arrayD1[i]+arrayD2[i]
    return sumA
    
#multiply component by component two arrays of the same size
cpdef multArrays (cnp.ndarray[dtype_t, ndim=1] arrayD1,cnp.ndarray[dtype_t, ndim=1] arrayD2):
    cdef cnp.ndarray[dtype_t, ndim=1] multA
    cdef int i=0
    cdef int sizeA = arrayD1.shape[0]
    
    multA=np.empty((sizeA), dtype=np.double)
    
    for i in range (sizeA):
            multA[i]=arrayD1[i]*arrayD2[i]
    return multA
    
#divide array by a double number and returns an array
cpdef divideArray (cnp.ndarray[dtype_t, ndim=1] arrayD,double dividingNumber):
    cdef cnp.ndarray[dtype_t, ndim=1] divA
    cdef int i=0
    cdef int sizeA = arrayD.shape[0]
    
    divA=np.empty((sizeA), dtype=np.double)
    
    for i in range (sizeA):
        divA[i]=arrayD[i]/dividingNumber
    return divA
    
cpdef expArray (cnp.ndarray[dtype_t, ndim=1] arrayD):
    cdef cnp.ndarray[dtype_t, ndim=1] expA
    cdef int i=0
    cdef int sizeA = arrayD.shape[0]
    
    expA=np.empty((sizeA), dtype=np.double)
    
    for i in range (sizeA):
        expA[i]=exp(arrayD[i])

    return expA
    
cpdef pythonArray_to_C_array(pythonArray):
    pass

    return 0

cpdef GetisZeroOrString(arrayD,str error1Label,str error2Label):
    Errors=np.array([])

    cdef double auxValue
    cdef int i=0,j=0,k=0,aux=0
    cdef int sizeA = arrayD.shape[0]
    cdef bint error = False
    
    for i in range(sizeA):
        try:
            auxValue=arrayD[i]
            if auxValue==0:
                Errors=np.append(Errors,[error1Label])
                error=True
                break
                
        except TypeError:
            Errors=np.append(Errors,[error2Label])
            error=True
            break
        
        #else:
            #try:
                #aux=int(auxValue)
            #except ValueError:
                
                #Errors=np.append(Errors,[error2Label])
                #error=True
    #print(arrayD)
    return error, Errors
    

