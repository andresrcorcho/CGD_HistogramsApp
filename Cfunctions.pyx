import numpy as np
cimport numpy as cnp
import cython
from libc.math cimport sqrt, exp,pi
from cython.parallel import prange

ctypedef double dtype_t
#ctypedef int dtype_I

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)

#Gaussian of a single value

cdef inline double singleGauss(double x, double sigma) nogil:
    cdef double produc1=0,divide1=0,exp1=0,divide2=0,gaussian=0
    cdef double term1= sigma*sqrt(2.*pi)
    cdef double term2= 2. * sigma * sigma
    
    product1=-x*x
    divide1=product1/term2
    exp1=exp(divide1)
    divide2=exp1/term1
    gaussian=divide2
    return gaussian
    
cdef inline double substract(double a, double b) nogil:
    cdef double subs=0
    subs=a-b
    return subs
    
cdef inline double division(double a, double b) nogil:
    cdef double div=0
    div=a/b
    return div
    
cdef inline double suma(double a, double b) nogil:
    cdef double sum=0
    sum=a+b
    return sum

#optimized calculation of gaussian - almost pure C
cpdef solve_gaussianC(cnp.ndarray[dtype_t, ndim=1] x,double sigma):
    cdef cnp.ndarray[dtype_t, ndim=1] gaussian
    cdef int sizeA = x.shape[0]
    cdef int i=0
    cdef double produc1,divide1,exp1,divide2
    cdef double term1= sigma*sqrt(2.*pi)
    cdef double term2= 2. * sigma * sigma
    
    gaussian=np.zeros((sizeA), dtype=np.double)
    
    for i in prange(sizeA,nogil=True):
        gaussian[i]=singleGauss(x[i], sigma)
    return gaussian


#Calculate gaussian fuction from a dataset and sigma (standar deviation -error)
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
    
    for i in prange(sizeA,nogil=True):
        subs[i]=substract(arrayD[i],dat)
    return subs
    
    
#sum component by component two arrays of the same size
cpdef sumArrays (cnp.ndarray[dtype_t, ndim=1] arrayD1,cnp.ndarray[dtype_t, ndim=1] arrayD2):
    cdef cnp.ndarray[dtype_t, ndim=1] sumA
    cdef int i=0
    cdef int sizeA = arrayD1.shape[0]
    
    sumA=np.empty((sizeA), dtype=np.double)
    
    for i in prange (sizeA,nogil=True):
        sumA[i]=suma(arrayD1[i],arrayD2[i])
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
    
    for i in prange (sizeA,nogil=True):
        divA[i]=division(arrayD[i],dividingNumber)
    return divA
    
#calculate the exp function of an entire array - element by element
cpdef expArray (cnp.ndarray[dtype_t, ndim=1] arrayD):
    cdef cnp.ndarray[dtype_t, ndim=1] expA
    cdef int i=0
    cdef int sizeA = arrayD.shape[0]
    
    expA=np.empty((sizeA), dtype=np.double)
    
    for i in range (sizeA):
        expA[i]=exp(arrayD[i])

    return expA
    
#function that detects 0 values or strings in input .txt files
cpdef GetisZeroOrString(arrayD,str error1Label,str error2Label):
    Errors=np.array([])

    cdef double auxValue
    cdef int i=0
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
    return error, Errors
    
cpdef PDP_C(cnp.ndarray[dtype_t, ndim=1] x1_grid,cnp.ndarray[dtype_t, ndim=1] data_array,cnp.ndarray[dtype_t, ndim=1] sigma_array):
    cdef int sizeA = data_array.shape[0]
    cdef int sizeB = x1_grid.shape[0]
    cdef int i=0,j=0
    cdef double auxF=0
    cdef cnp.ndarray[dtype_t, ndim=1] x
    cdef cnp.ndarray[dtype_t, ndim=1] PDF
    x=np.empty((sizeB), dtype=np.double)
    PDF=np.empty((sizeB), dtype=np.double)
 
    for i in range(sizeA):
        #for j in range(sizeB):
            #x[j]=x1_grid[j]-data_array[i]
        x=substractArray(x1_grid,data_array[i])
            #auxF=solve_gaussian_s(x[j],sigma_array[i])
            #PDF[j]=PDF[j]+auxF
        if i==0:
            PDF=solve_gaussianC(x,sigma_array[i])
        else:
            PDF=sumArrays(PDF,solve_gaussianC(x,sigma_array[i]))
    return divideArray(PDF,sizeA)

#KDE- Fixed -Method 1
#Implemented from the equations presented in:
#Vermeesch, P. (2012). On the visualisation of detrital age distributions.
#Chemical Geology, 312, 190-194.

cpdef KDEp_C(cnp.ndarray[dtype_t, ndim=1] x1_grid,cnp.ndarray[dtype_t, ndim=1] data_array,double bandwidth):
    cdef int sizeA = data_array.shape[0]
    cdef int sizeB = x1_grid.shape[0]
    cdef int i=0,j=0
    cdef cnp.ndarray[dtype_t, ndim=1] x
    cdef cnp.ndarray[dtype_t, ndim=1] KDE
    x=np.empty((sizeB), dtype=np.double)
    KDE=np.empty((sizeB), dtype=np.double)
        
    for i in range(sizeA):
        #for j in prange(sizeB):
            #x[j]=x1_grid[j]-data_array[i]
        x=substractArray(x1_grid,data_array[i])
            #auxF=solve_gaussian_s(x[j],bandwidth)
            #KDE[j]=KDE[j]+auxF
        if i==0:
            KDE=solve_gaussianC(x,bandwidth)
        else:
            KDE=sumArrays(KDE,solve_gaussianC(x,bandwidth))
            
    return divideArray(KDE,sizeA)
    
#cut ages in a defined interval
cpdef cutAgesHistogramC(cnp.ndarray[dtype_t, ndim=1] AgesArray, double Min, double Max):
    DataAgesHist=np.array([])
    cdef int i=0
    cdef int sizeA = AgesArray.shape[0]
    
    for i in range(sizeA):
        age=AgesArray[i]
        if InInterval(age,Min,Max):
            DataAgesHist=np.append(DataAgesHist,[age])
            
    return DataAgesHist

#count number of samples in the defined Min and Max interval
cpdef getLocalNsamplesC(cnp.ndarray[dtype_t, ndim=1] AgesArray,double Min, double Max):
        cdef int i=0, NumAges=0
        cdef int sizeA = AgesArray.shape[0]
        
        for i in range(sizeA):
            age=(AgesArray[i])
            if InInterval(age,Min,Max):
                NumAges=NumAges+1
        return NumAges
        
cdef inline bint InInterval(double age, double Min, double Max):
    cdef bint isIn=False
    if age>=Min and age<=Max:
        isIn=True
    return isIn


