import capi_demo1 as cp1
# To use this run in the terminal: python setup.py build_ext --inplace
# And on Nova: python3.8.5 setup.py build_ext --inplace
import sys
from random import choices
import numpy as np
import pandas as pd
from enum import Enum
from decimal import Decimal


class Goal (Enum):
    '''
    convert class from string to number 
    '''
    spk=1
    wam=2
    ddg=3
    lnorm=4
    jacobi=5

def convert_array_to_pandas(arr):
    '''
    convert list to data frame
    '''
    n = len(arr)
    pd1 = pd.DataFrame(arr)
    return pd1


def fromFile(path):
    '''
    function that reads the vectors data 
    into matrix
    '''
    output = []
    c = open(path, 'r')
    lines = c.readlines()
    for line in lines:
        arr = line.split(',')
        for i in range(len(arr)):
            arr[i] = float(arr[i])
        output.append(arr)
    return (output)

def fromFilecsv(path):
    '''
    function that reads the vectors data (csv file)
    into matrix
    '''
    tmp= pd.read_csv(path,header=None)
    names= tmp.columns
    tmp= tmp.to_numpy().tolist()
    return np.array(tmp)

def choose_the_next_center(file, vec_of_mins, n):
    '''
    function to get the first centers 
    according to the weigthed values
    '''
    cnt = 0
    total = np.sum(vec_of_mins)
    choice = [i for i in range(len(vec_of_mins))]
    v = vec_of_mins / total
    val = np.random.uniform(low=0, high=total, size=(1,))[0]
    for i in range(n):
        cnt = cnt + vec_of_mins[i]
        if cnt >= val:
            return i
    return 0


def min_distance(vec, val, k, i):
    '''
    function that calcultes the minimum dustance
    for vector to specific center
    '''
    distance = float('inf')
    for index in range(i):
        d = vec - val[index]
        d = np.sum(d * d)
        if (d < distance):
            distance = d
    return distance


def check_distance(file, val, n, k, i):
    '''
    function that calculates
    the distance of specific vector
    from all of the centers
    '''
    vec_of_mins = np.zeros((n))
    for index in range(n):
        vec = file[index]
        vec_of_mins[index] = min_distance(vec, val, k, i)
    return vec_of_mins


def intalize(k, file,indices):
    ''''
    function that initialize the first 
    k cnters according k-means++ method
    '''
    return_val = np.zeros((k, file.shape[1]))
    cnt = 0
    n= file.shape[0]
    np.random.seed(0)
    index=file.shape[0]
    index = np.random.choice(index)
    indices[0]=index
    if k > 0:
        cnt = cnt + 1
        return_val[0] = file[index]
    '''
    adding first k-1 centers
    '''
    while cnt < k:
        vec_of_mins = check_distance(file, return_val, n, k, cnt)
        indexx= np.arange(vec_of_mins.size)
        #index = choose_the_next_center(file, vec_of_mins, n)
        total=(np.sum(vec_of_mins))
        prob_values = (vec_of_mins) / total
        indexx = np.random.choice(indexx, p = prob_values)
        indices[cnt] = indexx
        return_val[cnt] = file[indexx]
        cnt = cnt + 1
    return return_val

def mainFunc(k, goal, path1):
    '''
    user function 
    gets the input from the user
    and decides what is the output according to this
    '''
    max_iter=300 
    n= len(path1)
    if path1[n-3:]=="txt":
        final_data = fromFile(path1)
        final_data=np.array(final_data)
    else:
        final_data= fromFilecsv(path1)
    if (goal==1): #spk
        datapointsAndK=cp1.geo(makeOne(final_data),final_data.shape[0],final_data.shape[1],goal,k) #datapoints=T matrix
        datapoints=datapointsAndK[0]
        k=int(datapointsAndK[1])
        if (k==0):
            return
        datapoints=np.array(datapoints)
        indices=[0 for i in range(k)]
        centers= intalize(k,datapoints,indices)
        for i in range(0,len(indices)):
            print(indices[i],end="")
            if (i!=len(indices)-1):
                print(",",end="")
            else:
                print("")
        rows= datapoints.shape[0]
        cols= datapoints.shape[1]
        cp1.fit(k,max_iter,makeOne(centers), makeOne(datapoints),makeOne(final_data), rows,cols,final_data.shape[1])
    else: #all the rest
        cp1.geo(makeOne(final_data),final_data.shape[0],final_data.shape[1],goal,k)

def makeOne(mat):
    ''''
    function that converts two dims
     matrix to one dim matrix'''
    arr1 = np.array(mat[0])
    for i in range (1,len(mat)):
        arr2=np.array(mat[i]);
        arr1=np.concatenate((arr1,arr2));
    return arr1.tolist();


def printSpecial(mat,lines,cols):
    ''''
    function that prints matrix
    according the instructions in the project
    '''
    for i in range(lines):
        for j in range(cols):
            print(round(mat[i*cols+j],3), end="")
            if (j != cols - 1):
                print(",", end="")
        print("\n", end="")

arr= sys.argv
if (len(arr) < 3):
    assert (3 < 0)
if (Decimal(arr[1])-int(Decimal(arr[1]))!=0 or Decimal(arr[1])<0):
    assert(1<0)
K=int(Decimal(arr[1]))
if (arr[2]!="spk"and arr[2]!="wam"and arr[2]!="jacobi" and arr[2]!="ddg" and arr[2]!="lnorm"):
    assert(1<0)
goal= Goal[arr[2]].value
file= arr[3]
temp=mainFunc(K, goal, file)
