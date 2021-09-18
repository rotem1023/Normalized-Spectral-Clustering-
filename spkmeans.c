#define PY_SSIZE_T_CLEAN                                 
#include <math.h>        
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#define EPSILON 0.000000000000001
#define MAX_LINE_LENGTH 10000

/*declerations*/
double get_distance(double* vec,double* center);
int get_center_index(double* vec, double** centers,int k);
double** get_final_centers(double** T, double** centers,double** data,int k,int n1,int n);
double ** oneDimToTwoDim(double *oneDim, int dimRow, int dimCol);
double * twoDimToOneDim(double **twoDim, int dimRow, int dimCol);
int goalEnum(char* goal);
static double distanceFromMean(int i1, int j1, int n,double vec [], double mean [] );
static void updateMean (int j, int place, int n,double vec [], double mean [], int numOfVec);
double** mainFuncV2 (int K,  int max_iter, double* arrInitialCentroids,double* arrFinalData, double** rawDataMat ,int lines, int dim,int dNew);
int AreSame(double a, double b);
double l2norm(double* point1,double* point2, int dim);
double calcWeight (double* point1, double* point2, int dim);
double** adjacencyMatrix(double*** setOfPoints, int n, int d);
double calc_index_mat( double** m_one,int i, double** m_two,int j, int row, int col);
double*** multiply(double *** mat1, double*** mat2, int N);
double*** transpose(double***A, int row);
double** copy_mat (double** matrix, int row, int col);
double*** get_matrix_d (double*** matrix, int row, int col, int isDDG);
double*** get_matrix_L (double*** d, double*** w, int row, int col);
int* largestOffDiagElem(double***mat, int row, int col);
double*** unitMatrix (int row, int col);
double frobenius(double***mat, int dim);
double diagSum(double***mat, int dim);
double* eigenVals (double*** mat, int dim);
void BubbleSort(double a[], int indices[], int array_size);
double* calcDifs(double* arr, int dim);
int getK(double* difs, int dim);
void putColumn(double*** targetMat,double*** sourceMat,int targetInd,int sourceInd, int dim);
double*** getKeigenVectors(double*** mat,int* indices,int k, int dim);
double calc_mat_i (double** matrix, int k, int i );
double*** get_matrix_T(double*** matrix, int row, int k);
double*** geo_c(double** a, int dimRow, int dimCol, int goal,int k);
void updateCheck (int n,double vec [], double mean []);
int toStop (int n,double vec [], double mean []);
void printMat(double** matrix, int rows, int columns);

/*END- declerations*/

double get_distance(double* vec,double* center){
    /*
    function that returns the distance 
    from specific vector to specific center
    */
    int n;
    double val;
    int i;
    int tmpOne;
    int tmpTwo;
    double tmp;
    val=0;
    tmpOne= sizeof(vec);
    tmpTwo=sizeof(double);
    n =tmpOne/tmpTwo;
    for (i=0;i<n;i++){
        tmp= vec[i]-center[i];
        val=val+(tmp*tmp);
    }
    return (val);
}

int get_center_index(double* vec, double** centers,int k){
    /*
    function that returns the index of a center in the 
    cnetres matrix
    */
    double min_val;
    int index;
    double* center;
    double val;
    int i;
    min_val=0;
    index=0;
    for (i=0;i<k ;i++){
        center=centers[i];
        val= get_distance(vec,center);
        if (i==0){
            min_val=val;
        }
        if (min_val>val){
            min_val=val;
            index=i;
        }
    }
    return (index);
}

double** get_final_centers(double** T, double** centers,double** data,int k,int n1,int n){
    /*
    function that calc from the final centers of T matrix 
    the final centers from  the initial and "full" vectors
    */
    int* num_of_vecs;
    double** final_centers;
    int i;
    int j;
    int index;
    double* vec;
    int constant;
    double* center;
    double* add_vec;
    num_of_vecs= (int*) calloc(k,sizeof(int));
    assert(num_of_vecs!=NULL&& "An Error Has Occured");
    final_centers= (double**) malloc(k*sizeof(double*));
    assert(final_centers!=NULL&& "An Error Has Occured");
    for (i=0; i<k; i++){
        final_centers[i]= (double*) calloc(n1,sizeof(double));
        assert(final_centers[i]!=NULL&& "An Error Has Occured");
    }
    for (i=0;i<n;i++){
        vec= T[i];
        index=get_center_index(vec,centers,k);
        constant= num_of_vecs[index];
        center= final_centers[index];
        add_vec=data[i];
        for (j=0; j<n1;j++){
            center[j]=(center[j]*constant+add_vec[j])/(constant+1);
        }
        num_of_vecs[index]=num_of_vecs[index]+1;
    }
    return (final_centers);
}

double ** oneDimToTwoDim(double *oneDim, int dimRow, int dimCol){
    /*
    function that convert one dim matrix (in the memory) to two 
    dims matrix
    */
    double **output;
    int i,j;
    output=(double**)malloc(sizeof(double*)*dimRow);
    assert(output!=NULL&& "An Error Has Occured");
    for (i=0;i<dimRow;i++){
        output[i]=(double*)malloc(sizeof(double)*dimCol);
        assert(output[i]!=NULL&& "An Error Has Occured");
    }
    for ( i = 0; i < dimRow; i++)
        {
            for (j=0;j<dimCol;j++){
                output[i][j]=oneDim[j+dimCol*i];
            }
        }
    return output;
}

double * twoDimToOneDim(double **twoDim, int dimRow, int dimCol){
    /*
    function that convert two dims matrix (in the memory) to one 
    dim matrix
    */
    double *output;
    int i,j;
    output=(double*)malloc(sizeof(double)*((dimRow+1)*(dimCol+1)));
    assert(output!=NULL&& "An Error Has Occured");
    for ( i = 0; i < dimRow; i++)
        {
            for (j=0;j<dimCol;j++){
                output[j+dimCol*i]=twoDim[i][j];
            }
        }
    return output;
}
int goalEnum(char* goal){
    /*
    convert func from string to number 
    */
    if (strcmp(goal,"spk")==0){
        return 1;
    }
    if (strcmp(goal,"wam")==0){
        return 2;
    }
    if (strcmp(goal,"ddg")==0){
        return 3;
    }
    if (strcmp(goal,"lnorm")==0){
        return 4;
    }
    if (strcmp(goal,"jacobi")==0){
        return 5;
    }
    printf("Invalid Input!\n");
    return 0;
}


/*---KMEANS PLUS PLUS---*/

static double distanceFromMean(int i1, int j1, int n,double vec [], double mean [] ){
    /*
    function that calcultes the  distance
    for specific vector to specific center
    */
    double d=0;
    double meanI;
    double x;
    double distance;
    int i;
    for ( i = 0; i < n; i++)
    {
        meanI= mean[n*j1+i];
        x=vec[n*i1+i];
        distance= (x-meanI)*(x-meanI);
        d=d+distance;         
    }
    return (d);
}

static void updateMean (int j, int place, int n,double vec [], double mean [], int numOfVec){
    /*
    function that updates specific cnter
    according to vector that we added to 
    the center
    */
    int i;
    for ( i = 0; i < n; i++)
    {
        double meanI= mean[i+place*n];
        double sum;
        if (numOfVec>0){
            sum= meanI*numOfVec+vec[i+j*n];
        }
        else{
            sum= vec[i+j*n];
        }
        mean[i+place*n]= sum/(numOfVec+1); 
    }
}

double** mainFuncV2 (int K,  int max_iter, double* arrInitialCentroids,double* arrFinalData, double** rawDataMat ,int lines, int dim,int dNew) {
    /*
    k means implemntion after initialize the first centers 
    */
    int i,j,cc,b1,k1,k2,iter,place,d;
    double distance,curDistance;
    double *vecOfMeansOld;
    double *vecOfMeansNew;
    double *vector;
    double ** vecOfMeansOldMat;
    double ** finalCenters;
    int *numberOfVecPerGroup;
    i=0;
    numberOfVecPerGroup=(int *)malloc(sizeof(int)*K);
    assert(numberOfVecPerGroup!=NULL&& "An Error Has Occured");
    i=0;
    iter=0;
    for ( iter = 0; iter < K; iter++)
    {
        numberOfVecPerGroup[iter]=0;
    }
    iter=0;
    vecOfMeansOld = (double *)malloc(sizeof(double)*(K+1)*(dim+1));
    assert(vecOfMeansOld!=NULL&& "An Error Has Occured");
    vecOfMeansNew = (double *)malloc(sizeof(double)*(K+1)*(dim+1));
    assert(vecOfMeansNew!=NULL&& "An Error Has Occured");
    vector= (double *)malloc(sizeof(double)*(lines+1)*(dim+1));
    assert(vector!=NULL&&"An Error Has Occured");
    for ( i = 0; i < lines; i++)
    {
        for ( cc = 0; cc < dim; cc++)
        {
        vector[i*dim+cc]= arrFinalData[i*dim+cc];    
        }   
    }
    /* get the first centers*/
    for ( i = 0; i < K; i++)
    {
        for ( cc = 0; cc < dim; cc++)
        {
        vecOfMeansOld[i*dim+cc]= arrInitialCentroids[i*dim+cc];    
        }   
    }
    d=dim;
    /* update the centers*/
    while(iter<max_iter){
            for ( i = 0; i < lines ; i++)
            {
                distance=0;
                place=0;
                
                for (  j = 0; j < K; j++)
                {
                    curDistance= distanceFromMean(i,j, d, vector, vecOfMeansOld );
                    if (j==0){
                            distance=curDistance;
                    }
                    if(curDistance< distance){
                        distance=curDistance;
                        place=j;
                    } 
                }
                updateMean (i,place,d, vector, vecOfMeansNew , numberOfVecPerGroup[place]);
                numberOfVecPerGroup[place]=numberOfVecPerGroup[place]+1;
            }

            b1=1;
            for (k1=0; k1<K; k1++)
            {
                for (k2 = 0; k2 < dim; k2++)
                {
                    if(vecOfMeansOld[k1*dim+k2]!=vecOfMeansNew[k1*dim+k2]){
                        b1=0;
                    }
                }    
            }
            if (b1==1){
                break;
            }
            for (k1=0; k1<K; k1++)
            {
                for (k2 = 0; k2 < dim; k2++)
                {
                    vecOfMeansOld[k1*dim+k2]=vecOfMeansNew[k1*dim+k2];
                    vecOfMeansNew[k1*dim+k2]=0;
                }
                
            }
            for ( i = 0; i < K; i++)
            {
                numberOfVecPerGroup[i]=0;
            }
            iter++;
        }
    
    vecOfMeansOldMat=oneDimToTwoDim(vecOfMeansOld,K,d);
    finalCenters=get_final_centers(oneDimToTwoDim(arrFinalData,lines,dim),vecOfMeansOldMat,rawDataMat,K,dNew,lines);
    for (i=0;i<K;i++){
        for (j=0;j<d;j++){
            printf("%.4f",vecOfMeansOldMat[i][j]);
            if(j!=d-1){
                printf(",");
            }
        }
        printf("\n");
    }
    free(vecOfMeansOld);
    free(numberOfVecPerGroup);
    free(vecOfMeansNew);
    free(arrInitialCentroids);
    free(arrFinalData);
    free(vector);
    /*ret=MatToPyObj(finalCenters,K,dNew);*/
    return finalCenters;
}




/*---NEW PROJECT---*/

int AreSame(double a, double b)
{
    return fabs(a - b) < EPSILON;
}

double l2norm(double* point1,double* point2, int dim){
    /* function to calc the l2 distance between two vectors*/
    int i;
    double sum;
    sum=0;
    for (i=0;i<dim;i++){
            sum+=pow(point1[i]-point2[i],2);
    }
    return sqrt(sum);
} 

double calcWeight (double* point1, double* point2, int dim){
    double res;
    res=exp((-1)*l2norm(point1,point2, dim)/2);
    return res;
}
/* every line is a vector*/
double** adjacencyMatrix(double*** setOfPoints, int n, int d){
    /* func that creates the adjacency Matrix from
     the raw data*/
    double **mat;
    int i,j,k;
    mat=(double**)malloc((n)*sizeof(double*));
    assert(mat!=NULL&& "An Error Has Occured");
    for (i=0;i<n;i++){
        mat[i]=(double*)malloc((n)*sizeof(double));
        assert(mat[i]!=NULL&& "An Error Has Occured");
    }
    for (j=0;j<n;j++){
        for (k=0;k<n;k++){
                    if (j==k){
                        mat[j][k]=0;
                    }
                    else{
                    mat[j][k]=calcWeight((*setOfPoints)[j],(*setOfPoints)[k],d);
                    }  
        }
    }
    return mat;
}

double calc_index_mat( double** m_one,int i, double** m_two,int j, int row, int col){
    /*
    help func in matrix mul
    calc the valuue A[i][j] in matrix A
    */
    double* row_vec;
    double* col_vec;
    int k;
    double sum;
    sum=0;
    row_vec= (double*) malloc(row* sizeof(double));
    assert(row_vec!= NULL&& "An Error Has Occured");
    col_vec= (double*) malloc(col* sizeof(double));
    assert(col_vec!=NULL&& "An Error Has Occured");
    for (k = 0; k < row; k++)
    {
       row_vec[k] = m_one[i][k];
       col_vec[k]= m_two[k][j];
    }
    /*
    Cross product
    */
    for (k = 0; k < row; k++)
    {
        sum= sum + (row_vec[k]*col_vec[k]);
    }
    return (sum); 
}


double*** multiply(double *** mat1, double*** mat2, int N)
{
    /*matrix mul func*/
    int i, j, k;
    double***res;
    res= (double***) malloc (1*sizeof(double**) );
    assert(res!= NULL&& "An Error Has Occured");
    res[0]= (double**) malloc (N *sizeof(double*) );
    assert (res[0]!=NULL&& "An Error Has Occured");
    for (i=0; i<N;i++){
        res[0][i]= (double*) calloc(N, sizeof(double));
        assert(res[0][i]!= NULL&& "An Error Has Occured");
    }
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            res[0][i][j] = 0;
            for (k = 0; k < N; k++)
                res[0][i][j] += mat1[0][i][k] * mat2[0][k][j];
        }
    }
    return res;
}


double*** transpose(double***A, int row){
    /* transpose squre matrix*/
    double*** output;
    int i;
    int j;
    output= (double***) malloc (1*sizeof(double**) );
    assert(output!= NULL&& "An Error Has Occured");
    output[0]= (double**) malloc (row *sizeof(double*) );
    assert (output[0]!=NULL&& "An Error Has Occured");
    for (i=0; i<row;i++){
        output[0][i]= (double*) calloc(row, sizeof(double));
        assert(output[0][i]!= NULL&& "An Error Has Occured");
    }
    for (i=0;i<row;i++){
        for (j=0;j<row;j++){
            output[0][i][j]=A[0][j][i];
        }
    }
    return output;
}


double** copy_mat (double** matrix, int row, int col){
    /*
    take matrix A and clone this matrix to matrix B ( output )
    */
    double** output;
    int i;
    int j;
    output= (double**) malloc (row*sizeof(double*));
    assert( output !=NULL&& "An Error Has Occured");
    for ( i = 0; i < row; i++)
    {
        output[i]= (double*) malloc(col*sizeof(double));
        assert( output[i] !=NULL&& "An Error Has Occured");
    }
    for ( i = 0; i<row ; i++)
    {
        for ( j = 0; j < col; j++)
        {
            output[i][j]= matrix[i][j];
        }
        
    }
    return (output);
}


double*** get_matrix_d (double*** matrix, int row, int col, int isDDG){
    /*
    calc matrix D^-0.5 (sub section 2)
    */
    double** m;
    double*** output;
    double ** printedMat;
    int i;
    int j;
    double summ;
    m=*matrix;
    /*
    creating the output matrix
    */
    output= (double***) malloc (1*sizeof(double**) );
    assert(output!= NULL&& "An Error Has Occured");
    output[0]= (double**) malloc (row *sizeof(double*) );
    printedMat= (double**) malloc (row *sizeof(double*) );
    assert (output[0]!=NULL&&"An Error Has Occured");
    assert (printedMat!=NULL&&"An Error Has Occured");
    for (i=0; i<row;i++){
        output[0][i]= (double*) calloc(col, sizeof(double));
        printedMat[i]= (double*) calloc(col, sizeof(double));
        assert(output[0][i]!= NULL&& "An Error Has Occured");
        assert(printedMat[i]!= NULL&& "An Error Has Occured");
    }
    /*
    calc matrix D^0.5
    */
    for (i=0; i<row;i++){
        summ=0;
        for (j = 0; j < col; j++)
        {
            summ=summ+m[i][j];
        }
        output[0][i][i]=1/(pow(summ,0.5));
        printedMat[i][i]=summ;
    }
    if (isDDG){
        printMat(printedMat,row,col);
    }
    free(printedMat);
    return(output);
}


double*** get_matrix_L (double*** d, double*** w, int row, int col){
    /*
    calc matrix L-norm (section 2)
    */
    double** tmp;
    double** tmp_d;
    int i;
    int j;
    tmp= copy_mat (*d, row,  col);
    tmp_d= copy_mat(*d, row, col);
    /*
    calc (D^-0.5)*W
    */
    for ( i = 0; i < row; i++)
    {
       for (j = 0; j < col; j++)
       {
           tmp[i][j]= calc_index_mat(*d,i, *w,j, row, col);
       }  
    }
    /*
    calc ((D^-0.5)*W)*D^-0.5
    */
    for ( i = 0; i < row; i++)
    {
       for (j = 0; j < col; j++)
       {
           (*d)[i][j]= calc_index_mat(tmp,i, tmp_d,j, row, col);
       }  
    }
    /*
    calc I- (((D^-0.5)*W)*D^-0.5)
    */
    for (i=0;i<row;i++){
        for (j=0;j<col;j++){
            if (i==j){
                (*d)[i][j]= 1-((*d)[i][j]);
            }
            else{
                (*d)[i][j]= (-1)*((*d)[i][j]);
            }
        }
    }
    return(d);
    
}


int* largestOffDiagElem(double***mat, int row, int col){
    /* calc the maximum differance between  eigen values*/
    double maxi;
    int* maxIndices;
    int i;
    int j;
    maxi=0;
    maxIndices=(int*)malloc(2*sizeof(int));
    assert(maxIndices!=NULL&& "An Error Has Occured");
    maxIndices[0]=-1;
    maxIndices[1]=-1;
    for (i=0;i<row;i++){
        for (j=i+1;j<col;j++){
            if (fabs((*mat)[i][j])>maxi){
                maxIndices[0]=i;
                maxIndices[1]=j;
                maxi=fabs((*mat)[i][j]);
            }
        }
    }
    return maxIndices;
}


double*** unitMatrix (int row, int col){
    /* create unit matrix with dimensions row*col*/
    double*** output;
    int i;
    output= (double***) malloc (1*sizeof(double**) );
    assert(output!= NULL&& "An Error Has Occured");
    output[0]= (double**) malloc (row *sizeof(double*) );
    assert (output[0]!=NULL&&"An Error Has Occured");
    for (i=0; i<row;i++){
        output[0][i]= (double*) calloc(col, sizeof(double));
        assert(output[0][i]!= NULL&& "An Error Has Occured");
        output[0][i][i]=1;
    }
    return output;
}


/*
Rotem's code- section 2
*/


double frobenius(double***mat, int dim){
    /* calc frobenius for specific matrix*/
    int i;
    int j;
    double sum;
    sum=0;
    for(i=0;i<dim;i++){
        for(j=0;j<dim;j++){
            sum=sum+mat[0][i][j]*mat[0][i][j];
        }
    }
    return sum;
}


double diagSum(double***mat, int dim){
    /*calc the sum of squres of the main cross in a matrix*/
    int i;
    double sum;
    sum=0;
    for(i=0;i<dim;i++){
        sum=sum+mat[0][i][i]*mat[0][i][i];
        }
    return sum;
}

double* eigenVals (double*** mat, int dim){
    /* get eigen valus from a matrix*/
    double*output;
    int i;
    output=(double*)malloc(dim*sizeof(double));
    assert(output!=NULL&& "An Error Has Occured");
    for (i=0;i<dim;i++){
        output[i]=mat[0][i][i];
    }
    return output;
}


void BubbleSort(double a[], int indices[], int array_size)
{
    /* sorting array*/
    int i, j;
    double temp;
    int temp2;
    for (i = 0; i < (array_size - 1); ++i)
    {
        for (j = 0; j < array_size - 1 - i; ++j )
        {
            if (a[j] > a[j+1])
            {
                temp = a[j+1];
                temp2=indices[j+1];
                a[j+1] = a[j];
                indices[j+1]=indices[j];
                a[j] = temp;
                indices[j]=temp2;
            }
        }
    }
}


double* calcDifs(double* arr, int dim){
    /*calc the distance between eigen values*/
    int i;
    double* output;
    output=(double*)malloc((dim-1)*sizeof(double));
    assert(output!=NULL&& "An Error Has Occured");
    for (i=1;i<dim;i++){
        output[i-1]=arr[i]-arr[i-1];
    }
    return output;
}


int getK(double* difs, int dim){
    /* calc the size of k*/
    double maxi;
    int maxIndex;
    int i;
    maxi=0;
    maxIndex=0;
    for (i=0;i<dim;i++){
        if (difs[i]>maxi){
            maxi=difs[i];
            maxIndex=i;
        }
    }
    return maxIndex+1;
}


void putColumn(double*** targetMat,double*** sourceMat,int targetInd,int sourceInd, int dim){
    /*copy only the first k eigen vectors*/
    int j;
    for (j=0;j<dim;j++){
        targetMat[0][j][targetInd]=sourceMat[0][j][sourceInd];
    }
    return;
}


double*** getKeigenVectors(double*** mat,int* indices,int k, int dim){
    /*get only the first k eigen vectors*/
    double*** output;
    int i;
    int getIndex;
    output=(double***)malloc(sizeof(double**));
    assert(output!=NULL&& "An Error Has Occured");
    output[0]=(double**)malloc(dim*sizeof(double*));
    assert(output[0]!=NULL&& "An Error Has Occured");
    for (i=0;i<dim;i++){
        output[0][i]=(double*)malloc(k*sizeof(double));
        assert(output[0][i]!=NULL&& "An Error Has Occured");
    }
    for (i=0;i<k;i++){
        getIndex=indices[i];
        putColumn(output,mat,i,getIndex,dim);
    }
    return output;
}


double calc_mat_i (double** matrix, int k, int i ){
    /*calc root of( sum of squres of a row in a matrix)*/
    double tmp;
    double output;
    int r;
    output= 0;
    for (r = 0; r < k; r++)
    {
        tmp= matrix[i][r];
        output= output+(tmp*tmp);
    }
    output= pow(output,0.5);
    return(output);

}



double*** get_matrix_T(double*** matrix, int row, int k){
    /*
    calc T matrix (section 5)
    */
    int i;
    int j;
    double tmp;
    double*** output;
    output= (double***) malloc (1*sizeof(double**) );
    assert(output!= NULL&& "An Error Has Occured");
    output[0]= (double**) malloc (row *sizeof(double*) );
    assert (output[0]!=NULL&& "An Error Has Occured");
    for (i=0; i<row;i++){
        output[0][i]= (double*) calloc(k, sizeof(double));
        assert(output[0][i]!= NULL&& "An Error Has Occured");
    }
    for ( i = 0; i < row; i++)
    {
        for (j = 0; j < k; j++)
        {
            tmp= calc_mat_i (*matrix,  k,  i );
            if(tmp==0){
                tmp=1;
            }
            (*output)[i][j]= (*matrix)[i][j]/tmp;
        }
        
    }
    return (output);
}
void printMat (double** matrix, int rows, int columns){
    int i,j;
    for (i=0;i<rows;i++){
        for(j=0;j<columns;j++){
            if (fabs(matrix[i][j])<0.00005){
                printf("0.0000");}
            else{
                printf("%.4f",matrix[i][j]);}
            if(j!=columns-1){
                printf(",");
            }
        }
        printf("\n");}
}


double*** geo_c(double** a, int dimRow, int dimCol, int goal,int k){
    /* main funct
     calcs the output according to the goal
     */
    int DIM;
    double***L;
    double*** d;
    double **aa;
    double **b;
    int i;
    int* maxIndices;
    double***V;
    double***P;
    double theta;
    int signTheta;
    double t;
    double c;
    double s;
    double***Ptranspose;
    double***temp;
    double***Ltag;
    double* eigenValues;
    double offL;
    double offLtag;
    int ind;
    double* eigenDifs;
    int* eigenValuesIndices;
    double*** eigenVectors;
    double***T;
    double***retForSpk;
    double ***emptyMat;
    ind=0;
    aa= a;
    b=adjacencyMatrix(&aa,dimRow,dimCol);

    DIM=dimRow;
    emptyMat= (double***) malloc(1*sizeof(double**));
    assert(emptyMat!=NULL&& "An Error Has Occured");
    emptyMat[0]= (double**) malloc(1*sizeof(double*));
    assert(emptyMat[0]!=NULL&& "An Error Has Occured");
    emptyMat[0][0]= (double*) calloc(1,sizeof(double));
    assert(emptyMat[0][0]!=NULL&& "An Error Has Occured");
    if (k>dimRow){
        printf("Invalid input!\n");
        return emptyMat;
    }
    if (goal==2){ /*wam*/
        printMat(b,DIM,DIM);
        free(aa);
        free(b);
        return emptyMat;
    }
    d= get_matrix_d (&b, DIM, DIM,goal==3);
    if(goal==3){ /*DDG*/
        free(aa);
        free(b);
        free(d);
        return emptyMat;
    }
    L= get_matrix_L(d,&b, DIM,DIM);
    if(goal==4){ /* lnorm */
        printMat(L[0],DIM,DIM);
        free(aa);
        free(b);
        free(L);
        return emptyMat;
    }
    maxIndices=largestOffDiagElem(L,DIM,DIM);
    V=unitMatrix(DIM,DIM);
    /* JACOBI algo*/
    while ((maxIndices[0]!=-1&&maxIndices[1]!=-1)&&(ind<100)){
        P=unitMatrix(DIM,DIM);
        theta=(L[0][maxIndices[1]][maxIndices[1]]-L[0][maxIndices[0]][maxIndices[0]])/(2*L[0][maxIndices[0]][maxIndices[1]]);
        signTheta=theta>=0? 1:-1;
        theta=fabs(theta);
        t=signTheta/(theta+sqrt(theta*theta+1));
        c=1/sqrt(t*t+1);
        s=c*t;
        P[0][maxIndices[0]][maxIndices[0]]=c;
        P[0][maxIndices[1]][maxIndices[1]]=c;
        P[0][maxIndices[0]][maxIndices[1]]=s;
        P[0][maxIndices[1]][maxIndices[0]]=s*(-1);
        temp=multiply(V,P,DIM);
        V[0]=copy_mat(temp[0],DIM,DIM);
        Ptranspose=transpose(P,DIM);
        Ltag=multiply(multiply(Ptranspose,L,DIM),P,DIM);
        offL=frobenius(L,DIM)-diagSum(L,DIM);
        offLtag=frobenius(Ltag,DIM)-diagSum(Ltag,DIM);
        if(fabs(offL-offLtag)<=EPSILON){
            L[0]=copy_mat(Ltag[0],DIM,DIM);
            break;
        }
        L[0]=copy_mat(Ltag[0],DIM,DIM);
        maxIndices=largestOffDiagElem(L,DIM,DIM);
        ind++;
    }

    if (goal==5){ /*Jacobi*/
        for(i=0;i<DIM;i++){
                if (fabs(L[0][i][i])<0.00005){
                         printf("0.0000");}
                    else{
                        printf("%.4f",L[0][i][i]);}
                    if(i!=DIM-1){
                        printf(",");
                    }
        }
            printf("\n");
        printMat(V[0],DIM,DIM);
        free(aa);
        free(b);
        free(L);
        free(V);
        return emptyMat;
    }
    /* finding k */
    eigenValuesIndices=(int*)malloc(DIM*sizeof(int));
    assert(eigenValuesIndices!=NULL&& "An Error Has Occured");
    for(i=0;i<DIM;i++){
        eigenValuesIndices[i]=i;
    }
    eigenValues=eigenVals(L,DIM);
    BubbleSort(eigenValues,eigenValuesIndices,DIM);
    eigenDifs=calcDifs(eigenValues,DIM);
    if (k==0){
        k=getK(eigenDifs,(int)floor(DIM/2));
    }
    eigenVectors=getKeigenVectors(V,eigenValuesIndices,k,DIM); /*U*/
    T=get_matrix_T(eigenVectors,DIM,k);
    retForSpk= (double***) malloc (2*sizeof(double**) );
    assert(retForSpk!= NULL&& "An Error Has Occured");
    retForSpk[0]=T[0];
    retForSpk[1]= (double**) malloc (1 *sizeof(double*) );
    assert (retForSpk[1]!=NULL&& "An Error Has Occured");
    retForSpk[1][0]= (double*) malloc(1*sizeof(double));
    assert(retForSpk[0][0]!= NULL&& "An Error Has Occured");
    retForSpk[1][0][0]=k*1.0;
    free(aa);
    free(b);
    free(L);
    free(V);
    free(eigenValues);
    free(eigenVectors);
    free(eigenDifs);
    free(eigenValuesIndices);
    free(emptyMat);
    return retForSpk;
}

/* KMEANS - FIRST PROJECT*/ 


void updateCheck (int n,double vec [], double mean []){
    int i;
    for ( i = 0; i < n; i++)
    {
        vec[i]=mean[i];
    }   
}

int toStop (int n,double vec [], double mean []){
    int i;
    for ( i = 0; i < n; i++){
        if(vec[i]!=mean[i]){
            return (0);
        }
    }
    return (1);
}


int main (int argc, char *argv[]) {
    /*main func for C implemntion*/
    int i,j,cc,row,columns,b1,k1,k2,lines,iter,place,d,dim,i1,j1,K, max_iter,goal;
    FILE *fp;
    double distance,curDistance,cur;
    char * token;
    double *vecOfMeansOld;
    double *vecOfMeansNew;
    double *vector;
    char *singleLine;
    int *numberOfVecPerGroup;
    double *** resultSpk;
    double ** T;
    double integerCheckK, fractionalCheckK;
    /*get raw data*/
    char line[MAX_LINE_LENGTH];
    max_iter=300;
    fractionalCheckK=modf(atof(argv[1]),&integerCheckK);
    if ((integerCheckK<0) || (argc!=4) || (fractionalCheckK>0) ){
        printf( "Invalid Input!\n");
        return 0;
    }
    K=atoi(argv[1]);
    lines=0;
    i=0;
    dim=1;
    fp = fopen(argv[3],"r");
    assert(fp!=NULL && "An Error Has Occured");
    while (fgets(line, MAX_LINE_LENGTH, fp)){
        if (dim==1){
            while (line[i]!='\0'){
                if (line[i]==','){
    
                    dim++;
                }
                i++;
            }
        }
            lines++;
        }
    if (K>lines){
        printf( "Invalid Input!\n");
        return 0;
    }
    fclose(fp);
    fp = fopen(argv[3],"r");
    vector = (double *)malloc(sizeof(double)*(lines+1)*(dim+1));
    assert(vector!=NULL && "An Error Has Occured");
    i1=0;
    j1=0;
    singleLine=(char *)malloc(sizeof(char)*dim*64);
    assert(singleLine!=NULL && "An Error Has Occured");
    while(fgets(singleLine, dim*64, fp)){
        assert(singleLine!=NULL && "An Error Has Occured");
        token = strtok(singleLine, ",");
        while( token != NULL ) {
            cur=atof(token);
            vector[dim*i1+j1]=cur;
            token = strtok(NULL, ",");
            j1++;
        }
        i1++;
        j1=0;
    }
    free(singleLine);
    goal=goalEnum(argv[2]);
    if (goal==0){
        return 0;
    }
    if (goal!=1){
        geo_c(oneDimToTwoDim(vector,lines,dim),lines,dim,goal,K);
        return 0;
    }
    numberOfVecPerGroup=(int *)malloc(sizeof(int)*K);
    assert(numberOfVecPerGroup!=NULL && "An Error Has Occured");
        for ( iter = 0; iter < K; iter++)
    {
        numberOfVecPerGroup[iter]=0;

    }
    /*k means algo on modified data*/
    iter=0;
    resultSpk=geo_c(oneDimToTwoDim(vector,lines,dim),lines,dim,goal,K);
    T=resultSpk[0];
    K=resultSpk[1][0][0];
    vector=twoDimToOneDim(T,lines,K);
    dim=K;
    d=K;
    vecOfMeansOld = (double *)malloc(sizeof(double)*(K+1)*(dim+1));
    assert(vecOfMeansOld!=NULL&& "An Error Has Occured");
    vecOfMeansNew = (double *)malloc(sizeof(double)*(K+1)*(dim+1));
    assert(vecOfMeansNew!=NULL&& "An Error Has Occured");
    for ( i = 0; i < K; i++)
    {
        for ( cc = 0; cc < dim; cc++)
        {
            vecOfMeansOld[i*dim+cc]= vector[i*dim+cc];
            
        }
        
    }
    d=dim;
    /*update centers*/
    while(iter<max_iter){
            
            for ( i = 0; i < lines ; i++)
            {
                distance=0;
                place=0;
                
                for (  j = 0; j < K; j++)
                {
                    curDistance= distanceFromMean(i,j, d, vector, vecOfMeansOld );
                    if (j==0){
                            distance=curDistance;
                    }
                    if(curDistance< distance){
                        distance=curDistance;
                        place=j;
                    } 
                }
                updateMean (i,place,d, vector, vecOfMeansNew , numberOfVecPerGroup[place]);
                numberOfVecPerGroup[place]=numberOfVecPerGroup[place]+1;
            }

            b1=1;
            for (k1=0; k1<K; k1++)
            {
                for (k2 = 0; k2 < dim; k2++)
                {
                    if(vecOfMeansOld[k1*dim+k2]!=vecOfMeansNew[k1*dim+k2]){
                        b1=0;
                    }
                }    
            }
            if (b1==1){
                break;
            }
            for (k1=0; k1<K; k1++)
            {
                for (k2 = 0; k2 < dim; k2++)
                {
                    vecOfMeansOld[k1*dim+k2]=vecOfMeansNew[k1*dim+k2];
                    vecOfMeansNew[k1*dim+k2]=0;
                }
                
            }
            for ( i = 0; i < K; i++)
            {
                numberOfVecPerGroup[i]=0;
            }
            iter++;
        }
        for (row=0; row<K; row++)
        {
            for(columns=0; columns<dim; columns++)
            {
            printf("%.4f", vecOfMeansOld[row*dim+columns]);
            if (columns!=(dim-1)){
                printf(",");
            }
            }
            printf("\n");
        }
        free(vector);
        free(vecOfMeansOld);
        free(numberOfVecPerGroup);
        free(vecOfMeansNew);
        fclose(fp);
        return (0);
}