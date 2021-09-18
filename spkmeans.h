#ifndef SPKMEANS_H_
#define SPKMEANS_H_

double*** geo_c(double** a, int dimRow, int dimCol, int goal, int k);
double** mainFuncV2 (int K,  int max_iter, double* arrInitialCentroids, double* arrFinalData, double** rawDataMat, int lines, int dim, int dNew);

#endif
