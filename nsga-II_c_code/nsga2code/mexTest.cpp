//
// Created by 啦啦啦 on 2023/3/8.
//
#include "mex.h";
double add(double x, double y){
    return x+y;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    double a, b;
    a = mxGetScalar(prhs[0]);
    b = mxGetScalar(prhs[1]);
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    double *c;
    c = mxGetPr(plhs[0]);
    *c = add(a, b);
}

