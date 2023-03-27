//
// Created by 啦啦啦 on 2023/2/16.
//

#include "mathUtil.h"

#ifndef NSGA2CODE_MARKOV_H
#define NSGA2CODE_MARKOV_H

void estimateAndSampleMRF(population *spop, int g, double *pv, double *alpha, int solLength);

void normalaiseFitnessToAllPositive(population *spop);

//void copy(population *to, population *from);
//
int toInt(int i);

void sample(double *pv, individual *ind);

double rand_finite_double();

void deumD(double *pv, double *alpha, double g);

void deumPv(double *pv, double *alpha, double lr, int solLength);

// 该函数的作用是修改pv的值，进行使得变异之间产生相互作用
void estimateAndSampleMRF(population *spop, int g, double *pv, double *alpha, int solLength) {
//    int popSize = sizeof(spop->ind) / sizeof(individual);  //chrom是基因的长度
    int rows = solLength, cols = solLength;
    double tc = 0.5;
    double lr = 1;
//    population *spopOrg;
//
//    copy(spopOrg, spop);

//    for(int ti = 0; ti < popsize; ti++){
//        for(int tj = 0; tj < nfunc; tj++){
//            printf("markov_f[%d] = %f\n", tj, spop->ind[ti].fitness[tj]);
//        }
//    }

    normalaiseFitnessToAllPositive(spop);

//    for (int ti = 0; ti < popsize; ti++) {
//        for (int tj = 0; tj < nfunc; tj++) {
//            printf("markov_f[%d] = %f\n", tj, spop->ind[ti].fitness[tj]);
//        }
//    }

    double a[rows][cols];
//    printf("rows = %d, cols = %d\n", rows, cols);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols - 1; j++) {
//            printf("num = %f\n", toInt(spop->ind[i].genes[j]));
            a[i][j] = toInt(spop->ind[i].genes[j]);
//              a[i][j] = spop->ind[i].genes[j];
        }
        a[i][cols - 1] = 1;  //for constant C
    }

//    for (int ti = 0; ti < rows; ti++){
//        for (int tj = 0; tj < cols; tj++){
//            printf("a[%d][%d] = %f\n", ti, tj, *((double *) a + cols * ti + tj));
//        }
//    }

//    for (int i = 0; i < rows; i++) {
//        for (int j = 0; j < cols; j++) {
//            printf("a[%d][%d] = %f\n", i, j, a[i][j]);
//        }
//    }

    double f[rows];
    for (int i = 0; i < rows; i++) {
        f[i] = -logf(spop->ind[i].fitness[0]);
//        printf("fitness = %f\n", spop->ind[i].fitness[0]);
//        printf("f[%d] = %f\n", i, f[i]);
    }

    solveWithSVD((double **) a, rows, cols, rows, cols, f, alpha);

    for(int fi = 0; fi < rows; fi++){
        printf("f[%d] = %f\n", fi, f[fi]);
    }

//    printf("solveWithSVD_flag");

    for (int i = 0; i < cols; ++i) {
//        printf("before, pv[%d] = %f\n", i, pv[i]);
        pv[i] = 1.0 / (1.0 + exp(alpha[i] * tc * g));
//        printf("alpha[%d] = %f\n", i, alpha[i]);
//        printf("x = %f\n", alpha[i]*tc*g);
//        printf("exp = %f\n", exp(alpha[i]*tc*g));
//        printf("expression = %f\n", 1.0 / (1.0 + exp(alpha[i]*tc*g)));
//        printf("after, pv[%d] = %f\n", i, pv[i]);
    }
//    deumD((double *) pv, (double *) alpha, g);
//    deumPv((double *) pv, (double *) alpha, lr, solLength);

    //popSize == popsize ?
    for (int i = 0; i < popsize; ++i) {
        sample((double *) pv, &spop->ind[i]);
    }
}

void normalaiseFitnessToAllPositive(population *spop) {
//    int popSize = sizeof(spop->ind) / sizeof(individual);
//    int fitnessSize = sizeof (spop->ind->fitness) / sizeof (float);
//    for (int i = 0; i < fitnessSize; ++i){
    float diff = 0;
    if (spop->ind[popsize - 1].fitness[0] <= 0 && spop->ind[0].fitness[0] <= 0) {
        diff = spop->ind[0].fitness[0] + spop->ind[popsize - 1].fitness[0];
        for (int i = popsize - 1; i >= 0; i--) {
            spop->ind[i].fitness[0] -= diff;
        }
    } else if (spop->ind[popsize - 1].fitness[0] <= 0 && spop->ind[0].fitness[0] > 0) {
        diff = spop->ind[0].fitness[0] - spop->ind[popsize - 1].fitness[0];
        for (int i = popsize - 1; i >= 0; i--) {
            spop->ind[i].fitness[0] += diff;
        }
    } else if (spop->ind[popsize - 1].fitness[0] > 0 && spop->ind[0].fitness[0] <= 0) {
        diff = spop->ind[popsize - 1].fitness[0] - spop->ind[0].fitness[0];
        for (int i = popsize - 1; i >= 0; i--) {
            spop->ind[i].fitness[0] += diff;
        }
    }
//    }
}

void deumD(double *pv, double *alpha, double g) {
    int cols = 51;
    double tc = 1;
    for (int i = 0; i < cols; ++i) {
        pv[i] = 1.0 / (1.0 + exp(alpha[i] * tc * g));
    }
}

void deumPv(double *pv, double *alpha, double lr, int solLength) {
    for (int i = 0; i < solLength; i++) {
        if (alpha[i] < 0) pv[i] = (pv[i] * (1 - lr)) + lr;
        else if (alpha[i] > 0) pv[i] = pv[i] * (1 - lr);
    }
}

void sample(double *pv, individual *ind) {
    // int length = sizeof(ind->genes) / sizeof(int);  //个体中记录基因的数组的长度，其实就是chrom
    for (int i = 0; i < chrom; ++i) {
        double random = randomperc();
//        printf("pv[%d] = %f\n", i, pv[i]);
        if (random < pv[i]) {
            ind->genes[i] = 1;
        } else {
//            ind->genes[i] = -1;
            ind->genes[i] = 0;
        }
        nmut++;
    }
}

//void copy(population *to, population *from) {
//    memcpy(to, from, sizeof(population));
//    to->ind_ptr = malloc(sizeof(individual));;
//    strcpy(to->ind_ptr, from->ind_ptr);
//}

int toInt(int i) {
    return i == 0 ? -1 : i;
}

double rand_finite_double() {
    union {
        double d;
        unsigned char uc[sizeof(double)];
    } u;
    do {
        for (unsigned i = 0; i < sizeof u.uc; i++) {
            u.uc[i] = (unsigned char) rand();
        }
    } while (!isfinite(u.d));
    return u.d;
}

#endif //NSGA2CODE_MARKOV_H

