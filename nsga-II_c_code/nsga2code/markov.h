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
    double tc = 0.55;
    double lr = 1;

    normalaiseFitnessToAllPositive(spop);

    double a[rows][cols];
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols - 1; j++) {
            a[i][j] = toInt(spop->ind[i].genes[j]);
        }
        a[i][cols - 1] = 1;  //for constant C
    }

    double f[rows];
    for (int i = 0; i < rows; i++) {
        float fitness = 0.0;
        for (int j = 0; j < nfunc; ++j) {
            fitness += spop->ind[i].fitness[j];
        }
        f[i] = -logf(fitness);
    }

    solveWithSVD((double **) a, rows, cols, rows, cols, f, alpha);

    //d
    for (int i = 0; i < cols; ++i) {
        pv[i] = 1.0 / (1.0 + exp(alpha[i] * tc * g));
    }

    //pv
//    for (int i = 0; i < solLength; i++) {
//        if (alpha[i] < 0) pv[i] = (pv[i] * (1 - lr)) + lr;
//        else if (alpha[i] > 0) pv[i] = pv[i] * (1 - lr);
//    }

    for (int i = 0; i < popsize; ++i) {
        sample((double *) pv, &spop->ind[i]);
    }
}

void normalaiseFitnessToAllPositive(population *spop) {
    for (int i = 0; i < nfunc; ++i) {
        float diff = 0;
        if (spop->ind[popsize - 1].fitness[i] <= 0 && spop->ind[0].fitness[i] <= 0) {
            diff = spop->ind[0].fitness[i] + spop->ind[popsize - 1].fitness[i];
            for (int j = popsize - 1; j >= 0; j--) {
                spop->ind[j].fitness[i] -= diff;
            }
        } else if (spop->ind[popsize - 1].fitness[i] <= 0 && spop->ind[0].fitness[i] > 0) {
            diff = spop->ind[0].fitness[i] - spop->ind[popsize - 1].fitness[i];
            for (int j = popsize - 1; j >= 0; j--) {
                spop->ind[j].fitness[i] += diff;
            }
        } else if (spop->ind[popsize - 1].fitness[i] > 0 && spop->ind[0].fitness[i] <= 0) {
            diff = spop->ind[popsize - 1].fitness[i] - spop->ind[0].fitness[i];
            for (int j = popsize - 1; j >= 0; j--) {
                spop->ind[j].fitness[i] += diff;
            }
        }
    }
}

void sample(double *pv, individual *ind) {
    //个体中基因的数组的长度就是chrom
    for (int i = 0; i < chrom; ++i) {
        double random = randomperc();
        if (random < pv[i]) {
            ind->genes[i] = 1;
        } else {
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
    return i == 0 ? -1 : 1;
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

