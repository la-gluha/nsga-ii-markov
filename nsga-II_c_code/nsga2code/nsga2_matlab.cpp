//
// Created by 啦啦啦 on 2023/3/8.
//
#define square(x) ((x)*(x))
#define maxpop   500  /*Max population */
#define maxchrom 200  /*Max chromosome length*/
#define maxvar    20  /*Max no. of variables*/
#define maxfun    10  /*Max no. of functions */
#define maxcons   20  /*Max no. of Constraints*/

int gener,       /*No of generations*/
nvar,nchrom,          /*No of variables*/
ncons,         /*No of Constraints*/
vlen[maxvar],  /*Array to store no of bits for each variable*/
nmut,          /* No of Mutations */
ncross,        /*No of crossovers*/
ans;
float seed,      /*Random Seed*/
pcross,        /*Cross-over Probability*/
pmut_b, pmut_r,          /*Mutation Probability*/
lim_b[maxvar][2], lim_r[maxvar][2];/*Limits of variable in array*/
float di,        /*Distribution Index for the Cross-over*/
dim,           /*Distribution Index for the Mutation*/
delta_fit,     /* variables required forfitness for fitness sharing */
min_fit,
        front_ratio;
int optype,      /*Cross-over type*/
nfunc,         /*No of functions*/
sharespace;    /*Sharing space (either parameter or fitness)*/

double coef[maxvar]; /*Variable used for decoding*/

static int popsize,  /*Population Size*/
chrom;             /*Chromosome size*/

typedef struct       /*individual properties*/
{
    int genes[maxchrom], /*bianry chromosome*/
    rank,              /*Rank of the individual*/
    flag;              /*Flag for ranking*/
    float xreal[maxvar], /*list of real variables*/
    xbin[maxvar];      /*list of decoded value of the chromosome */
    float fitness[maxfun],/*Fitness values */
    constr[maxcons],     /*Constraints values*/
    cub_len,             /*crowding distance of the individual*/
    error;              /* overall constraint violation for the individual*/
}individual;        /*Structure defining individual*/


typedef struct
{
    int maxrank;            /*Maximum rank present in the population*/
    float rankrat[maxpop];  /*Rank Ratio*/
    int rankno[maxpop];     /*Individual at different ranks*/
    individual ind[maxpop], /*Different Individuals*/
    *ind_ptr;
}population ;             /*Popuation Structure*/

#include "random.h"       /*Random Number Generator*/

#include "input.h"        /*File Takes Input from user*/

#include "realinit.h"     /*Random Initialization of the populaiton*/

#include "init.h"         /*Random Initialization of the population*/

#include "decode.h"       /*File decoding the binary dtrings*/

#include "ranking.h"      /*File Creating the Pareto Fronts*/

#include "rancon.h"       /*File Creating the Pareto Fronts when
			    Constraints are specified*/

#include "func-con.h"     /*File Having the Function*/

#include "select.h"       /*File for Tournament Selection*/


#include "crossover.h"    /*Binary Cross-over*/

#include "uniformxr.h"    /*Uniform Cross-over*/

#include "realcross2.h"   /*Real Cross-over*/

#include "mut.h"          /*Binary Mutation*/

#include "realmut1.h"     /*Real Mutation*/

#include "keepaliven.h"   /*File For Elitism and Sharing Scheme*/

#include "report.h"       /*Printing the report*/
#include "markov.h"
#include "mex.h"


population oldpop,
        newpop,
        matepop,
        *old_pop_ptr,
        *new_pop_ptr,
        *mate_pop_ptr;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    int i,j,l,f,maxrank1;
    float *ptr,tot;
    FILE
            *rep_ptr,
            *gen_ptr,
            *rep2_ptr,
            *end_ptr,
            *g_var,
            *lastit;
    /*File Pointers*/
    rep_ptr = fopen("output.out","w");
    gen_ptr =fopen("all_fitness.out","w");
    rep2_ptr = fopen("ranks.out","w");
    end_ptr = fopen("final_fitness.out","w");
    g_var = fopen("final_var.out","w");
    lastit = fopen("plot.out","w");
    /*Opening the files*/

    old_pop_ptr = &(oldpop);

    nmut = 0;
    ncross = 0;

    /*Get the input from the file input.h*/
//    input(rep_ptr);

    //输入的参数，这里直接写死用于对马尔可夫网络进行debug，算法正常后再使用input确定参数
    nvar = 0, nchrom = 4, nfunc = 2, ncons = 0, popsize = 100, gener = 100;
    pcross = 0.9, optype = 1;

    for (i = 0; i < nchrom; i++)
    {
        vlen[i] = 15;
        chrom += vlen[i];

        lim_b[i][0] = 0.0, lim_b[i][1] = 10.0;
    }

    pmut_b = 0.015, seed = 0.5;


    //马尔可夫使用到的参数
    int solLength = 15; //使用时需要对该变量进行赋值，该变量指个体数量，gens的长度对应的是目标的个数
    /* pv、alpha是用来保存马尔可夫网络相关参数的全局数组，通过它们可以使得每次变异之间产生联系 */
    double pv[solLength];
    double alpha[solLength+1];

    fprintf(rep_ptr,"Results in a file\n");
    fprintf(end_ptr,"# Last generation population (Feasible and non-dominated)\n");
    fprintf(end_ptr,"# Fitness_vector (first %d)  Constraint_violation (next %d)  Overall_penalty\n",nfunc,ncons);
    fprintf(g_var,"#Feasible Variable_vectors for non-dominated solutions at last generation\n");
    fprintf(g_var,"# Real (first %d)  Binary (next %d)\n",nvar,nchrom);
    fprintf(lastit,"# Feasible and Non-dominated Objective Vector\n");

    /*Initialize the random no generator*/
    warmup_random(seed);

    /*Binary Initializaton*/
    if (nchrom > 0)
        init(old_pop_ptr);
    if (nvar > 0)
        realinit(old_pop_ptr);

    old_pop_ptr = &(oldpop);

    // decode binary strings
    decode(old_pop_ptr);

    old_pop_ptr = &(oldpop);
    new_pop_ptr = &(newpop);

    for(j = 0;j < popsize;j++)
    {
        /*Initializing the Rank array having different individuals
      at a particular  rank to zero*/
        old_pop_ptr->rankno[j] = 0;
        new_pop_ptr->rankno[j] = 0;
    }

    old_pop_ptr = &(oldpop);

    func(old_pop_ptr);
    /*Function Calculaiton*/

    fprintf(rep_ptr,"----------------------------------------------------\n");
    fprintf(rep_ptr,"Statistics at Generation 0 ->\n");
    fprintf(rep_ptr,"--------------------------------------------------\n");

    /********************************************************************/
    /*----------------------GENERATION STARTS HERE----------------------*/
    for (i = 0;i < gener;i++)
    {
        printf("Generation = %d\n",i+1);
        old_pop_ptr = &(oldpop);
        mate_pop_ptr = &(matepop);
        fprintf(rep_ptr,"Population at generation no. -->%d\n",i+1);
        fprintf(gen_ptr,"#Generation No. -->%d\n",i+1);
        fprintf(gen_ptr,"#Variable_vector  Fitness_vector Constraint_violation Overall_penalty\n");

        /*--------SELECT----------------*/
        //选择
        // old_pop_ptr 父代P
        // mate_pop_ptr 从old_pop_ptr中选择出来的种群D，在这上面进行操作，old_pop_ptr保留
        // new_pop_ptr 根据D以及MRF生成的新种群
        nselect(old_pop_ptr ,mate_pop_ptr );

        new_pop_ptr = &(newpop);
        mate_pop_ptr = &(matepop);

        /*CROSSOVER----------------------------*/
        //交叉
        if (nchrom > 0)
        {

            if(optype == 1)
            {
                crossover(new_pop_ptr ,mate_pop_ptr );
                /*Binary Cross-over*/
            }

            if(optype == 2)
            {
                unicross(new_pop_ptr ,mate_pop_ptr );
                /*Binary Uniform Cross-over*/
            }
        }
        if (nvar > 0)
            realcross(new_pop_ptr ,mate_pop_ptr );
        /*Real Cross-over*/


        // Markov用于在变异中
        /*------MUTATION-------------------*/
        //变异
        new_pop_ptr = &(newpop);
//        mate_pop_ptr = &(matepop);

        if (nchrom > 0)
            mutate(new_pop_ptr );
//            markovMut(new_pop_ptr, i, pv, alpha, solLength);
        /*Binary Mutation */

        if (nvar > 0)
            real_mutate(new_pop_ptr );
        /*Real Mutation*/

        new_pop_ptr = &(newpop);

        /*-------DECODING----------*/
        //译码为二进制串，对二进制串进行操作
        if(nchrom > 0)
            decode(new_pop_ptr );
        /*Decoding for binary strings*/

        // 上面已经得到了新的子代，这里在进行精英策略，选择出新的子代
        // 如果不进行精英策略则是直接把生成的子代作为子代进行后面的操作
        // 函数值评估
        // 这里包括keepalive部分都是精英策略，这个部分是为其做准备，而keepalive部分则是执行精英策略
        /*----------FUNCTION EVALUATION-----------*/
        new_pop_ptr = &(newpop);
        func(new_pop_ptr );

        /*-------------------SELECTION KEEPING FRONTS ALIVE--------------*/
        old_pop_ptr = &(oldpop);
        new_pop_ptr = &(newpop);
        mate_pop_ptr = &(matepop);

        /*Elitism And Sharing Implemented*/
        // 父代中选出子代放到mate_pop_ptr作为新的父代，然后产生new_pop_ptr
        // 可以将mate_pop_ptr认为是父代，而old_pop_ptr则是作为备份
        keepalive(old_pop_ptr ,new_pop_ptr ,mate_pop_ptr,i+1);

        mate_pop_ptr = &(matepop);
        if(nchrom > 0)
            decode(mate_pop_ptr );

        mate_pop_ptr = &(matepop);
        /*------------------REPORT PRINTING--------------------------------*/
        report(i ,old_pop_ptr ,mate_pop_ptr ,rep_ptr ,gen_ptr, lastit );

        /*==================================================================*/

        /*----------------Rank Ratio Calculation------------------------*/
        new_pop_ptr = &(matepop);
        old_pop_ptr = &(oldpop);

        /*Finding the greater maxrank among the two populations*/

        if(old_pop_ptr->maxrank > new_pop_ptr->maxrank)
            maxrank1 = old_pop_ptr->maxrank;
        else
            maxrank1 = new_pop_ptr->maxrank;

        fprintf(rep2_ptr,"--------RANK AT GENERATION %d--------------\n",i+1);
        fprintf(rep2_ptr,"Rank old ranks   new ranks     rankratio\n");

        for(j = 0;j < maxrank1 ; j++)
        {
            /*Sum of the no of individuals at any rank in old population
              and the new populaion*/

            tot = (old_pop_ptr->rankno[j])+ (new_pop_ptr->rankno[j]);

            /*Finding the rank ratio for new population at this rank*/

            new_pop_ptr->rankrat[j] = (new_pop_ptr->rankno[j])/tot;

            /*Printing this rank ratio to a file called ranks.dat*/

            fprintf(rep2_ptr," %d\t  %d\t\t %d\t %f\n",j+1,old_pop_ptr->rankno[j],new_pop_ptr->rankno[j],new_pop_ptr->rankrat[j]);

        }

        fprintf(rep2_ptr,"-----------------Rank Ratio-------------------\n");
        /*==================================================================*/

        /*=======Copying the new population to old population======*/

        old_pop_ptr = &(oldpop);
        new_pop_ptr = &(matepop);

        for(j = 0;j < popsize;j++)
        {
            old_pop_ptr->ind_ptr = &(old_pop_ptr->ind[j]);
            new_pop_ptr->ind_ptr = &(new_pop_ptr->ind[j]);
            if(nchrom > 0)
            {
                /*For Binary GA copying of the chromosome*/

                for(l = 0;l < chrom;l++)
                    old_pop_ptr->ind_ptr->genes[l]=new_pop_ptr->ind_ptr->genes[l];

                for(l = 0;l < nchrom;l++)
                    old_pop_ptr->ind_ptr->xbin[l] = new_pop_ptr->ind_ptr->xbin[l];
            }
            if(nvar > 0)
            {
                /*For Real Coded GA copying of the chromosomes*/
                for(l = 0;l < nvar;l++)
                    old_pop_ptr->ind_ptr->xreal[l] = new_pop_ptr->ind_ptr->xreal[l];
            }

            /*Copying the fitness vector */
            for(l = 0 ; l < nfunc ;l++)
                old_pop_ptr->ind_ptr->fitness[l] = new_pop_ptr->ind_ptr->fitness[l];

            /*Copying the dummy fitness*/
            old_pop_ptr->ind_ptr->cub_len = new_pop_ptr->ind_ptr->cub_len;

            /*Copying the rank of the individuals*/
            old_pop_ptr->ind_ptr->rank = new_pop_ptr->ind_ptr->rank;

            /*Copying the error and constraints of the individual*/

            old_pop_ptr->ind_ptr->error = new_pop_ptr->ind_ptr->error;
            for(l = 0;l < ncons;l++)
            {
                old_pop_ptr->ind_ptr->constr[l] = new_pop_ptr->ind_ptr->constr[l];
            }

            /*Copying the flag of the individuals*/
            old_pop_ptr->ind_ptr->flag = new_pop_ptr->ind_ptr->flag;
        }   // end of j

        maxrank1 = new_pop_ptr->maxrank ;

        /*Copying the array having the record of the individual
      at different ranks */
        for(l = 0;l < popsize;l++)
        {
            old_pop_ptr->rankno[l] = new_pop_ptr->rankno[l];
        }

        /*Copying the maxrank */
        old_pop_ptr->maxrank = new_pop_ptr->maxrank;

        /*Printing the fitness record for last generation in a file last*/
        if(i == gener-1)
        {  // for the last generation
            old_pop_ptr = &(matepop);
            for(f = 0;f < popsize ; f++) // for printing
            {
                old_pop_ptr->ind_ptr = &(old_pop_ptr->ind[f]);

                if ((old_pop_ptr->ind_ptr->error <= 0.0) && (old_pop_ptr->ind_ptr->rank == 1))  // for all feasible solutions and non-dominated solutions
                {
                    for(l = 0;l < nfunc;l++)
                        fprintf(end_ptr,"%f\t",old_pop_ptr->ind_ptr->fitness[l]);
                    for(l = 0;l < ncons;l++)
                    {
                        fprintf(end_ptr,"%f\t",old_pop_ptr->ind_ptr->constr[l]);
                    }
                    if (ncons > 0)
                        fprintf(end_ptr,"%f\t",old_pop_ptr->ind_ptr->error);
                    fprintf(end_ptr,"\n");

                    if (nvar > 0)
                    {
                        for(l = 0;l < nvar ;l++)
                        {
                            fprintf(g_var,"%f\t",old_pop_ptr->ind_ptr->xreal[l]);
                        }
                        fprintf(g_var,"  ");
                    }

                    if(nchrom > 0)
                    {
                        for(l = 0;l < nchrom;l++)
                        {
                            fprintf(g_var,"%f\t",old_pop_ptr->ind_ptr->xbin[l]);
                        }
                    }
                    fprintf(g_var,"\n");
                }  // feasibility check
            } // end of f (printing)

        } // for the last generation
    }  // end of i

    /*                   Generation Loop Ends                                */
    /************************************************************************/

    fprintf(rep_ptr,"NO. OF CROSSOVER = %d\n",ncross);
    fprintf(rep_ptr,"NO. OF MUTATION = %d\n",nmut);
    fprintf(rep_ptr,"------------------------------------------------------------\n");
    fprintf(rep_ptr,"---------------------------------Thanks---------------------\n");
    fprintf(rep_ptr,"-------------------------------------------------------------\n");
    printf("NOW YOU CAN LOOK IN THE FILE OUTPUT2.DAT\n");

    /*Closing the files*/
    fclose(rep_ptr);
    fclose(gen_ptr);
    fclose(rep2_ptr);
    fclose(end_ptr);
    fclose(g_var);
    fclose(lastit);
}
