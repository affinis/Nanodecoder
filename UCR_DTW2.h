/***********************************************************************/
/************************* DISCLAIMER **********************************/
/***********************************************************************/
/** This UCR Suite software is copyright protected � 2012 by          **/
/** Thanawin Rakthanmanon, Bilson Campana, Abdullah Mueen,            **/
/** Gustavo Batista and Eamonn Keogh.                                 **/
/**                                                                   **/
/** Unless stated otherwise, all software is provided free of charge. **/
/** As well, all software is provided on an "as is" basis without     **/
/** warranty of any kind, express or implied. Under no circumstances  **/
/** and under no legal theory, whether in tort, contract,or otherwise,**/
/** shall Thanawin Rakthanmanon, Bilson Campana, Abdullah Mueen,      **/
/** Gustavo Batista, or Eamonn Keogh be liable to you or to any other **/
/** person for any indirect, special, incidental, or consequential    **/
/** damages of any character including, without limitation, damages   **/
/** for loss of goodwill, work stoppage, computer failure or          **/
/** malfunction, or for any and all other damages or losses.          **/
/**                                                                   **/
/** If you do not agree with these terms, then you you are advised to **/
/** not use this software.                                            **/
/***********************************************************************/
/***********************************************************************/
#ifndef UCR_DTW_H
#define UCR_DTW_H


#ifndef level_t
#define level_t float
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
//#include <stdbool.h>

#define min(x,y) ((x)<(y)?(x):(y))
#define max(x,y) ((x)>(y)?(x):(y))
//#define dist(x,y) ((x-y)*(x-y))
#define dist(x,y) fabs(x-y)
#define INF 1000000       //Pseudo Infitinte number for this code

#ifndef GAP_PENALTY
#define GAP_PENALTY 5
#endif

/// Print function for debugging                                                                                                                                             
void printArray(level_t *x, int len)
{   for(int i=0; i<len; i++)
        printf("%6.2lf,",(float) x[i]);
    printf("\n");
}
//define trace_t
typedef struct{
	char step;
	level_t a,b;
    char qual, qual2;
}trace_t;
void printTrace(trace_t *t){
    int i=0;
    level_t dis=0;
    int v=0;
    int h=0;
    while(t[i].step=='d' || t[i].step=='v' || t[i].step=='h' || t[i].step=='i' || t[i].step=='w'){
        printf("TRACEstep,A,B,qual,qual2\t%c\t%.2f\t%.2f\t%i\t%i\t%i\n", t[i].step, t[i].a, t[i].b, t[i].qual, t[i].qual2,i);
        dis+=dist(t[i].a,t[i].b);
        if(t[i].step=='h'){
            dis+=GAP_PENALTY; h++;
        }else if(t[i].step=='v'){
            dis+=GAP_PENALTY; v++;
        }
        i++;
    }
    printf("TRACEstep,A,B,qual,qual2\t%c\t%.2f\t%.2f\t%i\t%i\t%i\n", t[i].step, t[i].a, t[i].b, t[i].qual, t[i].qual2,i);
    dis+=dist(t[i].a, t[i].b);
    char comment[16]="shit";
    if(dis==0)
        strcpy(comment,"perfect");
    else if(dis<20)
        strcpy(comment,"good");
    else if(dis<40)
        strcpy(comment,"mid");
    else if(dis<60)
        strcpy(comment,"bad");
    printf("TRACEend,dis,v,h,v+h,comment\t%.4f\t%i\t%i\t%i\t%s\n",dis,v,h,v+h,comment);
    ///
    int j;
    printf("TRACEsummary,");
    for(j=0; j<i; j++){
    if(t[j].step=='d'){
    	printf("d");
    }else if(t[j].step=='h'){
    	printf("h");
    }else if(t[j].step=='v'){
    	printf("dv");
    }
    }
    printf("\n");
}

//define result queue
typedef struct result_t
{   long long location;
    level_t distance;
    int flag;
    int net_indel;
}result_t;

void reset_resultq(result_t *res){
    res[0].location=res[1].location=-1;
    res[0].distance=res[1].distance=INF;
    res[0].flag=res[1].flag=0;
    res[0].net_indel=res[1].net_indel=0;
}
//add result to result queue, return the end of queue distance
//current version does not decrease number of result, may want to change this later!!
level_t add_result_cap2(result_t *res, long long loc, level_t dis, int f, int indel, int loc_mask){
    if(res[0].location == -1){
        res[0].location =loc; res[0].distance = dis; res[0].flag = f; res[0].net_indel=indel;
        return dis+30;  //tolerance at 20
    }else if(res[0].distance > dis) {
        if(abs(res[0].location-loc)> loc_mask){   //check for consecutive match
            long long tmploc=res[0].location;
            level_t tmpdis=res[0].distance;
            int tmpf=res[0].flag;
            int tmpi=res[0].net_indel;
            res[0].location =loc; res[0].distance = dis; res[0].flag = f;res[0].net_indel=indel;
            res[1].location =tmploc; res[1].distance = tmpdis; res[1].flag = tmpf;res[1].net_indel=tmpi;
            return tmpdis;
        }else{                           //consecutive match
            res[0].location =loc; res[0].distance = dis; res[0].flag = f; res[0].net_indel=indel;
            if(res[1].location != -1)
                return res[1].distance;
            return dis+30;
        }
    }else if(abs(res[0].location-loc)>loc_mask){
        if(res[1].location==-1 || res[1].distance> dis){
            res[1].location =loc; res[1].distance = dis; res[1].flag = f; res[1].net_indel=indel;
            return dis;
        }
    }
    if(res[1].location == -1)
        return res[0].distance+30;
    else
        return res[1].distance;
}

void printResult(result_t *res){
    printf("loc,dis,flag,indel,%lld,%.6f,%i,%i", res[0].location, res[0].distance, res[0].flag, res[0].net_indel);
}

/// Data structure for sorting the query
typedef struct index_t
    {   level_t value;
        int    index;
    } index_t;

/// Data structure (circular array) for finding minimum and maximum for LB_Keogh envolop
typedef struct deque
{   int *dq;
    int size,capacity;
    int f,r;
} deque;


/// Sorting function for the query, sort by abs(z_norm(q[i])) from high to low
int index_comp(const void *a, const void* b)
{   index_t* x = (index_t*)a;
    index_t* y = (index_t*)b;
    return abs(y->value) - abs(x->value);   // high to low
}

/// Initial the queue at the begining step of envelop calculation
void init(deque *d, int capacity)
{
    d->capacity = capacity;
    d->size = 0;
    d->dq = (int *) malloc(sizeof(int)*d->capacity);
    d->f = 0;
    d->r = d->capacity-1;
}

/// Destroy the queue
void destroy(deque *d)
{
    free(d->dq);
}

/// Insert to the queue at the back
void push_back(struct deque *d, int v)
{
    d->dq[d->r] = v;
    d->r--;
    if (d->r < 0)
        d->r = d->capacity-1;
    d->size++;
}

/// Delete the current (front) element from queue
void pop_front(struct deque *d)
{
    d->f--;
    if (d->f < 0)
        d->f = d->capacity-1;
    d->size--;
}

/// Delete the last element from queue
void pop_back(struct deque *d)
{
    d->r = (d->r+1)%d->capacity;
    d->size--;
}

/// Get the value at the current position of the circular queue
int front(struct deque *d)
{
    int aux = d->f - 1;

    if (aux < 0)
        aux = d->capacity-1;
    return d->dq[aux];
}

/// Get the value at the last position of the circular queueint back(struct deque *d)
int back(struct deque *d)
{
    int aux = (d->r+1)%d->capacity;
    return d->dq[aux];
}

/// Check whether or not the queue is empty
int empty(struct deque *d)
{
    return d->size == 0;
}

/// Finding the envelop of min and max value for LB_Keogh
/// Implementation idea is intoruduced by Danial Lemire in his paper
/// "Faster Retrieval with a Two-Pass Dynamic-Time-Warping Lower Bound", Pattern Recognition 42(9), 2009.
void lower_upper_lemire(level_t *t, int len, int r, level_t *l, level_t *u)
{
    struct deque du, dl;

    init(&du, 2*r+2);
    init(&dl, 2*r+2);

    push_back(&du, 0);
    push_back(&dl, 0);

    for (int i = 1; i < len; i++)
    {
        if (i > r)
        {
            u[i-r-1] = t[front(&du)];
            l[i-r-1] = t[front(&dl)];
        }
        if (t[i] > t[i-1])
        {
            pop_back(&du);
            while (!empty(&du) && t[i] > t[back(&du)])
                pop_back(&du);
        }
        else
        {
            pop_back(&dl);
            while (!empty(&dl) && t[i] < t[back(&dl)])
                pop_back(&dl);
        }
        push_back(&du, i);
        push_back(&dl, i);
        if (i == 2 * r + 1 + front(&du))
            pop_front(&du);
        else if (i == 2 * r + 1 + front(&dl))
            pop_front(&dl);
    }
    for (int i = len; i < len+r+1; i++)
    {
        u[i-r-1] = t[front(&du)];
        l[i-r-1] = t[front(&dl)];
        if (i-front(&du) >= 2 * r + 1)
            pop_front(&du);
        if (i-front(&dl) >= 2 * r + 1)
            pop_front(&dl);
    }
    destroy(&du);
    destroy(&dl);
    //printArray(u, len);
    //printArray(l, len);
}

/// Calculate quick lower bound
/// Usually, LB_Kim take time O(m) for finding top,bottom,fist and last.
/// However, because of z-normalization the top and bottom cannot give siginifant benefits.
/// And using the first and last points can be computed in constant time.
/// The prunning power of LB_Kim is non-trivial, especially when the query is not long, say in length 128.
level_t lb_kim_hierarchy_meanless(level_t *t, level_t *q, int j, int len, level_t best_so_far)
{
    /// 1 point at front and back
    level_t d, lb;
    level_t x0 = t[j];
    //level_t y0 = (t[(len-1+j)] - mean);
    lb = dist(x0,q[0]);
    if (lb >= best_so_far)   return lb;

    /// 2 points at front
    level_t x1 = t[(j+1)];
    d = min(dist(x1,q[0]), dist(x0,q[1]));
    d = min(d, dist(x1,q[1]));
    lb += d;
    if (lb >= best_so_far)   return lb;

    /// 2 points at back


    /// 3 points at front
    level_t x2 = t[(j+2)];
    d = min(dist(x0,q[2]), dist(x1, q[2]));
    d = min(d, dist(x2,q[2]));
    d = min(d, dist(x2,q[1]));
    d = min(d, dist(x2,q[0]));
    lb += d;
    //if (lb >= best_so_far)   return lb;

    /// 3 points at back

    return lb;
}

/// LB_Keogh 1: Create Envelop for the query
/// Note that because the query is known, envelop can be created once at the begenining.
///
/// Variable Explanation,
/// order : sorted indices for the query.
/// uo, lo: upper and lower envelops for the query, which already sorted.
/// t     : a circular array keeping the current data.
/// j     : index of the starting location in t
/// cb    : (output) current bound at each position. It will be used later for early abandoning in DTW.
level_t lb_keogh_cumulative_meanless(int* order, level_t *t, level_t *uo, level_t *lo, level_t *cb, int j, int len, level_t best_so_far)
{
    level_t lb = 0;
    level_t x, d;

    for (int i = 0; i < len && lb < best_so_far; i++)
    {
        x = t[(order[i]+j)];
        d = 0;
        if (x > uo[i])
            d = dist(x,uo[i]);
        else if(x < lo[i])
            d = dist(x,lo[i]);
        lb += d;
        cb[order[i]] = d;
    }
    return lb;
}

/// LB_Keogh 2: Create Envelop for the data
/// Note that the envelops have been created (in main function) when each data point has been read.
///
/// Variable Explanation,
/// tz: Z-normalized data
/// qo: sorted query
/// cb: (output) current bound at each position. Used later for early abandoning in DTW.
/// l,u: lower and upper envelop of the current data
level_t lb_keogh_data_cumulative_meanless(int* order, level_t *tz, level_t *qo, level_t *cb, level_t *l, level_t *u, int len, level_t best_so_far)
{
    level_t lb = 0;
    level_t uu,ll,d;

    for (int i = 0; i < len && lb < best_so_far; i++)
    {
        uu = u[order[i]];
        ll = l[order[i]];
        d = 0;
        if (qo[i] > uu)
            d = dist(qo[i], uu);
        else
        {   if(qo[i] < ll)
            d = dist(qo[i], ll);
        }
        lb += d;
        cb[order[i]] = d;
    }
    return lb;
}

/// Calculate Dynamic Time Wrapping distance
/// A,B: data and query, respectively
/// cb : cummulative bound used for early abandoning
/// r  : size of Sakoe-Chiba warpping band



//a INDEL_TOLERANCE seperate from r maybe too much hussle
//#ifndef INDEL_TOLERANCE
//#define INDEL_TOLERANCE 2
//#endif
///my dtw
///to do search two more rows...!!
level_t dtw2(level_t* B, level_t* A, level_t *cb, int m, int r, level_t best_so_far, int *indel) // A B switched from original ucr_suite
{

    level_t *cost;
    level_t *cost_prev;
    level_t *cost_prev2;
    level_t *cost_tmp;
    int i,j,k;  //for A B cost respectively, A query, B data
    level_t x,y,z,min_cost;
    level_t final_dtw=INF;

    /// Instead of using matrix of size O(m^2) or O(mr), we will reuse 3 array of size O(r).
    cost = (level_t*)malloc(sizeof(level_t)*(2*r+1));
    for(k=0; k<2*r+1; k++)    cost[k]=INF;
    cost_prev = (level_t*)malloc(sizeof(level_t)*(2*r+1));
    for(k=0; k<2*r+1; k++)    cost_prev[k]=INF;
    cost_prev2 = (level_t*)malloc(sizeof(level_t)*(2*r+1));
    for(k=0; k<2*r+1; k++)    cost_prev2[k]=INF;
    
    for (i=0; i<m; i++)
    {
        k = max(0,r-i);
        min_cost = INF;

        for(j=max(0,i-r); j<=min(m-1+r,i+r); j++, k++) //major mod! data has to be m+r long, function itself has no way to know, so be careful!
        {
            /// Initialize all row and column
            if ((i==0)&&(j==0))
            {
                cost[k]=dist(A[0],B[0]);
                min_cost = cost[k];
                continue;
            }

            //if ((j-1<0)||(k-1<0))     y = INF;
            if(j==0 || k==0) y=INF;
            else                      y = cost_prev[k-1]+dist(A[i],B[j-1]); //horizontal
            //if ((i-1<0)||(k+1>2*r))   x = INF;
            if(i==0 || k==2*r) x=INF;
            else                      x = cost_prev2[k+1]+dist(A[i-1],B[j]); //vertical
            //if ((i-1<0)||(j-1<0))     z = INF;
            if(i==0 || j==0) z=INF;
            else                      z = cost_prev[k]; //diag

            /// Classic DTW calculation
            cost[k] = min( min( x, y) + GAP_PENALTY, z) + dist(A[i],B[j]);
            /// Find minimum cost in row for early abandoning (possibly to use column instead of row).
            if (cost[k] < min_cost)
            {   min_cost = cost[k];
            	*indel=k-r;
            }
        }
        /// We can abandon early if the current cummulative distace with lower bound together are larger than best_so_far
        if (i+r < m-1 && min_cost + cb[i+r+1] >= best_so_far)
        {   free(cost);
            free(cost_prev);
	        free(cost_prev2);
            return min_cost + cb[i+r+1];
        }
        ///To allow variable suffix
        if(i==m-1){
            final_dtw=min_cost;
        }
        /// Move current array to previous array.
        cost_tmp = cost;
        cost = cost_prev2;
        cost_prev2=cost_prev;
        cost_prev = cost_tmp;
    }
    k--;

    //final_dtw = min(min(cost_prev[k], cost_prev2[k+1]), cost[k+2]);
    free(cost);
    free(cost_prev);
    free(cost_prev2);
    return final_dtw;
}
///rewrite dtw, clean up and allow trace back
level_t dtw3(level_t* B, level_t* A, level_t *cb, int m, int r, level_t best_so_far, int *indel) // A B switched from original ucr_suite
{
    int i,j,k;  //for A B cost respectively, A query, B data
    level_t x,y,z,min_cost;
    level_t final_dtw=INF;

    level_t **cost=malloc(m * sizeof *cost);
    for (i = 0; i < m; i++){
        cost[i] = malloc((2*r+1) * sizeof *cost[i]);
        if(cost[i]==NULL)
            printf("ERROR: Memory can't be allocated!\n");
        for(k=0; k<2*r+1;k++) cost[i][k]=INF;
    }

    
    for (i=0; i<m; i++)
    {
        k = max(0,r-i);
        min_cost = INF;

        for(j=max(0,i-r); j<=min(m-1+r,i+r); j++, k++) //major mod! data has to be m+r long, function itself has no way to know, so be careful!
        {
            /// Initialize all row and column
            if ((i==0)&&(j==0))
            {
                cost[i][k]=dist(A[0],B[0]);
                min_cost = cost[i][k];
                continue;
            }

            //if ((j-1<0)||(k-1<0))     y = INF;
            if(j==0 || k==0 || i==0) y=INF;
            else                      y = cost[i-1][k-1]+dist(A[i],B[j-1]); //horizontal
            //if ((i-1<0)||(k+1>2*r))   x = INF;
            if(i<2 || k==2*r) x=INF;
            else                      x = cost[i-2][k+1]+dist(A[i-1],B[j]); //vertical
            //if ((i-1<0)||(j-1<0))     z = INF;
            if(i==0 || j==0) z=INF;
            else                      z = cost[i-1][k]; //diag

            /// Classic DTW calculation
            cost[i][k] = min( min( x, y) + GAP_PENALTY, z) + dist(A[i],B[j]);
            /// Find minimum cost in row for early abandoning (possibly to use column instead of row).
            if (cost[i][k] < min_cost)
            {   min_cost = cost[i][k];
            	*indel=k-r;
            }
        }
        /// We can abandon early if the current cummulative distace with lower bound together are larger than best_so_far
        if (i+r < m-1 && min_cost + cb[i+r+1] >= best_so_far)
        {   for(j=0; j<m; j++)
                free(cost[j]);
            free(cost);
            return min_cost + cb[i+r+1];
        }
        ///To allow variable suffix
        if(i==m-1){
            final_dtw=min_cost;
        }
        /// Move current array to previous array.
    }
    //k--;

    //final_dtw = min(min(cost_prev[k], cost_prev2[k+1]), cost[k+2]);

    for(j=0; j<m; j++)
        free(cost[j]);
    free(cost);
    return final_dtw;
}
/*
        \w
      d \
 i\   \ |v
  \___\|
    h

*/
level_t dtw3_bt(level_t* B, level_t* A, char* Q, char* Q2, int m, int r) // Q same size as B, store char qual
{
    int i,j,k;  //for A B cost respectively, A query, B data
    int indel;
    level_t x,y,z,min_cost;
    level_t final_dtw=INF;

    level_t **cost=malloc(m * sizeof *cost);

    trace_t *t=malloc(2*m * sizeof *t);
    for (i = 0; i < m; i++){
        cost[i] = malloc((2*r+1) * sizeof *cost[i]);
        if(cost[i]==NULL)
            printf("ERROR: Memory can't be allocated!\n");
        for(k=0; k<2*r+1;k++) cost[i][k]=INF;
    }

    
    for (i=0; i<m; i++)
    {
        k = max(0,r-i);
        min_cost = INF;

        for(j=max(0,i-r); j<=min(m-1+r,i+r); j++, k++) //major mod! data has to be m+r long, function itself has no way to know, so be careful!
        {
            /// Initialize all row and column
            if ((i==0)&&(j==0))
            {
                cost[i][k]=dist(A[0],B[0]);
                min_cost = cost[i][k];
                continue;
            }

            //if ((j-1<0)||(k-1<0))     y = INF;
            if(j==0 || k==0 || i==0) y=INF;
            else                      y = cost[i-1][k-1]+dist(A[i],B[j-1]); //horizontal
            //if ((i-1<0)||(k+1>2*r))   x = INF;
            if(i<2 || k==2*r) x=INF;
            else                      x = cost[i-2][k+1]+dist(A[i-1],B[j]); //vertical
            //if ((i-1<0)||(j-1<0))     z = INF;
            if(i==0 || j==0) z=INF;
            else                      z = cost[i-1][k]; //diag

            /// Classic DTW calculation
            cost[i][k] = min( min( x, y) + GAP_PENALTY, z) + dist(A[i],B[j]);
            /// Find minimum cost in row for early abandoning (possibly to use column instead of row).
            if (cost[i][k] < min_cost)
            {   min_cost = cost[i][k];
            	indel=k-r;
            }
        }
        ///To allow variable suffix
        if(i==m-1){
            final_dtw=min_cost;
        }
        /// Move current array to previous array.
    }
    //k--;

    //final_dtw = min(min(cost_prev[k], cost_prev2[k+1]), cost[k+2]);
    ///back tracing
    i=m-1;
    k=indel+r;
    j=m-1+indel;
    int ii=0;
    while(i>0){
    	z=cost[i-1][k]-GAP_PENALTY;
    	if(k==0 || j==0) y=INF;
    	else y=cost[i-1][k-1]+dist(A[i],B[j-1]);
    	if(i<2 || k==2*r) x=INF; 
    	else x=cost[i-2][k+1]+dist(A[i-1],B[j]);
    	if(z<=y && z<=x){
    	  //trace diag
          t[ii].step='d';t[ii].a=A[i];t[ii].b=B[j];t[ii].qual=Q[j];t[ii].qual2=Q2[j];
          ii++;
          //
          i--;
          j--;
    	}else if(x<=y){
    	  //trace vert
          t[ii].step='v';t[ii].a=A[i];t[ii].b=B[j];t[ii].qual=Q[j];t[ii].qual2=Q2[j];
          ii++;
          t[ii].step='w';t[ii].a=A[i-1];t[ii].b=B[j];t[ii].qual=Q[j];t[ii].qual2=Q2[j];
          ii++;
          //
          i-=2;
          j--;
          k++;
    	}else{
    	  //trace hori
          t[ii].step='h';t[ii].a=A[i];t[ii].b=B[j];t[ii].qual=Q[j];t[ii].qual2=Q2[j];
          ii++;
          t[ii].step='i';t[ii].a=A[i];t[ii].b=B[j-1];t[ii].qual=Q[j-1];t[ii].qual2=Q2[j-1];
          ii++;
          //
          i--;
          j-=2;
          k--;
    	}
    }
    t[ii].step='s';t[ii].a=A[0];t[ii].b=B[0];t[ii].qual=Q[0];t[ii].qual2=Q2[0];
    //t[ii].step='s';t[ii].a=A[i];t[ii].b=B[j];t[ii].qual=Q[j];
    if(i!=0){
    	printf("DTW_BT: impossible!\n");
    }
    
    for(j=0; j<m; j++)
        free(cost[j]);
    free(cost);
    printTrace(t);
    return final_dtw;
}




/// If expected error happens, teminated the program.
void error(int id)
{
    if(id==1)
        printf("ERROR : Memory can't be allocated!!!\n\n");
    else if ( id == 2 )
        printf("ERROR : File not Found!!!\n\n");
    else if ( id == 3 )
        printf("ERROR : Can't create Output File!!!\n\n");
    else if ( id == 4 )
    {
        printf("ERROR : Invalid Number of Arguments!!!\n");
        printf("Command Usage:  UCR_DTW.exe  data-file  query-file   m   R\n\n");
        printf("For example  :  UCR_DTW.exe  data.txt   query.txt   128  0.05\n");
    }
    exit(1);
}

/// Calculate the nearest neighbor of a times series in a larger time series expressed as location and distance,
/// using the UCR suite optimizations. from github/klon
int ucrdtw(level_t* data, long long data_size, level_t* query, long query_size, level_t warp_width, int verbose, result_t* res, level_t bsf, int flag) {
    long m = query_size;
    int r = warp_width <= 1 ? floor(warp_width * m) : floor(warp_width);
    //r=max(r, INDEL_TOLERANCE);
    long M= m+r;
    //best_so_far = INF;
    level_t best_so_far=bsf; /// best-so-far
    level_t *q; /// data array
    int *order; ///new order of the query
    level_t *u, *l, *qo, *uo, *lo, *cb, *cb1, *cb2;

    //level_t d = 0.0;
    long long i;
    //level_t mean=90;   //use a fixed mean value!! 90 is just quick and dirty

    long long loc = -1; //my mod
    double t1, t2;
    int kim = 0, keogh = 0, keogh2 = 0;
    level_t dist = 0, lb_kim = 0, lb_k = 0, lb_k2 = 0;
    level_t *u_buff, *l_buff;
    index_t *q_tmp;

    if (verbose) {
        t1 = clock();
    }
    //switch which best_so_far to use, supplied or INF

    /// calloc everything here
    q = (level_t*) calloc(m, sizeof(level_t));
    if (q == NULL) {
        printf("ERROR: Memory can't be allocated!\n");
        return -1;
    }
    memcpy((void*)q, (void*)query, m * sizeof(level_t));

    qo = (level_t*) calloc(m, sizeof(level_t));
    if (qo == NULL) {
        printf("ERROR: Memory can't be allocated!\n");
        return -1;
    }

    uo = (level_t*) calloc(m, sizeof(level_t));
    if (uo == NULL) {
        printf("ERROR: Memory can't be allocated!\n");
        return -1;
    }

    lo = (level_t*) calloc(m, sizeof(level_t));
    if (lo == NULL) {
        printf("ERROR: Memory can't be allocated!\n");
        return -1;
    }

    order = (int *) calloc(m, sizeof(int));
    if (order == NULL) {
        printf("ERROR: Memory can't be allocated!\n");
        return -1;
    }

    q_tmp = (index_t *) calloc(m, sizeof(index_t));
    if (q_tmp == NULL) {
        printf("ERROR: Memory can't be allocated!\n");
        return -1;
    }

    u = (level_t*) calloc(m, sizeof(level_t));
    if (u == NULL) {
        printf("ERROR: Memory can't be allocated!\n");
        return -1;
    }

    l = (level_t*) calloc(m, sizeof(level_t));
    if (l == NULL) {
        printf("ERROR: Memory can't be allocated!\n");
        return -1;
    }

    cb = (level_t*) calloc(m, sizeof(level_t));
    if (cb == NULL) {
        printf("ERROR: Memory can't be allocated!\n");
    }

    cb1 = (level_t*) calloc(m, sizeof(level_t));
    if (cb1 == NULL) {
        printf("ERROR: Memory can't be allocated!\n");
        return -1;
    }

    cb2 = (level_t*) calloc(m, sizeof(level_t));
    if (cb2 == NULL) {
        printf("ERROR: Memory can't be allocated!\n");
        return -1;
    }


    u_buff = (level_t*) calloc(data_size, sizeof(level_t));
    if (u_buff == NULL) {
        printf("ERROR: Memory can't be allocated!\n");
        return -1;
    }

    l_buff = (level_t*) calloc(data_size, sizeof(level_t));
    if (l_buff == NULL) {
        printf("ERROR: Memory can't be allocated!\n");
        return -1;
    }


    /// Create envelope of the query: lower envelope, l, and upper envelope, u
    lower_upper_lemire(q, m, r, l, u);

    /// Sort the query one time by abs(z-norm(q[i]))
    for (i = 0; i < m; i++) {
        q_tmp[i].value = q[i];
        q_tmp[i].index = i;
    }
    qsort(q_tmp, m, sizeof(index_t), index_comp);

    /// also create another arrays for keeping sorted envelope
    for (i = 0; i < m; i++) {
        int o = q_tmp[i].index;
        order[i] = o;
        qo[i] = q[o];
        uo[i] = u[o];
        lo[i] = l[o];
    }

    /// Initial the cummulative lower bound
    for (i = 0; i < m; i++) {
        cb[i] = 0;
        cb1[i] = 0;
        cb2[i] = 0;
    }

    i = 0;          /// current index of the data in current chunk of size EPOCH
    //j = 0;          /// the starting index of the data in the circular array, t
    int done = 0;
    //int it = 0, ep = 0, k = 0;
    int k = 0;
    //long long I; /// the starting index of the data in current chunk of size EPOCH
    //long long data_index =0; //??
    int indel; //this is for storing the net_indel inside dtw.
    while (!done) {

        if (data_size < M-1) {
            done = 1;
        } else {
            // for(i=0; i< data_size; i++){
            //     buffer[i]=data[i]-mean;
            // }
            lower_upper_lemire(data, data_size, r, l_buff, u_buff);
            /// Do main task here..
            for (i = M-1; i < data_size; i++) { 
                /// A bunch of data has been read and pick one of them at a time to use
                //d = data[i];

                /// Start the task when there are more than m-1 points in the current chunk


                // if (i >= M - 1) {
                //     for (k = 0; k < M; k++) {
			    //         tz[k] = data[(i-M+1+ k)];
                //     }

                    /// Use a constant lower bound to prune the obvious subsequence
                    lb_kim = lb_kim_hierarchy_meanless(data, q, i-M+1, m, best_so_far); //!!!

                    if (lb_kim < best_so_far) {
                        /// Use a linear time lower bound to prune; z_normalization of t will be computed on the fly.
                        /// uo, lo are envelope of the query.
                        //lb_k = lb_keogh_cumulative(order, data, uo, lo, cb1, j, m, mean, best_so_far);
                        lb_k = lb_keogh_cumulative_meanless(order, data, uo, lo, cb1, i-M+1, m, best_so_far); ///!!!
                        if (lb_k < best_so_far) {
                            /// Take another linear time to compute z_normalization of t.
                            /// Note that for better optimization, this can merge to the previous function.
                
                    // for (k = 0; k < M; k++) {
			        //     tz[k] = buffer[(i-M+1+ k)];
                    // }

                            /// Use another lb_keogh to prune
                            /// qo is the sorted query. tz is unsorted z_normalized data.
                            /// l_buff, u_buff are big envelope for all data in this chunk
                            lb_k2 = lb_keogh_data_cumulative_meanless(order, data+i-M+1, qo, cb2, l_buff + i-M+1, u_buff + i-M+1, m, best_so_far); ///!!!
                            if (lb_k2 < best_so_far) {
                                /// Choose better lower bound between lb_keogh and lb_keogh2 to be used in early abandoning DTW
                                /// Note that cb and cb2 will be cumulative summed here.
                                /// Actually this might not be air tight, 
                                if (lb_k > lb_k2) {
                                    cb[m - 1] = cb1[m - 1];
                                    for (k = m - 2; k >= 0; k--)
                                        cb[k] = cb[k + 1] + cb1[k];
                                } else {
                                    cb[m - 1] = cb2[m - 1];
                                    for (k = m - 2; k >= 0; k--)
                                        cb[k] = cb[k + 1] + cb2[k];
                                }

                                /// Compute DTW and early abandoning if possible
                                dist = dtw3(data+i-M+1, q, cb, m, r, best_so_far, &indel);
//printf("tmpdebug i: %lli, dist: %.4f\n", i, dist);
                                if (dist < best_so_far) {   /// Update best_so_far
                                                    /// loc is the real starting location of the nearest neighbor in the file                                   
                                    //loc = (it) * (EPOCH - M + 1) + i - M + 1;
                                    loc=i-M+1;
                                    best_so_far = add_result_cap2(res, loc, dist, flag, indel, 5);
                                }
                            } else
                                keogh2++;
                        } else
                            keogh++;
                    } else
                        kim++;

                    /// Reduce absolute points from sum and sum square

                
            }
            done=1;
            /// If the size of last chunk is less then EPOCH, then no more data and terminate.
            // if (ep < EPOCH)
            //     done = 1;
            // else
            //     it++;
        }
    }

   // i = (it) * (EPOCH - M + 1) + ep;

    free(q);
    free(qo);
    free(uo);
    free(lo);
    free(order);
    free(q_tmp);
    free(u);
    free(l);
    free(cb);
    free(cb1);
    free(cb2);
    //free(t);
    //free(tz);
    //free(buffer);
    free(u_buff);
    free(l_buff);

    if (verbose) {
        t2 = clock();
        printResult(res);
        printf("\n");
        printf("Data Scanned : %lld\n", i);
        printf("Total Execution Time : %.4f secs\n", (t2 - t1) / CLOCKS_PER_SEC);
        printf("\n");
        printf("Pruned by LB_Kim    : %6.2f%%\n", ((double) kim / i) * 100);
        printf("Pruned by LB_Keogh  : %6.2f%%\n", ((double) keogh / i) * 100);
        printf("Pruned by LB_Keogh2 : %6.2f%%\n", ((double) keogh2 / i) * 100);
        printf("DTW Calculation     : %6.2f%%\n", 100 - (((double) kim + keogh + keogh2) / i * 100));
    }

    return 0;
}


#endif
