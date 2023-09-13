#ifndef DECODER_H
#define DECODER_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "UCR_DTW2.h"
#include "pore_model.h"

using std::cout;
using std::endl;

void reverse_levels(level_t *q, long len){
	level_t tmp;
	for(long i=0; i<floor(len/2); i++){
		tmp=q[i];
		q[i]=q[len-1-i];
		q[len-1-i]=tmp;
	}
}

///for now, this is a hack
result_t pick_result_cap2_alt(result_t *res, long qsize){
    result_t ret;
    level_t dis_cut=80;
    ret.location=-1;
    ret.distance=INF;
    ret.flag=1;
    if(res[0].location == -1){
        ret=res[0];
        return ret;
    }
    if(res[0].flag==1){
        if(res[0].distance < dis_cut-40 && res[0].location > 200)
            ret=res[0];
        return ret;
    }
    return ret;

}



result_t pick_result_cap2(result_t *res, long qsize){
    result_t ret;
    level_t dis_cut=80;
    ret.location=-1;
    ret.distance=INF;
    ret.flag=1;
    if(res[0].location == -1){
        ret=res[0];
        return ret;
    }
    if(res[0].flag==1){
        if(res[0].distance < dis_cut && res[0].location < 40)
            ret=res[0];
        return ret;
    }else{
        if(res[1].flag==1 && res[1].distance < dis_cut-20 && res[1].location < 40){
            ret=res[1];
        }else if(res[0].distance < dis_cut && res[0].location> qsize-60){
            ret=res[0];
        }else if(res[1].flag==-1 && res[1].distance < dis_cut-20 && res[1].location > qsize-60){
            ret=res[1];
        }
        
    }    
    return ret;

}

typedef struct {
	char **seq;
	level_t **lvl, **lvl_r, **u, **l, **ur, **lr;
	int **order, **order_r;
	int size, max_size, length, max_length, warp_width;
}WL_DB;

const char* getSeq(WL_DB *db, long loc,int offset){
	if(loc<0 || loc> db->size -1 ){
		return "NONE";	
	}else{
		return db->seq[loc]+offset;	
	}
}

WL_DB *init_WL_DB(int db_size, int len){
	WL_DB *ret=malloc(sizeof *ret);
	if(ret==NULL)
		printf("Error: unable to malloc\n");
	(ret->seq = malloc(db_size * sizeof *ret->seq)) || printf("Error: unable to malloc\n");
	(ret->lvl=malloc(db_size * sizeof *ret->lvl)) || printf("Error: unable to malloc\n");
	(ret->lvl_r=malloc(db_size * sizeof *ret->lvl_r)) || printf("Error: unable to malloc\n");
	(ret->u=malloc(db_size * sizeof *ret->u)) || printf("Error: unable to malloc\n");
	(ret->l=malloc(db_size * sizeof *ret->l)) || printf("Error: unable to malloc\n");
	(ret->ur=malloc(db_size * sizeof *ret->ur)) || printf("Error: unable to malloc\n");
	(ret->lr=malloc(db_size * sizeof *ret->lr)) || printf("Error: unable to malloc\n");
	(ret->order=malloc(db_size * sizeof *ret->order)) || printf("Error: unable to malloc\n");
	(ret->order_r=malloc(db_size * sizeof *ret->order_r)) || printf("Error: unable to malloc\n");
	for(int i=0; i<db_size; i++){
		(ret->seq[i] = malloc(len * sizeof *ret->seq[i])) || printf("Error: unable to malloc\n");
		(ret->lvl[i]=malloc(len * sizeof *ret->lvl[i])) || printf("Error: unable to malloc\n");
		(ret->lvl_r[i]=malloc(len * sizeof *ret->lvl_r[i])) || printf("Error: unable to malloc\n");
		(ret->u[i]=malloc(len * sizeof *ret->u[i])) || printf("Error: unable to malloc\n");
		(ret->l[i]=malloc(len * sizeof *ret->l[i])) || printf("Error: unable to malloc\n");
		(ret->ur[i]=malloc(len * sizeof *ret->ur[i])) || printf("Error: unable to malloc\n");
		(ret->lr[i]=malloc(len * sizeof *ret->lr[i])) || printf("Error: unable to malloc\n");
		(ret->order[i]=malloc(len * sizeof *ret->order[i])) || printf("Error: unable to malloc\n");
		(ret->order_r[i]=malloc(len * sizeof *ret->order_r[i])) || printf("Error: unable to malloc\n");
	}
	ret->max_size=db_size;
	ret->max_length=len;
	return ret;
}
void free_WL_DB(WL_DB *db){
	for(int i=0; i<db->max_size; i++){
		free(db->seq[i]);
		free(db->lvl[i]);
		free(db->lvl_r[i]);
		free(db->u[i]);
		free(db->l[i]);
		free(db->ur[i]);
		free(db->lr[i]);
		free(db->order[i]);
		free(db->order_r[i]);
	}
	free(db->seq);
	free(db->lvl);
	free(db->lvl_r);
	free(db->u);
	free(db->l);
	free(db->ur);
	free(db->lr);
	free(db->order);
	free(db->order_r);
	free(db);
}
//long read_WL(const char *f, const entry_t *m, const char *prefix, long max_size, long append_inf, char **whitelist, level_t **whitelist_level, level_t **whitelist_level_r){
long read_WL(const char *f, const entry_t *m, const char *prefix, long append_inf, int r, WL_DB *db){ //return number of entries read, may want to remove append_inf
    FILE *file;
    char line[1024];
    char tmp[1024];
    int i=0;
    int len;
    index_t q_tmp[1024];
    //long plen=strlen(prefix);
    file = fopen(f, "r");
    if (file == NULL) {
        fprintf(stderr, "read_whitelist error: Unable to open the file %s\n", f);
        return -1;
    } 
    while(fgets(line, sizeof(line), file)){
        line[strcspn(line, "\n")]=0; // remove newline
        //strcpy(whitelist[i],prefix);
        //strcat(whitelist[i],line);
        //seq2level(m, whitelist[i], whitelist_level[i], max_size, append_inf);
        //strcpy(tmp,whitelist[i]);
        //revcom(tmp);
        //seq2level(m, tmp, whitelist_level_r[i],max_size, append_inf);
        strcpy(db->seq[i],prefix);
        strcat(db->seq[i],line);
        seq2level(m, db->seq[i], db->lvl[i],db->max_length, append_inf);
        strcpy(tmp,db->seq[i]);
        revcom(tmp);
        seq2level(m,tmp, db->lvl_r[i], db->max_length, append_inf);
        
        len=strlen(db->seq[i])+append_inf;
        lower_upper_lemire(db->lvl[i], len, r, db->l[i], db->u[i]);
        lower_upper_lemire(db->lvl_r[i],len, r, db->lr[i],db->ur[i]);
        ///fill in order
        for (int j = 0; j < len; j++) {
        	q_tmp[j].value = db->lvl[i][j];
        	q_tmp[j].index = j;
    	}
    	qsort(q_tmp, len, sizeof(index_t), index_comp);


    	for (int j = 0; j < len; j++) {
        	db->order[i][j] = q_tmp[j].index;
    	}
    	///fill in order_r
        for (int j = 0; j < len; j++) {
        	q_tmp[j].value = db->lvl_r[i][j];
        	q_tmp[j].index = j;
    	}
    	qsort(q_tmp, len, sizeof(index_t), index_comp);


    	for (int j = 0; j < len; j++) {
        	db->order_r[i][j] = q_tmp[j].index;
    	}
        i+=1;
    }
    db->size=i;
    db->warp_width=r;
    db->length=strlen(db->seq[0])+append_inf-MER_LENGTH+1;   //record level length
    fprintf(stderr,"read_whitelist: %i entries read.\n", db->size);
    return i;
}

level_t dtw3_lite(level_t* B, level_t* A, int m, int r, level_t best_so_far) // A B switched from original ucr_suite
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

        for(j=max(0,i-r); j<=min(m-1,i+r); j++, k++) 
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
            }
        }
//early abandon
        if (min_cost >= best_so_far)
        {   for(j=0; j<m; j++)
                free(cost[j]);
            free(cost);
            return min_cost;
        }
        ///To allow variable suffix
        if(i==m-1){
            //final_dtw=min_cost;
            final_dtw=cost[i][r];
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

///dtw wrapper for whitelist check
//int ucrdtw_checkWL(level_t* data, long long data_size, level_t* query, long query_size, level_t warp_width, int verbose, result_t* res, level_t bsf, int flag) {
int ucrdtw_checkWL(level_t* data, long long data_size, level_t* u_d, level_t* l_d, WL_DB* db,long wl_i ,int reverse,int verbose, result_t* res, level_t bsf, int flag){ //reverse 0 forward otherwise reverse
//    cout << "ucrdtw_checkWL start" << endl;
//    cout << wl_i << endl; //0
    long m = db->length;
//    cout << m << endl; //3420
    int r = db->warp_width;
//    cout << r << endl; //3828
    long M= m+r;
    //switch which best_so_far to use, supplied or INF
    level_t best_so_far=bsf; /// best-so-far
    level_t *q, *u, *l;
    int *order;
    
    if(reverse==0){
    	q=db->lvl[wl_i];
    	u=db->u[wl_i];
    	l=db->l[wl_i];
    	order=db->order[wl_i];
    }else{
    	q=db->lvl_r[wl_i];
    	u=db->ur[wl_i];
    	l=db->lr[wl_i];
    	order=db->order_r[wl_i];    
    }

    level_t *qo, *uo, *lo, *cb, *cb1, *cb2;

    //level_t d = 0.0;
    long long i;
    //level_t mean=90;   //use a fixed mean value!! 90 is just quick and dirty

    long long loc = -1; //my mod
    double t1, t2;
    int kim = 0, keogh = 0, keogh2 = 0;
    level_t dist = 0, lb_kim = 0, lb_k = 0, lb_k2 = 0;
    //level_t *u_buff, *l_buff;
    index_t *q_tmp;


    if (verbose) {
        t1 = clock();
    }
    

    /// calloc everything here


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
    
    /// also create another arrays for keeping sorted envelope
    for (int i = 0; i < m; i++) {
        int o=order[i];
        qo[i] = q[o];
        uo[i] = u[o];
        lo[i] = l[o];
    }
    
    //cout << "ucrdtw_checkWL: Initial the cummulative lower bound" << endl;
    /// Initial the cummulative lower bound
    for (int i = 0; i < m; i++) {
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
//    cout << "ucrdtw_checkWL: while" << endl;
    while (!done) {

        if (data_size < M-1) {
            done = 1;
        } else {
            // for(i=0; i< data_size; i++){
            //     buffer[i]=data[i]-mean;
            // }
            //lower_upper_lemire(data, data_size, r, l_buff, u_buff);
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
                            lb_k2 = lb_keogh_data_cumulative_meanless(order, data+i-M+1, qo, cb2, l_d + i-M+1, u_d + i-M+1, m, best_so_far); ///!!!
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

    free(qo);
    free(uo);
    free(lo);
    free(cb);
    free(cb1);
    free(cb2);
    //free(t);
    //free(tz);
    //free(buffer);

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
