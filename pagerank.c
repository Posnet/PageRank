#include <stdlib.h>
#include <unistd.h>
#include <pthread.h>

#include "pagerank.h"

long double * granks;
long double * gprevranks;
long double * gnorms_t;
long double gdamp;
long double gnorm;
long double dampratio;
int gnpages;
int nthreads;
int twidth;
page** gplist;
pthread_t * threads;


void * get_new_rank(void * pin){

    long double sum;
    int lpin = (int)pin;
    int start = lpin * twidth;
    int end = lpin + twidth;
    page * p;
    if (end+twidth > gnpages){
        end = gnpages;
    }
    for (int i = start; i < end; i++){
        sum = 0;
        gnorms_t[i] = 0;
        p = gplist[i];
        if(p->inlinks){
            node * curr = p->inlinks->head;
            node * prev = NULL;
            while(curr){
                sum += (gprevranks[curr->page->index]/curr->page->noutlinks);
                prev = curr;
                curr = prev->next;
            }
            sum *= gdamp;
        }
        granks[i] = dampratio + sum;
        sum = (granks[i] - gprevranks[i]);
        gnorms_t[i] += (sum*sum);
    }
    return NULL;
}

void pagerank(list* plist, int ncores, int npages, int nedges, long double dampener)
{
    nthreads = npages;
    granks = (long double *)malloc(sizeof(long double)*npages);
    gprevranks = (long double *)malloc(sizeof(long double)*npages);
    threads = (pthread_t *)malloc(sizeof(pthread_t)*nthreads);
    gplist = (page **)malloc(sizeof(page *)*npages);
    gnorms_t = (long double *)malloc(sizeof(long double)*nthreads);
    dampratio = ((1 - dampener)/npages);
    twidth = (int)(npages/nthreads);
    gnpages = npages;
    gdamp = dampener;
    node* curr;
    node* prev;
    curr = plist->head;
    for (int i = 0; i<npages; i++){
        granks[i] = 1.0/npages;
        gprevranks[i] = 1.0/npages;
        gplist[i] = curr->page;
        prev = curr;
        curr = prev->next;
    }

    long double * temp;
    do{
        gnorm = 0;
        for (int i = 0; i < nthreads; i++){
            int * index = &i;
            pthread_create(&threads[i], NULL, get_new_rank, (void *) index);
        }
        for (int i = 0; i < nthreads; i++){
            pthread_join(threads[i], NULL);
        }
        for (int i = 0; i<npages; i++){
            gnorm += gnorms_t[i];
        }
        temp = gprevranks;
        gprevranks = granks;
        granks = temp;
    } while(gnorm > (EPSILON*EPSILON));

    curr = plist->head;
    prev = NULL;
     while(curr){
         printf("%s %.4Lf\n",curr->page->name, gprevranks[curr->page->index]);
         prev = curr;
         curr = prev->next;
    }
    free(temp);
    free(gprevranks);
    free(threads);
    free(plist);
    free(gnorms_t);
}

/* DO NOT MODIFY BELOW THIS POINT */
int main(void)
{
    /*
     * DO NOT MODIFY THE MAIN FUNCTION OR HEADER FILE
     */

    double damping_factor;
    list* plist = NULL;
    int ncores, npages, nedges;

    /* read the input into the appropriate variables and populate the list of pages */
    read_input(&plist, &ncores, &npages, &nedges, &damping_factor);

    /* run pagerank and print the results */
    pagerank(plist, ncores, npages, nedges, damping_factor);

    /* clean up the memory used by the list of pages */
    page_list_destroy(plist);

    return 0;
}
