#include <stdlib.h>
#include <unistd.h>
#include <pthread.h>

#include "pagerank.h"


double * granks;
double * gprevranks;
double * gnorms_t;
double gdamp;
double gnorm;
double dampratio;
int gnpages;
int nthreads;
page** gplist;
pthread_t * threads;


void * get_new_rank(void * pin){

    double sum;
    int lpin = (int) pin;
    page * p;
    gnorms_t[lpin] = 0;
    double prevrank;
    for (int i = lpin; i < gnpages; i += nthreads){
        sum = 0;
        p = gplist[i];
        if(p->inlinks){
            node * curr = p->inlinks->head;
            node * prev = NULL;
            while(1){
                prevrank = dampratio;
                if(curr->page->inlinks){
                    prevrank = gprevranks[curr->page->index];
                }
                sum += (prevrank/curr->page->noutlinks);
                if (p->inlinks->tail == curr){
                    break;
                }
                prev = curr;
                curr = prev->next;
            }
            sum *= gdamp;
        }
        granks[i] = dampratio + sum;
        sum = (granks[i] - gprevranks[i]);
        gnorms_t[lpin] += (sum*sum);
    }
    return NULL;
}

void pagerank(list* plist, int ncores, int npages, int nedges, double dampener)
{
    nthreads = npages;
    if (nthreads > ncores){
        nthreads = ncores;
    }
    granks = (double *)malloc(sizeof(double)*npages);
    gprevranks = (double *)malloc(sizeof(double)*npages);
    threads = (pthread_t *)malloc(sizeof(pthread_t)*nthreads);
    gplist = (page **)malloc(sizeof(page *)*npages);
    gnorms_t = (double *)malloc(sizeof(double)*nthreads);
    dampratio = ((1 - dampener)/npages);
    gnpages = npages;
    gdamp = dampener;
    node* curr;
    node* prev;
    curr = plist->head;
    int count = 0;
    for (int i = 0; i<npages; i++){
        granks[i] = 1.0/npages;
        gprevranks[i] = 1.0/npages;
        if(curr->page->inlinks){
            gplist[count] = curr->page;
            curr->page->index = count;
            count += 1;
        }
        prev = curr;
        curr = prev->next;
    }
    nthreads = count;
    gnpages = count;
    if (nthreads > ncores){
        nthreads = ncores;
    }

    double * temp;
    do{
        gnorm = 0;
        for (int i = 0; i < nthreads; i++){
            pthread_create(&threads[i], NULL, get_new_rank, (void *) i);
        }
        for (int i = 0; i < nthreads; i++){
            pthread_join(threads[i], NULL);
        }
        for (int i = 0; i < nthreads; i++){
            gnorm += gnorms_t[i];
        }
        temp = gprevranks;
        gprevranks = granks;
        granks = temp;
    } while(gnorm > (EPSILON*EPSILON));

    curr = plist->head;
    prev = NULL;
    double rank;
    while (1) {
        rank = dampratio;
        if(curr->page->inlinks){
           rank = gprevranks[curr->page->index];
        }
        printf("%s %.4f\n",curr->page->name, rank);
        if (plist->tail == curr){
           break;
        }
        prev = curr;
        curr = prev->next;
     }
    free(granks);
    free(gprevranks);
    free(gnorms_t);
    free(gplist);
    free(threads);

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
