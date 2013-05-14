#include <stdlib.h>
#include <unistd.h>
#include <pthread.h>

#include "pagerank.h"

long double * granks;
long double * gprevranks;
int gnpages;
long double gdamp;

long double find_norm2(long double * ranks, long double * prevranks, int npages){
    long double sum = 0;
    long double x = 0;
    for(int i = 0; i<npages; i++){
        x = ranks[i] - prevranks[i];
        sum += (x*x);
    }
    return sum;
}

long double get_new_rank2(page * p, long double * prevranks, int npages, long double damp){
    long double a = (1-damp)/npages;
    long double sum = 0;
    if(p->inlinks){
        node * curr = p->inlinks->head;
        node * prev = NULL;
        while(curr){
            sum += (prevranks[curr->page->index]/curr->page->noutlinks);
            prev = curr;
            curr = prev->next;
        }
        sum *= damp;
    }
    return a + sum;
}

void pagerank2(list* plist, int ncores, int npages, int nedges, long double dampener)
{
    long double * ranks = (long double *)malloc(sizeof(long double)*npages);
    long double * prevranks = (long double *)malloc(sizeof(long double)*npages);
    for (int i = 0; i<npages; i++){
        ranks[i] = 1.0/npages;
        prevranks[i] = 1.0/npages;
    }
    node* curr;
    node* prev;
    int count = 1;
    long double norm = 0;
    do{
        for(int i = 0; i<npages; i++){
            prevranks[i] = ranks[i];
        }
        curr = plist->head;
        prev = NULL;
        while(curr){
            ranks[curr->page->index] = get_new_rank2(curr->page, prevranks, npages, dampener);
            prev = curr;
            curr = prev->next;
        }
        norm = find_norm2(ranks, prevranks, npages);
        //printf("%d: %Lf\n", count, norm);
        count += 1;
    } while(norm > (EPSILON*EPSILON));

    curr = plist->head;
    prev = NULL;
    while(curr){
        printf("%s %.4Lf\n",curr->page->name, ranks[curr->page->index]);
        prev = curr;
        curr = prev->next;
    }
    free(ranks);
    free(prevranks);
}

long double find_norm(){
    long double sum = 0;
    long double x = 0;
    for(int i = 0; i<gnpages; i++){
        x = granks[i] - gprevranks[i];
        sum += (x*x);
    }
    return sum;
}

void * get_new_rank(void * pin){
    long double a = (1-gdamp)/gnpages;
    long double sum = 0;
    page * p = (page *)pin;
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
    granks[p->index] = a + sum;
    return NULL;
}

void pagerank(list* plist, int ncores, int npages, int nedges, long double dampener)
{
    if (npages > 200){
        pagerank2(plist, ncores, npages, nedges, dampener);
        return;
    }
    granks = (long double *)malloc(sizeof(long double)*npages);
    gprevranks = (long double *)malloc(sizeof(long double)*npages);
    gnpages = npages;
    gdamp = dampener;
    pthread_t * threads = (pthread_t *)malloc(sizeof(pthread_t)*npages);
    for (int i = 0; i<npages; i++){
        granks[i] = 1.0/npages;
        gprevranks[i] = 1.0/npages;
    }

    node* curr;
    node* prev;
    int count = 1;
    long double norm;
    do{
        for(int i = 0; i<npages; i++){
            gprevranks[i] = granks[i];
        }
        curr = plist->head;
        prev = NULL;
        for (int i = 0; i < npages; i++){
            pthread_create(&threads[i], NULL, get_new_rank, curr->page);
            prev = curr;
            curr = prev->next;
        }
        for (int i = 0; i < npages; i++){
            pthread_join(threads[i], NULL);
        }
        norm = find_norm();
        //printf("%d: %Lf\n", count, norm);
        count += 1;
    } while(norm > (EPSILON*EPSILON));

    curr = plist->head;
    prev = NULL;
    while(curr){
        printf("%s %.4Lf\n",curr->page->name, granks[curr->page->index]);
        prev = curr;
        curr = prev->next;
    }
    free(granks);
    free(gprevranks);
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
