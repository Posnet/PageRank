/////////////
//Includes //
/////////////

#include <stdlib.h>
#include <unistd.h>
#include <pthread.h>

#include "pagerank.h"

////////////
//Globals //
////////////

int ***nodes;
double *PageRank;
double *PrevRank;
double *outlinks;
double *TempRank = NULL;
double epsilon   = EPSILON * EPSILON;
double jumpProb;
double norm = 1;
double next_norm = 1;
int * gnpages;
int gnthreads;
int thread_done = 0;
double *local_norms;

pthread_barrier_t threadSync;
pthread_mutex_t threadLock = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t normLock = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t thread_cond = PTHREAD_COND_INITIALIZER;


///////////////////////
//MAIN PAGERANK WORK //
///////////////////////

/**
 * pagerank worker function
 * @param id worker id
 */
// extern inline void * worker(void * id){
//     int threadID = (int )id;
//     int lnpages = gnpages[threadID];
//     register int ** local_nodes = nodes[threadID];
//     int *links;
//     int did_not_wait;
//     int pn, cp, il, index, nlinks, limit;
//     register double rank;
//     register double localnorm = 1;
//     pthread_barrier_wait(&threadSync);
//     while (*(long int*)&epsilon < *(long int*)&norm){
//     // while (norm > epsilon){
//         localnorm = 0;
//         for(pn = 0; pn < lnpages; pn++){
//             links = local_nodes[pn];
//             index = links[0];
//             nlinks = links[1];
//             limit = nlinks + 2;
//             rank = 0;
//             for(il = 2; il < limit; il++){
//                 cp = links[il];
//                 rank += (PrevRank[cp] * outlinks[cp]);
//             }
//             rank += jumpProb;
//             PageRank[index] = rank;
//             rank = rank - PrevRank[index];
//             localnorm += rank * rank;
//         }

//         pthread_mutex_lock(&normLock);
//         next_norm += localnorm;
//         pthread_mutex_unlock(&normLock);
//         did_not_wait = pthread_barrier_wait(&threadSync);
//         pthread_mutex_lock(&threadLock);
//         if(did_not_wait){
//         TempRank = PrevRank;
//         PrevRank = PageRank;
//         PageRank = TempRank;
//         norm = next_norm;
//         next_norm = 0;
//         }
//         pthread_mutex_unlock(&threadLock);
//         pthread_barrier_wait(&threadSync);
//     }
//     return NULL;
// }

extern inline void * worker(void * id){
    int threadID = (int )id;

    int lnpages = gnpages[threadID];
    register int ** local_nodes = nodes[threadID];
    int *links;
    int pn, cp, il, index, nlinks, limit, i;
    register double rank;
    register double localnorm = 1;
    while (*(long int*)&epsilon < *(long int*)&norm){
       localnorm = 0;
       for(pn = 0; pn < lnpages; pn++){
           links = local_nodes[pn];
           index = links[0];
           nlinks = links[1];
           limit = nlinks + 2;
           rank = 0;
           for(il = 2; il < limit; il++){
               cp = links[il];
               rank += (PrevRank[cp] * outlinks[cp]);
           }
           rank += jumpProb;
           PageRank[index] = rank;
           rank = rank - PrevRank[index];
           localnorm += rank * rank;
       }
        local_norms[threadID] = localnorm;
        pthread_mutex_lock(&threadLock);
        thread_done++;
        if(thread_done != gnthreads){
            pthread_cond_wait(&thread_cond, &threadLock);
        }else{
            TempRank = PrevRank;
            PrevRank = PageRank;
            PageRank = TempRank;
            thread_done = 0;
            norm = 0;
            for (i = 0; i < gnthreads; i++)
            norm += local_norms[i];
            pthread_cond_broadcast(&thread_cond);
        }
        pthread_mutex_unlock(&threadLock);
   }
   return NULL;
}

static inline void pagerank1(list *plist, int ncores, int npages, int nedges, double dampener)
{
    int **nodes = (int **)malloc(sizeof(int *) * npages);
    double *PageRank = (double *)malloc(sizeof(double) * npages);
    register double *PrevRank = (double *)malloc(sizeof(double) * npages);
    register double *outlinks = (double *)malloc(sizeof(double) * npages);
    double *TempRank = NULL;
    double epsilon   = EPSILON * EPSILON;
    double jumpProb  = ((1.0 - dampener) / npages);
    double baseProb  = 1.0 / npages;
    double norm = 0;
    register double noutlinks, rank;
    register int index, il, pn, count, nlinks, limit, cp;
    int p = 0;
    register int * links;

    register node *curr = plist->head;
    node *c;

    while (curr)
    {
        index = curr->page->index;
        noutlinks = dampener/curr->page->noutlinks;
        outlinks[index] = noutlinks;
        if (curr->page->inlinks)  // Not dangling page
        {
            c = curr->page->inlinks->head;
            nlinks = curr->page->inlinks->length;
            links = (int *)malloc(sizeof(int)*(nlinks+2));
            links[0] = index;
            links[1] = nlinks;
            rank = 0;
            count = 2;
            while (c)
            {
                links[count] = c->page->index;
                rank += (baseProb * (1.0/c->page->noutlinks));
                c = c->next;
                count++;
            }
            nodes[p] = links;
            rank = jumpProb + (dampener * rank);
            PrevRank[index] = rank;
            PageRank[index] = rank;
            p++;
        }
        else   // Dangling page
        {
            PrevRank[index] = jumpProb;
            PageRank[index] = jumpProb;
        }
        curr = curr->next;
    }
    npages = p;
    // printf("%d\n", npages);

    // print_nodes(plist);6
    do{
        TempRank = PrevRank;
        PrevRank = PageRank;
        PageRank = TempRank;
        norm = 0;
        for(pn = 0; pn < npages; pn++){
            links = nodes[pn];
            index = links[0];
            nlinks = links[1];
            limit = nlinks + 2;
            rank = 0;
            for(il = 2; il < limit; il++){
                cp = links[il];
                rank += (PrevRank[cp] * outlinks[cp]);
            }
            rank += jumpProb;
            PageRank[index] = rank;
            rank = rank - PrevRank[index];
            norm += rank * rank;
        }
    } while (*(long int*)&epsilon < *(long int*)&norm);



    curr = plist->head;
    // printf("n: %f\n", norm);
    while (curr)
    {
        printf("%s %.4f\n", curr->page->name, PageRank[curr->page->index]);
        curr = curr->next;
    }


    if (PageRank){
        free(PageRank);
    }
    if (PrevRank){
        free(PrevRank);
    }
    if (outlinks){
        free(outlinks);
    }
    if(nodes){
        for(pn = 0; pn < npages; pn++){
            if(nodes[pn]){
                free(nodes[pn]);
            }
        }
        free(nodes);
    }
}

/**
 * Main pagerank function
 * @param plist    list of pages
 * @param ncores   number of cores
 * @param npages   number of pages
 * @param nedges   number of edges
 * @param dampener dampener
 */
static inline void pagerank2(list *plist, int ncores, int npages, int nedges, double dampener)
{
    double baseProb  = 1.0 / npages;
    int nthreads = 1;
    // if (){
    nthreads = ncores;
    // }
    // nthreads = 1;

    gnthreads = nthreads;
    // nodes = (int **)malloc(sizeof(int *) * npages);
    // PageRank = (double *)malloc(sizeof(double) * npages);
    // PrevRank = (double *)malloc(sizeof(double) * npages);
    // outlinks = (double *)malloc(sizeof(double) * npages);
    int ** lnodes[nthreads];
    nodes = &lnodes[0];
    double storage_array[3*npages];
    int edges[nthreads];
    PageRank = &storage_array[0];
    PrevRank = &storage_array[npages];
    outlinks = &storage_array[2*npages];
    jumpProb  = ((1.0 - dampener) / npages);
    // pthread_t * threads = (pthread_t *)malloc(sizeof(pthread_t)*nthreads);
    pthread_t threads[nthreads];
    double tln[nthreads];
    local_norms = &tln[0];
    int lnpages[nthreads];
    double noutlinks, rank;
    int index, pn, count, nlinks, i, mi, min;
    int * links;
    for (i = 0; i < nthreads; i++){
        nodes[i] = (int **)malloc(sizeof(int *)*npages);
        edges[i] = 0;
        lnpages[i] = 0;
    }
    pthread_barrier_init(&threadSync, NULL, nthreads);

    node *curr = plist->head;
    node *c;

    while (curr)
    {
        index = curr->page->index;
        noutlinks = dampener/curr->page->noutlinks;
        outlinks[index] = noutlinks;
        if (curr->page->inlinks)  // Not dangling page
        {
            c = curr->page->inlinks->head;
            nlinks = curr->page->inlinks->length;
            // int links[nlinks+2];
            // links = &lt[0];
            links = malloc(sizeof(int)*(nlinks+2));
            links[0] = index;
            links[1] = nlinks;
            rank = 0;
            count = 2;
            while (c)
            {
                links[count] = c->page->index;
                rank += (baseProb * (1.0/c->page->noutlinks));
                c = c->next;
                count++;
            }

            min = edges[0];
            mi = 0;
            for(i = 1; i < nthreads; i++){
                if (edges[i] < min){
                    min = edges[i];
                    mi = i;
                }
            }
            nodes[mi][lnpages[mi]] = links;
            edges[mi] += nlinks;
            lnpages[mi] += 1;
            rank = jumpProb + (dampener * rank);
            PrevRank[index] = rank;
            PageRank[index] = rank;
        }
        else   // Dangling page
        {
            PrevRank[index] = jumpProb;
            PageRank[index] = jumpProb;
        }
        curr = curr->next;
    }
    gnpages = &lnpages[0];
    // printf("%d\n", npages);

    // print_nodes(plist);
    // printf("next norm %f \n", next_norm);

    for (i = 0; i<nthreads; i++){
        // printf("thread: %d has %d pages and %d edges\n", i, lnpages[i], edges[i]);
        pthread_create(&threads[i], NULL, worker, (void *)i);
    }
    for (i = 0; i<nthreads; i++){
        pthread_join(threads[i], NULL);
    }

    curr = plist->head;
    // printf("n: %f\n", norm);
    while (curr)
    {
        printf("%s %.4f\n", curr->page->name, PrevRank[curr->page->index]);
        curr = curr->next;
    }

    pthread_mutex_destroy(&threadLock);
    pthread_barrier_destroy(&threadSync);
    pthread_cond_destroy(&thread_cond);
    // if (PageRank){
    //     free(PageRank);
    // }
    // if (PrevRank){
    //     free(PrevRank);
    // }
    // if (outlinks){
    //     free(outlinks);
    // }
    // if (threads){
    //     free(threads);
    // }
    if(nodes){
        for (i = 0; i < nthreads; i++){
            if(nodes[i]){
                for(pn = 0; pn < lnpages[i]; pn++){
                    if(nodes[i][pn]){
                        free(nodes[i][pn]);
                    }

                }
            free(nodes[i]);
            }
        }
    }
}

static inline void pagerank(list *plist, int ncores, int npages, int nedges, double dampener){
    if(ncores == 1 || nedges < 1000){
        pagerank1(plist, ncores, npages, nedges, dampener);
    }else{
        pagerank2(plist, ncores, npages, nedges, dampener);
    }
    // pagerank2(plist, ncores, npages, nedges, dampener);
}

////////////////
//NOT MY CODE //
////////////////

/* DO NOT MODIFY BELOW THIS POINT */
int main(void)
{
    /*
     * DO NOT MODIFY THE MAIN FUNCTION OR HEADER FILE
     */

    double damping_factor;
    list *plist = NULL;
    int ncores, npages, nedges;

    /* read the input into the appropriate variables and populate the list of pages */
    read_input(&plist, &ncores, &npages, &nedges, &damping_factor);

    /* run pagerank and print the results */
    pagerank(plist, ncores, npages, nedges, damping_factor);

    /* clean up the memory used by the list of pages */
    page_list_destroy(plist);

    return 0;
}
