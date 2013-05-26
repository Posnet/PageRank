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

int **nodes;
double *PageRank;
double *PrevRank;
double *outlinks;
double *TempRank = NULL;
double epsilon   = EPSILON * EPSILON;
double jumpProb;
double norm = 1;
int gnpages;
int gnthreads;
pthread_barrier_t threadSync;
pthread_mutex_t threadLock = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t normLock = PTHREAD_MUTEX_INITIALIZER;


///////////////////////
//MAIN PAGERANK WORK //
///////////////////////

/**
 * pagerank worker function
 * @param id worker id
 */
extern inline void * worker(void * id){
    int pn, cp, il, index, nlinks, limit;
    int threadID = (int )id;
    int lnpages = gnpages;
    int nthreads = gnthreads;
    double localnorm = 1;
    register int *links;
    register double rank;
    pthread_barrier_wait(&threadSync);
    while (*(long int*)&epsilon < *(long int*)&norm){
        pthread_barrier_wait(&threadSync);
        pthread_mutex_lock(&threadLock);
        if(norm != 0){
        TempRank = PrevRank;
        PrevRank = PageRank;
        PageRank = TempRank;
        norm = 0;
        }
        pthread_mutex_unlock(&threadLock);
        pthread_barrier_wait(&threadSync);
        localnorm = 0;
        for(pn = threadID; pn < lnpages; pn += nthreads){
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
            localnorm += rank * rank;
        }

        pthread_mutex_lock(&normLock);
        norm += localnorm;
        pthread_mutex_unlock(&normLock);
        pthread_barrier_wait(&threadSync);
    }
    return NULL;
}

// extern inline void * worker(void * id)
// {
//   int pn, cp;
//   int nlinks;
//   int *links;
//   int limit;
//   int lid = (int)id;
//   signed int il;
//   int lnthreads;
//   register double rank;
//   double localnorm = 0;
//   int lnpages;
//   int index;

//   lnthreads = gnthreads;
//   lnpages = gnpages;
//   while ( 1 )
//   {
//     pthread_barrier_wait(&threadSync);
//     if (*(long int*)&norm < *(long int*)&epsilon)
//       break;
//     pthread_barrier_wait(&threadSync);
//     pthread_mutex_lock(&threadLock);
//     if(!norm)
//     {
//       TempRank = PrevRank;
//       PrevRank = PageRank;
//       PageRank = TempRank;
//       norm = 0;
//     }
//     pthread_mutex_unlock(&threadLock);
//     pthread_barrier_wait(&threadSync);
//     if ( lid < lnpages )
//     {
//       local_outlinks = outlinks;
//       local_prevrank = PrevRank;
//       localnodes = nodes;
//       pn = lid;
//       do
//       {
//         links = localnodes[pn];
//         index = links[0];
//         nlinks = links[1];
//         if ( nlinks > 0 )
//         {
//           il = 2;
//           limit = nlinks + 2;
//           do
//           {
//             cp = links[il];
//             rank += local_prevrank[cp] * local_outlinks[cp];
//             ++il;
//           }
//           while ( il < limit );
//         }
//         rank += local_jumprob;
//         PageRank[index] = rank;
//         rank = rank - local_prevrank[index];
//         localnorm += rank * rank;
//         pn += lnthreads;
//       }
//       while ( pn < lnpages );
//     }
//     pthread_mutex_lock(&normLock);
//     norm += localnorm;
//     pthread_mutex_unlock(&normLock);
//     // pthread_barrier_wait(&threadSync);
//   }
//   return 0;
// }

/**
 * Main pagerank function
 * @param plist    list of pages
 * @param ncores   number of cores
 * @param npages   number of pages
 * @param nedges   number of edges
 * @param dampener dampener
 */
extern inline void pagerank(list *plist, int ncores, int npages, int nedges, double dampener)
{
    double baseProb  = 1.0 / npages;
    int p = 0;
    int nthreads = 1;
    // if (){
    nthreads = ncores;
    // }
    // nthreads = 1;
    int edgeAvg = nedges;
    gnthreads = nthreads;
    // nodes = (int **)malloc(sizeof(int *) * npages);
    // PageRank = (double *)malloc(sizeof(double) * npages);
    // PrevRank = (double *)malloc(sizeof(double) * npages);
    // outlinks = (double *)malloc(sizeof(double) * npages);
    int * lnodes[npages];
    nodes = &lnodes[0];
    double storage_array[3*npages];
    PageRank = &storage_array[0];
    PrevRank = &storage_array[npages];
    outlinks = &storage_array[2*npages];
    jumpProb  = ((1.0 - dampener) / npages);
    // pthread_t * threads = (pthread_t *)malloc(sizeof(pthread_t)*nthreads);
    pthread_t threads[nthreads];
    register double noutlinks, rank;
    register int index, pn, count, nlinks, i;
    int * links;
    pthread_barrier_init(&threadSync, NULL, nthreads);

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
    gnpages = p;
    // printf("%d\n", npages);

    // print_nodes(plist);

    for (i = 0; i<nthreads; i++){
        pthread_create(&threads[i], NULL, worker, (void *)i);
    }
    for (i = 0; i<nthreads; i++){
        pthread_join(threads[i], NULL);
    }

    curr = plist->head;
    // printf("n: %f\n", norm);
    while (curr)
    {
        printf("%s %.4f\n", curr->page->name, PageRank[curr->page->index]);
        curr = curr->next;
    }

    pthread_mutex_destroy(&threadLock);
    pthread_barrier_destroy(&threadSync);
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
    // if(nodes){
    //     for(pn = 0; pn < npages; pn++){
    //         if(nodes[pn]){
    //             free(nodes[pn]);
    //         }
    //     }
    //     // free(nodes);
    // }
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
