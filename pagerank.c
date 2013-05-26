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
const double epsilon   = EPSILON * EPSILON;
double jumpProb;
double norm = 1;
int *gnpages;
int gnthreads;
int thread_done = 0;
double *local_norms;

pthread_barrier_t threadSync;
pthread_mutex_t threadLock = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t normLock = PTHREAD_MUTEX_INITIALIZER;


///////////////////////
//MAIN PAGERANK WORK //
///////////////////////

static inline void *worker(void *id)
{
    int threadID = (int )id;
    int lnpages = gnpages[threadID];
    int **local_nodes = nodes[threadID];
    register int *links;
    int pn, cp, il, index, nlinks, limit;
    int localgnthreads = gnthreads;
    register const double ljumpprob = jumpProb;
    // register double lepsilon = epsilon;
    register double rank;
    register double localnorm;

    while (*(long int *)&epsilon < * (long int *)&norm)
    // while (norm > lepsilon)
    {
        localnorm = 0;
        for (pn = 0; pn < lnpages; pn++)
        {
            links = local_nodes[pn];
            index = links[0];
            nlinks = links[1];
            limit = nlinks + 2;
            rank = 0;
            for (il = 2; il < limit; il++)
            {
                cp = links[il];
                rank += (PrevRank[cp] * outlinks[cp]);
            }
            rank += ljumpprob;
            PageRank[index] = rank;
            rank = rank - PrevRank[index];
            localnorm += rank * rank;
        }

        local_norms[threadID] = localnorm;
        pthread_mutex_lock(&threadLock);
        thread_done++;

        if (thread_done == localgnthreads)
        {
            TempRank = PrevRank;
            PrevRank = PageRank;
            PageRank = TempRank;
            thread_done = 0;
            norm = 0;
            for (il = 0; il < localgnthreads; il++)
                norm += local_norms[il];
        }

        pthread_mutex_unlock(&threadLock);
        pthread_barrier_wait(&threadSync);
    }
    return NULL;
}


/**
 * Main pagerank function
 * @param plist    list of pages
 * @param ncores   number of cores
 * @param npages   number of pages
 * @param nedges   number of edges
 * @param dampener dampener
 */
static inline void pagerank(list *plist, int ncores, int npages, int nedges, double dampener)
{
    const double baseProb  = 1.0 / npages;
    int nthreads = 1;
    nthreads = ncores;

    gnthreads = nthreads;
    int **lnodes[nthreads];
    nodes = &lnodes[0];
    double storage_array[3 * npages];
    int edges[nthreads];
    PageRank = &storage_array[0];
    PrevRank = &storage_array[npages];
    outlinks = &storage_array[2 * npages];
    jumpProb  = ((1.0 - dampener) / npages);
    pthread_t threads[nthreads];
    double tln[nthreads];
    local_norms = &tln[0];
    int lnpages[nthreads];
    register double noutlinks, rank;
    register int index, pn, count, nlinks, i, mi, min;
    register int *links;
    for (i = 0; i < nthreads; i++)
    {
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
        noutlinks = dampener / curr->page->noutlinks;
        outlinks[index] = noutlinks;
        if (curr->page->inlinks)  // Not dangling page
        {
            c = curr->page->inlinks->head;
            nlinks = curr->page->inlinks->length;
            // int links[nlinks+2];
            // links = &lt[0];
            links = malloc(sizeof(int) * (nlinks + 2));
            links[0] = index;
            links[1] = nlinks;
            rank = 0;
            count = 2;
            while (c)
            {
                links[count] = c->page->index;
                rank += (baseProb / c->page->noutlinks);
                c = c->next;
                count++;
            }

            min = edges[0];
            mi = 0;
            for (i = 1; i < nthreads; i++)
            {
                if (edges[i] < min)
                {
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

    for (i = 0; i < nthreads; i++)
    {
        nodes[i] = (int**)realloc(nodes[i], (sizeof(int *)*lnpages[i]));
        pthread_create(&threads[i], NULL, worker, (void *)i);
    }
    for (i = 0; i < nthreads; i++)
    {
        pthread_join(threads[i], NULL);
    }

    curr = plist->head;
    while (curr)
    {
        printf("%s %.4f\n", curr->page->name, PrevRank[curr->page->index]);
        curr = curr->next;
    }

    pthread_mutex_destroy(&threadLock);
    pthread_barrier_destroy(&threadSync);

    if (nodes)
    {
        for (i = 0; i < nthreads; i++)
        {
            if (nodes[i])
            {
                for (pn = 0; pn < lnpages[i]; pn++)
                {
                    if (nodes[i][pn])
                    {
                        free(nodes[i][pn]);
                    }

                }
                free(nodes[i]);
            }
        }
    }
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
