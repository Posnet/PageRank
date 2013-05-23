/////////////
//Includes //
/////////////

#include <stdlib.h>
#include <unistd.h>
#include <pthread.h>

#include "pagerank.h"

///////////////////////
//MAIN PAGERANK WORK //
///////////////////////

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

    // print_nodes(plist);
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
