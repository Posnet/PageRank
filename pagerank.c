/////////////
//Includes //
/////////////

#include <stdlib.h>
#include <unistd.h>
#include <pthread.h>

#include "pagerank.h"

///////////////////////////////
//PAGE RANK GLOBAL VARIABLES //
///////////////////////////////

double *PageRank;
double *PrevRank;
double *TempRank;
double *hasConverged;
double *partialSum;
double *outlinks;
int **nodes;
int realPages;
double damp;
double pages;
double cores;
double edges;
double jumpProb;
double baseProb;
double epsilon = EPSILON * EPSILON;
double norm;


///////////////////////////////
//PAGE RANK HELPER FUNCTIONS //
///////////////////////////////

/**
 * A global get next node method for the main
 * linked list of pagerank nodes. Stores itself as
 * the next node of the one being returned to prevent
 * deletion.
 * @param  void
 * @return      Next Node
 */
// Node *get_next(void)
// {
//     Node *prev = globalCurr;
//     if(globalCurr){
//         globalCurr = globalCurr->next;
//     }
//     return prev;
// }

/**
 * Resets the global node pointer to the head
 * of the list
 */
// void reset_next(void)
// {
//     globalCurr = Head->next;
// }

/**
 * Declare and allocate required global variables
 * @param ncores   number of cores
 * @param npages   number of pages
 * @param nedges   number of edges
 * @param dampener dampener ratio
 */
void init_globals(int ncores, int npages, int nedges, double dampener)
{
    //Init Globals
    damp = dampener;
    pages = npages;
    cores = ncores;
    edges = nedges;
    jumpProb = ((1.0 - damp) / pages);
    baseProb = 1.0 / npages;

    PageRank = (double *)malloc(sizeof(double) * pages);
    PrevRank = (double *)malloc(sizeof(double) * pages);
    TempRank = NULL;
    hasConverged = (double *)malloc(sizeof(double) * pages);
    partialSum = (double *)malloc(sizeof(double) * pages);
    outlinks = (double *)malloc(sizeof(double) * pages);
    nodes = (int **)malloc(sizeof(int *) * pages);
    norm = 0;
}

/**
 * Initialises the global PageRank linked list
 * Initialises the base values of PageRank
 * Does not add constant nodes to the linked list
 * generates the initial normal value
 * @param plist list of pages to be processed
 */
void process_data(list *plist)
{
    node *curr = plist->head;
    node *c;
    int index;
    int p = 0;
    double noutlinks;
    int count;
    int nlinks;
    int *links;
    double rank;
    while (curr)
    {
        index = curr->page->index;
        noutlinks = damp/curr->page->noutlinks;
        outlinks[index] = noutlinks;
        if (curr->page->inlinks)  // Not dangling page
        {
            c = curr->page->inlinks->head;
            nlinks = curr->page->inlinks->length;
            // printf("%f\n", outlinks[index]);
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
            rank = jumpProb + (damp * rank);
            PrevRank[index] = rank;
            PageRank[index] = rank;
            partialSum[index] = 0;
            hasConverged[index] = 0;
            p++;
        }
        else   // Dangling page
        {
            PrevRank[index] = jumpProb;
            PageRank[index] = jumpProb;
            hasConverged[index] = (jumpProb * noutlinks);
        }
        curr = curr->next;
    }
    realPages = p;
}

/**
 * Process and individual node.
 * If the node has only constant inlinks delete from main list
 * and store value in hasConverges array.
 * If individual inlink is constant, add its value to partial
 * Rank and delete in from the inlinks list.
 * Update global PageRank and norm.
 * @param n node to be processed
 */
void process_node(int p)
{
    int * links = nodes[p];
    int index = links[0];
    int nlinks = links[1];
    int c;
    double rank = 0;
    double ps = partialSum[index];
    for(int i = 2; i < nlinks+2; i++){
        c = links[i];
        if((hasConverged[c])){
            // rank += hasConverged[c];
            ps += hasConverged[c];
            links[i] = links[nlinks+1];
            links[1]--;
        }else{
            rank += (PrevRank[c] * outlinks[c]);
        }
    }
    rank = jumpProb + (rank + ps);
    partialSum[index] = ps;
    PageRank[index] = rank;
    if(links[1] == 0){
        free(nodes[p]);
        nodes[p] = nodes[realPages-1];
        realPages--;
        hasConverged[index] = rank * outlinks[index];
        PrevRank[index] = rank;
    }
    rank = rank - PrevRank[index];
    norm += rank * rank;
}

/**
 * Single parse through global node list.
 * Process every node.
 * rotate current and previous PageRanks.
 */
void tick(void)
{
    TempRank = PrevRank;
    PrevRank = PageRank;
    PageRank = TempRank;
    norm = 0;
    int tempPages = realPages;
    for(int i = 0; i < tempPages; i++){
        process_node(i);
    }
}

/**
 * Free all allocated data.
 */
void free_all(void)
{
    if (PageRank){
        free(PageRank);
    }
    if (PrevRank){
        free(PrevRank);
    }
    if (hasConverged){
        free(hasConverged);
    }
    if (outlinks){
        free(outlinks);
    }
    if(partialSum){
        free(partialSum);
    }
    if(nodes){
        for(int i = 0; i < realPages; i++){
            if(nodes[i]){
                free(nodes[i]);
            }
        }
        free(nodes);
    }
}

/**
 * helper function to print nodes in plist
 * uses global pagerank to find result values.
 * @param plist page_list to be printed
 */
void print_nodes(list *plist)
{
    node *curr = plist->head;
    // printf("n: %f\n", norm);
    while (curr)
    {
        printf("%s %.4f\n", curr->page->name, PageRank[curr->page->index]);
        curr = curr->next;
    }
    // printf("\n");
}

/**
 * accounts for simple edge cases to avoid unnecessary
 * calculations.
 * No Nodes
 * No edges
 * Single Node
 * @param plist  list of pages to be printed
 * @param npages number of pages in graph
 * @param nedges number of edges in graph
 */
// void edge_cases(list *plist, int npages, int nedges)
// {
//     //TODO
// }


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
void pagerank(list *plist, int ncores, int npages, int nedges, double dampener)
{
    //edge_cases(plist, ncores, npages);
    init_globals(ncores, npages, nedges, dampener);
    process_data(plist);
    // print_nodes(plist);
    do{
        tick();
        // print_nodes(plist);
    } while ((*(long int*)&epsilon < *(long int*)&norm) && realPages);
    print_nodes(plist);
    free_all();
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
