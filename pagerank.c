/////////////
//Includes //
/////////////

#include <stdlib.h>
#include <unistd.h>
#include <pthread.h>

#include "pagerank.h"

/////////////////////////////
//GENERIC HELPER FUNCTIONS //
/////////////////////////////

/**
 * A memory safe version of calloc
 * Throws error if allocation not successful
 * @param num
 * @param size
 */
void *scalloc(size_t num, size_t size)
{
    void *res;
    if ((res = calloc(num, size)) == (void *)0)
    {
        exit(EXIT_FAILURE);
    }
    return res;
}

/**
 * A memory safe version of malloc
 * Throws error if allocation not successful
 * @param size
 */
void *smalloc(size_t size)
{
    void *res;
    if ((res = malloc(size)) == (void *)0)
    {
        exit(EXIT_FAILURE);
    }
    return res;
}

///////////////////////////////
//PAGE RANK GLOBAL VARIABLES //
///////////////////////////////

double *PageRank;
double *PrevRank;
double *TempRank;
double *hasConverged;
int *outlinks;
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

    PageRank = (double *)smalloc(sizeof(double) * pages);
    PrevRank = (double *)smalloc(sizeof(double) * pages);
    TempRank = NULL;
    hasConverged = (double *)smalloc(sizeof(double) * pages);
    outlinks = (int *)smalloc(sizeof(int) * pages);
    nodes = (int **)smalloc(sizeof(int *) * pages);
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
    int noutlinks;
    int count;
    int nlinks;
    int *links;
    double rank;
    while (curr)
    {
        index = curr->page->index;
        noutlinks = c->page->noutlinks;
        if (curr->page->inlinks)  // Not dangling page
        {
            c = curr->page->inlinks->head;
            nlinks = curr->page->inlinks->length;
            outlinks[index] = noutlinks;
            links = (int *)smalloc(sizeof(int)*(nlinks+2));
            links[0] = index;
            links[1] = nlinks;
            rank = 0;
            count = 2;
            while (c)
            {
                links[count] = c->page->index; //Possible remove has converged here
                rank += (baseProb / noutlinks);
                c = c->next;
                count++;
            }
            nodes[p] = links;
            rank = jumpProb + (damp * rank);
            PrevRank[index] = rank;
            hasConverged[index] = 0;
            p++;
        }
        else   // Dangling page
        {
            PrevRank[index] = jumpProb;
            hasConverged[index] = jumpProb/noutlinks;
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
double process_node(int page)
{
    int * links = nodes[page];
    int index = links[0];
    int nlinks = links[1];
    int c;
    double t;
    double rank = 0;
    for(int i = 0; i < nlinks; i++){
        c = links[i];
        if((t = hasConverged[c])){
            rank += t;
        }else{
            rank += (PrevRank[c]/outlinks[c]);
        }
    }
    rank *= damp;
    rank += jumpProb;
    rank = PrevRank[index] - rank;
    return rank*rank;
}

/**
 * Single parse through global node list.
 * Process every node.
 * rotate current and previous PageRanks.
 */
void tick(void)
{
    norm = 0;
    for(int i = 0; i < realPages; i++){
        norm += process_node(i);
    }
}

/**
 * Free all allocated data.
 */
void free_all(void)
{
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
    } while ((*(int*)&epsilon < *(int*)&norm));
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
