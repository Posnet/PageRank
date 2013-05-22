/////////////
//Includes //
/////////////

#include <stdlib.h>
#include <unistd.h>
#include <pthread.h>

#include "pagerank.h"

/////////////////////
//Type Definitions //
/////////////////////

typedef struct SuperPage SuperPage;
typedef struct Node Node;

////////////////////////
//Function Prototypes //
////////////////////////

/* pagerank.c */



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

/**
 * Simple function to calculate the square of
 * the difference of 2 number
 * @param  curr current rank
 * @param  prev previous rank
 * @return      difference squared
 */
double diff_squared(double curr, double prev)
{
    double temp = curr - prev;
    return temp * temp;
}





///////////////////////////////
//PAGE RANK GLOBAL VARIABLES //
///////////////////////////////

double *PageRank;
double *PrevRank;
double *TempRank;
int *hasConverged;
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
Node *get_next(void)
{
    Node *prev = globalCurr;
    if(globalCurr){
        globalCurr = globalCurr->next;
    }
    return prev;
}

/**
 * Resets the global node pointer to the head
 * of the list
 */
void reset_next(void)
{
    globalCurr = Head->next;
}

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
    hasConverged = (int *)smalloc(sizeof(int) * pages);
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
void process_node(Node *n)
{

}

/**
 * Single parse through global node list.
 * Process every node.
 * rotate current and previous PageRanks.
 */
void tick(void)
{
    Node *curr;
    norm = 0;
    TempRank = PrevRank;
    PrevRank = PageRank;
    PageRank = TempRank;
    reset_next();
    while ((curr = get_next()))
    {
        process_node(curr);
    }
}

/**
 * Free all allocated data.
 */
void free_all(void)
{
    if (PageRank)
    {
        free(PageRank);
    }
    if (PrevRank)
    {
        free(PrevRank);
    }
    if (hasConverged)
    {
        free(hasConverged);
    }
    if(Head){
        Free_List(Head);
    }
    if (superPages)
    {
        for (int i = 0; i < pages; i++)
        {
            if (superPages[i])
            {
                SuperPage_free(superPages[i]);
            }
        }
        free(superPages);
    }
}

/**
 * Check the convergence of the ranks, by
 * comparing the epsilon value to the normal
 * square difference sum
 * @return return 1 if not converges
 * otherwise 0
 */
int check_convergence()
{
    return (norm > epsilon);
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
void edge_cases(list *plist, int npages, int nedges)
{
    //TODO
}


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
    while (check_convergence())
    {
        tick();
        // print_nodes(plist);
    }
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
