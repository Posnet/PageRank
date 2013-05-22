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
SuperPage *SuperPage_create(page *p);
Node *Node_create(void);
Node *List_create(void);
Node *insert_after(Node *n, int index);
Node *delete_node(Node *n);
Node *get_next(void);
void *scalloc(size_t num, size_t size);
void *smalloc(size_t size);
void SuperPage_free(SuperPage *p);
void Free_List(Node *head);
void reset_next(void);
void init_globals(int ncores, int npages, int nedges, double dampener);
void process_data(list *plist);
double process_node(Node *n);
void tick(void);
void free_all(void);
void print_nodes(list *plist);
void edge_cases(list *plist, int ncores, int nedges);
void pagerank(list *plist, int ncores, int npages, int nedges, double dampener);
double diff_squared(double curr, double prev);
int main(void);
int check_convergance(void);


//////////////////////
//Stuct Definitions //
//////////////////////

struct SuperPage
{
    int index;
    int noutlinks;
    double partialRank;
    Node *inlinks;
};

struct Node
{
    Node *prev;
    Node *next;
    int type; //  tail <  0  = normal head > 0
    int elem;
};

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

/**
 * Initialises and Allocated memory for a new
 * Superpage
 * @param  p the page being copied
 * @return   new super page pointer
 */
SuperPage *SuperPage_create(page *p)
{
    SuperPage *newPage = (SuperPage *)smalloc(sizeof(SuperPage));
    newPage->index = p->index;
    newPage->noutlinks = p->noutlinks;
    newPage->partialRank = 0;
    newPage->inlinks = NULL;
    return newPage;
}

/**
 * frees a given SuperPage's memory
 * @param p the page to be freed
 */
void SuperPage_free(SuperPage *p)
{
    if (p)
    {
        Free_List(p->inlinks);
        free(p);
    }
}


//////////////////////////////
//Linked List Access Methods//
//////////////////////////////

/**
 * Initialises and allocates
 * a new Node
 * @param  void
 * @return      new Node pointer
 */
Node *Node_create(void)
{
    Node *n = (Node *)smalloc(sizeof(Node));
    n->elem = -1;
    n->next = NULL;
    n->prev = NULL;
    n->type = 0;
    return n;
}

/**
 * Creates a new doubly linked list
 * @param  void
 * @return      returns the new virtual head node
 */
Node *List_create(void)
{
    Node *head = Node_create();
    head->type = 1;
    return head;
}

/**
 * Inserts a node after a given node
 * @param  n     Node to be inserted after
 * @param  index the elem value of the new Node
 * @return       returns a pointer to the new Node
 */
Node *insert_after(Node *n, int index)
{
    Node *newNode = Node_create();
    newNode->elem = index;
    if (n)
    {
        if (n->next)
        {
            n->next->prev = newNode;
            newNode->next = n->next;
            n->next = newNode;
            newNode->prev = n;
        }
        else
        {
            n->next = newNode;
            newNode->prev = n;
        }
    }
    return newNode;
}

/**
 * Deletes a given Node and frees it
 * @param  n Node to be deleted
 * @return   returns the Node after the one being
 * deleted, if at the end of the list returns NULL
 * If the node is the head of the list return NULL
 */
Node *delete_node(Node *n)
{
    if (n->type)
    {
        return NULL;
    }
    Node *res = NULL;
    n->prev->next = n->next;
    if (n->next)
    {
        res = n->next;
        n->next->prev = n->prev;
    }
    free(n);
    return res;
}

/**
 * Frees all the nodes and the virtual head
 * of a given list
 * @param head head of list to be freed
 */
void Free_List(Node *head)
{
    Node *n = head->next;
    Node *t = NULL;
    while (n)
    {
        t = n->next;
        free(n);
        n = t;
    }
    free(head);
}

///////////////////////////////
//PAGE RANK GLOBAL VARIABLES //
///////////////////////////////

double *PageRank;
double *PrevRank;
double *TempRank;
double * thread_norm;
double norm;
int *hasConverged;
int nthreads;
double damp;
double pages;
double cores;
double edges;
double jumpProb;
double baseProb;
double epsilon = EPSILON * EPSILON;

pthread_mutex_t nextLock = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t delnLock = PTHREAD_MUTEX_INITIALIZER;
pthread_t * threads;
SuperPage **superPages;
Node *Head;
Node *globalCurr;


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
    pthread_mutex_lock(&nextLock);
    Node *prev = globalCurr;
    if(globalCurr){
        globalCurr = globalCurr->next;
    }
    pthread_mutex_unlock(&nextLock);
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
    norm = 0;
    nthreads = ncores;
    Head = List_create();
    jumpProb = ((1.0 - damp) / pages);
    baseProb = 1.0 / npages;

    threads = (pthread_t *)smalloc(sizeof(pthread_t)*nthreads);
    PageRank = (double *)smalloc(sizeof(double) * pages);
    PrevRank = (double *)smalloc(sizeof(double) * pages);
    TempRank = NULL;
    hasConverged = (int *)smalloc(sizeof(int) * pages);
    superPages = (SuperPage **)smalloc(sizeof(SuperPage *)*pages);
    thread_norm = (double *)smalloc(sizeof(double)*nthreads);
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
    norm = 0;
    Node *current = Head;
    Node *ptr;
    node *curr = plist->head;
    node *c;
    int index;
    double rank;
    while (curr)
    {
        index = curr->page->index;
        if (curr->page->inlinks)  // Not dangling page
        {
            SuperPage *newPage = SuperPage_create(curr->page);
            superPages[index] = newPage;
            current = insert_after(current, index);
            newPage->inlinks = List_create();
            ptr = newPage->inlinks;
            c = curr->page->inlinks->head;
            rank = 0;
            while (c)
            {
                insert_after(ptr, c->page->index); //Possible remove has converged here
                rank += (baseProb / c->page->noutlinks);
                c = c->next;
            }
            rank = jumpProb + (damp * rank);
            PageRank[index] = rank;
            PrevRank[index] = rank;
            hasConverged[index] = 0;
            norm += diff_squared(baseProb, rank);
        }
        else   // Dangling page
        {
            PageRank[index] = jumpProb;
            PrevRank[index] = jumpProb;
            superPages[index] = NULL;
            hasConverged[index] = curr->page->noutlinks;
            norm += diff_squared(baseProb, jumpProb);
        }
        curr = curr->next;
    }
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
double process_node(Node *n)
{
    SuperPage *p = superPages[n->elem];
    double rank = 0;
    int index;
    Node *curr = p->inlinks->next;
    while (curr)
    {
        index = curr->elem;
        if (hasConverged[index])
        {
            p->partialRank += PrevRank[index]/hasConverged[index];
            curr = delete_node(curr);
        }
        else
        {
            rank += (PrevRank[index] / superPages[index]->noutlinks);
            curr = curr->next;
        }
    }
    index = n->elem;
    rank += p->partialRank;
    rank = jumpProb + (damp * rank);
    PageRank[index] = rank;
    double result = diff_squared(rank, PrevRank[index]);
    if (!(p->inlinks->next))
    {
        pthread_mutex_lock(&delnLock);
        hasConverged[index] = (rank/p->noutlinks);
        PrevRank[index] = rank;
        delete_node(n);
        pthread_mutex_unlock(&delnLock);
    }
    return result;
}

void * worker(void * args){
    Node * curr = NULL;
    int threadID = (int)args;
    thread_norm[threadID] = 0;
    while ((curr = get_next()))
    {
        thread_norm[threadID] += process_node(curr);
    }
    return NULL;
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
    reset_next();
    double total = 0;
    for(int i = 0; i < nthreads; i++){
        total += thread_norm[i];
    }
    norm = total;
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
    reset_next();
    while (check_convergence())
    {
        for(int i = 0; i < nthreads; i++){
        pthread_create(&threads[i], NULL, worker, (void *)i);
        }
        for(int i = 0; i < nthreads; i++){
        pthread_join(threads[i], NULL);
        }
        tick();
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
