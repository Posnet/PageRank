#include <stdlib.h>
#include <unistd.h>
#include <pthread.h>

#include "pagerank.h"

typedef struct SuperPage SuperPage;
typedef struct Node Node;

struct SuperPage
{
    int index;
    int noutlinks;
    double partialRank;
    Node * inlinks;
};

struct Node
{
    Node * prev;
    Node * next;
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
static void *scalloc(size_t num, size_t size)
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
static void *smalloc(size_t size)
{
    void *res;
    if ((res = malloc(size)) == (void *)0)
    {
        exit(EXIT_FAILURE);
    }
    return res;
}

static double diff_squared(double curr, double prev){
    double temp = curr - prev;
    return temp*temp;
}

static SuperPage * SuperPage_create(page* p){
    SuperPage * newPage = (SuperPage *)smalloc(sizeof(SuperPage));
    newPage->index = p->index;
    newPage->noutlinks = p->noutlinks;
    newPage->partialRank = 0;
    Node * inlinks;
}

static void SuperPage_free(void){
    //TODO
}


//////////////////////////////
//Linked List Access Methods//
//////////////////////////////

static Node* Node_create(void){
    Node * n = (Node *)smalloc(sizeof(Node));
    n->elem = NULL;
    n->next = NULL:
    n->prev = NULL;
    int type = 0;
    return n;
}

static Node* List_create(void){
    node * head = Node_create();
    head->type = 1;
    return head;
}

static Node * insert_after(Node * n, int index){
    Node * newNode = Node_create();
    newNode->elem = index;
    if(n){
        if (n->next){
            n->next->prev = newNode;
            newNode->next = n->next;
            n->next = newNode;
            newNode->prev = n;
        }else{
            n->next = newNode;
            newNode->prev = n;
        }
    }
    return newNode;
}

static SuperPage * delete_node(Node * n){
    if (n->type){
        return NULL;
    }
    n->prev = n->next;
    n->next = n->prev;
    SuperPage * p = n->elem;
    free(n);
    return p;
}

static void Free_List(Node * head){
    Node * n = head->next;
    Node * t = NULL;
    while (n){
        t = n->next;
        free(n);
        n = t;
    }
    free(head);
}

///////////////////////////////
//PAGE RANK GLOBAL VARIABLES //
///////////////////////////////

double * PageRank;
double * PrevRank;
double * TempRank;
double * hasConverged;
SuperPage * superPages;
double damp;
double pages;
double cores;
double edges;
double jumpProb;
double baseProb;
double epsilon = EPSILON * EPSILON;
double norm;

Node * head;


///////////////////////////////
//PAGE RANK HELPER FUNCTIONS //
///////////////////////////////

void init_globals(int ncores, int npages, int nedges, double dampener){
    //Init Globals
    damp = dampener;
    pages = npages;
    cores = ncores;
    edges = nedges;
    head = List_create();
    jumpProb = ((1.0 - damp)/pages);
    baseProb = 1.0/npages;

    PageRank = (double *)smalloc(sizeof(double)*pages);
    PrevRank = (double *)smalloc(sizeof(double)*pages);
    TempRank = NULL;
    hasConverged = (double *)smalloc(sizeof(double)*pages);
    superPages = (SuperPage *)smalloc(sizeof(SuperPage)*pages);
    norm = 0;
}

void process_data(list* plist){
    norm = 0;
    Node * current = head;
    node * curr = plist->head;
    node * temp;
    node * c;
    node * t;
    Node * ptr;
    int index;
    while(curr){
        index = curr->page->index;
        if (curr->page->inlinks){ // Not dangling page
            SuperPage * newPage = SuperPage_create(curr->page);
            superPages = newPage->index;
            newPage->inlinks = List_create();
            ptr = newPage->inlinks;
            c = curr->page->inlinks->head;
            double rank =
            while(c){
                insert_after(ptr, c->index);
                t = c->next;
                c = t;
            }
        }else{ // Dangling Page
            PageRank[index] = jumpProb;
            PrevRank[index] = jumpProb;
            hasConverged[index] = (jumpProb/curr->page->noutlinks);
            norm += diff_squared(jumpProb, baseProb);
        }
        temp = curr;
        curr = temp->next;
    }
}

Node * process_node(Node * n){
    //TODO
    return n;
}

void tick(void){
    //TODO
}

void load_balance(void){
    //TODO
}

void free_all(void){
    //TODO
}

bool check_convergance(){
    //TODO
    return false;
}

void edge_cases(list* plist, int ncores, int nedges){
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
void pagerank(list* plist, int ncores, int npages, int nedges, double dampener)
{
    //edge_cases(plist, ncores, npages);
    init_globals(ncores, npages, nedges, dampener);
    process_data(plist);
    //load_balance();
    while (check_convergance()){
        tick();
    }
    print_nodes();
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
