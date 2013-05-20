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
    //Sanity Checks
    // assert(plist->head != null);
    // assert(npages > 0);
    // assert(ncores > 0);
    // assert(dampener > 0 && dampener < 1);

    //Initial inits
    double jump_prob = ((1.0 - dampener)/npages);
    double base_conv = 1.0/npages;
    double norm = 0;
    double trank = 0;
    double temp_store = 0;
    int has_non_constants = 0;

    //saves time if only one node
    // if (npages == 1)
    // {
    //     printf("%s %.4f\n", plist->head->page->name, jump_prob);
    //     return;
    // }

    //more inits
    double * has_converged = (double *)smalloc((sizeof(double *))*npages);
    double * page_rank = (double *)smalloc((sizeof(double *))*npages);
    double * prev_rank = (double *)smalloc((sizeof(double *))*npages);
    double * temp_rank;
    list * page_list = page_list_create();
    node * curr = plist->head;

    node * prev = NULL;
    page * cpage = NULL;
    node * c = NULL;
    node * p = NULL;

    //Initial loop, explicit for speed
    while(curr)
    {
        cpage = curr->page;
        if (cpage->inlinks) //If the page isn't a dangling page.
        {
            has_converged[cpage->index] = 0;
            page_list_add_end(page_list, cpage);
            c = cpage->inlinks->head;
            trank = 0;
            while(c)
            {
                trank += base_conv/c->page->noutlinks;
                //printf("%s, %f, %i\n", c->page->name, trank, c->page->noutlinks);
                p = c;
                c = p->next;
            }
            trank = jump_prob + trank * dampener;
            prev_rank[cpage->index] = trank;
            page_rank[cpage->index] = trank;
            temp_store = base_conv - trank;
            norm += temp_store * temp_store;

        }
        else // if the page is a dangling page.
        {
            prev_rank[cpage->index] = jump_prob;
            page_rank[cpage->index] = jump_prob;
            has_converged[cpage->index] = (jump_prob/cpage->noutlinks);
            temp_store = base_conv - jump_prob;
            norm += temp_store*temp_store;
        }
        prev = curr;
        curr = prev->next;
    }

    //Convergence loop
    while(norm > EPSILON * EPSILON)
    {
    //         //Printing loop
    // printf("norm: %f, eps: %f\n", norm, EPSILON*EPSILON);
    // curr = plist->head;
    // while(curr){
    //     printf("%s %.4f\n", curr->page->name, jump_prob + prev_rank[curr->page->index]);
    //     prev = curr;
    //     curr = prev->next;
    // }

        curr = page_list->head;
        prev = NULL;
        norm = 0;
        while(curr)
        {
            cpage = curr->page;
            has_non_constants = 0;
            c = cpage->inlinks->head;
            trank = 0;
            while(c)
            {
                if (has_converged[c->page->index] > 0)
                {
                    trank += has_converged[c->page->index];
                }else
                {
                    trank += (prev_rank[c->page->index]/c->page->noutlinks);
                    has_non_constants = 1;
                }
                p = c;
                c = p->next;
            }
            trank = jump_prob + trank * dampener;
            if (has_non_constants == 0 ){ //|| prev_rank[cpage->index] == trank
                has_converged[cpage->index] = (trank)/cpage->noutlinks;
                temp_store = (trank - prev_rank[cpage->index]);
                norm += temp_store * temp_store;
                prev_rank[cpage->index] = trank;
                page_rank[cpage->index] = trank;
                if (page_list->head == curr)
                {
                    page_list->head = curr->next;
                    free(curr);
                    if(page_list->head)
                    {
                        curr = page_list->head->next;
                    }else
                    {
                        curr = NULL;
                    }
                }else
                {
                    prev->next = curr->next;
                    free(curr);
                    curr = prev->next;
                }
            }else
            {
                page_rank[cpage->index] = trank;
                temp_store = (trank - prev_rank[cpage->index]);
                norm += temp_store * temp_store;
                prev = curr;
                curr = prev->next;
            }
        }
        temp_rank = prev_rank;
        prev_rank = page_rank;
        page_rank = temp_rank;

    }

    //Printing loop
    curr = plist->head;
    while(curr){
        printf("%s %.4f\n", curr->page->name, prev_rank[curr->page->index]);
        prev = curr;
        curr = prev->next;
    }


    // Free allocated memory
    if(has_converged){
        free(has_converged);
    }
    if(prev_rank){
        free(prev_rank);
    }
    if(page_rank){
        free(page_rank);
    }

    node* next;
    node* current;
    if (page_list == NULL) /* null page list */
      return;

    current = page_list->head;
    while (current != NULL)
    {
      next = current->next;
      //page_destroy(current->page);
      // if (current->page){
        // free(current->page);
      // }
      free(current);
      current = next;
    }
    free(page_list);
    // page_list_destroy(page_list);
    return;
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
