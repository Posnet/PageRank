#include <stdlib.h>
#include <unistd.h>
#include <pthread.h>

#include "pagerank.h"

void * get_new_rank(void * pin){
}

void pagerank(list* plist, int ncores, int npages, int nedges, long double dampener)
{

}

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
