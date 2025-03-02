#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char *argv[]) {
    /* Variable declarations must be first */
    char *goal;
    char *filename;

    /* Ensure correct number of arguments */
    if (argc != 3) {
        printf("An Error Has Occurred\n");
        return 1;
    }

    /* Assign values AFTER all declarations */
    goal = argv[1];
    filename = argv[2];
}
