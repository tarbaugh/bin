#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv) {
    FILE *fp;

	fp = fopen(argv[1], "r");


    while(1) {
            while(fscanf(fp, "%lf%lf%lf", &t, &theta, &w) != EOF) {
        }
    }
}
