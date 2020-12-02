#include <time.h>

/* retourne le temps CPU depuis le début de l'execution du programme */
double mytimer(void){
    return (double) clock() / CLOCKS_PER_SEC;
}
