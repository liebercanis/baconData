#include <stdio.h>  /* printf, fopen */
#include <stdlib.h> /* exit, EXIT_FAILURE */

int main(int argc, char* argv[])
{
  printf(" executing %s \n" , argv[1]);
  if(argc<1) {
    printf(" usage: test <tag>  \n ");
    exit(0);
   }
  printf(" test done for  %s \n" , argv[1]);
  exit(0);
}
