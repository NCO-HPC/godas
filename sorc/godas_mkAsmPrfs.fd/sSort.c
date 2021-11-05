#include <stdlib.h>

int ssort__(slist)
char *slist;
{
  char cmd[30];

  sprintf(cmd, "/bin/sort -o %s %s",slist,slist);
  system(cmd);
  return(0);
}
