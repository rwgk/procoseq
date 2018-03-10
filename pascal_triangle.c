#include <stdio.h>
#include <stdlib.h>


#include "longer.h"


static const char *progn = "pt";


#define progerror(message) p_rogerror(message, __LINE__)

static void p_rogerror(char *message, int source_code_line)
{
  fflush(stdout);
  fprintf(stderr, "%s(%d): %s\n", progn, source_code_line, message);
  exit(1);
}


#define NotEnoughCore() N_otEnoughCore(__LINE__)

static void N_otEnoughCore(int source_code_line)
{
  p_rogerror("Not enough core", source_code_line);
}


#define OutOfIntegerRange() O_utOfIntegerRange(__LINE__)

static void O_utOfIntegerRange(int source_code_line)
{
  p_rogerror("Out of integer range", source_code_line);
}


static void usage(void)
{
  fprintf(stderr, "usage: %s max_m\n", progn);
  exit(1);
}


int main(int argc, char *argv[])
{
  int    max_m, m, n;
  llint  *prevln, *nextln;


  if (argc < 2)
    max_m = 10;
  else
  {
    if (argc != 2) usage();
    if (argv[1][0] == '-') usage();
    if (sscanf(argv[1], "%d", &max_m) != 1) usage();
    if (max_m < 0) usage();
  }

  prevln = NULL;

  for (m = 0; m <= max_m; m++)
  {
        nextln = malloc((m / 2 + 1) * sizeof (*nextln));
    if (nextln == NULL)
      NotEnoughCore();

    nextln[0] = 1;

    for (n = 1; n <= (m - 1) / 2; n++)
    {
      nextln[n] = prevln[n - 1] + prevln[n];

      if (   nextln[n] - prevln[n] != prevln[n - 1]
          || nextln[n] < prevln[n - 1]
          || nextln[n] < prevln[n])
        OutOfIntegerRange();
    }

    if (m && m % 2 == 0)
    {
      nextln[n] = prevln[n - 1] + prevln[n - 1];

      if (   nextln[n] - prevln[n - 1] != prevln[n - 1]
          || nextln[n] < prevln[n - 1])
        OutOfIntegerRange();
    }

    if (prevln) free(prevln);

    prevln = nextln;
    nextln = NULL;

    fprintf(stdout, "%d:", m);

    for (n = 0; n <= m / 2; n++)
      fprintf(stdout, " %s", llint_to_str(prevln[n], NULL, 0));

    putc('\n', stdout);
  }

  exit(0);
  return 0;
}
