#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>


#include "longer.h"
#include "inumber.h"
#include "fperiod.h"


#define fprintf (void) fprintf
#define fflush  (void) fflush
#define sprintf (void) sprintf


/* reconstruction of coordination sequences */

static char *progn = "abc";


#define progerror(message) p_rogerror(message, __LINE__)

static void p_rogerror(char *message, int source_code_line)
{
  fflush(stdout);
  fprintf(stderr, "%s(%d): %s\n", progn, source_code_line, message);
  exit(1);
}


#define InternalError() I_nternalError(__LINE__)

static void I_nternalError(const int source_code_line)
{
  p_rogerror("Internal Error", source_code_line);
}


#define NotEnoughCore() N_otEnoughCore(__LINE__)

static void N_otEnoughCore(int source_code_line)
{
  p_rogerror("Not enough core", source_code_line);
}


#define IllegalLine(fnin, lcount) I_llegalLine(fnin, lcount, __LINE__)

static void I_llegalLine(const char *fnin, int lcount, int source_code_line)
{
  fflush(stdout);
  fprintf(stderr, "%s(%d): Illegal line #%d",
    progn, source_code_line, lcount);
  if (fnin != NULL) fprintf(stderr, " in %s", fnin);
  putc('\n', stderr);
  exit(1);
}


static int fgetline(FILE *fpin, char s[], int size_s)
{
  int         last_s, c, i;


  last_s = size_s - 1;

  i = 0;

  while ((c = getc(fpin)) != EOF && c != '\n')
  {
#ifdef __MSDOS__
    if (c == 0x1A /* CtrlZ */) {
      ungetc(c, fpin);
      c = EOF;
      break;
    }
#endif

    if (i < last_s) s[i++] = (char) c;
  }

  s[i] = '\0';

  if (i == 0 && c == EOF)
    return 0;

  return 1;
}


static char firstnonblank(const char *s)
{
  while (isspace(*s)) s++;
  return *s;
}


static int str_icmp(const char *s, const char *t)
{
  char  cs, ct;


  while (*s || *t)
  {
    cs = toupper(*s++);
    ct = toupper(*t++);
    if (cs < ct) return -1;
    if (cs > ct) return  1;
  }

  return 0;
}


static void llint_PrintFraction(FILE *fpout, llint n, llint d)
{
  if (d < 0) {
    n *= -1;
    d *= -1;
  }

  fprintf(fpout, "%s",  llint_to_str(n, NULL, 0)); if (d != 1)
  fprintf(fpout, "/%s", llint_to_str(d, NULL, 0));
}


#ifdef FLOAT_ABC

static void abc(int k1, llint Nk1,
                int k2, llint Nk2,
                int k3, llint Nk3)
{
  int      Len;
  ldouble  a, b, c;


  Len = k2 - k1;

  a = (ldouble)(Nk3 - Nk2 - Nk2 + Nk1) / ((llint) 2 * Len * Len);
  b = (ldouble)(Nk3 - Nk2) / Len - a * (k3 + k2);
  c = (ldouble) Nk3 - (b + a * k3) * k3;

  fprintf(stdout, "%d %lld %d %lld %d %lld\n",
    k1, Nk1, k2, Nk2, k3, Nk3);

  fprintf(stdout, "abc(%d) ", k3 % Len);
  fprintf(stdout, "%.9f %.9f %.9f\n", (double) a, (double) b, (double) c);
}

#endif


typedef struct
  {
    int     k;
    llint  Nk;
  }
  T_k_Nk;


static int k_Nk_sort_low_k_first(const T_k_Nk *a, const T_k_Nk *b)
{
  if (a->k < b->k) return -1;
  if (a->k > b->k) return  1;
  return 0;
}


static void usage(void)
{
  fprintf(stderr, "usage: %s [file]\n", progn);
  exit(1);
}


#define BUFLEN  255
#define MaxBuf_ab  4096


int main(int argc, char *argv[])
{
  int   i, n;
  int   lcount;
  char  buf[BUFLEN + 1];
  char  *fnin;
  FILE  *fpin;

  int   F_a;

  int                 nSeq, mSeq, nSeq3, iSeq;
  T_k_Nk              *Seq;
  int                 k, k1, k2, k3;
  llint               Nk, Nk1, Nk2, Nk3;
  T_llintFrac         a, b, c;
  int                 n_a, n_b;
  static T_llintFrac  ai[MaxBuf_ab], bi[MaxBuf_ab];


  F_a  = 0;
  fnin = NULL;

  for (i = 1; i < argc; i++)
  {
    if (argv[i][0] == '-')
    {
      if (str_icmp(argv[i], "-a") == 0)
        F_a = 1;
      else
        usage();
    }
    else
    {
      if (fnin)
        usage();

      fnin = argv[i];
    }
  }

  if (fnin)
  {
        fpin = fopen(fnin, "r");
    if (fpin == NULL)
    {
      fprintf(stderr, "%s: Can't open %s\n", progn, fnin);
      exit(1);
    }
  }
  else
    fpin = stdin;

  mSeq = 1000;

      Seq = malloc(mSeq * sizeof (*Seq));
  if (Seq == NULL)
    NotEnoughCore();

  nSeq = 0;

  lcount = 0;

  while (fgetline(fpin, buf, sizeof buf / sizeof (*buf)))
  {
    lcount++;

        n = firstnonblank(buf);
    if (n && n != '#')
    {
          n = sscanf(buf, "%ld%lld", &k, &Nk);
      if (n != 2)
        IllegalLine(fnin, lcount);

      if (nSeq == mSeq)
      {
        mSeq += 1000;

            Seq = realloc(Seq, mSeq * sizeof (*Seq));
        if (Seq == NULL)
          NotEnoughCore();
      }

      Seq[nSeq].k  = k;
      Seq[nSeq].Nk = Nk;
          nSeq++;
    }
  }

  mSeq = nSeq;

      Seq = realloc(Seq, mSeq * sizeof (*Seq));
  if (Seq == NULL)
    InternalError();

  if (nSeq % 3)
    progerror("Number of elements is not a multiple of 3");

  qsort(Seq, nSeq, sizeof (*Seq),
        (int (*)(const void *, const void *)) k_Nk_sort_low_k_first);

  n_a = 0;
  n_b = 0;

  nSeq3 = nSeq / 3;

  for (iSeq = 0; iSeq < nSeq3; iSeq++)
  {
    k1 = Seq[iSeq            ].k; Nk1 = Seq[iSeq            ].Nk;
    k2 = Seq[iSeq +     nSeq3].k; Nk2 = Seq[iSeq +     nSeq3].Nk;
    k3 = Seq[iSeq + 2 * nSeq3].k; Nk3 = Seq[iSeq + 2 * nSeq3].Nk;

    if (k2 - k1 != k3 - k2)
      progerror("Corrupt k");

#ifdef FLOAT_ABC
         abc(k1, Nk1, k2, Nk2, k3, Nk3);
#endif

    if (F_a)
    {
      a.n = Nk3 - Nk2 - Nk2 + Nk1;
      a.d = (llint) 2 * (k2 - k1) * (k2 - k1);
      llint_RemoveGCD2(&a.n, &a.d);

      fprintf(stdout, "a(%d) ", k3 % (k2 - k1));
      llint_PrintFraction(stdout, a.n, a.d); putc('\n',  stdout);

      if (n_a >= MaxBuf_ab) InternalError();

      ai[n_a].n = a.n;
      ai[n_a].d = a.d;
         n_a++;
    }
    else
    {
      if (iabc(k1, Nk1, k2, Nk2, k3, Nk3, &a, &b, &c) != 0)
      {
        fprintf(stdout, "abc(%d) Out of range\n", k3 % (k2 - k1));
        a.n = b.n = c.n = 0;
        a.d = b.d = c.d = 1;
      }
      else
      {
        fprintf(stdout, "abc(%d) ", k3 % (k2 - k1));
        llint_PrintFraction(stdout, a.n, a.d); putc(' ',  stdout);
        llint_PrintFraction(stdout, b.n, b.d); putc(' ',  stdout);
        llint_PrintFraction(stdout, c.n, c.d); putc('\n', stdout);
      }

      if (n_a >= MaxBuf_ab) InternalError();
      if (n_b >= MaxBuf_ab) InternalError();

      ai[n_a].n = a.n;
      ai[n_a].d = a.d;
         n_a++;
      bi[n_b].n = b.n;
      bi[n_b].d = b.d;
         n_b++;
    }
  }

  if (n_a)
  {
    T_FP_Info    FPI;
    T_llintFrac  Mean_a;

    (void) FindPeriod0Frac(NULL, ai, n_a, 1, &FPI);
    llintMeanFrac(ai, n_a, FPI.BestStart, FPI.BestLen, &Mean_a);

    fprintf(stdout, "# 3Dai{%d %d %d} a{",
      FPI.BestStart, FPI.BestLen, FPI.BestMatches);
    llint_PrintFraction(stdout, Mean_a.n, Mean_a.d);
    fprintf(stdout, "}\n");
  }

  if (n_b)
    (void) FindPeriod0Frac("# bi", bi, n_b, 1, NULL);

  return 0;
}
