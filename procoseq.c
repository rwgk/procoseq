#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>


#include "gic_vs_zeoatlas.h"
#include "gic_vs_icsd.h"
#include "gic_vs_typix.h"


#include "longer.h"
#include "fperiod.h"


#define Fprintf (void) fprintf
#define Fflush  (void) fflush
#define Sprintf (void) sprintf


static char  *progn = "procoseq";


#define progerror(message) p_rogerror(message, __LINE__)

static void p_rogerror(char *message, int source_code_line)
{
  Fflush(stdout);
  Fprintf(stderr, "%s(%d): %s\n", progn, source_code_line, message);
  exit(1);
}


#define NotEnoughCore() N_otEnoughCore(__LINE__)

static void N_otEnoughCore(int source_code_line)
{
  p_rogerror("Not enough core", source_code_line);
}


#define IllegalLine(fnin, lcount) I_llegalLine(fnin, lcount, __LINE__)

static void I_llegalLine(const char *fnin, int lcount, int source_code_line)
{
  Fflush(stdout);
  Fprintf(stderr, "%s(%d): Illegal line #%d",
    progn, source_code_line, lcount);
  if (fnin != NULL) Fprintf(stderr, " in %s", fnin);
  putc('\n', stderr);
  exit(1);
}


#undef VERBOSE_MALLOC
#undef VERY_VERBOSE

#ifdef VERBOSE_MALLOC

#define mActiveMallocs 100
static void *ActiveMallocs[mActiveMallocs];
static int  nActiveMallocs = 0;

static void CheckMallocs(void)
{
  char buf[16];
  char *cptr;
  int i;
  for (i = 0; i < nActiveMallocs; i++) {
    cptr = ActiveMallocs[i];
    if (cptr != NULL) {
      size_t s = *((size_t *) (cptr - 16));
      strcpy(buf, "ABCDEFGHIJKGLMN");
      *((size_t *) buf) = s;
      if (strcmp(cptr - 16, buf) != 0) {
        printf("CorruptFront: s=%ld ptr=%ld\n", (long) s, (long) (cptr - 16));
      }
      else if (strcmp(cptr + s, "ABCDEFGHIJKGLMN") != 0) {
        printf("CorruptBack: s=%ld ptr=%ld\n", (long) s, (long) (cptr - 16));
      }
#ifdef VERY_VERBOSE
      else {
        printf("Good: s=%ld ptr=%ld\n", (long) s, (long) (cptr - 16));
      }
#endif
    }
  }
}

static void *verbose_malloc(size_t s, const long line)
{
  char *cptr;
  Fprintf(stdout, "#(%ld) malloc(%ld)\n", line, (long) s);
  CheckMallocs();
  cptr = malloc(s + 32);
  if (cptr == NULL) return NULL;
  strcpy(cptr, "ABCDEFGHIJKGLMN");
  *((size_t *) cptr) = s;
  strcpy(&cptr[16 + s], "ABCDEFGHIJKGLMN");
  if (nActiveMallocs >= mActiveMallocs) progerror("Internal Error");
  ActiveMallocs[nActiveMallocs++] = cptr + 16;
  CheckMallocs();
  return cptr + 16;
}

static void *verbose_realloc(void *vptr, size_t s, const long line)
{
  char *cptr = (char *) vptr;
  int iAM;
  Fprintf(stdout, "#(%ld) realloc(%ld)\n", line, (long) s);
  CheckMallocs();
  for (iAM = 0; iAM < nActiveMallocs; iAM++) {
    if (ActiveMallocs[iAM] == vptr) {
      ActiveMallocs[iAM] = NULL;
      break;
    }
  }
  if (iAM == nActiveMallocs) progerror("Internal Error");
  cptr = realloc(cptr - 16, s + 32);
  if (cptr == NULL) return NULL;
  strcpy(cptr, "ABCDEFGHIJKGLMN");
  *((size_t *) cptr) = s;
  strcpy(&cptr[16 + s], "ABCDEFGHIJKGLMN");
  ActiveMallocs[iAM] = cptr + 16;
  CheckMallocs();
  return cptr + 16;
}

static void verbose_free(void *vptr, const long line)
{
  char *cptr = (char *) vptr;
  size_t s = *((size_t *) (cptr - 16));
  int i;
  Fprintf(stdout, "#(%ld) free %ld\n", line, (long) s);
  CheckMallocs();
  for (i = 0; i < nActiveMallocs; i++) {
    if (ActiveMallocs[i] == vptr) {
      ActiveMallocs[i] = NULL;
      break;
    }
  }
  if (i == nActiveMallocs) progerror("Internal Error");
  free(cptr - 16);
  CheckMallocs();
}

#define malloc(s) verbose_malloc(s, __LINE__)
#define free(p) verbose_free(p, __LINE__)
#define realloc(p, s) verbose_realloc(p, s, __LINE__)

#endif


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


static char *str_dup(char *s)
{
  char *c = (char *) malloc((strlen(s) + 1) * sizeof (*s));
  if (c != NULL) (void) strcpy(c, s);
  return c;
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


static int str_ibegin(char *s1, char *s2) /* string ignore-case begin */
{
  char     u1, u2;

  while (*s1 && *s2)
  {
    u1 = toupper(*s1++);
    u2 = toupper(*s2++);
    if      (u1 < u2) return -1;
    else if (u1 > u2) return  1;
  }
  if (*s2) return -1;
  return 0;
}


static int getnum(char *str, int **num, int max_nnum)
{
  int     i, n, nnum;
  char    buf[128], xtrac, *s;
  double  fn;


  nnum = 0;
  i = 0;

  for (s = str; *s; s++)
  {
    if (*s == ',')
    {
      nnum++;
      i = -1;
    }
    else if (isspace(*s))
    {
      if (i == 1)
      {
        nnum++;
        i = 0;
      }
    }
    else
      i = 1;
  }

  if (i)
    nnum++;

  if (max_nnum && nnum > max_nnum)
    return -1;

  if (nnum)
  {
        *num = malloc(nnum * sizeof (*(*num)));
    if (*num == NULL)
      progerror("Not Enough Core");
  }

  s = str;

  for (i = 0; i < nnum; i++)
  {
    while (*s && isspace(*s)) s++;

    n = 0;

    while (n < sizeof buf / sizeof (*buf))
    {
      if (*s == '\0' || *s == ',' || isspace(*s))
        break;

      buf[n++] = *s++;
    }

    if (n == sizeof buf / sizeof (*buf))
    {
      free(*num);
      *num = NULL;
      return -(i + 1);
    }

    buf[n] = '\0';

    if (n == 0)
      fn = 0.;
    else
    {
          n = sscanf(buf, "%lf%c", &fn, &xtrac);
      if (n != 1)
      {
        free(*num);
        *num = NULL;
        return -(i + 1);
      }

      if (fn < 0.)
          fn -= .5;
      else
          fn += .5;
    }

    (*num)[i] = (int) fn;

    if (*s == ',')
      s++;
  }

  return nnum;
}


static int str_arg(char *arg, char *s, int n)  /* - string argument -       */
{                                         /* return n-th argument of string */
  char *a = arg;                          /* s in arg.                      */
                                          /* n = 0 for the first argument.  */
  while (*s == ' ' || *s == '\t') s++;    /* Blank and Tab are seperators.  */

  while (n--)
  { while (*s != '\0' && *s != ' ' && *s != '\t') s++;
    while (*s == ' ' || *s == '\t') s++;
  }

  while (*s != '\0' && *s != ' ' && *s != '\t') *a++ = *s++;
  *a = '\0';

  return (a != arg);
}


static char firstnonblank(char *s)
{
  while (isspace(*s)) s++;
  return *s;
}


static int int_sort(const int *a, const int *b)
{
  if (*a > *b) return -1;
  if (*a < *b) return  1;
  return 0;
}


static double gammln(double xx)
{
  double x,tmp,ser;
  static double cof[6]={76.18009173,-86.50532033,24.01409822,
    -1.231739516,0.120858003e-2,-0.536382e-5};
  int j;

  x=xx-1.0;
  tmp=x+5.5;
  tmp -= (x+0.5)*log(tmp);
  ser=1.0;
  for (j=0;j<=5;j++) {
    x += 1.0;
    ser += cof[j]/x;
  }
  return -tmp+log(2.50662827465*ser);
}

#define ITMAX 100
#define EPS 3.0e-7

static void gcf(double *gammcf, double a, double x, double *gln)
{
  int n;
  double gold=0.0,g,fac=1.0,b1=1.0;
  double b0=0.0,anf,ana,an,a1,a0=1.0;

  *gln=gammln(a);
  a1=x;
  for (n=1;n<=ITMAX;n++) {
    an=(double) n;
    ana=an-a;
    a0=(a1+a0*ana)*fac;
    b0=(b1+b0*ana)*fac;
    anf=an*fac;
    a1=x*a0+anf*a1;
    b1=x*b0+anf*b1;
    if (a1) {
      fac=1.0/a1;
      g=b1*fac;
      if (fabs((g-gold)/g) < EPS) {
        *gammcf=exp(-x+a*log(x)-(*gln))*g;
        return;
      }
      gold=g;
    }
  }
  progerror("a too large, ITMAX too small in routine GCF");
}

#undef ITMAX
#undef EPS

#define ITMAX 100
#define EPS 3.0e-7

static void gser(double *gamser, double a, double x, double *gln)
{
  int n;
  double sum,del,ap;

  *gln=gammln(a);
  if (x <= 0.0) {
    if (x < 0.0) progerror("x less than 0 in routine GSER");
    *gamser=0.0;
    return;
  } else {
    ap=a;
    del=sum=1.0/a;
    for (n=1;n<=ITMAX;n++) {
      ap += 1.0;
      del *= x/ap;
      sum += del;
      if (fabs(del) < fabs(sum)*EPS) {
        *gamser=sum*exp(-x+a*log(x)-(*gln));
        return;
      }
    }
    progerror("a too large, ITMAX too small in routine GSER");
    return;
  }
}

#undef ITMAX
#undef EPS

static double gammq(double a, double x)
{
  double gamser,gammcf,gln;

  if (x < 0.0 || a <= 0.0)
    progerror("Invalid arguments in routine GAMMQ");

  if (x < (a+1.0)) {
    gser(&gamser,a,x,&gln);
    return 1.0-gamser;
  }
  gcf(&gammcf,a,x,&gln);
  return gammcf;
}

static double sqrarg;
#define SQR(a) (sqrarg=(a),sqrarg*sqrarg)

static void fit(double x[], double y[], int ndata,
                double *sig, int mwt,
                double *a, double *b,
                double *siga, double *sigb,
                double *chi2, double *q)
{
  int i;
  double wt,t,sxoss,sx=0.0,sy=0.0,st2=0.0,ss,sigdat;

  *b=0.0;
  if (mwt) {
    ss=0.0;
    for (i=0;i<ndata;i++) {
      wt=1.0/SQR(sig[i]);
      ss += wt;
      sx += x[i]*wt;
      sy += y[i]*wt;
    }
  } else {
    for (i=0;i<ndata;i++) {
      sx += x[i];
      sy += y[i];
    }
    ss=ndata;
  }
  sxoss=sx/ss;
  if (mwt) {
    for (i=0;i<ndata;i++) {
      t=(x[i]-sxoss)/sig[i];
      st2 += t*t;
      *b += t*y[i]/sig[i];
    }
  } else {
    for (i=0;i<ndata;i++) {
      t=x[i]-sxoss;
      st2 += t*t;
      *b += t*y[i];
    }
  }
  *b /= st2;
  *a=(sy-sx*(*b))/ss;
  *siga=sqrt((1.0+sx*sx/(ss*st2))/ss);
  *sigb=sqrt(1.0/st2);
  *chi2=0.0;
  if (mwt == 0) {
    for (i=0;i<ndata;i++)
      *chi2 += SQR(y[i]-(*a)-(*b)*x[i]);
    *q=1.0;
    sigdat=sqrt((*chi2)/(ndata-2));
    *siga *= sigdat;
    *sigb *= sigdat;
  } else {
    for (i=0;i<ndata;i++)
      *chi2 += SQR((y[i]-(*a)-(*b)*x[i])/sig[i]);
    *q=gammq(0.5*(ndata-2),0.5*(*chi2));
  }
}

#undef SQR


static void CopySeq(const llint *Seq, int nSeq, llint *N, int nN)
{
  /* ATTENTION: Sequence manipulated
   */

  if (nN > nSeq)
    *N++ = 1;

  while (nSeq--)
    *N++ = *Seq++;
}


static void Admin_itpl(char *fnout, FILE **fpout, int BeginEnd)
{
  if (*fpout == NULL)
  {
    if (fnout)
    {
          *fpout = fopen(fnout, "a");
      if (*fpout == NULL)
        *fpout = stdout;
    }
    else
      *fpout = stdout;

    if (BeginEnd)
      Fprintf(*fpout, ">Begin itpl\n");
  }
  else
  {
    if (BeginEnd)
      Fprintf(*fpout, ">End itpl\n");

    if (*fpout != stdout)
      (void) fclose(*fpout);
    else if (fnout)
      progerror("Can't append to itpl file");
    else
      Fflush(stdout);

    *fpout = NULL;
  }
}


static void Print_itpl(char *Label, llint *it, int nit, int *pls, int npls,
                       char *fnout)
{
  int   i;
  FILE  *fpout;

                     fpout = NULL;
  Admin_itpl(fnout, &fpout, 1);

  for (i = 0; i < nit; i++)
  {
    if (i)
    {
      if (i % 10) putc(' ',  fpout);
      else        putc('\n', fpout);
    }

    Fprintf(fpout, "%s", llint_to_str(it[i], NULL, 0));
  }

  Fprintf(fpout, "\n:");

  for (i = 0; i < npls; i++)
    Fprintf(fpout, " %d", pls[i]);

  Fprintf(fpout, " (%d) %s\n", nit, Label);

  Admin_itpl(fnout, &fpout, 1);
}


static void DCS(char *Label, const llint *Seq, int nSeq, int *pls, int npls,
                char *fnout)
{
  int        i, j, k, n;
  llint      *Na;
  int        nNa, nZero;
  FILE       *fpout;
  T_FP_Info  FPI;


  /* ATTENTION: Sequence manipulated
   */
                                nNa = nSeq;
  if (nSeq == 0 || Seq[0] != 1) nNa++;

      Na = malloc(nNa * sizeof (*Na));
  if (Na == NULL)
    NotEnoughCore();

  CopySeq(Seq, nSeq, Na, nNa);

  for (i = 0; i < npls; i++)
    for (j = nNa - 1, k = j - pls[i]; k >= 0; j--, k--)
      Na[j] -= Na[k];

  nZero = 0;
  n = 0;

  for (i = 0; i < nNa; i++)
  {
    if (Na[i] == 0)
    {
      if (nZero)
        nZero++;
      else
        nZero = 1;
    }
    else
    {
      if (n < nZero)
          n = nZero;

      nZero = 0;
    }
  }

  if (nZero >= n + 5)
    Print_itpl(Label, Na, nNa - nZero, pls, npls, fnout);
  else
  {
    n = FindPeriod0(NULL, Na, nNa, pls[npls - 1], &FPI);

                       fpout = NULL;
    Admin_itpl(fnout, &fpout, 1);

    if (n)
      Fprintf(fpout, ": DCS Period Length %d\n", n);

    for (i = 0; i < nNa; i++)
      Fprintf(fpout, "%s\n", llint_to_str(Na[i], NULL, 0));

    Admin_itpl(fnout, &fpout, 1);
  }

  free(Na);
}


static void Crack(char *Label, llint *Seq, int nSeq, int *maxpls, int npls,
                  char *fnout, int ReportDepth)
{
  int        i, j, k, n;
  llint      *Na;
  int        nNa;
  llint      *best_it;
  int        best_nit, *best_pls, best_hpls = -1;
  int        *pls, ipls, hpls, nZero, FirstCrack;
  int        *maxplscp;
  FILE       *fpout;
  T_FP_Info  FPI;
  int        MaxMatches;


  /* ATTENTION: Sequence manipulated
   */
                                nNa = nSeq;
  if (nSeq == 0 || Seq[0] != 1) nNa++;

      Na = malloc(nNa * sizeof (*Na));
  if (Na == NULL)
    NotEnoughCore();

      pls = malloc(npls * sizeof (*pls));
  if (pls == NULL)
    NotEnoughCore();

      best_it = malloc(nNa * sizeof (*best_it));
  if (best_it == NULL)
    NotEnoughCore();

      best_pls = malloc(npls * sizeof (*best_pls));
  if (best_pls == NULL)
    NotEnoughCore();

      maxplscp = malloc(npls * sizeof (*maxplscp));
  if (maxplscp == NULL)
    NotEnoughCore();

  for (i = 0; i < npls; i++)
  {
    maxplscp[i] = maxpls[i];
    pls[i] = 0;
  }

  FirstCrack = 0;
  best_nit = 0;
  MaxMatches = 0;

  ipls = 0;
  hpls = 0;

  for (;;)
  {
        pls[ipls]++;
    if (pls[ipls] > maxplscp[ipls])
    {
          ipls++;
      if (ipls == npls)
        break;

      if (hpls < ipls)
      {
          hpls = ipls;

                           fpout = NULL;
        Admin_itpl(fnout, &fpout, 0);

        Fprintf(fpout, "#     Crack depth %d\n", hpls + 1);

        Admin_itpl(fnout, &fpout, 0);
      }

      continue;
    }

    if (ipls)
    {
      for (i = 0; i < ipls; i++)
        pls[i] = pls[ipls];

      if (ipls >= ReportDepth)
      {
                           fpout = NULL;
        Admin_itpl(fnout, &fpout, 0);

        putc('@', fpout);

        for (i = 0; i <= hpls; i++)
          Fprintf(fpout, " %d", pls[i]);

        putc('\n', fpout);

        Admin_itpl(fnout, &fpout, 0);
      }

      ipls = 0;
    }

    CopySeq(Seq, nSeq, Na, nNa);

    for (i = 0; i <= hpls; i++)
      for (j = nNa - 1, k = j - pls[i]; k >= 0; j--, k--)
        Na[j] -= Na[k];

    nZero = 0;
    n = 0;

    for (i = 0; i < nNa; i++)
    {
      if (Na[i] == 0)
      {
        if (nZero)
          nZero++;
        else
          nZero = 1;
      }
      else
      {
        if (n < nZero)
            n = nZero;

        nZero = 0;
      }
    }

    if (nZero >= n + 5)
    {
      n = nNa - nZero;

      if (FirstCrack == 0)
      {
                           fpout = NULL;
        Admin_itpl(fnout, &fpout, 0);

        Fprintf(fpout, "#     First Crack");

        for (i = 0; i <= hpls; i++)
          Fprintf(fpout, " %d", pls[i]);

        Fprintf(fpout, " (%d)\n", n);

        Admin_itpl(fnout, &fpout, 0);

        FirstCrack = 1;
      }

      if (best_nit > n || best_nit == 0)
      {
          best_nit  = n;
          best_hpls = hpls;

        for (i = 0; i < best_nit; i++)
          best_it[i] = Na[i];

        for (i = 0; i <= best_hpls; i++)
          best_pls[i] = pls[i];
      }

            n = pls[0];
      for (i = 1; i <= hpls; i++)
        if (n < pls[i])
            n = pls[i];

      for (i = 0; i <  npls; i++)
        if (maxplscp[i] > n)
            maxplscp[i] = n;
    }
    else if ((n = FindPeriod0(NULL, Na, nNa, maxplscp[hpls], &FPI)) != 0)
    {
      if (MaxMatches < FPI.BestMatches)
      {
                           fpout = NULL;
        Admin_itpl(fnout, &fpout, 0);

        Fprintf(fpout, "#     %s Crack Status", Label);

        for (i = 0; i <= hpls; i++)
          Fprintf(fpout, " %d", pls[i]);

        Fprintf(fpout, " Period S,L,M %d %d %d\n",
          FPI.BestStart,
          FPI.BestLen,
          FPI.BestMatches);

        Admin_itpl(fnout, &fpout, 0);

        MaxMatches = FPI.BestMatches;
      }
    }
  }

  if (best_nit)
    Print_itpl(Label, best_it, best_nit, best_pls, best_hpls + 1, fnout);
  else
  {
                       fpout = NULL;
    Admin_itpl(fnout, &fpout, 1);

    Fprintf(fpout, "#     %s Crack", Label);

    for (i = 0; i < npls; i++)
      Fprintf(fpout, " %d", maxpls[i]);

    Fprintf(fpout, " FAIL\n");

    Admin_itpl(fnout, &fpout, 1);
  }

  free(maxplscp);
  free(best_pls);
  free(best_it);
  free(pls);
  free(Na);
}


static int CodeCmp(char *cza, char *cgic)
{
  char *cgic0 = cgic;

  do
  {
    while (*cza && *cza != '=')
    {
      if (*cza != *cgic) break;
           cza++;  cgic++;
    }

    if (*cgic == '\0' && (*cza == '\0' || *cza == '='))
      return 0;

    while (*cza && *cza != '=') cza++;
    if (*cza == '=') cza++;

    cgic = cgic0;
  }
  while (*cza);

  return -1;
}


static void PrintExpanded(FILE *fpout, char *s, char c, char *e)
{
  char  *e0 = e;

  while (*s)
  {
    if (*s != c)
      putc(*s++, fpout);
    else
    {
      for (e = e0; *e;)
        putc(*e++, fpout);

      s++;
    }
  }
}


static char *Find_TYPIX_Label(T_TYPIX_LblTab *TLT, char *gic_Label)
{
  char  buf[128];
  char  *TYPIX_Label;


  for (TYPIX_Label = NULL; TLT->TYPIX_Label; TLT++)
  {
    if (strcmp(TLT->gic_Label, gic_Label) == 0)
    {
      if (TYPIX_Label)
      {
        Sprintf(buf,
          "Internal Error: Duplicate gic Label %s (TYPIX)", gic_Label);
        progerror(buf);
      }

      TYPIX_Label = TLT->TYPIX_Label;
    }
  }

  return TYPIX_Label;
}


typedef struct
  {
    int                ID;
    char               *Label;
    int                Length;
    llint              *N;
    int                Same;
    T_gic_vs_zeoatlas  *ZA;
    T_gic_vs_icsd      *ICSD;
    T_gic_vs_typix     *TYPIX;
  }
  T_Coseq;


static T_Coseq *CopyCoseq(T_Coseq *Coseq, int *nCoseq, int *MaxCoseq,
                          T_Coseq *BufCoseq)
{
  T_Coseq  *Cs;


  if (*nCoseq == *MaxCoseq)
  {
                               *MaxCoseq += 10;
        Coseq = realloc(Coseq, *MaxCoseq * sizeof (*Coseq));
    if (Coseq == NULL) NotEnoughCore();
  }

  Cs = &Coseq[*nCoseq];

  Cs->ID     = *nCoseq;
  Cs->Label  = BufCoseq->Label;
  Cs->Length = BufCoseq->Length;
  Cs->Same   = BufCoseq->Same;
  Cs->ZA     = BufCoseq->ZA;
  Cs->ICSD   = BufCoseq->ICSD;
  Cs->TYPIX  = BufCoseq->TYPIX;

      Cs->N = malloc(Cs->Length * sizeof (*Cs->N));
  if (Cs->N == NULL) NotEnoughCore();

  (void) memcpy(Cs->N, BufCoseq->N, Cs->Length * sizeof (*(Cs->N)));

  BufCoseq->Label  = NULL;
  BufCoseq->Length = 0;
  BufCoseq->Same   = -1;
  BufCoseq->ZA     = NULL;
  BufCoseq->ICSD   = NULL;
  BufCoseq->TYPIX  = NULL;

  (*nCoseq)++;

  return Coseq;
}


static int CoseqSF(const T_Coseq *a, const T_Coseq *b)
{
  int  i;


  if (  a->ZA && ! b->ZA) return -1;
  if (! a->ZA &&   b->ZA) return  1;
  if (  a->ZA) return strcmp(a->ZA->ZA_Label, b->ZA->ZA_Label);

  for (i = 0; i < a->Length && i < b->Length; i++)
  {
    if (a->N[i] < b->N[i]) return -1;
    if (a->N[i] > b->N[i]) return  1;
  }

  if (a->Length < b->Length) return -1;
  if (a->Length > b->Length) return  1;

  return 0;
}


static void usage(void)
{
  Fprintf(stderr,
    "usage: %s"
    " [-Code] [-Maple] [-ACEgr] [-EoIS=#] [-Select=Label]\n"
    " [-NoPeriods] [-1.Diff] [-3.Diff] [-zoe2] [-Crack[=#,...]] [-DCS=#,...]\n"
    " [-itpl=file_name] [-ReportDepth=#] [-Show_abc[=file]] [-iFP]\n"
    " [focus-log-section-coseq]\n",
    progn);
  exit(1);
}


#define BufLen  255


int main(int argc, char *argv[])
{
  int   i, j, k, n, c;
  int   lcount;
  char  *fnin, *fn_abc;
  FILE  *fpin, *fp_abc;
  char  buf[BufLen + 1], *cp, *lbl, xtrac;
  char  CombLbl[2 * (BufLen + 1)];

  int   iarg, Must_be_equal_those, iOut;
  char  arg[BufLen + 1];

  int  MaxCoseq, Max_k;
  int    nCoseq, nDoubleSet, nUniq;
  int    iCoseq, jCoseq;

  T_Coseq    *Coseq, *Cs;
  T_Coseq  BufCoseq;

  llint   *FirstDiff, *SecondDiff, *ThirdDiff;
  llint   Sum_k2, Sum_Nk;
  double       *k_over_Sum_k2;
  double  *Sum_Nk_over_Sum_k2;
  double  a0, b;
  double  sig_a0, sig_b;
  double  chi2, q;

  char  *Code;
  int   F_Maple, F_ACEgr, F_EoIS_No, F_NoPeriods;
  int   F_1Diff, F_3Diff, F_zoe2;
  int   F_Crack, F_DCS, F_ReportDepth, F_Show_abc, F_iFP;
  int   *maxpls, npls, *DCSpls, DCSnpls;
  char  *F_itpl;
  char  *F_Select;

  T_FindPeriod2  xFindPeriod2;


  MaxCoseq = 10;
  Max_k = 50;

  putc('#', stdout);
  for (i = 0; i < argc;)
    Fprintf(stdout, " %s", argv[i++]);
  putc('\n', stdout);

  Fprintf(stdout, "# sizeof (llint)   = %d\n", sizeof (llint));
  Fprintf(stdout, "# sizeof (ldouble) = %d\n", sizeof (ldouble));

  Code          = NULL;
  F_Maple       = 0;
  F_ACEgr       = 0;
  F_EoIS_No     = 0;
  F_NoPeriods   = 0;
  F_1Diff       = 0;
  F_3Diff       = 0;
  F_zoe2        = 0;
  F_Crack       = 0;
  maxpls        = NULL;
  npls          = 0;
  F_DCS         = 0;
  DCSpls        = NULL;
  DCSnpls       = 0;
  F_itpl        = NULL;
  F_ReportDepth = 0;
  F_Select      = NULL;
  F_iFP         = 0;
  F_Show_abc    = 0;
  fn_abc        = NULL;
  fnin          = NULL;

  for (i = 1; i < argc; i++)
  {
    if (argv[i][0] == '-')
    {
      if   (   str_icmp(&argv[i][1], "h") == 0
            || str_icmp(&argv[i][1], "help") == 0)
        usage();
      else if (str_icmp(&argv[i][1], "Maple") == 0)
      {
        if (F_Maple) usage();
        F_Maple = 1;
      }
      else if (str_icmp(&argv[i][1], "ACEgr") == 0)
      {
        if (F_ACEgr) usage();
        F_ACEgr = 1;
      }
      else if (str_ibegin(&argv[i][1], "EoIS=") == 0)
      {
        if (F_EoIS_No) usage();
            n = sscanf(&argv[i][6], "%d%c", &F_EoIS_No, &xtrac);
        if (n != 1 || F_EoIS_No < 1) usage();
      }
      else if (str_icmp(&argv[i][1], "NoPeriods") == 0)
      {
        if (F_NoPeriods) usage();
        F_NoPeriods = 1;
      }
      else if (str_icmp(&argv[i][1], "1.Diff") == 0)
      {
        if (F_1Diff) usage();
        F_1Diff = 1;
      }
      else if (str_icmp(&argv[i][1], "3.Diff") == 0)
      {
        if (F_3Diff) usage();
        F_3Diff = 1;
      }
      else if (str_icmp(&argv[i][1], "zoe2") == 0)
      {
        if (F_zoe2) usage();
        F_zoe2 = 1;
      }
      else if (str_ibegin(&argv[i][1], "Crack") == 0)
      {
        if (F_Crack) usage();

             cp = &argv[i][1 + 5];
        if (*cp == '=')
          cp++;
        else if (*cp == '\0')
          cp = "10,10,10";
        else
          usage();

            npls = getnum(cp, &maxpls, 0);
        if (npls < 1)
          usage();

        for (n = 0; n < npls; n++)
          if (maxpls[n] < 1) usage();

        qsort(maxpls, npls, sizeof (*maxpls),
              (int (*)(const void *, const void *)) int_sort);

        F_Crack = 1;
      }
      else if (str_ibegin(&argv[i][1], "DCS=") == 0)
      {
        if (F_DCS) usage();

            DCSnpls = getnum(&argv[i][1 + 4], &DCSpls, 0);
        if (DCSnpls < 1)
          usage();

        for (n = 0; n < DCSnpls; n++)
          if (DCSpls[n] < 1) usage();

        qsort(DCSpls, DCSnpls, sizeof (*DCSpls),
              (int (*)(const void *, const void *)) int_sort);

        F_DCS = 1;
      }
      else if (str_ibegin(&argv[i][1], "itpl=") == 0)
      {
        if (F_itpl)
          usage();

                   cp = &argv[i][1 + 5];
        if (strlen(cp) == 0)
          usage();

            F_itpl = str_dup(cp);
        if (F_itpl == NULL)
          NotEnoughCore();

            fpin = fopen(F_itpl, "r");
        if (fpin)
        {
          Fprintf(stderr,
            "%s: itpl file %s already existing: remove first\n",
            progn, F_itpl);
          exit(1);
        }
      }
      else if (str_ibegin(&argv[i][1], "ReportDepth=") == 0)
      {
        if (F_ReportDepth) usage();

            n = sscanf(&argv[i][13], "%d%c", &F_ReportDepth, &xtrac);
        if (n != 1 || F_ReportDepth < 1) usage();
      }
      else if (str_ibegin(&argv[i][1], "Select=") == 0)
      {
        if (F_Select) usage();

        F_Select = &argv[i][8];
      }
      else if (str_icmp(&argv[i][1], "iFP") == 0)
      {
        if (F_iFP) usage();
        F_iFP = 1;
      }
      else if (str_ibegin(&argv[i][1], "Show_abc") == 0)
      {
        if (F_Show_abc) usage();
            F_Show_abc = 1;

             cp = &argv[i][9];
        if (*cp)
        {
          if (*cp++ != '=') usage();
          fn_abc = cp;
        }
      }
      else
      {
        if (Code) usage();
        Code = &argv[i][1];
      }
    }
    else if (fnin == NULL)
      fnin = argv[i];
    else
      usage();
  }

  if (F_iFP == 0) xFindPeriod2 =  FindPeriod2;
  else            xFindPeriod2 = iFindPeriod2;

  if (Code == NULL)
  {
    if (fnin)
      Code = fnin;
    else
      Code = progn;
  }

  if (F_ReportDepth == 0)
      F_ReportDepth = npls + 1;

  if (fnin)
  {
        fpin = fopen(fnin, "r");
    if (fpin == NULL)
    {
      Fprintf(stderr, "%s: Can't open %s\n", progn, fnin);
      exit(1);
    }
  }
  else
    fpin = stdin;

  if (fn_abc)
  {
        fp_abc = fopen(fn_abc, "w");
    if (fp_abc == NULL)
    {
      Fprintf(stderr, "%s: Can't open %s\n", progn, fn_abc);
      exit(1);
    }
  }
  else if (F_Show_abc)
    fp_abc = stdout;
  else
    fp_abc = NULL;

      Coseq = malloc(MaxCoseq * sizeof (*Coseq));
  if (Coseq == NULL) NotEnoughCore();

      BufCoseq.N = malloc(Max_k * sizeof (*BufCoseq.N));
  if (BufCoseq.N == NULL) NotEnoughCore();

  BufCoseq.ID     = -1;
  BufCoseq.Label  = NULL;
  BufCoseq.Length = 0;
  BufCoseq.Same   = -1;
  BufCoseq.ZA     = NULL;
  BufCoseq.ICSD   = NULL;
  BufCoseq.TYPIX  = NULL;

  nUniq = 0;
  nCoseq = 0;
  Must_be_equal_those = 0;

  lcount = 0;

  while (fgetline(fpin, buf, sizeof buf / sizeof (*buf)))
  {
    lcount++;

        c = firstnonblank(buf);
    if (c && c != '#')
    {
          cp = strstr(buf, " equal those of ");
      if (cp != NULL)
      {
        if (Must_be_equal_those != 1)
          IllegalLine(fnin, lcount);

        Must_be_equal_those = 0;

        if (BufCoseq.Length == 0)
          IllegalLine(fnin, lcount);

        cp += 16;
        while (isspace(*cp)) cp++;
        lbl = cp;
        while (*cp && isspace(*cp) == 0) cp++;
        *cp = '\0';

        for (iCoseq = 0; iCoseq < nCoseq; iCoseq++)
        {
          if (strcmp(Coseq[iCoseq].Label, lbl) == 0)
          {
            if (memcmp(Coseq[iCoseq].N, BufCoseq.N,
                       BufCoseq.Length * sizeof (*(BufCoseq.N))) != 0)
              IllegalLine(fnin, lcount);

            BufCoseq.Same = Coseq[iCoseq].ID;
            break;
          }
        }

        if (iCoseq == nCoseq)
          IllegalLine(fnin, lcount);

        Coseq = CopyCoseq(Coseq, &nCoseq, &MaxCoseq, &BufCoseq);
      }
      else
      {
        for (iarg = 0; ; iarg++)
        {
          if (str_arg(arg, buf, iarg) == 0)
            break;

          if (Must_be_equal_those != 0)
            IllegalLine(fnin, lcount);

          if (strcmp(arg, "...") == 0)
          {
            Must_be_equal_those = 1;
          }
          else if (isalpha(arg[0]))
          {
            if (BufCoseq.Length != 0)
            {
              Coseq = CopyCoseq(Coseq, &nCoseq, &MaxCoseq, &BufCoseq);
              nUniq++;
            }

                BufCoseq.Label = str_dup(arg);
            if (BufCoseq.Label == NULL) NotEnoughCore();
          }
          else
          {
            if (BufCoseq.Length == Max_k)
            {
                                          Max_k += 100;
                  BufCoseq.N
                    = realloc(BufCoseq.N, Max_k * sizeof (*BufCoseq.N));
              if (BufCoseq.N == NULL) NotEnoughCore();
            }

            if (read_llint(arg, &BufCoseq.N[BufCoseq.Length]) != 0)
              IllegalLine(fnin, lcount);

            BufCoseq.Length++;
          }
        }
      }
    }
  }

  if (BufCoseq.Length != 0)
  {
    if (Must_be_equal_those)
      IllegalLine(fnin, lcount);

    Coseq = CopyCoseq(Coseq, &nCoseq, &MaxCoseq, &BufCoseq);
    nUniq++;
  }
  else if (BufCoseq.Label != NULL)
    IllegalLine(fnin, lcount);

  (void) fclose(fpin);

  free(BufCoseq.N);

  Cs = Coseq;

  for (iCoseq = 0; iCoseq < nCoseq; iCoseq++, Cs++)
  {
    if (Cs->Label == NULL)
    {
                      Sprintf(buf, "Q%3.3d", nCoseq);
          Cs->Label = str_dup(buf);
      if (Cs->Label == NULL)
        NotEnoughCore();
    }

#ifdef USE_GIC_VS
    if (Cs->Same == -1)
    {
      int NeedOne = 0;

      for (i = 0; gic_vs_zeoatlas[i].Code; i++)
      {
        if (CodeCmp(gic_vs_zeoatlas[i].Code, Code) == 0)
        {
          NeedOne = 1;

          if (strcmp(gic_vs_zeoatlas[i].gic_Label, Cs->Label) == 0)
          {
            if (Cs->ZA)
            {
              Sprintf(buf,
                "Internal Error: Duplicate gic_Label %s %s", Code, Cs->Label);
              progerror(buf);
            }

            Cs->ZA = &gic_vs_zeoatlas[i];
          }
        }
      }

      if (NeedOne && Cs->ZA == NULL)
      {
        Sprintf(buf, "No gic_Label %s for Code %s", Cs->Label, Code);
        progerror(buf);
      }

      for (i = 0; gic_vs_icsd[i].Code; i++)
      {
        if (strcmp(gic_vs_icsd[i].Code, Code) == 0)
        {
          if (Cs->ICSD)
          {
            Sprintf(buf,
              "Internal Error: Duplicate gic Code (ICSD) %s", Code);
            progerror(buf);
          }

          Cs->ICSD = &gic_vs_icsd[i];
        }
      }

      for (i = 0; gic_vs_typix[i].Code; i++)
      {
        if (strcmp(gic_vs_typix[i].Code, Code) == 0)
        {
          if (Cs->TYPIX)
          {
            Sprintf(buf,
              "Internal Error: Duplicate gic Code (TYPIX) %s", Code);
            progerror(buf);
          }

          Cs->TYPIX = &gic_vs_typix[i];

          if (Find_TYPIX_Label(Cs->TYPIX->TLT, Cs->Label) == NULL)
          {
            Sprintf(buf, "Error: No TYPIX Label for %s %s", Code, Cs->Label);
            progerror(buf);
          }
        }
      }
    }
#endif
  }

  if (nCoseq > 1)
    qsort(Coseq, nCoseq, sizeof (*Coseq),
          (int (*)(const void *, const void *)) CoseqSF);

  if (F_Maple)
  {
    Cs = Coseq;

    for (iCoseq = 0; iCoseq < nCoseq; iCoseq++, Cs++)
    {
      if (Cs->Same == -1)
      {
        Fprintf(stdout, "print(`%s %s`);\n", Code, Cs->Label);
        Fprintf(stdout, "guessgf([");

        for (i = 0;;)
        {
          if (i && i % 10 == 0)
            Fprintf(stdout, "\\\n");

          Fprintf(stdout, "%s", llint_to_str(Cs->N[i], NULL, 0));

          if (++i == Cs->Length)
            break;

          putc(',', stdout);
        }

        Fprintf(stdout, "], x);\n");
      }
    }

    return 0;
  }

  if (F_EoIS_No)
  {
    int   iSeq, lLine, iLine, NeedComma;
    char  Letter_iLine[3] = { 'S', 'T', 'U' };

    iOut = 1;

    Cs = Coseq;

    for (iCoseq = 0; iCoseq < nCoseq; iCoseq++, Cs++)
    {
      if (Cs->Same == -1)
      {
        Fprintf(stdout, "%%I A%6.6d\n", F_EoIS_No);

        iSeq = 0;
        lLine = 0;

        for (iLine = 0; iLine < 3; iLine++)
        {
          if (iSeq == Cs->Length)
            break;

          if (iLine)
            putc('\n', stdout);

          lLine = fprintf(stdout, "%%%c A%6.6d ",
                          Letter_iLine[iLine], F_EoIS_No);
          NeedComma = 0;

          if (iSeq == 0)
          {
            putc('1', stdout);
            lLine++;
            NeedComma = 1;
          }

          for (;;)
          {
            lLine += NeedComma;
            lLine += sprintf(buf, "%s", llint_to_str(Cs->N[iSeq], NULL, 0));

            if (NeedComma && lLine > 79)
            {
              if (iLine < 2)
              {
                putc(',', stdout);
                break;
              }

              if (lLine > 80)
                break;
            }

            if (NeedComma) putc(',', stdout);
            Fprintf(stdout, "%s", buf);
            NeedComma = 1;

            if (++iSeq == Cs->Length)
              break;
          }
        }

        if (iLine || lLine)
          putc('\n', stdout);

        Fprintf(stdout, "%%N A%6.6d Coordination sequence ", F_EoIS_No);

        if      (Cs->ZA)
        {
          PrintExpanded(stdout, Cs->ZA->ZA_Label, (char) '=', " and ");
          Fprintf(stdout, " for Zeolite Code ");
          PrintExpanded(stdout, Cs->ZA->Code, (char) '=', " and ");
          putc('\n', stdout);
        }
        else if (Cs->ICSD)
        {
          if (nUniq > 1)
            Fprintf(stdout, "T%d ", iOut);

          Fprintf(stdout, "for %s\n", Cs->ICSD->Name);
        }
        else if (Cs->TYPIX)
        {
          char  *TL;

          Fprintf(stdout, "for %s", Cs->TYPIX->Name);

          TL = Find_TYPIX_Label(Cs->TYPIX->TLT, Cs->Label);

          if (TL[0])
            Fprintf(stdout, ", %s", TL);

          putc('\n', stdout);
        }
        else
          Fprintf(stdout, "for %s\n", Code);

        Fprintf(stdout, "%%R A%6.6d\n",      F_EoIS_No);
        Fprintf(stdout, "%%O A%6.6d 0,2\n",  F_EoIS_No);
        Fprintf(stdout, "%%K A%6.6d nonn\n", F_EoIS_No);
        Fprintf(stdout, "%%A A%6.6d %s\n",   F_EoIS_No,
          "rwgk@cci.lbl.gov (R.W. Grosse-Kunstleve)");

        if      (Cs->ICSD) {
          Fprintf(stdout, "%%D A%6.6d",        F_EoIS_No);
          Fprintf(stdout, " %s %d\n",
            "Inorganic Crystal Structure Database: Collection Code",
            Cs->ICSD->CollectionCode);
        }
        else if (Cs->TYPIX) {
          Fprintf(stdout, "%%D A%6.6d",        F_EoIS_No);
          Fprintf(stdout, " Gmelin Handbook of Inorg. and Organomet. Chem.,");
          Fprintf(stdout, " 8th Ed., 1994, TYPIX search code");
          Fprintf(stdout, " %s\n", Cs->TYPIX->SearchCode);
        }

        if      (Cs->ZA) {
          Fprintf(stdout, "%%H A%6.6d %s%s\n", F_EoIS_No,
            "<a href=\"http://www.iza-structure.org/\">",
            "Intl. Zeolite Assn. (Structure Commission)</a>");
        }

        Fprintf(stdout, "%%H A%6.6d %s%s\n", F_EoIS_No,
          "<a href=\"http://cci.lbl.gov/~rwgk/EIS/CS.html\">",
          "Coordination Sequences & Encyclopedia of Integer Sequences</a>");

        Fprintf(stdout, "%%p A%6.6d\n",      F_EoIS_No);
        putc('\n', stdout);

        F_EoIS_No++;
        iOut++;
      }
    }

    Fprintf(stdout, "# NEXT %d\n", F_EoIS_No);

    return 0;
  }

  if (F_ACEgr)
  {
        Sum_Nk_over_Sum_k2 = malloc(Max_k * sizeof (*Sum_Nk_over_Sum_k2));
    if (Sum_Nk_over_Sum_k2 == NULL) NotEnoughCore();

        k_over_Sum_k2 = malloc(Max_k * sizeof (*k_over_Sum_k2));
    if (k_over_Sum_k2 == NULL) NotEnoughCore();
  }
  else
  {
    Sum_Nk_over_Sum_k2 = NULL;
         k_over_Sum_k2 = NULL;
  }

  if (F_ACEgr)
  {
    Fprintf(stdout, "#\n");
    Fprintf(stdout, "@flush\n");
    Fprintf(stdout, "@WITH G0\n");
    Fprintf(stdout, "@G0 on\n");

    Cs = Coseq;

    for (iCoseq = 0; iCoseq < nCoseq; iCoseq++, Cs++)
    {
      if      (Cs->ZA)
      {
        Fprintf(stdout, "@    title \"Zeo-Atlas ");
        PrintExpanded(stdout, Cs->ZA->Code, (char) '=', " and ");
        Fprintf(stdout, "\"\n");
        break;
      }
      else if (Cs->ICSD)
      {
        Fprintf(stdout, "@    title \"ICSD %s\"\n", Cs->ICSD->Name);
        break;
      }
      else if (Cs->TYPIX)
      {
        Fprintf(stdout, "@    title \"TYPIX %s\"\n", Cs->TYPIX->Name);
        break;
      }
    }

    if (iCoseq == nCoseq)
      Fprintf(stdout, "@    title \"%s\"\n", Code);

    Fprintf(stdout, "#\n");
  }

  if (F_Select)
  {
    Cs = Coseq;

    for (iCoseq = 0; iCoseq < nCoseq; iCoseq++, Cs++)
    {
      if (strcmp(Cs->Label, F_Select) == 0)
      {
        if (Cs->Same == -1)
          F_Select = Cs->Label;
        else
        {
          for (jCoseq = 0; jCoseq < nCoseq; jCoseq++)
          {
            if (Coseq[jCoseq].ID == Cs->Same)
            {
              F_Select = Coseq[jCoseq].Label;
              break;
            }
          }

          if (jCoseq == nCoseq)
            progerror("Internal Error: Corrupt Cs->Same");
        }

        break;
      }
    }
  }

  SecondDiff = NULL;
  nDoubleSet = 0;
  Cs = Coseq;

  for (iCoseq = 0; iCoseq < nCoseq; iCoseq++, Cs++)
  {
    if (Cs->Same == -1 && (! F_Select || F_Select == Cs->Label))
    {
      if (nDoubleSet != 0 && F_ACEgr == 0)
        putc('\f', stdout);

      CombLbl[2 * BufLen + 1] = '\0';

      Sprintf(CombLbl, "%s %s", Code, Cs->Label);

      for (jCoseq = 0; jCoseq < nCoseq; jCoseq++)
      {
        if (Coseq[jCoseq].Same == Cs->ID)
        {
          Sprintf(buf, " = %s", Coseq[jCoseq].Label);
          (void) strcat(CombLbl, buf);
        }
      }

      if      (Cs->ZA)
      {
        Sprintf(buf, " (Zeo-Atlas %s %s)", Cs->ZA->Code, Cs->ZA->ZA_Label);
        (void) strcat(CombLbl, buf);
      }
      else if (Cs->ICSD)
      {
        Sprintf(buf, " (ICSD %s)", Cs->ICSD->Name);
        (void) strcat(CombLbl, buf);
      }
      else if (Cs->TYPIX)
      {
        Sprintf(buf, " (TYPIX %s)", Cs->TYPIX->Name);
        (void) strcat(CombLbl, buf);
      }

      if (CombLbl[2 * BufLen + 1] || strlen(CombLbl) > 2 * BufLen + 1)
        progerror("Internal Error: CombLbl overflow");

      Fprintf(stdout, "# %s\n", CombLbl);

      Fprintf(stdout, "#   Nk (%d)", Cs->Length);

      for (i = 0; i < Cs->Length; i++)
      {
        if (i % 10 == 0)
          Fprintf(stdout, "\n#   ");

        Fprintf(stdout, "  %s", llint_to_str(Cs->N[i], NULL, 0));

        if ((i + 1) % 50 == 0)
          Fprintf(stdout, " <%d>", i + 1);
      }

      putc('\n', stdout);

      if (F_NoPeriods == 0)
      {
                    Sprintf(buf, "#     %s 0.D(a,b,c)", CombLbl);
           n = xFindPeriod2(buf, Cs->N, Cs->Length, 0, NULL);

        if (n && fp_abc)
        (void) xFindPeriod2(buf, Cs->N, Cs->Length, n, fp_abc);

                   Sprintf(buf, "#     %s 0.D(a,0,c)", CombLbl);
        (void) FindPeriod1(buf, Cs->N, Cs->Length, 0);
                   Sprintf(buf, "#     %s 0.D linear", CombLbl);
        (void) FindPeriod0(buf, Cs->N, Cs->Length, 2, NULL);

        if (F_DCS)
          DCS(CombLbl, Cs->N, Cs->Length, DCSpls, DCSnpls, F_itpl);

        if (F_Crack)
          Crack(CombLbl, Cs->N, Cs->Length, maxpls, npls,
                F_itpl, F_ReportDepth - 1);
      }

      if (F_1Diff)
      {
            FirstDiff  = malloc((Cs->Length - 1) * sizeof (*FirstDiff));
        if (FirstDiff  == NULL) NotEnoughCore();

        for (i = 0, j = 1; j < Cs->Length; i++, j++)
          FirstDiff[i] = Cs->N[j] - Cs->N[i];

        Fprintf(stdout, "#   First Difference");

        for (i = 0; i < Cs->Length - 1; i++)
        {
          if (i % 10 == 0)
            Fprintf(stdout, "\n#   ");

          Fprintf(stdout, "  %s", llint_to_str(FirstDiff[i], NULL, 0));

          if ((i + 1) % 50 == 0)
            Fprintf(stdout, " <%d>", i + 1);
        }

        putc('\n', stdout);

        free(FirstDiff); FirstDiff = NULL;
      }

      if (SecondDiff) progerror("Internal Error");
          SecondDiff = malloc((Cs->Length - 2) * sizeof (*SecondDiff));
      if (SecondDiff == NULL) NotEnoughCore();

      for (i = 0, j = 1, k = 2; k < Cs->Length; i++, j++, k++)
        SecondDiff[i] = Cs->N[k] - Cs->N[j] - Cs->N[j] + Cs->N[i];

      Fprintf(stdout, "#   Second Difference");

      for (i = 0; i < Cs->Length - 2; i++)
      {
        if (i % 10 == 0)
          Fprintf(stdout, "\n#   ");

        Fprintf(stdout, "  %s", llint_to_str(SecondDiff[i], NULL, 0));

        if ((i + 1) % 50 == 0)
          Fprintf(stdout, " <%d>", i + 1);
      }

      putc('\n', stdout);

      if (F_NoPeriods == 0)
      {
                   Sprintf(buf, "#     %s 2.D", CombLbl);
        (void) FindPeriod0(buf, SecondDiff, Cs->Length - 2, 2, NULL);
      }

      if (F_zoe2)
      {
        Fprintf(stdout, ">Begin zoe2\n");

        Fprintf(stdout, " %d\n", Cs->Length - 2);
        Fprintf(stdout, " 1\n");
        Fprintf(stdout, " 10\n");

        for (i = 0; i < Cs->Length - 2; i++)
          Fprintf(stdout, "%s\n", llint_to_str(SecondDiff[i], NULL, 0));

        Fprintf(stdout, ">End zoe2\n");
      }

      if (F_ACEgr)
      {
        Sum_k2 = 0;
        Sum_Nk = 0;

        for (i = 0, k = 1; i < Cs->Length; i = k++)
        {
          Sum_k2 += k * k;
          Sum_Nk += Cs->N[i];

               k_over_Sum_k2[i] = (double)      k / Sum_k2;
          Sum_Nk_over_Sum_k2[i] = (double) Sum_Nk / Sum_k2;
        }

        fit(k_over_Sum_k2, Sum_Nk_over_Sum_k2, Cs->Length,
            NULL, 0,
            &a0, &b,
            &sig_a0, &sig_b,
            &chi2, &q);

        Fprintf(stdout, "#         a0 = %.10g      b = %.10g\n",     a0,     b);
        Fprintf(stdout, "#     sig a0 = %.10g  sig b = %.10g\n", sig_a0, sig_b);
        Fprintf(stdout, "#       chi2 = %.10g\n",                chi2);
      }

      if (F_3Diff)
      {
            ThirdDiff = malloc((Cs->Length - 3) * sizeof (*ThirdDiff));
        if (ThirdDiff == NULL) NotEnoughCore();

        for (j = 0, i = 1; i < Cs->Length - 2; j = i++)
          ThirdDiff[j] = SecondDiff[i] - SecondDiff[j];

        Fprintf(stdout, "#   Third Difference");

        for (i = 0; i < Cs->Length - 3; i++)
        {
          if (i % 10 == 0)
            Fprintf(stdout, "\n#   ");

          Fprintf(stdout, "  %s", llint_to_str(ThirdDiff[i], NULL, 0));

          if ((i + 1) % 50 == 0)
            Fprintf(stdout, " <%d>", i + 1);
        }

        putc('\n', stdout);

        free(ThirdDiff);
      }

      if (F_ACEgr)
      {
        i = nDoubleSet * 2;
        k = nDoubleSet % 15 + 1;

        Fprintf(stdout, "#\n");
        Fprintf(stdout,
          "@    s%d symbol %d\n"
          "@    s%d symbol size 0.5\n"
          "@    s%d symbol color %d\n"
          "@    s%d linestyle 0\n"
          "@    s%d color %d\n",
          i, nDoubleSet + 2,
          i,
          i, k,
          i,
          i + 1, k);

        Fprintf(stdout, "#\n");
        Fprintf(stdout, "#   k_over_Sum_k2  Sum_Nk_over_Sum_k2\n");

        for (i = 0; i < Cs->Length; i++)
        {
          Fprintf(stdout, "    %.10g %.10g\n",
                 k_over_Sum_k2[i],
            Sum_Nk_over_Sum_k2[i]);
        }

        Fprintf(stdout, "&\n");

        Fprintf(stdout, "    %.10g %.10g\n", 0., a0);
        Fprintf(stdout, "    %.10g %.10g\n", 1., a0 + b);

        Fprintf(stdout, "&\n");
      }

      nDoubleSet++;
    }

    if (SecondDiff) { free(SecondDiff); SecondDiff = NULL; }
  }

  if (fp_abc && fp_abc != stdout)
    (void) fclose(fp_abc);

  if (F_ACEgr)
    Fprintf(stdout, "@autoscale\n");

  return 0;
}
