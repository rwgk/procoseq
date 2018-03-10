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


#define Max_PL 128


/* reconstruction of coordination sequences */

static char *progn = "rcs";


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


static void *Rmalloc(size_t s, int Verbose)
{
  if (Verbose) {
    fprintf(stdout, "# malloc(%ld)\n", (long) s);
    fflush(stdout);
  }

  return malloc(s);
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


static int fgetfld(FILE *fpin, char s[], int size_s)
{
  int  c, i, ok;


  ok = 1;
  i = 0;

  while ((c = getc(fpin)) != EOF)
  {
#ifdef __MSDOS__
    if (c == 0x1A /* CtrlZ */) {
      ungetc(c, fpin);
      c = EOF;
      break;
    }
#endif

    if (isspace(c))
    {
      if (i) break;
    }
    else
    {
      if (i < size_s - 1)
        s[i++] = (char) c;
      else
        ok = -1;
    }
  }

  s[i] = '\0';

  if (i == 0 && c == EOF)
    return 0;

  return ok;
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


static signed char *PrepareIsSumOf2(int *PL, int F_v)
{
  int  PLC[2], lowPL[3], low_lcm, i, j;

  int          LCM2, Limit_rk;
  int          rk, iPL0;
  signed char  *ISO2;


  low_lcm = -1;

  for (i = 0; i < 3; i++)
  {
    PLC[0] = PL[i];

    for (j = i + 1; j < 3; j++)
    {
      PLC[1] = PL[j];

      LCM2 = FindLCM(PLC, 2);

      if (LCM2 < low_lcm || low_lcm < 0)
      {
        low_lcm  = LCM2;
        lowPL[0] = i;
        lowPL[1] = j;
      }
    }
  }

  for (i = 0; i < 3; i++)
    if (i != lowPL[0] && i != lowPL[1]) {
      lowPL[2] = i;
      break;
    }

  for (i = 0; i < 3; i++)
    lowPL[i] = PL[lowPL[i]];

  for (i = 0; i < 3; i++)
    PL[i] = lowPL[i];

      LCM2 = FindLCM(PL, 2);
  if (LCM2 != low_lcm) InternalError();

      Limit_rk = PL[0] + PL[1];
  if (Limit_rk < LCM2 + 1 - Limit_rk)
      Limit_rk = LCM2 + 1 - Limit_rk;

      ISO2 = Rmalloc(Limit_rk * sizeof (*ISO2), F_v);
  if (ISO2 == NULL)
    return NULL;

  for (rk = 0; rk < Limit_rk; rk++)
  {
    ISO2[rk] = 0;

    for (iPL0 = 0; (rk - iPL0) >= 0; iPL0 += PL[0])
      if ((rk - iPL0) % PL[1] == 0) {
        ISO2[rk] = 1;
        break;
      }
  }

  return ISO2;
}


static llint Closed3Nk(int k, llint *IT, int nIT, int *PL, signed char *ISO2)
{
  int    GCD3, GCD2, LCM2, Limit_rk;
  int    iIT, kk, kkk, rk;
  llint  Nk, Nkk;


  GCD3 = FindGCD(PL, 3);
  GCD2 = FindGCD(PL, 2);
  LCM2 = FindLCM(PL, 2);

      Limit_rk = PL[0] + PL[1];
  if (Limit_rk < LCM2 + 1 - Limit_rk)
      Limit_rk = LCM2 + 1 - Limit_rk;

  Nk = 0;

  for (iIT = 0; iIT < nIT; iIT++)
  {
        kk = k - iIT;
    if (kk < 0) break;

    if (kk % GCD3 == 0)
    {
      Nkk = 0;

      for (kkk = kk % PL[2]; kkk <= kk; kkk += PL[2])
      {
        if (kkk % GCD2 == 0)
        {
          Nkk += kkk / LCM2;

              rk = kkk % LCM2;
          if (rk < Limit_rk) {
            if (ISO2[rk]) Nkk++;
          }
          else
            Nkk++;
        }
      }

      Nk += Nkk * IT[iIT];
    }
  }

  return Nk;
}


static void RCS(      llint *Nk, const int nNk,
                const llint *IT, const int nIT,
                const   int *PL, const int nPL)
{
  int  i, j, n, p;


      n = nIT;
  if (n > nNk)
      n = nNk;

  for (i = 0; i < n  ; i++) Nk[i] = IT[i];
  for (     ; i < nNk; i++) Nk[i] = 0;

  for (p = 0; p < nPL; p++)
    for (i = 0, j = PL[p]; j < nNk; i++, j++)
      Nk[j] += Nk[i];
}


static int DCS(llint *Nk, llint *NkC, const int nNk,
               const int *PL, const int nPL)
{
  int         i, j, k, n;
  int         nZero, S_PL;


            S_PL = 0; for (i = 0; i < nPL; i++) S_PL += PL[i];
  if (nNk < S_PL)
    return 0;

  if (NkC == NULL)
    NkC = Nk;
  else
    for (i = 0; i < nNk; i++)
      NkC[i] = Nk[i];

  for (i = 0; i < nPL; i++)
    for (j = nNk - 1, k = j - PL[i]; k >= 0; j--, k--)
      NkC[j] -= NkC[k];

  nZero = 0;
  n = 0;

  for (i = 0; i < nNk; i++)
  {
    if (NkC[i] == 0)
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
    return nNk - nZero;

  if (nNk >= S_PL + 100)
  {
        n = FindPeriod0(NULL, NkC + S_PL, nNk - S_PL, 100, NULL);
    if (n)
      return -n;
  }

  return 0;
}


static int int_sort_high_first(const int *a, const int *b)
{
  if (*a > *b) return -1;
  if (*a < *b) return  1;
  return 0;
}


static int int_sort_low_first(const int *a, const int *b)
{
  if (*a < *b) return -1;
  if (*a > *b) return  1;
  return 0;
}


static int ReduceIT(      llint *Nk, const int  nNk,
                    const llint *IT, const int  nIT,
                            int *PL,       int *nPL)
{
  int    i, n, ipl, jpl;
  int    New_nIT, Div, Extra1;
  int    nPLA;
  int     PLA[Max_PL];
  int    nPLB;
  int     PLB[Max_PL];


  New_nIT = nIT;

  nPLA = *nPL; for (i = 0; i < nPLA; i++) PLA[i] = PL[i];

  RestartTop:

  qsort(PLA, nPLA, sizeof (*PLA),
        (int (*)(const void *, const void *)) int_sort_high_first);

  for (ipl = 0; ipl < nPLA; )
  {
    nPLB = 0;

    for (jpl = 0; jpl < nPLA; jpl++)
      if (ipl != jpl)
        PLB[nPLB++] = PLA[jpl];

    RCS(Nk, nNk, IT, nIT, PL, *nPL);

             n = DCS(Nk, NULL, nNk, PLB, nPLB);
    if      (n > 0 && n < New_nIT)
    {
      nPLA = nPLB; for (i = 0; i < nPLB; i++) PLA[i] = PLB[i];
      New_nIT = n;
      continue;
    }
    else if (n < 0 && -n < PLA[ipl])
    {
      PLB[nPLB++] = -n;

      RCS(Nk, nNk, IT, nIT, PL, *nPL);

          n = DCS(Nk, NULL, nNk, PLB, nPLB);
      if (n > 0 && n < New_nIT)
      {
        nPLA = nPLB; for (i = 0; i < nPLB; i++) PLA[i] = PLB[i];
        New_nIT = n;
        goto RestartTop;
      }
    }

    ipl++;
  }

  for (ipl = 0; ipl < nPLA; ipl++)
  {
    for (i = 0; i < nPLA; i++) PLB[i] = PLA[i];

    for (Div = 2; Div < PLA[ipl]; Div++)
    {
          PLB[ipl] = PLA[ipl] / Div;
      if (PLB[ipl] * Div != PLA[ipl])
        continue;

      for (Extra1 = 0; Extra1 <= 1; Extra1++)
      {
        nPLB = nPLA;

        if (Extra1)
        {
          if (nPLB >= Max_PL) InternalError();
          PLB[nPLB++] = 1;
        }

        RCS(Nk, nNk, IT, nIT, PL, *nPL);

                 n = DCS(Nk, NULL, nNk, PLB, nPLB);
        if      (n > 0 && n < New_nIT)
        {
          PLA[ipl] = PLB[ipl];
          nPLA = nPLB;
          New_nIT = n;
          goto RestartTop;
        }
        else if (n < 0 && -n < PLA[ipl])
        {
          if (nPLB >= Max_PL) InternalError();
          PLB[nPLB++] = -n;

          RCS(Nk, nNk, IT, nIT, PL, *nPL);

              n = DCS(Nk, NULL, nNk, PLB, nPLB);
          if (n > 0 && n < New_nIT)
          {
            nPLA = nPLB; for (i = 0; i < nPLB; i++) PLA[i] = PLB[i];
            New_nIT = n;
            goto RestartTop;
          }
        }
      }
    }
  }

  RCS(Nk, nNk, IT, nIT, PL, *nPL);

  *nPL = nPLA; for (i = 0; i < nPLA; i++) PL[i] = PLA[i];

  return DCS(Nk, NULL, nNk, PL, *nPL);
}


static int ReducePL(      llint *Nk, const int  nNk,
                    const llint *IT, const int  nIT,
                            int *PL,       int *nPL)
{
  int    i, n, ipl, jpl;
  int    nPLA;
  int     PLA[Max_PL], S_PLA;
  int    nPLB;
  int     PLB[Max_PL], S_PLB;
  int    nPLC;
  int     PLC[Max_PL], S_PLC;
  int    PreviousLargest, lcm, Opt_lcm, Opt_jpl, Choice;


  nPLC = 0;
  PreviousLargest = -1;

  nPLA = *nPL; for (i = 0; i < nPLA; i++) PLA[i] = PL[i];

  RestartTop:

  qsort(PLA, nPLA, sizeof (*PLA),
        (int (*)(const void *, const void *)) int_sort_low_first);

  for (ipl = 0; ipl < nPLA; )
  {
    nPLB = 0;

    for (i = 0; i < nPLA; i++)
      if (i != ipl)
        PLB[nPLB++] = PLA[i];

    RCS(Nk, nNk, IT, nIT, PL, *nPL);

             n = DCS(Nk, NULL, nNk, PLB, nPLB);
    if      (n > 0)
    {
      nPLA = nPLB; for (i = 0; i < nPLA; i++) PLA[i] = PLB[i];
      continue;
    }
    else if (n < 0 && -n < PLA[ipl])
    {
      PLB[nPLB++] = -n;

      RCS(Nk, nNk, IT, nIT, PL, *nPL);

          n = DCS(Nk, NULL, nNk, PLB, nPLB);
      if (n > 0)
      {
        nPLA = nPLB; for (i = 0; i < nPLB; i++) PLA[i] = PLB[i];
        goto RestartTop;
      }
    }

    ipl++;
  }

  for (;;)
  {
    for (ipl = 0; ipl < nPLA; ipl++)
    {
      for (jpl = ipl + 1; jpl < nPLA; jpl++)
      {
        nPLB = 0;

        for (i = 0; i < nPLA; i++)
          if (i != ipl && i != jpl)
            PLB[nPLB++] = PLA[i];

        S_PLB = 0; for (i = 0; i < nPLB; i++) S_PLB += PLB[i];

        if (S_PLB < nNk)
        {
          RCS(Nk, nNk, IT, nIT, PL, *nPL);

                   n = DCS(Nk, NULL, nNk, PLB, nPLB);
          if      (n > 0)
          {
            nPLA = nPLB; for (i = 0; i < nPLA; i++) PLA[i] = PLB[i];
            goto RestartTop;
          }
          else if (n < 0)
          {
            PLB[nPLB++] = -n;

            RCS(Nk, nNk, IT, nIT, PL, *nPL);

                n = DCS(Nk, NULL, nNk, PLB, nPLB);
            if (n > 0)
            {
              nPLA = nPLB; for (i = 0; i < nPLB; i++) PLA[i] = PLB[i];
              goto RestartTop;
            }
          }
        }
      }
    }

    qsort(PLA, nPLA, sizeof (*PLA),
          (int (*)(const void *, const void *)) int_sort_high_first);

    if (nPLA == 0 || PreviousLargest >= PLA[0])
      break;

    if (nPLC == 0) {
      nPLC = nPLA; for (i = 0; i < nPLA; i++) PLC[i] = PLA[i];
    }

    Opt_lcm = 0;
    Opt_jpl = 0;

    for (ipl = 0; ipl < nPLA; ipl++)
    {
      PLB[0] = PLA[ipl];

      for (jpl = ipl + 1; jpl < nPLA; jpl++)
      {
        PLB[1] = PLA[jpl];

        lcm = FindLCM(PLB, 2);

        if (lcm > PLA[0])
        {
          if (Opt_lcm > lcm || Opt_lcm < PLA[0]) {
              Opt_lcm = lcm;
              Opt_jpl = jpl;
          }
        }
        else
        {
          if (Opt_lcm < lcm) {
              Opt_lcm = lcm;
              Opt_jpl = jpl;
          }
        }
      }
    }

    if (Opt_lcm == 0 || PLA[Opt_jpl] == Opt_lcm)
      break;

    PLA[Opt_jpl] = Opt_lcm;

    qsort(PLA, nPLA, sizeof (*PLA),
          (int (*)(const void *, const void *)) int_sort_low_first);

    PreviousLargest = PLA[nPLA - 1];
  }

  RCS(Nk, nNk, IT, nIT, PL, *nPL);

  Choice = 0;

      n = DCS(Nk, NULL, nNk, PLA, nPLA);
  if (n > 0)
  {
    if (nPLC != 0 && nPLC == nPLA)
    {
      S_PLA = 0; for (i = 0; i < nPLA; i++) S_PLA += PLA[i];
      S_PLC = 0; for (i = 0; i < nPLC; i++) S_PLC += PLC[i];

      if (S_PLC < S_PLA)
        Choice = 1;
    }
  }
  else if (nPLC)
    Choice = 1;

  if (Choice == 0)
  {
    *nPL = nPLA; for (i = 0; i < nPLA; i++) PL[i] = PLA[i];
    return n;
  }

  RCS(Nk, nNk, IT, nIT, PL, *nPL);

  *nPL = nPLC; for (i = 0; i < nPLC; i++) PL[i] = PLC[i];

  return DCS(Nk, NULL, nNk, PLC, nPLC);
}


static int FactorPL(      llint *Nk, const int  nNk,
                    const llint *IT, const int  nIT,
                            int *PL,       int *nPL,
                            int *TP,       int  nTP)
{
  int    i, n, ipl;


  if (nTP + 1 > Max_PL) InternalError();

  n = 0;

  for (ipl = 0; ipl < nTP + 1; ipl++)
  {
    RCS(Nk, nNk, IT, nIT, PL, *nPL);

        n = DCS(Nk, NULL, nNk, TP, nTP);
    if (n >= 0)
      break;

    TP[ipl] = -n;
  }

  if (n < 0)
  {
    nTP++;
    RCS(Nk, nNk, IT, nIT, PL, *nPL);
  }

  qsort(TP, nTP, sizeof (*TP),
        (int (*)(const void *, const void *)) int_sort_high_first);

  *nPL = nTP; for (i = 0; i < nTP; i++) PL[i] = TP[i];

  if (n < 0)
    return DCS(Nk, NULL, nNk, PL, *nPL);

  return n;
}


static void iFPabc(char *Label, llint *Nk, int nNk, int Start, int Len,
                   T_llintFrac *ai, int m_ai,
                   T_llintFrac *bi, int m_bi,
                   T_llintFrac *ci, int m_ci)
{
  int          OutOfRange;
  int          k1, k2, k3;
  int          n_ai, n_bi, n_ci;
  T_llintFrac  a, b, c;
  T_FP_Info    FPIa, FPIb, FPIc;


  OutOfRange = 0;

  n_ai = n_bi = n_ci = 0;

  k1 = Start;
  k2 = Start + Len;
  k3 = Start + Len + Len;

  for (; k3 < nNk; k1++, k2++, k3++)
  {
    if (iabc(k1, Nk[k1], k2, Nk[k2], k3, Nk[k3], &a, &b, &c) != 0) {
      OutOfRange = 1;
      break;
    }

    if (n_ai < m_ai) { ai[n_ai].n = a.n; ai[n_ai].d = a.d; n_ai++; }
    if (n_bi < m_bi) { bi[n_bi].n = b.n; bi[n_bi].d = b.d; n_bi++; }
    if (n_ci < m_ci) { ci[n_ci].n = c.n; ci[n_ci].d = c.d; n_ci++; }
  }

  (void) FindPeriod0Frac(NULL, ai, n_ai, 1, &FPIa);
  (void) FindPeriod0Frac(NULL, bi, n_bi, 1, &FPIb);
  (void) FindPeriod0Frac(NULL, ci, n_ci, 1, &FPIc);

  fprintf(stdout, "# %s: QE{%d %d} a{%d %d %d} b{%d %d %d} c{%d %d %d}",
    Label, Start, Len,
    FPIa.BestStart, FPIa.BestLen, FPIa.BestMatches,
    FPIb.BestStart, FPIb.BestLen, FPIb.BestMatches,
    FPIc.BestStart, FPIc.BestLen, FPIc.BestMatches);
  if (OutOfRange) fprintf(stdout, "OutOfRange");
  putc('\n', stdout);

  fprintf(stdout, ": %d %d %d %s",
    FPIc.BestLen, FPIb.BestLen, FPIa.BestLen, Label);
  if (OutOfRange) fprintf(stdout, "OutOfRange");
  putc('\n', stdout);
}


static void Print_itpl(char *Label,
                       const llint *IT, const int nIT,
                       const   int *PL, const int nPL)
{
  int   i;


  if (IT)
  {
    for (i = 0; i < nIT; i++)
    {
      if (i)
      {
        if (i % 10) putc(' ',  stdout);
        else        putc('\n', stdout);
      }

      fprintf(stdout, "%s", llint_to_str(IT[i], NULL, 0));
    }

    putc('\n', stdout);
  }

  putc(':', stdout);

  for (i = 0; i < nPL; i++)
    fprintf(stdout, " %d", PL[i]);

  fprintf(stdout, " (%d) %s", nIT, Label);

  if (IT) {
    putc('\n', stdout);
    putc('\n', stdout);
  }
}


static void CleanLabel(char *Label)
{
  char  *lbl, *l;


  lbl = Label;

  while (*lbl && isspace(*lbl)) lbl++;

  if (*lbl == '(')
  {
    l = lbl;
                              l++;
    while (*l && isspace(*l)) l++;
    while (      isdigit(*l)) l++;
    while (*l && isspace(*l)) l++;

    if (*l == ')')
    {                           l++;
      while (*l && isspace(*l)) l++;

      lbl = l;
    }
  }

  for (l = lbl; *l; )
    if (isspace(*l++) == 0)
      while (lbl != l) *Label++ = *lbl++;

  *Label = '\0';
}


static void PrintSignificantLabel(FILE *fpout, const char *Label)
{
  const char  *lbl, *lnb;


  while (*Label &&   isspace(*Label))       Label++;
  while (*Label && ! isspace(*Label)) putc(*Label++, fpout);

  for (lnb = lbl = Label; *lbl; lbl++)
  {
    if (! isspace(*lbl))
    {
      if (! isdigit(*lbl) && ! isalpha(*lbl))
        break;

      lnb = lbl + 1;
    }
  }

  while (Label != lnb) putc(*Label++, fpout);
}


static int SumInt(const int *N, int nN)
{
  int Sum = 0;

  while (nN--)
    Sum += *N++;

  return Sum;
}


static void usage(void)
{
  fprintf(stderr,
    "usage: %s [-s] [-t] [-v] [-iFP] [-Fix=#,#] [-ExtraPL=#[,...]] [-TD#]\n"
    "           [-ReduceIT] [-ReducePL] [-FactorPL[=#,...]] [-Select=Label]\n"
    "           [-OneLine] [-Show_abc] [-ACEgr] [-NoPeriods] [file] [k_max]\n"
    "           [-DCS=#[,...]] [-EchoITPL] [-CheckClosed3Nk] [-3Dai]\n"
    "           [-Nk=#,#] [-UseLCM] [-iFPabc] [-MapleMx]\n",
    progn);

  exit(1);
}


#define BUFLEN  255
#define MaxFvNk  32


int main(int argc, char *argv[])
{
  int   i, j, n, mode, iostat, n_itpl;
  int   F_s, F_t, F_v, F_iFP, *F_Fix, F_ExtraPL, *ExtraPL, F_TDn, F_NoPeriods;
  int   F_ReduceIT, F_ReducePL, F_FactorPL, *FvFactorPL, F_OneLine, F_Show_abc;
  int   F_ACEgr, F_CheckClosed3Nk, F_DCS, *DCS_PL, F_EchoITPL, k_max;
  int   F_3Dai, F_Nk, *FvNk[MaxFvNk], F_UseLCM, F_iFPabc, F_MapleMx;
  int   F_ITstat;
  char  *F_Select;
  char  buf[BUFLEN + 1], xtrac, *cp;
  char  lbl[BUFLEN + 1 + 20];
  char  *fnin;
  FILE  *fpin;

    int  nIT, mIT;
  llint  *IT;
    int  nPL, S_PL, gcd, lcm, QEs;
    int   PL[Max_PL];
    int  nNk, Auto_k_max;
  llint  *Nk, S_Nk;

  T_FindPeriod2  xFindPeriod2;
  signed char    *ISO2;


  F_s   = 0;
  F_t   = 0;
  F_v   = 0;
  F_iFP = 0;
  F_Fix = NULL;
  F_ExtraPL = 0;
    ExtraPL = NULL;
  F_DCS    = 0;
    DCS_PL = NULL;
  F_TDn = -1;
  F_Select = NULL;
  F_NoPeriods = 0;
  F_ReduceIT = 0;
  F_ReducePL = 0;
  F_FactorPL = 0;
  FvFactorPL = NULL;
  F_OneLine = 0;
  F_Show_abc = 0;
  F_ACEgr = 0;
  F_CheckClosed3Nk = 0;
  F_EchoITPL = 0;
  F_3Dai = 0;
  F_Nk = 0;
  F_UseLCM = 0;
  F_iFPabc = 0;
  F_MapleMx = 0;
  F_ITstat = 0;
  k_max = -1;
  fnin = NULL;

  for (i = 1; i < argc; i++)
  {
    if (argv[i][0] == '-')
    {
      if      (str_icmp(argv[i], "-s") == 0)
      {
        if (F_s) usage();
        F_s = 1;
      }
      else if (str_icmp(argv[i], "-t") == 0)
      {
        if (F_t) usage();
        F_t = 1;
      }
      else if (str_icmp(argv[i], "-v") == 0)
      {
        if (F_v) usage();
        F_v = 1;
      }
      else if (str_icmp(argv[i], "-iFP") == 0)
      {
        if (F_iFP) usage();
        F_iFP = 1;
      }
      else if (str_ibegin(argv[i], "-Fix=") == 0)
      {
        if (F_Fix) usage();

        if (getnum(argv[i] + 5, &F_Fix, 2) != 2)
          usage();

        if (F_Fix[0] != 1 && F_Fix[0] != 2 || F_Fix[1] < 0)
          usage();
      }
      else if (str_ibegin(argv[i], "-ExtraPL=") == 0)
      {
        if (F_ExtraPL) usage();

            F_ExtraPL = getnum(argv[i] + 9, &ExtraPL, 0);
        if (F_ExtraPL < 1)
          usage();

        for (j = 0; j < F_ExtraPL; j++)
          if (ExtraPL[j] < 1)
            usage();
      }
      else if (str_ibegin(argv[i], "-DCS=") == 0)
      {
        if (F_DCS) usage();

            F_DCS = getnum(argv[i] + 5, &DCS_PL, 0);
        if (F_DCS < 1)
          usage();

        for (j = 0; j < F_DCS; j++)
          if (DCS_PL[j] < 1)
            usage();
      }
      else if (str_ibegin(argv[i], "-TD") == 0)
      {
        if (F_TDn >= 0) usage();

            n = sscanf(argv[i] + 3, "%d%c", &F_TDn, &xtrac);
        if (n != 1 && (n != 2 || isspace(xtrac) == 0))
          usage();

        if (F_TDn < 0)
          usage();
      }
      else if (str_icmp(argv[i], "-NoPeriods") == 0)
      {
        if (F_NoPeriods) usage();
        F_NoPeriods = 1;
      }
      else if (str_icmp(argv[i], "-ReduceIT") == 0)
      {
        if (F_ReduceIT) usage();
        F_ReduceIT  = 1;
      }
      else if (str_icmp(argv[i], "-ReducePL") == 0)
      {
        if (F_ReducePL) usage();
        F_ReducePL = 1;
      }
      else if (str_ibegin(argv[i], "-FactorPL") == 0)
      {
        if (F_FactorPL) usage();

                  cp = argv[i] + 9;
        if      (*cp == '\0')
          F_FactorPL = -1;
        else if (*cp++ != '=')
          usage();
        else
        {
              F_FactorPL = getnum(cp, &FvFactorPL, 0);
          if (F_FactorPL < 1)
            usage();

          for (j = 0; j < F_FactorPL; j++)
            if (FvFactorPL[j] < 1)
              usage();

          if (F_FactorPL > 2)
            qsort(&FvFactorPL[1], F_FactorPL - 1, sizeof (*FvFactorPL),
                  (int (*)(const void *, const void *)) int_sort_high_first);
        }
      }
      else if (str_ibegin(argv[i], "-Select=") == 0)
      {
        if (F_Select) usage();

        F_Select = argv[i] + 8;

        if (F_Select[0] == '\0')
          usage();
      }
      else if (str_icmp(argv[i], "-OneLine") == 0)
      {
        if (F_OneLine)
          usage();

        F_OneLine  = 1;
      }
      else if (str_icmp(argv[i], "-Show_abc") == 0)
      {
        if (F_Show_abc) usage();
        F_Show_abc  = 1;
      }
      else if (str_icmp(argv[i], "-ACEgr") == 0)
      {
        if (F_ACEgr) usage();
        F_ACEgr  = 1;
      }
      else if (str_icmp(argv[i], "-CheckClosed3Nk") == 0)
      {
        if (F_CheckClosed3Nk) usage();
        F_CheckClosed3Nk  = 1;
      }
      else if (str_icmp(argv[i], "-EchoITPL") == 0)
      {
        if (F_EchoITPL) usage();
        F_EchoITPL  = 1;
      }
      else if (str_ibegin(argv[i], "-3Dai") == 0)
      {
        if (F_3Dai) usage();

                  cp = argv[i] + 5;
        if      (*cp == '\0')
          F_3Dai = -1;
        else if (*cp++ != '=')
          usage();
        else
        {
              n = sscanf(cp, "%d%c", &F_3Dai, &xtrac);
          if (n != 1 && (n != 2 || isspace(xtrac) == 0))
            usage();

          if (F_3Dai < 1)
            usage();
        }
      }
      else if (str_ibegin(argv[i], "-Nk=") == 0)
      {
        if (F_Nk >= MaxFvNk) InternalError();

            n = getnum(&argv[i][4], &FvNk[F_Nk], 2);
        if (n != 2)
          usage();

        if (FvNk[F_Nk][0] < 0 || FvNk[F_Nk][0] > FvNk[F_Nk][1])
          usage();

        F_Nk++;
      }
      else if (str_icmp(argv[i], "-UseLCM") == 0)
      {
        if (F_UseLCM) usage();
        F_UseLCM = 1;
      }
      else if (str_icmp(argv[i], "-iFPabc") == 0)
      {
        if (F_iFPabc) usage();
        F_iFPabc = 1;
      }
      else if (str_ibegin(argv[i], "-MapleMx=") == 0)
      {
        if (F_MapleMx) usage();

            n = sscanf(argv[i] + 9, "%d%c", &F_MapleMx, &xtrac);
        if (n != 1 && (n != 2 || isspace(xtrac) == 0))
          usage();

        if (F_MapleMx < 1)
          usage();
      }
      else if (str_icmp(argv[i], "-ITstat") == 0)
      {
        if (F_ITstat) usage();
        F_ITstat = 1;
      }
      else
        usage();
    }
    else
    {
          n = sscanf(argv[i], "%d%c", &j, &xtrac);
      if (n == 1)
      {
        if (j < 0 || k_max >= 0)
          usage();

        k_max = j;
      }
      else
      {
        if (fnin)
          usage();

        fnin = argv[i];
      }
    }
  }

  if (F_TDn < 0)
    F_TDn = 10;

  if (F_iFP == 0) xFindPeriod2 =  FindPeriod2;
  else            xFindPeriod2 = iFindPeriod2;

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

  if (F_Fix)
  {
    fprintf(stdout, "# Period Length for FindPeriod%d() is fixed at %d\n",
      F_Fix[0], F_Fix[1]);
  }

  IT     = NULL;
  iostat = 1;

  for (n_itpl = 0; ; n_itpl++)
  {
    if (IT) free(IT);

    if (iostat == 0)
      break;
                    mIT = 1000;
        IT = malloc(mIT * sizeof (*IT));
    if (IT == NULL)
      NotEnoughCore();

    nIT = 0;
    nPL = 0;

    mode = 0;

    while ((iostat = fgetfld(fpin, buf, BUFLEN + 1)) != 0)
    {
          n = sscanf(buf, "%d%c", &j, &xtrac);
      if (n != 1)
      {
        mode = 1;
        break;
      }

      if (nIT == mIT)
      {
            n = (int)(mIT * .2 + .5);
        if (n < 1000)
            n = 1000;

        mIT += n;

            IT = realloc(IT, mIT * sizeof (*IT));
        if (IT == NULL)
          NotEnoughCore();
      }

      IT[nIT++] = j;
    }

    if (iostat == 0 && nIT == 0)
      continue;

    if (mode == 1 && iostat != 0)
    {
      while ((iostat = fgetfld(fpin, buf, BUFLEN + 1)) != 0)
      {
            n = sscanf(buf, "%d%c", &j, &xtrac);
        if (n != 1 || j < 1)
        {
          mode = 2;
          break;
        }

        if (nPL >= Max_PL)
          progerror("nPL >= Max_PL");

        PL[nPL++] = j;
      }
    }

    if (mode == 2 && iostat)
    {
          i = strlen(buf);
      if (i < BUFLEN)
        buf[i++] = ' ';
      else
        i = BUFLEN;

      iostat = fgetline(fpin, buf + i, BUFLEN + 1 - i);
    }

    n = (F_Select && strstr(buf, F_Select) == NULL);

    CleanLabel(buf);

    if (buf[0] == '\0')
    {
      if (fnin)
        sprintf(buf, "%s:%d", fnin, n_itpl);
      else
        sprintf(buf, "itpl:%d", n_itpl);
    }

    if (n && F_Select && strstr(buf, F_Select) == NULL)
      continue;

    if (nPL + F_ExtraPL > Max_PL)
      InternalError();

    for (i = 0; i < F_ExtraPL; i++)
      PL[nPL++] = ExtraPL[i];

    S_PL = SumInt(PL, nPL);
    gcd = FindGCD(PL, nPL);
    lcm = FindLCM(PL, nPL);

        QEs = nIT - S_PL;
    if (QEs < 0)
        QEs = 0;

    Auto_k_max = k_max;

    if (F_UseLCM)
    {                   Auto_k_max = n = QEs + 4 * lcm + 9;
      if (k_max >= 0 && Auto_k_max > k_max)
        Auto_k_max = k_max;

      if (F_v)
        fprintf(stdout, "# %s: Auto_k_max %d + 4 * %d + 9 = %d (%d)\n",
          buf, QEs, lcm, n, Auto_k_max);
    }

    if (F_DCS)
    {
      nNk = SumInt(DCS_PL, F_DCS) + QEs + 100;

          Nk = Rmalloc(nNk * sizeof (*Nk), F_v);
      if (Nk == NULL)
        NotEnoughCore();

      RCS(Nk, nNk, IT, nIT, PL, nPL);

          n = DCS(Nk, NULL, nNk, DCS_PL, F_DCS);
      if (n <= 0)
      {
        fprintf(stdout, "# DCS FAILURE (%d): %s\n", n, buf);
        free(Nk);
        continue;
      }

      free(IT);

          Nk = realloc(Nk, n * sizeof (*Nk));
      if (Nk == NULL)
        InternalError();

      IT = Nk; nIT = mIT = n;

      if (Max_PL < F_DCS) InternalError();

      nPL = F_DCS; for (i = 0; i < F_DCS; i++) PL[i] = DCS_PL[i];
    }

    if (F_EchoITPL)
    {
      Print_itpl(buf, IT, nIT, PL, nPL);
      fflush(stdout);
    }

    if (F_ITstat)
    {
      j = 0;

      for (i = 0; i < nIT; i++)
        if (IT[i] < 0) j++;

      fprintf(stdout, "ITstat %d %d %d %s\n", nIT, nIT - j, j, buf);
    }

    if (F_ReduceIT)
    {
          nNk = nIT * 2 + 200;
      if (nNk < 300)
          nNk = 300;

          Nk = Rmalloc(nNk * sizeof (*Nk), F_v);
      if (Nk == NULL)
        NotEnoughCore();

      nIT = ReduceIT(Nk, nNk, IT, nIT, PL, &nPL);
      Print_itpl(buf, Nk, nIT, PL, nPL);
      fflush(stdout);

      free(Nk);
    }
    else if (F_ReducePL)
    {
          n = (nIT - S_PL) * 2;
      if (n < 0)
          n = nIT;
      if (n < 200)
          n = 200;

          nNk = lcm * 3 + n;
      if (nNk < 300)
          nNk = 300;

      if (nNk > k_max + 1 && k_max >= 0)
          nNk = k_max + 1;

          Nk = Rmalloc(nNk * sizeof (*Nk), F_v);
      if (Nk == NULL)
        NotEnoughCore();

      nIT = ReducePL(Nk, nNk, IT, nIT, PL, &nPL);

      if (nIT > 0 && F_s == 0)
        Print_itpl(buf, Nk,   nIT, PL, nPL);
      else
        Print_itpl(buf, NULL, nIT, PL, nPL);

      if (nIT <= 0)
        fprintf(stdout, " ReducePL FAILURE");

      if (nIT <= 0 || F_s != 0)
        putc('\n', stdout);

      fflush(stdout);

      free(Nk);
    }
    else if (F_FactorPL)
    {
      int  TP[Max_PL], nTP;

      if (F_FactorPL < 0)
      {
        nTP = 3 - 1;
        n   = 0;
      }
      else
      {
        if (FvFactorPL[0] > Max_PL) InternalError();

        nTP = FvFactorPL[0] - 1;
        n   = F_FactorPL - 1;
      }

      for (i = 0; i < nTP - n; i++)
        TP[i] = lcm;

          j = n - nTP;
      if (j < 0)
          j = 0;

      while (j < n)
        TP[i++] = FvFactorPL[++j];

          n = (nIT - S_PL) * 2;
      if (n < 0)
          n = nIT;
      if (n < 20)
          n = 20;

          nNk = lcm + SumInt(TP, nTP) + n;
      if (nNk < 50)
          nNk = 50;

      nNk += 150;

      if (nNk > k_max + 1 && k_max >= 0)
          nNk = k_max + 1;

          Nk = Rmalloc(nNk * sizeof (*Nk), F_v);
      if (Nk == NULL)
        NotEnoughCore();

      nIT = FactorPL(Nk, nNk, IT, nIT, PL, &nPL, TP, nTP);

      if (nIT > 0 && F_s == 0)
        Print_itpl(buf, Nk,   nIT, PL, nPL);
      else
        Print_itpl(buf, NULL, nIT, PL, nPL);

      if (nIT <= 0)
        fprintf(stdout, " FactorPL FAILURE");

      if (nIT <= 0 || F_s != 0)
        putc('\n', stdout);

      fflush(stdout);

      free(Nk);
    }
    else if (F_ACEgr)
    {
      for (i = 0; i < nIT; i++)
        fprintf(stdout, "%d %s\n", i, llint_to_str(IT[i], NULL, 0));

      fprintf(stdout, "&\n");
    }
    else if (F_MapleMx)
    {
      int    iMx, k;

      fprintf(stdout, "with(linalg):\n");
      fprintf(stdout, "interface(prettyprint=false);\n");
      fprintf(stdout, "interface(screenwidth=1024);\n");

      nNk = F_MapleMx * F_MapleMx + 1;

          Nk = Rmalloc(nNk * sizeof (*Nk), F_v);
      if (Nk == NULL)
        NotEnoughCore();

      RCS(Nk, nNk, IT, nIT, PL, nPL);

      for (iMx = 0; iMx < F_MapleMx; iMx++)
      {
        fprintf(stdout, "b%d:=matrix(%d,1,[",
          iMx, F_MapleMx);

        k = iMx + 1;

        for (i = 0; i < F_MapleMx; i++, k += F_MapleMx)
        {
          fprintf(stdout, "%s", llint_to_str(Nk[k], NULL, 0));

          if (i + 1 < F_MapleMx) putc(',', stdout);
        }

        fprintf(stdout, "]);\n");

        fprintf(stdout, "m%d:=matrix(%d,%2$d,[\n",
          iMx, F_MapleMx);

        k = iMx + 1;

        for (i = 0; i < F_MapleMx; i++, k += F_MapleMx)
        {
          for (j = 0; j < F_MapleMx; j++)
          {
            if (j == 0)
              putc('1', stdout);
            else
            {
              fprintf(stdout, "%d", k);
              if (j > 1) fprintf(stdout, "^%d", j);
            }

            if (i + 1 < F_MapleMx || j + 1 < F_MapleMx) putc(',', stdout);
          }

          if (i + 1 < F_MapleMx) putc('\n', stdout);
        }

        fprintf(stdout, "]);\n");

        fprintf(stdout, "inverse(m%d);\n", iMx);
        fprintf(stdout, "d%dc%d:=multiply(\", b%2$d);\n",
          F_MapleMx, iMx);
      }

      fprintf(stdout, "quit\n");
    }
    else if (F_3Dai && nPL == 3)
    {
      int          k, ia, na, ma;
      llint        NkBuf[3], TDgcd;
      T_llintFrac  *TDfrac;

      Nk = NkBuf;

      fprintf(stdout, ">Begin 3Dai\n");

          ISO2 = PrepareIsSumOf2(PL, F_v);
      if (ISO2 == NULL)
        NotEnoughCore();

      fprintf(stdout, "# QE{ %d %d } ((%d %d) %d) %s\n",
        QEs, lcm, PL[0], PL[1], PL[2], buf);

      fflush(stdout);

      if (F_3Dai < 0)
      {
        na = PL[0]; for (i = 1; i < nPL; i++) if (na > PL[i]) na = PL[i];
        ma = na;
      }
      else
      {
        na = F_3Dai + 1;
        ma = 1;
      }

      na += ma;

          TDfrac = Rmalloc(na * sizeof (*TDfrac), F_v);
      if (TDfrac == NULL)
        NotEnoughCore();

      for (ia = 0; ia < na; ia++)
      {
        k = QEs + ia;

        for (i = 0; i < 3; i++, k += lcm)
        {
          Nk[i] = Closed3Nk(k, IT, nIT, PL, ISO2);

          fprintf(stdout, "N%d %s\n", k, llint_to_str(Nk[i], NULL, 0));
        }

        TDfrac[ia].n = Nk[2] - Nk[1] - Nk[1] + Nk[0];
        TDfrac[ia].d = (llint) lcm * lcm;
        TDfrac[ia].d += TDfrac[ia].d;

        TDgcd = llint_FindGCD2(TDfrac[ia].n, TDfrac[ia].d);

        fprintf(stdout, "a%d=(%s", ia, llint_to_str(TDfrac[ia].n, NULL, 0));
        fprintf(stdout, "/%s)",        llint_to_str(TDfrac[ia].d, NULL, 0));
        fprintf(stdout, " / %s \n",    llint_to_str(TDgcd,        NULL, 0));

        if (TDgcd > 1) {
          TDfrac[ia].n /= TDgcd;
          TDfrac[ia].d /= TDgcd;
        }

        fprintf(stdout, "a%d=%s", ia, llint_to_str(TDfrac[ia].n, NULL, 0));
        fprintf(stdout, "/%s\n",      llint_to_str(TDfrac[ia].d, NULL, 0));

        fflush(stdout);
      }

      (void) FindPeriod3Dai(buf, TDfrac, na, ma);

      fprintf(stdout, ">End 3Dai\n");
       fflush(stdout);

      free(TDfrac);
      free(ISO2);
    }
    else if (F_iFPabc)
    {
      int          m_ai, m_bi, m_ci;
      T_llintFrac   *ai,  *bi,  *ci;

          nNk = QEs + 3 * lcm + 100;
      if (nNk > k_max && k_max >= 0)
          nNk = k_max + 1;

      m_ci = nNk - QEs - 2 * lcm; if (m_ci <    0) m_ci =    0;
      m_bi = 32768 * 4;           if (m_bi > m_ci) m_bi = m_ci;
      m_ai = 256;                 if (m_ai > m_ci) m_ai = m_ci;

      ci = bi = ai = NULL;

          Nk = Rmalloc(nNk * sizeof (*Nk), F_v);
      if (Nk == NULL)
        NotEnoughCore();

      if (m_ai) {
            ai = Rmalloc(m_ai * sizeof (*ai), F_v);
        if (ai == NULL)
          NotEnoughCore();
      }

      if (m_bi) {
            bi = Rmalloc(m_bi * sizeof (*bi), F_v);
        if (bi == NULL)
          NotEnoughCore();
      }

      if (m_ci) {
            ci = Rmalloc(m_ci * sizeof (*ci), F_v);
        if (ci == NULL)
          m_ci = 0;
      }

      RCS(Nk, nNk, IT, nIT, PL, nPL);

      iFPabc(buf, Nk, nNk, QEs, lcm, ai, m_ai, bi, m_bi, ci, m_ci);

              free(Nk);
      if (ci) free(ci);
      if (bi) free(bi);
      if (ai) free(ai);
    }
    else if (F_DCS == 0 || F_EchoITPL == 0)
    {
          nNk = F_TDn;
      if (nNk < Auto_k_max)
          nNk = Auto_k_max;

      nNk++;

          Nk = Rmalloc(nNk * sizeof (*Nk), F_v);
      if (Nk == NULL)
        NotEnoughCore();

      RCS(Nk, nNk, IT, nIT, PL, nPL);

      S_Nk = 0;

      for (i = 0; i <= F_TDn; i++)
        S_Nk += Nk[i];

      fprintf(stdout,
        "# O(IT) = %d O(PL) = %d S(pl) = %d TD%d = %s %s { %d %d }\n",
        nIT, nPL, S_PL, F_TDn,
        llint_to_str(S_Nk, NULL, 0),
        buf, gcd, lcm);

      if (! F_s && Auto_k_max > 0)
      {
        if (! F_t)
          i = 0;
        else if (Auto_k_max + 1 >= 10)
          i = Auto_k_max + 1 - 10;

        if (nPL == 3 && F_CheckClosed3Nk)
        {
              ISO2 = PrepareIsSumOf2(PL, F_v);
          if (ISO2 == NULL)
            NotEnoughCore();

          while (i < Auto_k_max + 1)
          {
            llint C3Nk = Closed3Nk(i, IT, nIT, PL, ISO2);

            fprintf(stdout, "%s ", llint_to_str(Nk[i], NULL, 0));
            fprintf(stdout, "%s",  llint_to_str(C3Nk,  NULL, 0));

            if (Nk[i] != C3Nk)
              fprintf(stdout, " Closed3Nk ERROR");

            putc('\n', stdout);

            i++;
          }

          free(ISO2);
        }
        else if (F_OneLine == 0)
        {
          while (i < Auto_k_max + 1)
            fprintf(stdout, "%s\n", llint_to_str(Nk[i++], NULL, 0));
        }
        else
        {
          PrintSignificantLabel(stdout, buf);

          while (i < Auto_k_max + 1)
            fprintf(stdout, " %s", llint_to_str(Nk[i++], NULL, 0));

          putc('\n', stdout);
        }

        fflush(stdout);
      }

      if (F_Nk && Auto_k_max > 0)
      {
        for (i = 0; i < F_Nk; i++)
          for (j = FvNk[i][0]; j <= FvNk[i][1]; j++)
            if (j < nNk)
              fprintf(stdout, "N%d %s\n", j, llint_to_str(Nk[j], NULL, 0));
      }

      if (F_NoPeriods == 0 && Auto_k_max > 0)
      {
        if (F_Fix == NULL)
        {
                       sprintf(lbl, "# %s (a,b,c)", buf);
              n = xFindPeriod2(lbl, Nk, nNk, 0, NULL);
          if (n == 0)
          {
                       sprintf(lbl, "# %s (a,0,c)", buf);
            (void) FindPeriod1(lbl, Nk, nNk, 0);
          }
          else if (F_Show_abc)
            (void) xFindPeriod2(lbl, Nk, nNk, n, stdout);
        }
        else if (F_Fix[0] == 2)
        {
                      sprintf(lbl, "# %s (a,b,c)", buf);
          (void) xFindPeriod2(lbl, Nk, nNk, F_Fix[1],
                              (F_Show_abc ? stdout : NULL));
        }
        else
        {
                     sprintf(lbl, "# %s (a,0,c)", buf);
          (void) FindPeriod1(lbl, Nk, nNk, F_Fix[1]);
        }
      }

      free(Nk);
    }
  }

  return 0;
}
