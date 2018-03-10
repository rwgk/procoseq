#include <stdio.h>
#include <stdlib.h>


#include "longer.h"
#include "inumber.h"
#include "fperiod.h"


#define Fprintf (void) fprintf
#define Sprintf (void) sprintf
#define Fflush  (void) fflush


int FindPeriod0(char *Label, llint *Seq, int nSeq, int MinMatches,
                T_FP_Info *FPI)
{
  int  BestStart, BestLen, BestMatches;
  int      Start,     Len,     Matches, nSeq_Len;
  int  i1, i2;


  if (nSeq && Seq[0] == 1)
  {
    /* ATTENTION: Sequence manipulated */

     Seq++;
    nSeq--;
  }

  BestLen     = 0;
  BestMatches = MinMatches;
  BestStart   = nSeq;

                nSeq_Len = nSeq - 1;
  for (Len = 1; nSeq_Len >= BestMatches; Len++, nSeq_Len--)
  {
    Start = 0;

    for (i2 = nSeq - 1; i2 >= nSeq_Len; i2--)
    {
      i1 = i2 - Len; if (i1 < Start) break;

      for (;;)
      {
        if (Seq[i1] != Seq[i2])
        {
          Start = i1 + 1;
          if (Start >= BestStart) goto NextLen;
          break;
        }

        i1 -= Len; if (i1 < Start) break;
      }
    }

    Matches = nSeq_Len - Start;

    if (    BestMatches <  Matches
        || (BestMatches == Matches && BestStart > Start))
    {
      BestLen     = Len;
      BestMatches = Matches;
      BestStart   = Start;
    }

    NextLen:;
  }

  if (FPI)
  {
    if (BestLen == 0) BestMatches = 0;

    FPI->BestLen     = BestLen;
    FPI->BestMatches = BestMatches;
    FPI->BestStart   = BestStart;
  }

  if (Label)
  {
    if (BestLen)
    {
      Fprintf(stdout, "%s: Period  Start = %d  Length = %d  Matches = %d",
        Label, BestStart + 1, BestLen, BestMatches);

      if (BestMatches < BestLen) Fprintf(stdout, "  Vague");

      putc('\n', stdout);
    }
    else
      Fprintf(stdout, "%s: No Period found\n", Label);

    Fflush(stdout);
  }

  if (BestMatches < BestLen)
    return 0;

  return BestLen;
}


int FindPeriod0Frac(char *Label, T_llintFrac *Seq, int nSeq, int MinMatches,
                    T_FP_Info *FPI)
{
  int  BestStart, BestLen, BestMatches;
  int      Start,     Len,     Matches, nSeq_Len;
  int  i1, i2;


  BestLen     = 0;
  BestMatches = MinMatches;
  BestStart   = nSeq;

                nSeq_Len = nSeq - 1;
  for (Len = 1; nSeq_Len >= BestMatches; Len++, nSeq_Len--)
  {
    Start = 0;

    for (i2 = nSeq - 1; i2 >= nSeq_Len; i2--)
    {
      i1 = i2 - Len; if (i1 < Start) break;

      for (;;)
      {
        if (Seq[i1].n != Seq[i2].n || Seq[i1].d != Seq[i2].d)
        {
          Start = i1 + 1;
          if (Start >= BestStart) goto NextLen;
          break;
        }

        i1 -= Len; if (i1 < Start) break;
      }
    }

    Matches = nSeq_Len - Start;

    if (    BestMatches <  Matches
        || (BestMatches == Matches && BestStart > Start))
    {
      BestLen     = Len;
      BestMatches = Matches;
      BestStart   = Start;
    }

    NextLen:;
  }

  if (FPI)
  {
    if (BestLen == 0) BestMatches = 0;

    FPI->BestLen     = BestLen;
    FPI->BestMatches = BestMatches;
    FPI->BestStart   = BestStart;
  }

  if (Label)
  {
    if (BestLen)
    {
      Fprintf(stdout, "# iFP0{%d %d %d} %s\n",
        BestStart, BestLen, BestMatches, Label);
    }
    else
      Fprintf(stdout, "# iFP0{No Period found} %s\n", Label);

    Fflush(stdout);
  }

  return BestLen;
}


int FindPeriod1(char *Label, llint *Seq, int nSeq, int FixedLen)
{
  int      iStart, Len, nDef, Best_iStart, BestLen, MinDef, Valid;
  int      i1, i2, i3;
  ldouble  x12, x22, x32, y2, y3_y2, a, c;
  int      InitSumMinMax;
  ldouble  Sum_a, BestSum_a;
  ldouble  Min_a, BestMin_a;
  ldouble  Max_a, BestMax_a;
  ldouble  Sum_c, BestSum_c;
  ldouble  Min_c, BestMin_c;
  ldouble  Max_c, BestMax_c;


  if (nSeq && Seq[0] == 1)
  {
    /* ATTENTION: Sequence manipulated */

     Seq++;
    nSeq--;
  }

  Best_iStart = 0;
  BestLen     = 0;
  MinDef      = nSeq - 3;

  if (FixedLen)
    Len = FixedLen;
  else
    Len = 1;

  while (Len * 2 < MinDef)
  {
    InitSumMinMax = 1;

    iStart = 0;

    for (i3 = nSeq - 1; i3 >= nSeq - Len; i3--)
    {
      i2 = i3 - Len;
      i1 = i2 - Len;

      if (FixedLen == 0 && i1 < iStart)
        break;

      x22 = (llint)(i2 + 1) * (i2 + 1);
      x32 = (llint)(i3 + 1) * (i3 + 1);
         y2 =           Seq[i2];
      y3_y2 = Seq[i3] - Seq[i2];

      a = y3_y2 / (x32 - x22);
      c = y2 - a * x22;

      if (InitSumMinMax)
      {
        Sum_a = a; Min_a = Max_a = a;
        Sum_c = c; Min_c = Max_c = c;
        InitSumMinMax = 0;
      }
      else
      {
        Sum_a += a; if (Min_a > a) Min_a = a; if (Max_a < a) Max_a = a;
        Sum_c += c; if (Min_c > c) Min_c = c; if (Max_c < c) Max_c = c;
      }

      c += .5; /* rounding */

      for ( ; i1 >= iStart; i1 -= Len)
      {
                                   x12 = (llint)(i1 + 1) * (i1 + 1);
        if (Seq[i1] != (llint)(a * x12 + c))
        {
          iStart = i1 + 1;
          break;
        }
      }
    }
                 nDef = iStart + Len * 2;
    if (MinDef > nDef)
    {
        Best_iStart = iStart;
        BestLen     = Len;
        MinDef      = nDef;
        BestSum_a   = Sum_a;
        BestMin_a   = Min_a;
        BestMax_a   = Max_a;
        BestSum_c   = Sum_c;
        BestMin_c   = Min_c;
        BestMax_c   = Max_c;
    }

    if (FixedLen)
      break;
    else
      Len++;
  }

  Valid = 0;

  if (BestLen)
  {
    Fprintf(stdout, "%s: P{%d %d %d}",
      Label, Best_iStart + 1, BestLen, nSeq - MinDef);

    if (nSeq - MinDef >= BestLen)
      Valid = 1;

    if (Valid || FixedLen)
      Fprintf(stdout, " a{%.9f %.9f %.9f} c{%.9f %.9f %.9f}",
        (double) BestMin_a,
        (double) BestMax_a,
        (double)(BestSum_a / BestLen),
        (double) BestMin_c,
        (double) BestMax_c,
        (double)(BestSum_c / BestLen));

    if (Valid == 0)
      Fprintf(stdout, " Vague");

    putc('\n', stdout);
  }
  else
    Fprintf(stdout, "%s: No Period found\n", Label);

  Fflush(stdout);

  return Valid;
}


static void SimplifyFraction(int nume, int deno, int *o_nume, int *o_deno)
{
  int gcd = FindGCD2(nume, deno);
  if (gcd)
  {
    *o_nume = nume / gcd;
    *o_deno = deno / gcd;

    if (*o_deno < 0) {
      *o_nume *= -1;
      *o_deno *= -1;
    }
  }
}


static const char *FormatFraction(int nume, int deno, int Decimal,
                                  char *Buffer, int SizeBuffer)
{
  int          n, d;
  char         *cp, *cpp;
  static char  StaticBuffer[40];


  if (NULL == Buffer) {
              Buffer =        StaticBuffer;
          SizeBuffer = sizeof StaticBuffer / sizeof (*StaticBuffer);
  }

  Buffer[SizeBuffer - 1] = '\0';

  if (nume == 0)
  {
    Buffer[0] = '0';
    Buffer[1] = '\0';
  }
  if (Decimal)
  {
    Sprintf(Buffer, "%.6g", (double) nume / deno);

         cp = Buffer;
    if (*cp == '-') cp++;
    if (*cp == '0') {
      cpp = cp + 1; while (*cp) *cp++ = *cpp++;
    }
  }
  else
  {
    SimplifyFraction(nume, deno, &n, &d);

    if (d == 1)
      Sprintf(Buffer, "%d", n);
    else
      Sprintf(Buffer, "%d/%d", n, d);
  }

  if (Buffer[SizeBuffer - 1] != '\0') {
      Buffer[SizeBuffer - 1] =  '\0';
    Fprintf(stderr, "Internal Error: FormatFraction(): Buffer too small");
    exit(1);
  }

  return Buffer;
}


void ShowCoeff(FILE *fpout, double origv, int ifac, int MaxFac)
{
  double   val, delta;
  int     ival, isgn, i;


  for (i = 1; i < MaxFac; i++)
  {
    val = origv * (i * ifac);

    if (val >= 0.)
      isgn =  1;
    else {
      isgn = -1;
       val = -val;
    }

    ival = (int)(val + .5);

        delta = val - (double) ival;
    if (delta < 0) delta = -delta;

    if (ival < 1)
      val = 1.;

    if (delta / val <= 1.e-6)
      break;
  }

  if (i < MaxFac)
    Fprintf(fpout, "%s", FormatFraction(isgn * ival, ifac * i, 0, NULL, 0));
  else
    Fprintf(fpout, "%.9f", origv);
}


int FindPeriod2(char *Label, llint *Seq, int nSeq, int FixedLen,
                FILE *fp_abc)
{
  int      BestStart, BestLen, BestMatches;
  int          Start,     Len,     Matches, nSeq_Len;
  int      i1, i2, i3, i4;
  ldouble  dx1, dx2, x1, x3, x4, y4, y4_y3, y4_y2, a, b, c;
  int      InitSumMinMax, Valid;
  ldouble  Sum_a, BestSum_a;
  ldouble  Min_a, BestMin_a;
  ldouble  Max_a, BestMax_a;
  ldouble  Sum_b, BestSum_b;
  ldouble  Min_b, BestMin_b;
  ldouble  Max_b, BestMax_b;
  ldouble  Sum_c, BestSum_c;
  ldouble  Min_c, BestMin_c;
  ldouble  Max_c, BestMax_c;


  if (nSeq && Seq[0] == 1)
  {
    /* ATTENTION: Sequence manipulated */

     Seq++;
    nSeq--;
  }

  if (FixedLen && fp_abc)
    Fprintf(fp_abc, ">Begin abc\n");

  BestLen     = 0;
  BestMatches = 4;
  BestStart   = nSeq;

  if (FixedLen)
    Len = FixedLen;
  else
    Len = 1;
                nSeq_Len = nSeq - Len;
  for (       ; nSeq_Len >= BestMatches; Len++, nSeq_Len--)
  {
    dx1 = Len;
    dx2 = Len + Len;

    InitSumMinMax = 1;

    Start = 0;

    for (i4 = nSeq - 1; i4 >= nSeq_Len; i4--)
    {
      i3 = i4 - Len; if (i3 < 0) goto NextLen;
      i2 = i3 - Len; if (i2 < 0) goto NextLen;
      i1 = i2 - Len; if (i1 < Start && FixedLen == 0) break;

      x3 = i3 + 1;
      x4 = i4 + 1;
      y4_y3 = Seq[i4] - Seq[i3];
      y4_y2 = Seq[i4] - Seq[i2];
      y4    = Seq[i4];

      a = (y4_y3 / dx1 - y4_y2 / dx2) / dx1;
      b =  y4_y3 / dx1 - a * (x4 + x3);
      c = y4 - (b + a * x4) * x4;

      if (FixedLen && fp_abc)
      {
        Fprintf(fp_abc, "abc(%d) ", i4 % Len + 1);
        ShowCoeff(fp_abc, (double) a, Len, 50); putc(' ',  fp_abc);
        ShowCoeff(fp_abc, (double) b, Len, 50); putc(' ',  fp_abc);
        ShowCoeff(fp_abc, (double) c, Len, 50); putc('\n', fp_abc);
      }

      if (InitSumMinMax)
      {
        Sum_a = a; Min_a = Max_a = a;
        Sum_b = b; Min_b = Max_b = b;
        Sum_c = c; Min_c = Max_c = c;
        InitSumMinMax = 0;
      }
      else
      {
        Sum_a += a; if (Min_a > a) Min_a = a; if (Max_a < a) Max_a = a;
        Sum_b += b; if (Min_b > b) Min_b = b; if (Max_b < b) Max_b = b;
        Sum_c += c; if (Min_c > c) Min_c = c; if (Max_c < c) Max_c = c;
      }

      c += .5; /* rounding */

      for ( ; i1 >= Start; i1 -= Len)
      {
                                    x1 = i1 + 1;
        if (Seq[i1] != (llint)((a * x1 + b) * x1 + c))
        {
          Start = i1 + 1;
          if (Start >= BestStart) goto NextLen;
          break;
        }
      }
    }

    Matches = nSeq - 3 * Len - Start;

    if (    BestMatches <  Matches
        || (BestMatches == Matches && BestStart > Start))
    {
      BestLen     = Len;
      BestMatches = Matches;
      BestStart   = Start;
      BestSum_a   = Sum_a;
      BestMin_a   = Min_a;
      BestMax_a   = Max_a;
      BestSum_b   = Sum_b;
      BestMin_b   = Min_b;
      BestMax_b   = Max_b;
      BestSum_c   = Sum_c;
      BestMin_c   = Min_c;
      BestMax_c   = Max_c;
    }

    NextLen:

    if (FixedLen) break;
  }

  Valid = 0;

  if (BestLen)
  {
    Fprintf(stdout, "%s: P{%d %d %d}",
      Label, BestStart + 1, BestLen, BestMatches);

    if (BestMatches >= BestLen)
      Valid = 1;

    if (Valid || FixedLen)
      Fprintf(stdout,
        " a{%.9f %.9f %.9f} b{%.9f %.9f %.9f} c{%.9f %.9f %.9f}",
        (double) BestMin_a,
        (double) BestMax_a,
        (double)(BestSum_a / BestLen),
        (double) BestMin_b,
        (double) BestMax_b,
        (double)(BestSum_b / BestLen),
        (double) BestMin_c,
        (double) BestMax_c,
        (double)(BestSum_c / BestLen));

    if (Valid == 0)
    {
      Fprintf(stdout, " Vague");
      Valid = -1;
    }

    putc('\n', stdout);
  }
  else
    Fprintf(stdout, "%s: No Period found\n", Label);

  if (FixedLen && fp_abc)
  {
    Fprintf(fp_abc, ">End abc\n");
    Fflush(fp_abc);
  }

  Fflush(stdout);

  if (Valid)
    return BestLen;

  return 0;
}


static void llint_PrintFraction(FILE *fpout, llint n, llint d)
{
  if (d < 0) {
    n *= -1;
    d *= -1;
  }

  Fprintf(fpout, "%s",  llint_to_str(n, NULL, 0)); if (d != 1)
  Fprintf(fpout, "/%s", llint_to_str(d, NULL, 0));
}


void llintMeanFrac(T_llintFrac *Seq, int nSeq, int Start, int Len,
                   T_llintFrac *Mean)
{
  int          iSeq;
  T_llintFrac  p, q;


  Mean->n = 0;
  Mean->d = 1;

  for (iSeq = Start; iSeq < nSeq && iSeq < Len; iSeq++)
  {
    p.n = Seq[iSeq].n;
    p.d = Len;
    llint_RemoveGCD2(&p.n, &p.d);
    p.d *= Seq[iSeq].d;

    q.d = llint_FindLCM2(p.d, Mean->d);
    if (q.d < p.d || q.d < Mean->d) { Mean->n = Mean->d = 0; break; }
    q.n = p.n * (q.d / p.d) + Mean->n * (q.d / Mean->d);
    llint_RemoveGCD2(&q.n, &q.d);

    Mean->n = q.n;
    Mean->d = q.d;
  }
}


int FindPeriod3Dai(char *Label, T_llintFrac *Seq, int nSeq, int MinMatches)
{
  T_FP_Info    FPI;
  T_llintFrac  Mean_a;


  (void) FindPeriod0Frac(NULL, Seq, nSeq, MinMatches, &FPI);
  llintMeanFrac(Seq, nSeq, FPI.BestStart, FPI.BestLen, &Mean_a);

  if (Label)
  {
    if (FPI.BestLen)
    {
      Fprintf(stdout, "# 3Dai{%d %d %d} a{",
        FPI.BestStart, FPI.BestLen, FPI.BestMatches);
      llint_PrintFraction(stdout, Mean_a.n, Mean_a.d);
      Fprintf(stdout, "} %s\n", Label);
    }
    else
      Fprintf(stdout, "# 3Dai{No Period found} %s\n", Label);

    Fflush(stdout);
  }

  return FPI.BestLen;
}


int iabc(int k1, llint Nk1,
         int k2, llint Nk2,
         int k3, llint Nk3,
         T_llintFrac *a, T_llintFrac *b, T_llintFrac *c)
{
  int          Len;
  T_llintFrac  p, q;


  Len = k2 - k1;

  /* a = (Nk3 - 2 * Nk2 + Nk1) / (2 * Len * Len)
   */

  a->n = Nk3 - Nk2 - Nk2 + Nk1;
  a->d = (llint) 2 * Len * Len;
  llint_RemoveGCD2(&a->n, &a->d);

  /* b =  (Nk3 - Nk2) / Len - a * (k3 + k2);
   */

  p.n = Nk3 - Nk2;
  p.d = Len;
  llint_RemoveGCD2(&p.n, &p.d);

  q.n = (llint) k3 + k2;
  q.d = a->d;
  llint_RemoveGCD2(&q.n, &q.d);

  q.n *= a->n;
  b->d = llint_FindLCM2(p.d, q.d);
  if (b->d < p.d || b->d < q.d) return -1;
  b->n = p.n * (b->d / p.d) - q.n * (b->d / q.d);
  llint_RemoveGCD2(&b->n, &b->d);

  /* c = Nk3 - (b + a * k3) * k3;
   */

  p.n = (llint) k3;
  p.d = a->d;
  llint_RemoveGCD2(&p.n, &p.d);
  p.n *= a->n;

  q.d = llint_FindLCM2(b->d, p.d);
  if (q.d < b->d || q.d < p.d) return -1;
  q.n = b->n * (q.d / b->d) + p.n * (q.d / p.d);
  llint_RemoveGCD2(&q.n, &q.d);

  p.n = (llint) k3;
  llint_RemoveGCD2(&p.n, &q.d);
  q.n *= p.n;

  c->n = Nk3 * q.d - q.n;;
  c->d = q.d;
  llint_RemoveGCD2(&c->n, &c->d);

  return 0;
}


#define Size_aiBuf  256


int iFindPeriod2(char *Label, llint *Seq, int nSeq, int FixedLen,
                 FILE *fp_abc)
{
  int      iStart, Len, nDef, Best_iStart, BestLen, MinDef, Valid;
  int      i1, i2, i3, i4;

  T_llintFrac  a, b, c, p, q, Mean_a, BestMean_a;
  int          OutOfRange;
  int          InitMean;

  int          i, j;
  int          n_aiBuf, Best_n_aiBuf;
  T_llintFrac       aiBuf[Size_aiBuf];
  T_llintFrac  Best_aiBuf[Size_aiBuf];


  if (nSeq && Seq[0] == 1)
  {
    /* ATTENTION: Sequence manipulated */

     Seq++;
    nSeq--;
  }

  OutOfRange = 0;

  if (FixedLen && fp_abc)
    Fprintf(fp_abc, ">Begin abc\n");

  Best_iStart = 0;
  BestLen     = 0;
  MinDef      = nSeq - 5;

  if (FixedLen)
    Len = FixedLen;
  else
    Len = 1;

  while (Len * 3 < MinDef)
  {
    InitMean = 1;
    n_aiBuf  = 0;

    iStart = 0;

    for (i4 = nSeq - 1; i4 >= nSeq - Len; i4--)
    {
      i3 = i4 - Len; if (i3 < 0) goto NextLen;
      i2 = i3 - Len; if (i2 < 0) goto NextLen;
      i1 = i2 - Len; if (i1 < iStart && FixedLen == 0) break;

      if (iabc(i2 + 1, Seq[i2],
               i3 + 1, Seq[i3],
               i4 + 1, Seq[i4],
               &a, &b, &c)      != 0)
      {
        OutOfRange++;
        goto NextLen;
      }

      if (FixedLen && fp_abc)
      {
        Fprintf(fp_abc, "abc(%d) ", i4 % Len + 1);
        llint_PrintFraction(fp_abc, a.n, a.d); putc(' ',  fp_abc);
        llint_PrintFraction(fp_abc, b.n, b.d); putc(' ',  fp_abc);
        llint_PrintFraction(fp_abc, c.n, c.d); putc('\n', fp_abc);
      }

      if (n_aiBuf < Size_aiBuf)
      {
        aiBuf[n_aiBuf].n = a.n;
        aiBuf[n_aiBuf].d = a.d;
              n_aiBuf++;
      }

      if (InitMean)
      {
        p.n = a.n;
        p.d = Len;
        llint_RemoveGCD2(&p.n, &p.d);
        Mean_a.n = p.n;
        Mean_a.d = p.d * a.d;

        InitMean = 0;
      }
      else
      {
        p.n = a.n;
        p.d = Len;
        llint_RemoveGCD2(&p.n, &p.d);
        p.d *= a.d;

        q.d = llint_FindLCM2(p.d, Mean_a.d);
        if (q.d < p.d || q.d < Mean_a.d) { OutOfRange++; goto NextLen; }
        q.n = p.n * (q.d / p.d) + Mean_a.n * (q.d / Mean_a.d);
        llint_RemoveGCD2(&q.n, &q.d);

        Mean_a.n = q.n;
        Mean_a.d = q.d;
      }

      for ( ; i1 >= iStart; i1 -= Len)
      {
        /* Seq[i1] != ((a * (i1 + 1) + b) * (i1 + 1) + c)
         */

        p.n = (llint) i1 + 1;
        p.d = a.d;
        llint_RemoveGCD2(&p.n, &p.d);
        p.n *= a.n;

        q.d = llint_FindLCM2(p.d, b.d);
        if (q.d < p.d || q.d < b.d) { OutOfRange++; goto NextLen; }
        q.n = p.n * (q.d / p.d) + b.n * (q.d / b.d);
        llint_RemoveGCD2(&q.n, &q.d);

        p.n = (llint) i1 + 1;
        llint_RemoveGCD2(&p.n, &q.d);
        q.n *= p.n;

        p.d = llint_FindLCM2(q.d, c.d);
        if (p.d < q.d || p.d < c.d) { OutOfRange++; goto NextLen; }
        p.n = q.n * (p.d / q.d) + c.n * (p.d / c.d);
        llint_RemoveGCD2(&p.n, &p.d);

        if (   (p.d !=  1 || Seq[i1] !=  p.n)
            && (p.d != -1 || Seq[i1] != -p.n))
        {
          iStart = i1 + 1;
          break;
        }
      }
    }
                 nDef = iStart + Len * 3;
    if (MinDef > nDef)
    {
        Best_iStart  = iStart;
        BestLen      = Len;
        MinDef       = nDef;
        BestMean_a.n = Mean_a.n;
        BestMean_a.d = Mean_a.d;
        Best_n_aiBuf = n_aiBuf;

        for (i = 0; i < n_aiBuf; i++)
        {
          Best_aiBuf[i].n = aiBuf[i].n;
          Best_aiBuf[i].d = aiBuf[i].d;
        }
    }

    NextLen:

    if (FixedLen)
      break;
    else
      Len++;
  }

  if (OutOfRange)
    Fprintf(stdout, "%s: Warning: %d time(s) out of integer range\n",
      Label, OutOfRange);

  Valid = 0;

  if (BestLen)
  {
    if (nSeq - MinDef >= BestLen)
      Valid = 1;

    if (Valid || FixedLen)
    {
          Len = 2 * BestLen;
      if (Len > Size_aiBuf)
          Len = Size_aiBuf;

      for (i = Best_n_aiBuf, j = 0; i < Len; i++, j++)
      {
        if (j == Best_n_aiBuf) j = 0;

        Best_aiBuf[i].n = Best_aiBuf[j].n;
        Best_aiBuf[i].d = Best_aiBuf[j].d;
      }

      i = 0; if (Label[0] == '#' && Label[1] == ' ') i = 2;

      Len = FindPeriod3Dai(&Label[i], Best_aiBuf, Len, 1);
    }

    Fprintf(stdout, "%s: P{%d %d %d}",
      Label, Best_iStart + 1, BestLen, nSeq - MinDef);

    if (Valid || FixedLen)
    {
      Fprintf(stdout, " <a>=");
      llint_PrintFraction(stdout, BestMean_a.n, BestMean_a.d);

      Fprintf(stdout, " P(a)=%d", Len);
    }

    if (Valid == 0)
    {
      Fprintf(stdout, " Vague");
      Valid = -1;
    }

    putc('\n', stdout);
  }
  else
    Fprintf(stdout, "%s: No Period found\n", Label);

  if (FixedLen && fp_abc)
  {
    Fprintf(fp_abc, ">End abc\n");
    Fflush(fp_abc);
  }

  Fflush(stdout);

  if (Valid)
    return BestLen;

  return 0;
}
