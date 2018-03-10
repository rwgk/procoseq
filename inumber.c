#include "longer.h"


int FindGCD(const int *S, int nS)
{
  int  ri, rj, rk;


  if (nS-- == 0) return 0;

      ri = *S++;
  if (ri < 0) ri = -ri;

  while (nS--)
  {
    if ((rj = *S++) != 0)
    {
      for (;;)
      {
        rk = ri % rj; if (rk == 0) { ri = rj; break; }
        ri = rj % rk; if (ri == 0) { ri = rk; break; }
        rj = rk % ri; if (rj == 0) {          break; }
      }

      if (ri < 0) ri = -ri;

      if (ri == 1) break;
    }
  }

  return ri;
}


llint llint_FindGCD(const llint *S, int nS)
{
  llint  ri, rj, rk;


  if (nS-- == 0) return 0;

      ri = *S++;
  if (ri < 0) ri = -ri;

  while (nS--)
  {
    if ((rj = *S++) != 0)
    {
      for (;;)
      {
        rk = ri % rj; if (rk == 0) { ri = rj; break; }
        ri = rj % rk; if (ri == 0) { ri = rk; break; }
        rj = rk % ri; if (rj == 0) {          break; }
      }

      if (ri < 0) ri = -ri;

      if (ri == 1) break;
    }
  }

  return ri;
}


int FindGCD2(int ri, int rj)
{
  int  rk;


  if (ri < 0) ri = -ri;

  if (rj)
  {
    for (;;)
    {
      rk = ri % rj; if (rk == 0) { ri = rj; break; }
      ri = rj % rk; if (ri == 0) { ri = rk; break; }
      rj = rk % ri; if (rj == 0) {          break; }
    }

    if (ri < 0) ri = -ri;
  }

  return ri;
}


llint llint_FindGCD2(llint ri, llint rj)
{
  llint  rk;


  if (ri < 0) ri = -ri;

  if (rj)
  {
    for (;;)
    {
      rk = ri % rj; if (rk == 0) { ri = rj; break; }
      ri = rj % rk; if (ri == 0) { ri = rk; break; }
      rj = rk % ri; if (rj == 0) {          break; }
    }

    if (ri < 0) ri = -ri;
  }

  return ri;
}


int FindLCM(const int *S, int nS)
{
  int  a, b, ri, rj, rk;


  if (nS-- == 0) return 1;

      ri = *S++;
  if (ri == 0) ri = 1;

  while (nS--)
  {
    if ((rj = *S++) != 0)
    {
      a = ri; b = rj;

      for (;;)
      {
        rk = ri % rj; if (rk == 0) { ri = rj; break; }
        ri = rj % rk; if (ri == 0) { ri = rk; break; }
        rj = rk % ri; if (rj == 0) {          break; }
      }

      ri = a / ri * b;
    }
  }

  if (ri < 0) return -ri;
              return  ri;
}


llint llint_FindLCM(const llint *S, llint nS)
{
  llint  a, b, ri, rj, rk;


  if (nS-- == 0) return 1;

      ri = *S++;
  if (ri == 0) ri = 1;

  while (nS--)
  {
    if ((rj = *S++) != 0)
    {
      a = ri; b = rj;

      for (;;)
      {
        rk = ri % rj; if (rk == 0) { ri = rj; break; }
        ri = rj % rk; if (ri == 0) { ri = rk; break; }
        rj = rk % ri; if (rj == 0) {          break; }
      }

      ri = a / ri * b;
    }
  }

  if (ri < 0) return -ri;
              return  ri;
}


int FindLCM2(int ri, int rj)
{
  int  a, b, rk;


  if (ri == 0) ri = 1;

  if (rj)
  {
    a = ri; b = rj;

    for (;;)
    {
      rk = ri % rj; if (rk == 0) { ri = rj; break; }
      ri = rj % rk; if (ri == 0) { ri = rk; break; }
      rj = rk % ri; if (rj == 0) {          break; }
    }

    ri = a / ri * b;
  }

  if (ri < 0) return -ri;
              return  ri;
}


llint llint_FindLCM2(llint ri, llint rj)
{
  llint  a, b, rk;


  if (ri == 0) ri = 1;

  if (rj)
  {
    a = ri; b = rj;

    for (;;)
    {
      rk = ri % rj; if (rk == 0) { ri = rj; break; }
      ri = rj % rk; if (ri == 0) { ri = rk; break; }
      rj = rk % ri; if (rj == 0) {          break; }
    }

    ri = a / ri * b;
  }

  if (ri < 0) return -ri;
              return  ri;
}


void llint_RemoveGCD2(llint *n, llint *d)
{
  llint  ri, rj, rk;


  ri = *n;
  rj = *d;

  if (rj)
  {
    for (;;)
    {
      rk = ri % rj; if (rk == 0) { ri = rj; break; }
      ri = rj % rk; if (ri == 0) { ri = rk; break; }
      rj = rk % ri; if (rj == 0) {          break; }
    }
  }

  if (ri < 0) ri = -ri;

  if (ri > 1)
  {
    *d /= ri;
    *n /= ri;
  }
}
