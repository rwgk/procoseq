#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>


#include "longer.h"


char *llint_to_str(llint lli, char *str, int sstr)
{
  int          neg, i, j;
#define            sbuf 64
  static char  buf[sbuf];


  if (str == NULL) {
      str =  buf;
     sstr = sbuf;
  }

  if (lli == 0)
  {
    str[0] = '0';
    str[1] = '\0';
  }
  else
  {
    neg = 0;

    if (lli < 0) {
        lli = -lli;
        neg = 1;
    }

    str[sstr - 1] = '\0';

    for (i = sstr - 2; i >= neg; i--)
    {
      if (lli == 0) break;

      switch ((int)(lli % 10)) {
        default: str[i] = '0'; break;
         case 1: str[i] = '1'; break;
         case 2: str[i] = '2'; break;
         case 3: str[i] = '3'; break;
         case 4: str[i] = '4'; break;
         case 5: str[i] = '5'; break;
         case 6: str[i] = '6'; break;
         case 7: str[i] = '7'; break;
         case 8: str[i] = '8'; break;
         case 9: str[i] = '9'; break;
      }

      lli /= 10;
    }

    if (lli == 0)
    {
      if (neg)
        str[i] = '-';
      else
        i++;

      if (i != 0) {
        j = 0; do str[j] = str[i++]; while (str[j++]);
      }
    }
    else
      for (i = 0; i < sstr - 1; i++) str[i] = '*';
  }

  return str;
}


int read_llint(char *str, llint *lli)
{
  int    sign, d;
  llint  val;


  while (*str && isspace(*str)) str++;

  if (*str == '-') {
    sign = -1;
    str++;
  }
  else {
    sign =  1;
    if (*str == '+') str++;
  }

  val = 0;

  for (;;)
  {
    switch (*str) {
      case '0': d = 0; break;
      case '1': d = 1; break;
      case '2': d = 2; break;
      case '3': d = 3; break;
      case '4': d = 4; break;
      case '5': d = 5; break;
      case '6': d = 6; break;
      case '7': d = 7; break;
      case '8': d = 8; break;
      case '9': d = 9; break;
      default: d = -1; break;
    }

    if (d == -1) break;

    val *= 10;
    val += d;

    str++;
  }
  
  while (*str && isspace(*str)) str++;

  *lli = sign * val;

  if (*str) return -1;

  return 0;
}
