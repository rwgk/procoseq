#ifndef LONGER_H__
#define LONGER_H__


#if   defined(NORMAL_TYPES)

#define llint     long int
#define ldouble   double

#elif defined(__ALPHA) && defined(__VMS)

#include <ints.h>
#define llint     int64
#define ldouble   long double

#elif defined(__sgi) || defined(_AIX)

#define llint     long long
#define ldouble   double

#elif defined(__x86_64__)

#define llint     long long
#define ldouble   long double

#else

#define llint     long int
#define ldouble   double

#endif


char *llint_to_str(llint lli, char *str, int sstr);
int read_llint(char *str, llint *lli);


#endif /* LONGER_H__ */
