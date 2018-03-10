#ifndef INUMBER_H__
#define INUMBER_H__


int FindGCD(const int *S, int nS);
llint llint_FindGCD(const llint *S, int nS);
int FindGCD2(int ri, int rj);
llint llint_FindGCD2(llint ri, llint rj);
int FindLCM(const int *S, int nS);
llint llint_FindLCM(const llint *S, llint nS);
int FindLCM2(int ri, int rj);
llint llint_FindLCM2(llint ri, llint rj);
void llint_RemoveGCD2(llint *n, llint *d);


#endif /* INUMBER_H__ */
