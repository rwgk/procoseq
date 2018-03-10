#ifndef FPERIOD_H__
#define FPERIOD_H__


typedef struct
  {
    int  BestLen;
    int  BestMatches;
    int  BestStart;
  }
  T_FP_Info;


typedef struct
  {
    llint  n;
    llint  d;
  }
  T_llintFrac;


int FindPeriod0(char *Label, llint *Seq, int nSeq, int MinMatches,
                T_FP_Info *FPI);
int FindPeriod0Frac(char *Label, T_llintFrac *Seq, int nSeq, int MinMatches,
                    T_FP_Info *FPI);
int FindPeriod1(char *Label, llint *Seq, int nSeq, int FixedLen);
void llintMeanFrac(T_llintFrac *Seq, int nSeq, int Start, int Len,
                   T_llintFrac *Mean);
int FindPeriod3Dai(char *Label, T_llintFrac *Seq, int nSeq, int MinMatches);
int iabc(int k1, llint Nk1,
         int k2, llint Nk2,
         int k3, llint Nk3,
         T_llintFrac *a, T_llintFrac *b, T_llintFrac *c);
int  FindPeriod2(char *Label, llint *Seq, int nSeq, int FixedLen,
                 FILE *fp_abc);
int iFindPeriod2(char *Label, llint *Seq, int nSeq, int FixedLen,
                 FILE *fp_abc);


typedef int (*T_FindPeriod2)(char *, llint *, int, int, FILE *);


#endif
