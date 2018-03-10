#ifndef GIC_VS_ICSD_H__
#define GIC_VS_ICSD_H__

typedef struct
  {
    char  *Code;
    char  *Name;
    int   CollectionCode;
  }
  T_gic_vs_icsd;


T_gic_vs_icsd gic_vs_icsd[] =
  {
    { "ban", "Banalsite",       4447 },
    { "coe", "Coesite",        18112 },
    { "cor", "Cordierite",     30947 },
    { "cri", "Cristobalite",    9327 },
    { "fel", "Feldspar",      100182 },
    { "kea", "Keatite",        34889 },
    { "mil", "Milarite",       71046 },
    { "mog", "Moganite",       67669 },
    { "par", "Paracelsian",    24690 },
    { "qua", "Quartz",         31048 },
    { "sca", "Scapolite",       9502 },
    { "tri", "Tridymite",      29343 },
    { NULL, NULL, 0 }
  };

#endif
