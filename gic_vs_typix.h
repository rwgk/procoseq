#ifndef GIC_VS_TYPIX_H__
#define GIC_VS_TYPIX_H__


typedef struct
  {
    char  *TYPIX_Label;
    char  *gic_Label;
  }
  T_TYPIX_LblTab;


typedef struct
  {
    char            *Code;
    char            *Name;
    char            *SearchCode;
    T_TYPIX_LblTab  *TLT;
  }
  T_gic_vs_typix;


T_TYPIX_LblTab TLT_CaF2[] =
  {
    { "Ca position", "T000" },
    { "F position",  "T001" },
    { NULL, NULL }
  };

T_TYPIX_LblTab TLT_NaCl[] =
  {
    { "Na & Cl position", "T000" },
    { NULL, NULL }
  };

T_TYPIX_LblTab TLT_FeS2[] =
  {
    { "Fe position", "T000" },
    { "S position",  "T001" },
    { NULL, NULL }
  };

T_TYPIX_LblTab TLT_NiAs[] =
  {
    { "Ni position", "T000" },
    { "As position", "T001" },
    { NULL, NULL }
  };

T_TYPIX_LblTab TLT_Me[] =
  {
    { "", "T000" },
    { NULL, NULL }
  };

T_TYPIX_LblTab TLT_aNd[] =
  {
    { "Position Nd2", "T000" },
    { "Position Nd1", "T001" },
    { NULL, NULL }
  };

T_TYPIX_LblTab TLT_Ni2In[] =
  {
    { "Position Ni2",      "T000" },
    { "Position Ni1 & In", "T001" },
    { NULL, NULL }
  };

T_TYPIX_LblTab TLT_aMn[] =
  {
    { "Position Mn1", "T003" },
    { "Position Mn2", "T002" },
    { "Position Mn3", "T001" },
    { "Position Mn4", "T000" },
    { NULL, NULL }
  };

T_TYPIX_LblTab TLT_Cr3Si[] =
  {
    { "Si position", "T000" },
    { "Cr position", "T001" },
    { NULL, NULL }
  };

T_TYPIX_LblTab TLT_sCrFe[] =
  {
    { "Position Xa", "T001" },
    { "Position Xb", "T004" },
    { "Position Xc", "T002" },
    { "Position Xd", "T003" },
    { "Position Xf", "T000" },
    { NULL, NULL }
  };

T_TYPIX_LblTab TLT_MgZn2[] =
  {
    { "Position Zn1", "T001" },
    { "Position Zn2", "T000" },
    { "Mg position",  "T002" },
    { NULL, NULL }
  };

T_TYPIX_LblTab TLT_MgCu2[] =
  {
    { "Mg position", "T000" },
    { "Cu position", "T001" },
    { NULL, NULL }
  };

T_TYPIX_LblTab TLT_MgNi2[] =
  {
    { "Position Mg1", "T001" },
    { "Position Mg2", "T000" },
    { "Position Ni1", "T004" },
    { "Position Ni2", "T003" },
    { "Position Ni3", "T002" },
    { NULL, NULL }
  };


T_gic_vs_typix gic_vs_typix[] =
  {
    { "CaF2_2", "CaF2(2)",        "(225) cF12", TLT_CaF2  },
    { "CaF2_1", "CaF2(1)",        "(225) cF12", TLT_CaF2  },
    { "NaCl",   "NaCl",           "(225) cF8",  TLT_NaCl  },
    { "mar",    "FeS2-Marcasite", "(58) oP6",   TLT_FeS2  },
    { "pyr",    "FeS2-Pyrite",    "(205) cP12", TLT_FeS2  },
    { "NiAs_2", "NiAs(2)",        "(194) hP4",  TLT_NiAs  },
    { "NiAs_1", "NiAs(1)",        "(194) hP4",  TLT_NiAs  },
    { "Cu",     "Cu",             "(225) cF4",  TLT_Me    },
    { "Mg",     "Mg",             "(194) hP2",  TLT_Me    },
    { "W_2",    "W(2)",           "(229) cI2",  TLT_Me    },
    { "W_1",    "W(1)",           "(229) cI2",  TLT_Me    },
    { "aLa",    "alpha-Nd",       "(194) hP4",  TLT_aNd   },
    { "Ni2In",  "Ni2In",          "(194) hP6",  TLT_Ni2In },
    { "aMn",    "alpha-Mn",       "(217) cI58", TLT_aMn   },
    { "Cr3Si",  "Cr3Si",          "(223) cP8",  TLT_Cr3Si },
    { "sCrFe",  "sigma-CrFe",     "(136) tP30", TLT_sCrFe },
    { "MgZn2",  "MgZn2",          "(194) hP12", TLT_MgZn2 },
    { "MgCu2",  "MgCu2",          "(227) cF24", TLT_MgCu2 },
    { "MgNi2",  "MgNi2",          "(194) hP24", TLT_MgNi2 },
    { NULL,     NULL,             NULL,         NULL      }
  };

#endif
