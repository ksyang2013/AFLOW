// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *                                                                         *
// ***************************************************************************
// Written by Stefano Curtarolo, David Hicks - 2016

#ifndef _AFLOW_ANRL_LIST_CPP
#define _AFLOW_ANRL_LIST_CPP
#include "aflow.h"


// ***************************************************************************
namespace anrl { 
  uint PrototypeANRL_LoadList(vector<string>& vproto,
			      vector<string>& vproto_label,
			      vector<uint>& vproto_nspecies,
			      vector<uint>& vproto_natoms,
			      vector<uint>& vproto_spacegroup,
			      vector<uint>& vproto_nunderscores,
			      vector<uint>& vproto_nparameters,
			      vector<string>& vproto_Pearson_symbol,
			      vector<string>& vproto_params,
			      vector<string>& vproto_Strukturbericht,
			      vector<string>& vproto_prototype,
			      vector<string>& vproto_dialect) {
    
    vproto.clear();
    vproto_label.clear();
    vproto_nspecies.clear();
    vproto_natoms.clear();
    vproto_spacegroup.clear();
    vproto_nunderscores.clear();
    vproto_nparameters.clear();
    vproto_Pearson_symbol.clear();
    vproto_params.clear();
    vproto_Strukturbericht.clear();
    vproto_prototype.clear();
    vproto_dialect.clear();

    
    //Label     # of Species    # atoms/primcell    #space_group_number   # Underscores     # Parameters    pearson_symbol    params    Strukturbericht     Prototype   Dialect        
    vproto.push_back("AB2_aP12_1_4a_8a;2;12;1;4;42;aP12;a,b/a,c/a,alpha,beta,gamma,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12;-;FeS2;FeS2");
    vproto.push_back("ABC2_aP16_1_4a_4a_8a;3;16;1;5;54;aP16;a,b/a,c/a,alpha,beta,gamma,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15,x16,y16,z16;-;AsKSe2;AsKSe2");
    vproto.push_back("A2B_aP6_2_2i_i;2;6;2;4;15;aP6;a,b/a,c/a,alpha,beta,gamma,x1,y1,z1,x2,y2,z2,x3,y3,z3;-;P2I4;P2I4");
    vproto.push_back("A_aP4_2_aci;1;4;2;3;9;aP4;a,b/a,c/a,alpha,beta,gamma,x3,y3,z3;-;Cf;Cf");
    vproto.push_back("A2B_mP12_3_bc3e_2e;2;12;3;4;21;mP12;a,b/a,c/a,beta,y1,y2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7;-;SiO2;SiO2");
    vproto.push_back("A_mP4_4_2a;1;4;4;3;10;mP4;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2;-;Te;High-Pressure Te");
    vproto.push_back("A_mC12_5_3c;1;6;5;3;13;mC12;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3;A19;Po;Po");
    vproto.push_back("A3BC_mC10_8_ab_a_a;3;5;8;5;13;mC10;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3,x4,y4,z4;-;Pb(Zr0.52Ti0.48)O3;Monoclinic PZT [PbO3]");
    vproto.push_back("A2B_mC144_9_24a_12a;2;72;9;4;112;mC144;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15,x16,y16,z16,x17,y17,z17,x18,y18,z18,x19,y19,z19,x20,y20,z20,x21,y21,z21,x22,y22,z22,x23,y23,z23,x24,y24,z24,x25,y25,z25,x26,y26,z26,x27,y27,z27,x28,y28,z28,x29,y29,z29,x30,y30,z30,x31,y31,z31,x32,y32,z32,x33,y33,z33,x34,y34,z34,x35,y35,z35,x36,y36,z36;-;SiO2;Monoclinic Low Tridymite");
    vproto.push_back("AB_mP4_11_e_e;2;4;11;4;8;mP4;a,b/a,c/a,beta,x1,z1,x2,z2;-;NiTi;NiTi2");
    vproto.push_back("ABC3_mP10_11_e_e_ef;3;10;11;5;13;mP10;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3,x4,y4,z4;G0_6;KClO3;KClO3");
    vproto.push_back("A_mP16_11_8e;1;16;11;3;20;mP16;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3,x4,z4,x5,z5,x6,z6,x7,z7,x8,z8;-;alphaPu;alphaPu");
    vproto.push_back("AB2_mC6_12_a_i;2;3;12;4;6;mC6;a,b/a,c/a,beta,x2,z2;C34;AuTe2;Calaverite");
    vproto.push_back("A_mC34_12_ah3i2j;1;17;12;3;17;mC34;a,b/a,c/a,beta,y2,x3,z3,x4,z4,x5,z5,x6,y6,z6,x7,y7,z7;-;betaPu;betaPu");
    vproto.push_back("AB3_mC16_12_g_ij;2;8;12;4;10;mC16;a,b/a,c/a,beta,y1,x2,z2,x3,y3,z3;D0_15;AlCl3;AlCl3");
    vproto.push_back("A5B2_mC14_12_a2i_i;2;7;12;4;10;mC14;a,b/a,c/a,beta,x2,z2,x3,z3,x4,z4;-;Au5Mn2;Au5Mn2");
    vproto.push_back("A_mC4_12_i;1;2;12;3;6;mC4;a,b/a,c/a,beta,x1,z1;-;O2;alphaO2");
    vproto.push_back("ABC4_mP12_13_e_a_2g;3;12;13;5;11;mP12;a,b/a,c/a,beta,y2,x3,y3,z3,x4,y4,z4;E1_b;AgAuTe4;Sylvanite");
    vproto.push_back("A_mP84_13_21g;1;84;13;3;67;mP84;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15,x16,y16,z16,x17,y17,z17,x18,y18,z18,x19,y19,z19,x20,y20,z20,x21,y21,z21;-;P;Monoclinic Phosphorus");
    vproto.push_back("A2B_mP12_14_2e_e;2;12;14;4;13;mP12;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3;C43;ZrO2;Baddeleyite");
    vproto.push_back("A_mP32_14_8e;1;32;14;3;28;mP32;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8;A_l;Se;betaSe");
    vproto.push_back("A_mP64_14_16e;1;64;14;3;52;mP64;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15,x16,y16,z16;A_k;Se;Se");
    vproto.push_back("A2B5_mC28_15_f_e2f;2;14;15;4;14;mC28;a,b/a,c/a,beta,y1,x2,y2,z2,x3,y3,z3,x4,y4,z4;-;B2Pd5;B2Pd5");
    vproto.push_back("AB_mC8_15_c_e;2;4;15;4;5;mC8;a,b/a,c/a,beta,y2;B26;CuO;Tenorite");
    vproto.push_back("A2B_mC48_15_ae3f_2f;2;24;15;4;20;mC48;a,b/a,c/a,beta,y2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7;-;SiO2;Coesite");
    vproto.push_back("ABC6D2_mC40_15_e_e_3f_f;4;20;15;6;18;mC40;a,b/a,c/a,beta,y1,y2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6;-;CaFeO6Si2;Esseneite");
    vproto.push_back("ABC4_oP12_16_ag_cd_2u;3;12;16;5;9;oP12;a,b/a,c/a,x5,y5,z5,x6,y6,z6;-;AlPS4;AlPS4");
    vproto.push_back("AB3_oP16_18_ab_3c;2;16;18;4;14;oP16;a,b/a,c/a,z1,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5;-;BaS3;BaS3");
    vproto.push_back("A2B_oP12_19_2a_a;2;12;19;4;12;oP12;a,b/a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3;-;Ag2Se;Naumannite");
    vproto.push_back("A2B_oC24_20_abc_c;2;12;20;4;11;oC24;a,b/a,c/a,x1,y2,x3,y3,z3,x4,y4,z4;-;SiO2;Orthorhombic Tridymite");
    vproto.push_back("AB_oP2_25_b_a;2;2;25;4;5;oP2;a,b/a,c/a,z1,z2;-;CdTe;High-Pressure CdTe");
    vproto.push_back("AB2_oP24_28_acd_2c3d;2;24;28;4;22;oP24;a,b/a,c/a,z1,y2,z2,y3,z3,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8;C46;AuTe2;Krennerite");
    vproto.push_back("AB3C4_oP16_31_a_ab_2ab;3;16;31;5;17;oP16;a,b/a,c/a,y1,z1,y2,z2,y3,z3,y4,z4,x5,y5,z5,x6,y6,z6;H2_5;AsCu3S4;Enargite");
    vproto.push_back("AB_oP8_33_a_a;2;8;33;4;9;oP8;a,b/a,c/a,x1,y1,z1,x2,y2,z2;-;CoAs;Modderite");
    vproto.push_back("AB3C4_oP32_33_a_3a_4a;3;32;33;5;27;oP32;a,b/a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8;-;AsK3S4;AsK3S4");
    vproto.push_back("A2B_oC12_36_2a_a;2;6;36;4;9;oC12;a,b/a,c/a,y1,z1,y2,z2,y3,z3;C24;HgBr2;HgBr2");
    vproto.push_back("A2BC_oC8_38_e_a_b;3;4;38;5;7;oC8;a,b/a,c/a,z1,z2,y3,z3;-;C2CeNi;C2CeNi");
    vproto.push_back("A2B_oC12_38_de_ab;2;6;38;4;9;oC12;a,b/a,c/a,z1,z2,y3,z3,y4,z4;-;Au2V;Au2V");
    vproto.push_back("AB4_oC20_41_a_2b;2;10;41;4;10;oC20;a,b/a,c/a,z1,x2,y2,z2,x3,y3,z3;D1_c;PtSn4;PtSn4");
    vproto.push_back("AB2_oC24_41_2a_2b;2;12;41;4;11;oC24;a,b/a,c/a,z1,z2,x3,y3,z3,x4,y4,z4;C_e;PdSn2;PdSn2");
    vproto.push_back("AB2_oF72_43_ab_3b;2;18;43;4;16;oF72;a,b/a,c/a,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5;C44;GeS2;GeS2");
    vproto.push_back("AB_oI4_44_a_b;2;2;44;4;5;oI4;a,b/a,c/a,z1,z2;-;AsGa;High-pressure GaAs");
    vproto.push_back("A2B3C7D_oP13_47_t_aq_eqrs_h;4;13;47;6;8;oP13;a,b/a,c/a,z4,z5,z6,z7,z8;-;YBa2Cu3O7-x;1212C [YBa2Cu3O7-x]");
    vproto.push_back("AB_oP4_51_e_f;2;4;51;4;5;oP4;a,b/a,c/a,z1,z2;B19;AuCd;beta'-AuCd");
    vproto.push_back("A3B2_oP20_56_ce_e;2;20;56;4;10;oP20;a,b/a,c/a,z1,x2,y2,z2,x3,y3,z3;D5_11;Sb2O3;Sb2O3");
    vproto.push_back("ABCD_oP16_57_d_c_d_d;4;16;57;6;10;oP16;a,b/a,c/a,x1,x2,y2,x3,y3,x4,y4;F5_9;KCNS;KCNS");
    vproto.push_back("AB_oP8_57_d_d;2;8;57;4;7;oP8;a,b/a,c/a,x1,y1,x2,y2;-;TlF;TlF-II, 57");
    vproto.push_back("AB2_oP6_58_a_g;2;6;58;4;5;oP6;a,b/a,c/a,x2,y2;C35/-/C18;CaCl2/etaFe2C/FeS2;Hydrophilite/eta-Fe2C/Marcasite");
    vproto.push_back("AB_oP4_59_a_b;2;4;59;4;5;oP4;a,b/a,c/a,z1,z2;-;CuTe;Vulcanite");
    vproto.push_back("ABC_oP6_59_a_a_a;3;6;59;5;6;oP6;a,b/a,c/a,z1,z2,z3;-;CNCl;CNCl");
    vproto.push_back("A3B_oP8_59_bf_a;2;8;59;4;7;oP8;a,b/a,c/a,z1,z2,x3,z3;D0_a;TiCu3;betaTiCu3");
    vproto.push_back("AB_oP16_61_c_c;2;16;61;4;9;oP16;a,b/a,c/a,x1,y1,z1,x2,y2,z2;B_e;CdSb;CdSb");
    vproto.push_back("A2B_oP24_61_2c_c;2;24;61;4;12;oP24;a,b/a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3;C21;TiO2;Brookite");
    vproto.push_back("A3B2_oP20_62_3c_2c;2;20;62;4;13;oP20;a,b/a,c/a,x1,z1,x2,z2,x3,z3,x4,z4,x5,z5;D5_8;Sb2S3;Stibnite");
    vproto.push_back("AB3C_oP20_62_c_cd_a;3;20;62;5;10;oP20;a,b/a,c/a,x2,z2,x3,z3,x4,y4,z4;-;CaTiO3;CaTiO3 Pnma Perovskite");
    vproto.push_back("A4B_oP20_62_2cd_c;2;20;62;4;12;oP20;a,b/a,c/a,x1,z1,x2,z2,x3,z3,x4,y4,z4;-;MgB4;MgB4");
    vproto.push_back("AB2C_oP16_62_c_2c_c;3;16;62;5;11;oP16;a,b/a,c/a,x1,z1,x2,z2,x3,z3,x4,z4;F5_6;CuSbS2;Chalcostibite");
    vproto.push_back("A2B_oP12_62_2c_c;2;12;62;4;9;oP12;a,b/a,c/a,x1,z1,x2,z2,x3,z3;C37/C25/C23;Co2Si/HgCl2/PbCl2;Co2Si/HgCl2/Cotunnite");
    vproto.push_back("AB_oP8_62_c_c;2;8;62;4;7;oP8;a,b/a,c/a,x1,z1,x2,z2;B16/B31/B27/B29;GeS/MnP/FeB/SnS;GeS/MnP/FeB/SnS");
    vproto.push_back("AB3_oP16_62_c_cd;2;16;62;4;10;oP16;a,b/a,c/a,x1,z1,x2,z2,x3,y3,z3;D0_11;Fe3C;Cementite");
    vproto.push_back("A3B7_oP40_62_cd_3c2d;2;40;62;4;20;oP40;a,b/a,c/a,x1,z1,x2,z2,x3,z3,x4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7;D10_1;C3Cr7;C3Cr7");
    vproto.push_back("A_oP8_62_2c;1;8;62;3;7;oP8;a,b/a,c/a,x1,z1,x2,z2;A_c;alphaNp;alphaNp");
    vproto.push_back("AB2C_oC16_63_c_2c_c;3;8;63;5;7;oC16;a,b/a,c/a,y1,y2,y3,y4;-;SrCuO2;SrCuO2");
    vproto.push_back("A2B_oC12_63_2c_c;2;6;63;4;6;oC12;a,b/a,c/a,y1,y2,y3;C49;ZrSi2;ZrSi2");
    vproto.push_back("AB_oC8_63_c_c;2;4;63;4;5;oC8;a,b/a,c/a,y1,y2;B33;CrB;CrB");
    vproto.push_back("A_oC4_63_c;1;2;63;3;4;oC4;a,b/a,c/a,y1;A20;alpha-U;alpha-Uranium");
    vproto.push_back("A_oC8_64_f;1;4;64;3;5;oC8;a,b/a,c/a,y1,z1;A11/A17/A14;alphaGa/P/I2;alphaGallium/Black Phosphorus/Molecular Iodine");
    vproto.push_back("A2B2C_oC80_64_efg_efg_df;3;40;64;5;18;oC80;a,b/a,c/a,x1,y2,y3,y4,z4,y5,z5,y6,z6,x7,y7,z7,x8,y8,z8;-;MgB2C2;MgB2C2");
    vproto.push_back("AB_oC8_65_j_g;2;4;65;4;5;oC8;a,b/a,c/a,x1,y2;-;alphaIrV;alphaIrV");
    vproto.push_back("A3B5_oC16_65_ah_bej;2;8;65;4;5;oC16;a,b/a,c/a,x4,y5;-;Ga3Pt5;Ga3Pt5");
    vproto.push_back("AB3_oC8_65_a_bf;2;4;65;4;3;oC8;a,b/a,c/a;L1_3;CdPt3;Predicted CdPt3");
    vproto.push_back("AB_oF8_69_a_b;2;2;69;4;3;oF8;a,b/a,c/a;B24;TlF;TlF");
    vproto.push_back("A_oF8_70_a;1;2;70;3;3;oF8;a,b/a,c/a;-;gammaPu;gamma-Pu");
    vproto.push_back("A2B_oF24_70_e_a;2;6;70;4;4;oF24;a,b/a,c/a,x2;C54;TiSi2;TiSi2");
    vproto.push_back("A_oF128_70_4h;1;32;70;3;15;oF128;a,b/a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4;A16;S;alpha-Sulfur");
    vproto.push_back("AB2_oI6_71_a_i;2;3;71;4;4;oI6;a,b/a,c/a,z2;-;ReSi2;ReSi2");
    vproto.push_back("AB2_oI6_71_a_g;2;3;71;4;4;oI6;a,b/a,c/a,y2;-;MoPt2;MoPt2");
    vproto.push_back("A2B_oI12_72_j_a;2;6;72;4;5;oI12;a,b/a,c/a,x2,y2;C42;SiS2;SiS2");
    vproto.push_back("AB4C_tI12_82_c_g_a;3;6;82;5;5;tI12;a,c/a,x3,y3,z3;H0_7;BPO4;BPO4");
    vproto.push_back("A2BC4_tI14_82_bc_a_g;3;7;82;5;5;tI14;a,c/a,x4,y4,z4;E3;CdAl2S4;CdAl2S4");
    vproto.push_back("AB_tP16_84_cej_k;2;16;84;4;7;tP16;a,c/a,x3,y3,x4,y4,z4;B34;PdS;PdS");
    vproto.push_back("A4B5_tI18_87_h_ah;2;9;87;4;6;tI18;a,c/a,x2,y2,x3,y3;-;Ti5Te4;Ti5Te4");
    vproto.push_back("AB4_tI10_87_a_h;2;5;87;4;4;tI10;a,c/a,x2,y2;D1_a;Ni4Mo;Ni4Mo");
    vproto.push_back("A2B_tP12_92_b_a;2;12;92;4;6;tP12;a,c/a,x1,x2,y2,z2;-;SiO2;alpha Cristobalite");
    vproto.push_back("A2B_tP36_96_3b_ab;2;36;96;4;15;tP36;a,c/a,x1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5;-;SiO2;Keatite");
    vproto.push_back("A_tP12_96_ab;1;12;96;3;6;tP12;a,c/a,x1,x2,y2,z2;-;Si;\"ST12'' of Si");
    vproto.push_back("A3BC_tP5_99_bc_a_b;3;5;99;5;6;tP5;a,c/a,z1,z2,z3,z4;-;Pb(Zr0.52Ti0.48)O3;Tetragonal PZT [PbO3]");
    vproto.push_back("AB3_tP8_113_a_ce;2;8;113;4;5;tP8;a,c/a,z2,x3,z3;D0_17;BaS3;BaS3");
    vproto.push_back("A2BC4D_tI16_121_d_a_i_b;4;8;121;6;4;tI16;a,c/a,x4,z4;H2_6;Cu2FeS4Sn;Stannite");
    vproto.push_back("ABC2_tI16_122_a_b_d;3;8;122;5;3;tI16;a,c/a,x3;E1_1;CuFeS2;Chalcopyrite");
    vproto.push_back("AB5C_tP7_123_b_ci_a;3;7;123;5;3;tP7;a,c/a,z4;-;HoCoGa5;HoCoGa5");
    vproto.push_back("AB3_tP4_123_a_ce;2;4;123;4;2;tP4;a,c/a;L6_0;CuTi3;CuTi3");
    vproto.push_back("AB_tP2_123_a_d;2;2;123;4;2;tP2;a,c/a;L1_0;CuAu;CuAu");
    vproto.push_back("ABC2_tP4_123_d_a_f;3;4;123;5;2;tP4;a,c/a;-;CaCuO2;CaCuO2");
    vproto.push_back("A2B3_tP10_127_g_ah;2;10;127;4;4;tP10;a,c/a,x2,x3;D5_a;Si2U3;Si2U3");
    vproto.push_back("ABCD_tP8_129_c_b_a_c;4;8;129;6;4;tP8;a,c/a,z3,z4;-;AsCuSiZr;AsCuSiZr");
    vproto.push_back("A_tP4_129_ac;1;4;129;3;3;tP4;a,c/a,z2;A_d;betaNp;betaNp");
    vproto.push_back("ABC_tP6_129_c_a_c;3;6;129;5;4;tP6;a,c/a,z2,z3;E0_1;PbFCl;Matlockite");
    vproto.push_back("A2B_tP6_129_ac_c;2;6;129;4;4;tP6;a,c/a,z2,z3;C38;Cu2Sb;Cu2Sb");
    vproto.push_back("AB_tP4_129_a_c;2;4;129;4;3;tP4;a,c/a,z2;B10;PbO;PbO");
    vproto.push_back("AB_tP4_129_c_c;2;4;129;4;4;tP4;a,c/a,z1,z2;B11;gammaCuTi;gammaCuTi");
    vproto.push_back("AB_tP4_131_c_e;2;4;131;4;2;tP4;a,c/a;B17;PtS;PtS");
    vproto.push_back("A_tP50_134_b2m2n;1;50;134;3;12;tP50;a,c/a,x2,z2,x3,z3,x4,y4,z4,x5,y5,z5;A_g;B;T-50 Boron");
    vproto.push_back("A_tP30_136_bf2ij;1;30;136;3;9;tP30;a,c/a,x2,x3,y3,x4,y4,x5,z5;A_b;betaU;betaU");
    vproto.push_back("AB_tP8_136_g_f;2;8;136;4;4;tP8;a,c/a,x1,x2;-;betaBeO;betaBeO");
    vproto.push_back("A2B_tP6_136_f_a;2;6;136;4;3;tP6;a,c/a,x2;C4;TiO2;CrSi2");
    //vproto.push_back("sigma_tP30_136_bf2ij;5;30;136;3;9;tP30;a,c/a,x2,x3,y3,x4,y4,x5,z5;D8_b;CrFe;sigmaCrFe");
    vproto.push_back("sigma_tP30_136_bf2ij;1;30;136;3;9;tP30;a,c/a,x2,x3,y3,x4,y4,x5,z5;D8_b;CrFe;sigmaCrFe");
    vproto.push_back("A_tP4_136_f;1;4;136;3;3;tP4;a,c/a,x1;-;N2;gammaN2");
    vproto.push_back("A_tP16_138_j;1;16;138;3;5;tP16;a,c/a,x1,y1,z1;A18;Cl2;Cl2");
    vproto.push_back("A3B_tI16_139_cde_e;2;8;139;4;4;tI16;a,c/a,z3,z4;D0_23;Al3Zr;Al3Zr");
    vproto.push_back("A_tI4_139_e;1;2;139;3;3;tI4;a,c/a,z1;-;Si;Hypothetical BCT5 Si");
    vproto.push_back("AB2C4_tI14_139_a_e_ce;3;7;139;5;4;tI14;a,c/a,z3,z4;-;(La,Ba)2CuO4;0201 [(La,Ba)2CuO4]");
    vproto.push_back("A12B_tI26_139_fij_a;2;13;139;4;4;tI26;a,c/a,x3,x4;D2_b;Mn12Th;Mn12Th");
    vproto.push_back("A_tI2_139_a;1;1;139;3;2;tI2;a,c/a;A6/A_a;In/Pa;In/alphaPa");
    vproto.push_back("A_tI8_139_h;1;4;139;3;3;tI8;a,c/a,x1;-;C;Hypothetical Tetrahedrally Bonded Carbon with 4-Member Rings");
    vproto.push_back("A3B_tI8_139_bd_a;2;4;139;4;2;tI8;a,c/a;D0_22;Al3Ti;Al3Ti");
    vproto.push_back("AB2_tI6_139_a_e;2;3;139;4;3;tI6;a,c/a,z2;C11_b;MoSi2;MoSi2");
    vproto.push_back("A4B5_tI18_139_i_ah;2;9;139;4;4;tI18;a,c/a,x2,x3;-;V4Zn5;V4Zn5");
    vproto.push_back("A4B_tI10_139_de_a;2;5;139;4;3;tI10;a,c/a,z3;D1_3;Al4Ba;Al4Ba");
    vproto.push_back("A8B_tI18_139_hi_a;2;9;139;4;4;tI18;a,c/a,x2,x3;-;Pt8Ti;Pt8Ti");
    vproto.push_back("A2B_tI6_139_d_a;2;3;139;4;2;tI6;a,c/a;L\'2;ThH2;ThH2");
    vproto.push_back("A2B_tI12_140_h_a;2;6;140;4;3;tI12;a,c/a,x2;C16;Al2Cu;Khatyrkite");
    vproto.push_back("AB3_tI16_140_b_ah;2;8;140;4;3;tI16;a,c/a,x3;D0_c;SiU3;SiU3");
    vproto.push_back("AB_tI16_140_ab_h;2;8;140;4;3;tI16;a,c/a,x3;B37;SeTl;SeTl");
    vproto.push_back("A4BC_tI24_141_h_b_a;3;12;141;5;4;tI24;a,c/a,y3,z3;-;ZrSiO4;Zircon");
    vproto.push_back("A_tI4_141_a;1;2;141;3;2;tI4;a,c/a;A5;Sn;betaSn");
    vproto.push_back("A3B4_tI28_141_ad_h;2;14;141;4;4;tI28;a,c/a,y3,z3;-;Mn3O4;Hausmannite");
    vproto.push_back("A2B_tI12_141_e_a;2;6;141;4;3;tI12;a,c/a,z2;C5;TiO2;Anatase");
    vproto.push_back("AB_tI16_141_e_e;2;8;141;4;4;tI16;a,c/a,z1,z2;B_g;MoB;MoB");
    vproto.push_back("A2B_tI24_141_2e_e;2;12;141;4;5;tI28;a,c/a,z1,z2,z3;-;Ga2Hf;Ga2Hf");
    vproto.push_back("AB_tI8_141_a_b;2;4;141;4;2;tI8;a,c/a;\"40\";NbP;NbP");
    vproto.push_back("A2B3_tI80_141_ceh_3h;2;40;141;4;11;tI80;a,c/a,z2,y3,z3,y4,z4,y5,z5,y6,z6;-;In2S3;betaIn2S3");
    vproto.push_back("ABC4_tI96_142_e_ab_2g;3;48;142;5;9;tI96;a,c/a,x3,x4,y4,z4,x5,y5,z5;-;PPrS4;PPrS4");
    vproto.push_back("A2B_hP9_147_g_ad;2;9;147;4;6;hP9;a,c/a,z2,x3,y3,z3;B_b;AgZn;zetaAgZn");
    vproto.push_back("AB_hR16_148_cf_cf;2;16;148;4;10;hR16;a,c/a,x1,x2,x3,y3,z3,x4,y4,z4;-;C8H8;Solid Cubane");
    vproto.push_back("AB3_hR8_148_c_f;2;8;148;4;6;hR8;a,c/a,x1,x2,y2,z2;D0_5;BiI3;BiI3");
    vproto.push_back("AB_hR26_148_b2f_a2f;2;26;148;4;14;hR26;a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4;-;PdAl;PdAl");
    vproto.push_back("AB3C_hR10_148_c_f_c;3;10;148;5;7;hR10;a,c/a,x1,x2,x3,y3,z3;-;FeTiO3;Ilmenite");
    vproto.push_back("A2B_hP9_150_ef_bd;2;9;150;4;5;hP9;a,c/a,z2,x3,x4;C22;Fe2P;Original Fe2P");
    vproto.push_back("A3B_hP24_151_3c_2a;2;24;151;4;13;hP24;a,c/a,x1,x2,x3,y3,z3,x4,y4,z4,x5,y5,z5;D0_4;CrCl3;CrCl3");
    vproto.push_back("A2B_hP9_152_c_a;2;9;152;4;6;hP9;a,c/a,x1,x2,y2,z2;-;SiO2;alphaQuartz");
    vproto.push_back("A_hP3_152_a;1;3;152;3;3;hP3;a,c/a,x1;A8;Se;gammaSe");
    vproto.push_back("AB_hP6_154_a_b;2;6;154;4;4;hP6;a,c/a,x1,x2;B9;HgS;Cinnabar");
    vproto.push_back("AB3_hR8_155_c_de;2;8;155;4;5;hR8;a,c/a,x1,y2,y3;D0_14;AlF3;AlF3");
    vproto.push_back("A3B2_hR5_155_e_c;2;5;155;4;4;hR5;a,c/a,x1,y2;D5_e;Ni3S2;Hazelwoodite");
    vproto.push_back("AB_hR6_160_b_b;2;6;160;4;6;hR6;a,c/a,x1,z1,x2,z2;B13;NiS;Millerite");
    vproto.push_back("AB_hR6_160_3a_3a;2;6;160;4;8;hR6;a,c/a,x1,x2,x3,x4,x5,x6;-;CSi;Moissanite 9R structure");
    vproto.push_back("ABC3_hR10_161_a_a_b;3;10;161;5;7;hR10;a,c/a,x1,x2,x3,y3,z3;-;LiNbO3;Ferroelectric LiNbO3");
    vproto.push_back("AB2_hP9_162_ad_k;2;9;162;4;4;hP9;a,c/a,x3,z3;-;NV2;betaV2N");
    vproto.push_back("AB2CD2_hP36_163_h_i_bf_i;4;36;163;6;10;hP36;a,c/a,z2,x3,x4,y4,z4,x5,y5,z5;F5_10;KAg(CN)2;KAg2");
    vproto.push_back("A3B2_hP5_164_ad_d;2;5;164;4;4;hP5;a,c/a,z2,z3;D5_13;Al3Ni2;Al3Ni2");
    vproto.push_back("AB2_hP3_164_a_d;2;3;164;4;3;hP3;a,c/a,z2;C6;CdI2;omega Phase");
    vproto.push_back("A3B_hP24_165_adg_f;2;24;165;4;7;hP24;a,c/a,z2,x3,x4,y4,z4;-;H3Ho;H3Ho");
    vproto.push_back("AB_hR2_166_a_b;2;2;166;4;2;hR2;a,c/a;L1_1;CuPt;CuPt");
    vproto.push_back("A_hR2_166_c;1;2;166;3;3;hR2;a,c/a,x1;A7/-/-;As/C/O2;alphaAs/Rhombohedral Graphite/betaO2");
    vproto.push_back("A_hR1_166_a;1;1;166;3;2;hR1;a,c/a;A_i/A10;Po/Hg;betaPo/alphaHg");
    vproto.push_back("A7B6_hR13_166_ah_3c;2;13;166;4;7;hR13;a,c/a,x2,x3,x4,x5,z5;D8_5;Fe7W6;Fe7W6 mu-phase");
    vproto.push_back("A_hR3_166_ac;1;3;166;3;3;hR3;a,c/a,x2;C19;Sm;alphaSm");
    vproto.push_back("A2B3_hR5_166_c_ac;2;5;166;4;4;hR5;a,c/a,x2,x3;C33;Bi2Te3;Bi2Te3");
    vproto.push_back("A5B2_hR7_166_a2c_c;2;7;166;4;5;hR7;a,c/a,x2,x3,x4;D8_i;B5Mo2;Mo2B5");
    vproto.push_back("A_hR12_166_2h;1;12;166;3;6;hR12;a,c/a,x1,z1,x2,z2;-;B;alphaBoron");
    vproto.push_back("ABC2_hR4_166_a_b_c;3;4;166;5;3;hR4;a,c/a,x3;F5_1;CrNiS2;Caswellsilverite");
    vproto.push_back("A_hR105_166_bc9h4i;1;105;166;3;33;hR105;a,c/a,z2,x3,z3,x4,z4,x5,z5,x6,z6,x7,z7,x8,z8,x9,z9,x10,z10,x11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15;-;B;betaBoron");
    vproto.push_back("A6B_hR7_166_g_a;2;7;166;4;3;hR7;a,c/a,x2;-;CaC6;CaC6");
    vproto.push_back("ABC3_hR10_167_a_b_e;3;10;167;5;3;hR10;a,c/a,x3;-/G0_1;LiNbO3/CaCO3;Paraelectric LiNbO3/Calcite");
    vproto.push_back("A2B3_hR10_167_c_e;2;10;167;4;4;hR10;a,c/a,x1,x2;D5_1;Al2O3;Corundum");
    vproto.push_back("A2B_hP18_180_fi_bd;2;18;180;4;4;hP18;a,c/a,z3,x4;C_a;Mg2Ni;Mg2Ni");
    vproto.push_back("AB2_hP9_180_d_j;2;9;180;4;3;hP9;a,c/a,x2;C40;CrSi2;CrSi2");
    vproto.push_back("A2B_hP9_180_j_c;2;9;180;4;3;hP9;a,c/a,x2;C8;SiO2;beta-Quartz");
    vproto.push_back("AB3_hP8_182_c_g;2;8;182;4;3;hP8;a,c/a,x2;-;Fe3C;Bainite");
    vproto.push_back("A_hP4_186_ab;1;4;186;3;4;hP4;a,c/a,z1,z2;-;C;Buckled Graphite");
    vproto.push_back("AB_hP8_186_ab_ab;2;8;186;4;6;hP8;a,c/a,z1,z2,z3,z4;B5;SiC;Moissanite-4H SiC");
    vproto.push_back("AB_hP4_186_b_b;2;4;186;4;4;hP4;a,c/a,z1,z2;B4;ZnS;Wurtzite");
    vproto.push_back("AB_hP12_186_a2b_a2b;2;12;186;4;8;hP12;a,c/a,z1,z2,z3,z4,z5,z6;B6;SiC;Moissanite-6H SiC");
    vproto.push_back("A5B3C_hP18_186_2a3b_2ab_b;3;18;186;5;11;hP18;a,c/a,z1,z2,z3,z4,z5,z6,z7,z8,z9;E9_4;Al5C3N;Al5C3N");
    vproto.push_back("AB_hP4_186_b_a;2;4;186;4;4;hP4;a,c/a,z1,z2;B12;BN;Original BN");
    vproto.push_back("ABC_hP3_187_a_d_f;3;3;187;5;2;hP3;a,c/a;-;BaPtSb;BaPtSb");
    vproto.push_back("AB_hP2_187_d_a;2;2;187;4;2;hP2;a,c/a;B_h;WC;Tungsten Carbide");
    vproto.push_back("A2B_hP9_189_fg_bc;2;9;189;4;4;hP9;a,c/a,x3,x4;C22;Fe2P;Revised Fe2P");
    vproto.push_back("AB4C_hP6_191_a_h_b;3;6;191;5;3;hP6;a,c/a,z3;-;AlB4Mg;AlB4Mg");
    vproto.push_back("AB5_hP6_191_a_cg;2;6;191;4;2;hP6;a,c/a;D2_d;CaCu5;CaCu5");
    vproto.push_back("A_hP1_191_a;1;1;191;3;2;hP1;a,c/a;A_f;gammaHgSn6-10;Simple Hexagonal Lattice");
    vproto.push_back("A3B_hP4_191_bc_a;2;4;191;4;2;hP4;a,c/a;-;Li3N;Li3 Ni");
    vproto.push_back("AB2_hP3_191_a_d;2;3;191;4;2;hP3;a,c/a;C32;AlB2;Hexagonal omega");
    vproto.push_back("A2B_hP6_191_h_e;2;6;191;4;4;hP6;a,c/a,z1,z2;C_h;Cu2Te;Cu2Te");
    vproto.push_back("AB_hP6_191_f_ad;2;6;191;4;2;hP6;a,c/a;B35;CoSn;CoSn");
    vproto.push_back("AB_hP8_194_ad_f;2;8;194;4;3;hP8;a,c/a,z3;B_i;AsTi;AsTi");
    vproto.push_back("A_hP6_194_h;1;6;194;3;3;hP6;a,c/a,x1;-;C;Hypothetical Tetrahedrally Bonded Carbon with 3-Member Rings");
    vproto.push_back("AB_hP12_194_af_bf;2;12;194;4;4;hP12;a,c/a,z3,z4;-;CMo;CMo structure");
    vproto.push_back("A_hP4_194_ac;1;4;194;3;2;hP4;a,c/a;A3\';La;alphaLa");
    vproto.push_back("AB3_hP8_194_c_bf;2;8;194;4;3;hP8;a,c/a,z3;D0_18;Na3As;Na3As");
    vproto.push_back("AB2_hP6_194_b_f;2;6;194;4;3;hP6;a,c/a,z2;-;CaIn2;CaIn2");
    vproto.push_back("AB_hP4_194_c_d;2;4;194;4;2;hP4;a,c/a;B_k;BN;BN");
    vproto.push_back("ABC2_hP8_194_d_a_f;3;8;194;5;3;hP8;a,c/a,z3;-;AlCCr2;AlCCr2");
    vproto.push_back("A3B_hP8_194_h_c;2;8;194;4;3;hP8;a,c/a,x2;D0_19;Ni3Sn;Ni3Sn");
    vproto.push_back("A_hP4_194_bc;1;4;194;3;2;hP4;a,c/a;A9;C;Hexagonal Graphite structure");
    vproto.push_back("AB2_hP6_194_c_f;2;6;194;4;3;hP6;a,c/a,z2;C7;MoS2;Molybdenite");
    vproto.push_back("A5B2_hP14_194_abdf_f;2;14;194;4;4;hP14;a,c/a,z4,z5;D8_h;W2B5;W2B5");
    vproto.push_back("AB2_hP12_194_f_ah;2;12;194;4;4;hP12;a,c/a,z2,x3;C14;MgZn2;MgZn2 Hexagonal Laves");
    vproto.push_back("ABC_hP6_194_c_d_a;3;6;194;5;2;hP6;a,c/a;-;LiBC;LiBC");
    vproto.push_back("A_hP4_194_f;1;4;194;3;3;hP4;a,c/a,z1;-;C;Lonsdaleite");
    vproto.push_back("AB2_hP6_194_c_ad;2;6;194;4;2;hP6;a,c/a;B8_2;Ni2In;Ni2In");
    vproto.push_back("AB3C4_hP16_194_c_af_ef;3;16;194;5;5;hP16;a,c/a,z3,z4,z5;-;AlN3Ti4;AlN3Ti4");
    vproto.push_back("A_hP2_194_c;1;2;194;3;2;hP2;a,c/a;A3;Mg;Hexagonal Close Packed");
    vproto.push_back("AB2_hP24_194_ef_fgh;2;24;194;4;6;hP24;a,c/a,z1,z2,z3,x5;C36;MgNi2;MgNi2 Hexagonal Laves");
    vproto.push_back("AB_hP12_194_df_ce;2;12;194;4;4;hP12;a,c/a,z3,z4;B18;CuS;Covellite");
    vproto.push_back("AB_hP4_194_c_a;2;4;194;4;2;hP4;a,c/a;B8_1;NiAs;NiAs");
    vproto.push_back("A2B_hP12_194_cg_f;2;12;194;4;3;hP12;a,c/a,z2;C10;SiO2;beta-Tridymite");
    vproto.push_back("A4B_cI40_197_cde_c;2;20;197;4;5;cI40;a,x1,x2,x3,x4;-;Ga4Ni;Ga4Ni3");
    vproto.push_back("ABC_cP12_198_a_a_a;3;12;198;5;4;cP12;a,x1,x2,x3;F0_1;NiSSb;Ullmanite");
    vproto.push_back("A3B_cP16_198_b_a;2;16;198;4;5;cP16;a,x1,x2,y2,z2;D1;NH3;Ammonia");
    vproto.push_back("A_cP8_198_2a;1;8;198;3;3;cP8;a,x1,x2;-;N2;alpha-N2");
    vproto.push_back("AB_cP8_198_a_a;2;8;198;4;3;cP8;a,x1,x2;B21/B20;CO/FeSi;alphaCO/FeSi");
    vproto.push_back("AB_cI16_199_a_a;2;8;199;4;3;cI16;a,x1,x2;B_a;CoU;CoU");
    vproto.push_back("AB32C48_cI162_204_a_2efg_2gh;3;81;204;5;13;cI162;a,x2,x3,x4,y5,z5,y6,z6,y7,z7,x8,y8,z8;-;Mg32(Al,Zn)49;Bergman [Mg32(Al,Zn)49]");
    vproto.push_back("A3B_cI32_204_g_c;2;16;204;4;3;cI32;a,y2,z2;D0_2;CoAs3;Skutterudite");
    vproto.push_back("A12B_cI26_204_g_a;2;13;204;4;3;cI26;a,y2,z2;-;Al12W;Al12W");
    vproto.push_back("A_cP8_205_c;1;8;205;3;2;cP8;a,x1;-;N2;alpha-N2");
    vproto.push_back("AB_cP16_205_c_c;2;16;205;4;3;cP16;a,x1,x2;-;CuCl;SC16");
    vproto.push_back("AB2_cP12_205_a_c;2;12;205;4;2;cP12;a,x2;C2;FeS2;Pyrite");
    vproto.push_back("AB3C6_cI80_206_a_d_e;3;40;206;5;5;cI80;a,x2,x3,y3,z3;D5_3;(Mn,Fe)2O3;Bixbyite");
    vproto.push_back("A_cI16_206_c;1;8;206;3;2;cI16;a,x1;-;Si;BC8");
    vproto.push_back("A_cP20_213_cd;1;20;213;3;3;cP20;a,x1,y2;A13;Mn;betaMn");
    vproto.push_back("A3B4C_cP8_215_d_e_a;3;8;215;5;2;cP8;a,x3;H2_4;Cu3S4V;Sulvanite");
    vproto.push_back("AB4_cP5_215_a_e;2;5;215;4;2;cP5;a,x2;-;Fe4C;Fe4C");
    vproto.push_back("AB3C4_cP8_215_a_c_e;3;8;215;5;2;cP8;a,x3;-;Cu3AsS4;Cubic Lazarevicite");
    vproto.push_back("AB5_cF24_216_a_ce;2;6;216;4;2;cF24;a,x3;C15_b;AuBe5;AuBe5");
    vproto.push_back("ABC_cF12_216_b_c_a;3;3;216;5;1;cF12;a;C1_b;AgAsMg;Half-Heusler");
    vproto.push_back("AB_cF8_216_c_a;2;2;216;4;1;cF8;a;B3;ZnS;NaTl");
    vproto.push_back("A4B_cI10_217_c_a;2;5;217;4;2;cI10;a,x2;-;SiF4;SiF4");
    vproto.push_back("A_cI58_217_ac2g;1;29;217;3;6;cI58;a,x2,x3,z3,x4,z4;A12;Mn;alphaMn");
    vproto.push_back("A5B8_cI52_217_ce_cg;2;26;217;4;6;cI52;a,x1,x2,x3,x4,z4;-;Cu5Zn8;gammaBrass");
    vproto.push_back("A_cI16_220_c;1;8;220;3;2;cI16;a,x1;-;Li;High Pressure cI16 Li");
    vproto.push_back("A3B2_cI40_220_d_c;2;20;220;4;3;cI40;a,x1,x2;D5_c;Pu2C3;Pu2C3");
    vproto.push_back("AB_cP2_221_b_a;2;2;221;4;1;cP2;a;B2;CsCl;CsCl");
    vproto.push_back("AB_cP6_221_c_d;2;6;221;4;1;cP6;a;-;NbO;NbO");
    vproto.push_back("AB3C_cP5_221_a_c_b;3;5;221;5;1;cP5;a;E2_1;CaTiO3;Cubic Perovskite");
    vproto.push_back("AB27CD3_cP32_221_a_dij_b_c;4;32;221;6;3;cP32;a,y5,y6;-;CrFe2525Ni6;Model of Austenite");
    vproto.push_back("AB3_cP4_221_a_c;2;4;221;4;1;cP4;a;L1_2;Cu3Au;Cu3Au");
    vproto.push_back("A_cP1_221_a;1;1;221;3;1;cP1;a;A_h;Po;alphaPo");
    vproto.push_back("AB11_cP36_221_c_agij;2;36;221;4;4;cP36;a,x3,y4,y5;D2_e;BaHg11;BaHg11");
    vproto.push_back("AB11CD3_cP16_221_a_dg_b_c;4;16;221;6;2;cP16;a,x5;-;CrFe11MoNi3;Model of Ferrite");
    vproto.push_back("A3B_cP4_221_d_a;2;4;221;4;1;cP4;a;D0_9;ReO3;alphaReO3");
    vproto.push_back("A6B_cP7_221_f_a;2;7;221;4;2;cP7;a,x2;D2_1;CaB6;CaB6");
    vproto.push_back("A3B_cP8_223_c_a;2;8;223;4;1;cP8;a;A15;Cr3Si;Cr3Si");
    vproto.push_back("A_cP46_223_dik;1;46;223;3;4;cP46;a,x2,y3,z3;-;Si;Si46 Clathrate");
    vproto.push_back("A2B_cP6_224_b_a;2;6;224;4;1;cP6;a;C3;Cu2O;Cuprite");
    vproto.push_back("A7B_cF32_225_bd_a;2;8;225;4;1;cF32;a;-;Ca7Ge;Ca7Ge");
    vproto.push_back("AB3_cF16_225_a_bc;2;4;225;4;1;cF16;a;D0_3;BiF3;BiF3");
    vproto.push_back("A9B16C7_cF128_225_acd_2f_be;3;32;225;5;4;cF128;a,x5,x6,x7;-;Cr9Fe16Ni7;Model of Ferrite");
    vproto.push_back("A12B_cF52_225_i_a;2;13;225;4;2;cF52;a,y2;D2_f;UB12;UB12");
    vproto.push_back("AB2_cF12_225_a_c;2;3;225;4;1;cF12;a;C1;CaF2;Cu2Mg Cubic Laves");
    vproto.push_back("A6B23_cF116_225_e_acfh;2;29;225;4;4;cF116;a,x3,x4,y5;D8_4;Cr23C6;Cr23C6");
    vproto.push_back("AB2C_cF16_225_a_c_b;3;4;225;5;1;cF16;a;L2_1;AlCu2Mn;Heusler");
    vproto.push_back("A_cF4_225_a;1;1;225;3;1;cF4;a;A1;Cu;Face-Centered Cubic");
    vproto.push_back("AB18C8_cF108_225_a_eh_f;3;27;225;5;4;cF108;a,x2,x3,y4;-;CrFe18Ni8;Model of Austenite");
    vproto.push_back("AB_cF8_225_a_b;2;2;225;4;1;cF8;a;B1;NaCl;Rock Salt");
    vproto.push_back("A2B_cF24_227_c_a;2;6;227;4;1;cF24;a;C9;SiO2;Ideal beta-Cristobalite");
    vproto.push_back("AB2_cF96_227_e_cf;2;24;227;4;3;cF96;a,x2,x3;-;NiTi2;NiTi2");
    vproto.push_back("AB_cF16_227_a_b;2;4;227;4;1;cF16;a;B32;NaTl;NaTl");
    vproto.push_back("A_cF136_227_aeg;1;34;227;3;4;cF136;a,x2,x3,z3;-;Si;Si34 Clathrate");
    vproto.push_back("A2B_cF24_227_d_a;2;6;227;4;1;cF24;a;C15;Cu2Mg;Cu2Mg Cubic Laves");
    vproto.push_back("A_cF8_227_a;1;2;227;3;1;cF8;a;A4;C;Diamond");
    vproto.push_back("A2BC4_cF56_227_d_a_e;3;14;227;5;2;cF56;a,x3;H1_1;Al2MgO4;Spinel");
    vproto.push_back("AB2_cF48_227_c_e;2;12;227;4;2;cF48;a,x2;-;CTi2;CTi2");
    vproto.push_back("AB3C3_cF112_227_c_de_f;3;28;227;5;3;cF112;a,x3,x4;E9_3;Fe3W3C;Fe3W3C");
    vproto.push_back("A_cI2_229_a;1;1;229;3;1;cI2;a;A2;W;Body-Centered Cubic");
    vproto.push_back("A3B_cI8_229_b_a;2;4;229;4;1;cI8;a;-;H3S;High Presssure H3S");
    vproto.push_back("A4B3_cI14_229_c_b;2;7;229;4;1;cI14;a;-;Pt3O4;Pt3O4");
    vproto.push_back("A2B7_cI54_229_e_afh;2;27;229;4;4;cI54;a,x2,x3,y4;L2_2;Sb2Tl7;L22");
    vproto.push_back("AB12C3_cI32_229_a_h_b;3;16;229;5;2;cI32;a,y3;-;CrFe12Ni3;Model of Austenite");
    vproto.push_back("AB4C3_cI16_229_a_c_b;3;8;229;5;1;cI16;a;-;CrFe4Ni3;Model of Ferrite");
    vproto.push_back("A4B3_cI112_230_af_g;2;56;230;4;3;cI112;a,x2,y3;-;Ga4Ni3;Ga4Ni3");

    // done now produce
    
    // FROM PROTO LIST
    for(uint i=0;i<vproto.size();i++) {
      vproto_label.push_back("");
      vproto_nspecies.push_back(0);
      vproto_natoms.push_back(0);
      vproto_spacegroup.push_back(0);
      vproto_nunderscores.push_back(0);
      vproto_nparameters.push_back(0);
      vproto_Pearson_symbol.push_back("");
      vproto_params.push_back("");
      vproto_Strukturbericht.push_back("");
      vproto_prototype.push_back("");
      vproto_dialect.push_back("");
      
      anrl::vproto2tokens(vproto.at(i),
			  vproto_label.at(i),vproto_nspecies.at(i),vproto_natoms.at(i),vproto_spacegroup.at(i),
			  vproto_nunderscores.at(i),vproto_nparameters.at(i),vproto_Pearson_symbol.at(i),
			  vproto_params.at(i),vproto_Strukturbericht.at(i),vproto_prototype.at(i),vproto_dialect.at(i));
    }
    
    return vproto.size();
  }
}

#endif // _AFLOW_ANRL_LIST_CPP

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *                                                                         *
// ***************************************************************************

