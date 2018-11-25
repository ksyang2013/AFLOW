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
    // -------------------------------------------------------------------------
    // Part 1
    // -------------------------------------------------------------------------
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
    vproto.push_back("AB3_oP16_18_ab_3c;2;16;18;4;14;oP16;a,b/a,c/a,z1,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5;D0_{17};BaS3;BaS3"); //DX 20180925 - moved Strukturberict designation from AB3_tP8_113_a_ce per M. Mehl's suggestion (historical reasons)
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
    vproto.push_back("A2B_oP12_62_2c_c;2;12;62;4;9;oP12;a,b/a,c/a,x1,z1,x2,z2,x3,z3;C37/C25/C23/C29;Co2Si/HgCl2/PbCl2/SrH2;Co2Si/HgCl2/Cotunnite/SrH2"); //added info from part 2 A2B_oP12_62_2c_c;2;12;62;4;9;oP12;a,b/a,c/a,x1,z1,x2,z2,x3,z3;C29;SrH2;SrH2
    vproto.push_back("AB_oP8_62_c_c;2;8;62;4;7;oP8;a,b/a,c/a,x1,z1,x2,z2;B16/B31/B27/B29/B14;GeS/MnP/FeB/SnS/FeAs;GeS/MnP/FeB/SnS/Westerveldite"); //added info from part 2 AB_oP8_62_c_c;2;8;62;4;7;oP8;a,b/a,c/a,x1,z1,x2,z2;B14;FeAs;Westerveldite
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
    vproto.push_back("AB3_tP8_113_a_ce;2;8;113;4;5;tP8;a,c/a,z2,x3,z3;-;BaS3;BaS3");//DX 20180925 - moved Strukturberict designation to AB3_oP16_18_ab_3c per M. Mehl's suggestion (historical reasons)
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
    vproto.push_back("A2B_tI12_141_e_a;2;6;141;4;3;tI12;a,c/a,z2;C5/C_{c};TiO2/alpha-ThSi2;Anatase/alpha-ThSi2"); // added info from part 2 - A2B_tI12_141_e_a;2;6;141;4;3;tI12;a,c/a,z2;C_{c};alpha-ThSi2;alpha-ThSi2
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
    vproto.push_back("A_hR105_166_bc9h4i;1;105;166;3;33;hR105;a,c/a,x2,x3,z3,x4,z4,x5,z5,x6,z6,x7,z7,x8,z8,x9,z9,x10,z10,x11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15;-;B;betaBoron");
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
    // -------------------------------------------------------------------------
    // Part 2
    // -------------------------------------------------------------------------
    vproto.push_back("A2B_aP6_2_aei_i;2;6;2;4;12;aP6;a,b/a,c/a,alpha,beta,gamma,x3,y3,z3,x4,y4,z4;-;H2S;H2S");
    vproto.push_back("A8B5_mP13_6_a7b_3a2b;2;13;6;4;30;mP13;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3,x4,z4,x5,z5,x6,z6,x7,z7,x8,z8,x9,z9,x10,z10,x11,z11,x12,z12,x13,z13;-;Mo8P5;Mo8P5");
    vproto.push_back("AB_mP4_6_2b_2a;2;4;6;4;12;mP4;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3,x4,z4;-;FeNi;FeNi");
    vproto.push_back("A2B_mP12_7_4a_2a;2;12;7;4;22;mP12;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6;-;H2S;H2S IV");
    vproto.push_back("A2B_mP18_7_6a_3a;2;18;7;4;31;mP18;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9;-;As2Ba;As2Ba");
    vproto.push_back("A3B_mP16_7_6a_2a;2;16;7;4;28;mP16;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8;-;epsilon-WO3;epsilon-WO3");
    vproto.push_back("A9B2_mP22_7_9a_2a;2;22;7;4;37;mP22;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11;-;Rh2Ga9;Rh2Ga9");
    vproto.push_back("A5B3_mC32_9_5a_3a;2;16;9;4;28;mC32;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8;-;alpha-P3N5;alpha-P3N5");
    vproto.push_back("AB3_mC16_9_a_3a;2;8;9;4;16;mC16;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4;-;H3Cl;H3Cl");
    vproto.push_back("A2B_mP6_10_mn_bg;2;6;10;4;8;mP6;a,b/a,c/a,beta,x3,z3,x4,z4;-;delta-PdCl2;delta-PdCl2");
    vproto.push_back("AB3_mP16_10_mn_3m3n;2;16;10;4;20;mP16;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3,x4,z4,x5,z5,x6,z6,x7,z7,x8,z8;-;H3Cl;H3Cl");
    vproto.push_back("ABC2_mP8_10_ac_eh_mn;3;8;10;5;8;mP8;a,b/a,c/a,beta,x5,z5,x6,z6;-;AuAgTe2;Muthmannite");
    vproto.push_back("AB_mP6_10_en_am;2;6;10;4;8;mP6;a,b/a,c/a,beta,x3,z3,x4,z4;-;LiSn;LiSn");
    vproto.push_back("A_mP8_10_2m2n;1;8;10;3;12;mP8;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3,x4,z4;-;C;S-carbon");
    vproto.push_back("A7B2C2_mC22_12_aij_h_i;3;11;12;5;12;mC22;a,b/a,c/a,beta,y2,x3,z3,x4,z4,x5,y5,z5;S2_{1};[Sc,Y]2Si2O7;Thortveitite");
    vproto.push_back("A_mC16_12_4i;1;8;12;3;12;mC16;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3,x4,z4;-;C;M-carbon");
    vproto.push_back("A2B_mP12_13_2g_ef;2;12;13;4;12;mP12;a,b/a,c/a,beta,y1,y2,x3,y3,z3,x4,y4,z4;-;H2S;H2S");
    vproto.push_back("A2B_mP6_14_e_a;2;6;14;4;7;mP6;a,b/a,c/a,beta,x2,y2,z2;-;gamma-PdCl2;gamma-PdCl2");
    vproto.push_back("A7B8_mP120_14_14e_16e;2;120;14;4;94;mP120;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15,x16,y16,z16,x17,y17,z17,x18,y18,z18,x19,y19,z19,x20,y20,z20,x21,y21,z21,x22,y22,z22,x23,y23,z23,x24,y24,z24,x25,y25,z25,x26,y26,z26,x27,y27,z27,x28,y28,z28,x29,y29,z29,x30,y30,z30;-;alpha-C7H8;alpha-Toluene");
    vproto.push_back("AB3_mC16_15_e_cf;2;8;15;4;8;mC16;a,b/a,c/a,beta,y2,x3,y3,z3;-;H3Cl;H3Cl");
    vproto.push_back("A_mC24_15_2e2f;1;12;15;3;12;mC24;a,b/a,c/a,beta,y1,y2,x3,y3,z3,x4,y4,z4;-;H;H-III");
    vproto.push_back("A2B_oP12_17_abe_e;2;12;17;4;11;oP12;a,b/a,c/a,x1,x2,x3,y3,z3,x4,y4,z4;-;alpha-Ag2Se;alpha-Naumannite");
    vproto.push_back("AB3_oP16_19_a_3a;2;16;19;4;15;oP16;a,b/a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4;-;H3Cl;H3Cl");
    vproto.push_back("AB2_oC6_21_a_k;2;3;21;4;4;oC6;a,b/a,c/a,z2;-;Ta2H;Ta2H");
    vproto.push_back("A2BC2_oF40_22_fi_ad_gh;3;10;22;5;7;oF40;a,b/a,c/a,y3,z4,z5,y6;-;CeRu2B2;CeRu2B2");
    vproto.push_back("AB_oF8_22_a_c;2;2;22;4;3;oF8;a,b/a,c/a;-;FeS;FeS");
    vproto.push_back("A3B_oI32_23_ij2k_k;2;16;23;4;14;oI32;a,b/a,c/a,z1,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5;-;H3S;H3S");
    vproto.push_back("A8B2C12D2E_oI50_23_bcfk_i_3k_j_a;5;25;23;7;18;oI50;a,b/a,c/a,x4,z5,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10;-;Cu8(Fe,Zn)3Sn2S12;Stannoidite");
    vproto.push_back("ABC2_oI16_23_ab_i_k;3;8;23;5;7;oI16;a,b/a,c/a,z3,x4,y4,z4;-;NaFeS2;NaFeS2");
    vproto.push_back("ABC4_oI12_23_a_b_k;3;6;23;5;6;oI12;a,b/a,c/a,x3,y3,z3;-;BPS4;BPS4");
    vproto.push_back("AB7CD2_oI44_24_a_b3d_c_ac;4;22;24;6;17;oI44;a,b/a,c/a,x1,x2,y3,z4,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8;-;Na2MgAlF7;Weberite");
    vproto.push_back("A2B_oP12_26_abc_ab;2;12;26;4;14;oP12;a,b/a,c/a,y1,z1,y2,z2,y3,z3,y4,z4,x5,y5,z5;-;H2S;H2S");
    vproto.push_back("A2B_oP12_26_abc_ab;2;12;26;4;14;oP12;a,b/a,c/a,y1,z1,y2,z2,y3,z3,y4,z4,x5,y5,z5;-;beta-SeO2;beta-SeO2");
    vproto.push_back("A5B_oP24_26_3a3b2c_ab;2;24;26;4;25;oP24;a,b/a,c/a,y1,z1,y2,z2,y3,z3,y4,z4,y5,z5,y6,z6,y7,z7,y8,z8,x9,y9,z9,x10,y10,z10;-;TlP5;TlP5");
    vproto.push_back("A6B4C16D_oP108_27_abcd4e_4e_16e_e;4;108;27;6;82;oP108;a,b/a,c/a,z1,z2,z3,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15,x16,y16,z16,x17,y17,z17,x18,y18,z18,x19,y19,z19,x20,y20,z20,x21,y21,z21,x22,y22,z22,x23,y23,z23,x24,y24,z24,x25,y25,z25,x26,y26,z26,x27,y27,z27,x28,y28,z28,x29,y29,z29;-;Ca4Al6O16S;Ca4Al6O16S");
    vproto.push_back("A2B_oP12_29_2a_a;2;12;29;4;12;oP12;a,b/a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3;-;ZrO2;ZrO2");
    vproto.push_back("AB2_oP12_29_a_2a;2;12;29;4;12;oP12;a,b/a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3;-;FeS2;Pyrite");
    vproto.push_back("ABC_oP12_29_a_a_a;3;12;29;5;12;oP12;a,b/a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3;-;CoAsS;Cobaltite");
    vproto.push_back("A5B3C15_oP46_30_a2c_bc_a7c;3;46;30;5;36;oP46;a,b/a,c/a,z1,z2,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13;-;Bi5Nb3O15;Bi5Nb3O15");
    vproto.push_back("ABC3_oP20_30_2a_c_3c;3;20;30;5;17;oP20;a,b/a,c/a,z1,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6;-;CuBrSe3;CuBrSe3");
    vproto.push_back("A13B2C2_oP34_32_a6c_c_c;3;34;32;5;28;oP34;a,b/a,c/a,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9;-;Re2O5[SO4]2;Re2O52");
    vproto.push_back("A2B3_oP40_33_4a_6a;2;40;33;4;33;oP40;a,b/a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10;-;kappa-Al2O3;kappa-alumina");
    vproto.push_back("A2B8C_oP22_34_c_4c_a;3;22;34;5;19;oP22;a,b/a,c/a,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6;-;TiAl2Br8;TiAl2Br8");
    vproto.push_back("AB2_oP6_34_a_c;2;6;34;4;7;oP6;a,b/a,c/a,z1,x2,y2,z2;-;FeSb2;FeSb2");
    vproto.push_back("AB8C2_oC22_35_a_ab3e_e;3;11;35;5;14;oC22;a,b/a,c/a,z1,z2,z3,y4,z4,y5,z5,y6,z6,y7,z7;-;V2MoO8;V2MoO8");
    vproto.push_back("AB_oC8_36_a_a;2;4;36;4;7;oC8;a,b/a,c/a,y1,z1,y2,z2;-;HCl;HCl");
    vproto.push_back("A2B5C2_oC36_37_d_c2d_d;3;18;37;5;16;oC36;a,b/a,c/a,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5;-;Li2Si2O5;Li2Si2O5");
    vproto.push_back("A2B3_oC40_39_2d_2c2d;2;20;39;4;19;oC40;a,b/a,c/a,x1,z1,x2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6;-;Ta3S2;Ta3S2");
    vproto.push_back("A9BC_oC44_39_3c3d_a_c;3;22;39;5;21;oC44;a,b/a,c/a,z1,x2,z2,x3,z3,x4,z4,x5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8;-;VPCl9;VPCl9");
    vproto.push_back("AB2C_oC16_40_a_2b_b;3;8;40;5;10;oC16;a,b/a,c/a,z1,y2,z2,y3,z3,y4,z4;-;K2CdPb;K2CdPb");
    vproto.push_back("AB3_oC16_40_b_3b;2;8;40;4;11;oC16;a,b/a,c/a,y1,z1,y2,z2,y3,z3,y4,z4;-;CeTe3;CeTe3");
    vproto.push_back("A10B3_oF52_42_2abce_ab;2;13;42;4;13;oF52;a,b/a,c/a,z1,z2,z3,z4,z5,y6,z6,x7,y7,z7;-;W3O10;W3O10");
    vproto.push_back("AB_oF8_42_a_a;2;2;42;4;5;oF8;a,b/a,c/a,z1,z2;-;BN;BN");
    vproto.push_back("A2BC2_oI20_45_c_b_c;3;10;45;5;10;oI20;a,b/a,c/a,z1,x2,y2,z2,x3,y3,z3;-;MnGa2Sb2;MnGa2Sb2");
    vproto.push_back("ABC_oI36_46_ac_bc_3b;3;18;46;5;18;oI36;a,b/a,c/a,z1,y2,z2,y3,z3,y4,z4,y5,z5,x6,y6,z6,x7,y7,z7;-;TiFeSi;TiFeSi");
    vproto.push_back("A2B8CD_oP24_48_k_2m_d_b;4;24;48;6;10;oP24;a,b/a,c/a,z3,x4,y4,z4,x5,y5,z5;-;alpha-RbPr[MoO4]2;alpha-RbPr2");
    vproto.push_back("A5B2_oP14_49_dehq_ab;2;14;49;4;5;oP14;a,b/a,c/a,x6,y6;-;beta-Ta2O5;beta-Ta2O5");
    vproto.push_back("AB2C8D_oP24_49_g_q_2qr_e;4;24;49;6;12;oP24;a,b/a,c/a,x3,y3,x4,y4,x5,y5,x6,y6,z6;-;CsPr[MoO4]2;CsPr2");
    vproto.push_back("A2BC4_oP28_50_ij_ac_ijm;3;28;50;5;10;oP28;a,b/a,c/a,y3,y4,y5,y6,x7,y7,z7;-;La2NiO4;La2NiO4");
    vproto.push_back("A3BC2_oP48_50_3m_m_2m;3;48;50;5;21;oP48;a,b/a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6;-;alpha-Tl2TeO3;alpha-Tl2TeO3");
    vproto.push_back("A2B_oP24_52_2e_cd;2;24;52;4;11;oP24;a,b/a,c/a,z1,x2,x3,y3,z3,x4,y4,z4;-;GaCl2;GaCl2");
    vproto.push_back("A3B2_oP20_52_de_cd;2;20;52;4;9;oP20;a,b/a,c/a,z1,x2,x3,x4,y4,z4;-;Sr2Bi3;Sr2Bi3");
    vproto.push_back("ABC2_oP16_53_h_e_gh;3;16;53;5;9;oP16;a,b/a,c/a,x1,y2,y3,z3,y4,z4;-;TaNiTe2;TaNiTe2");
    vproto.push_back("ABC3_oP20_53_e_g_hi;3;20;53;5;10;oP20;a,b/a,c/a,x1,y2,y3,z3,x4,y4,z4;-;CuBrSe3;CuBrSe3");
    vproto.push_back("ABC3_oP20_54_e_d_cf;3;20;54;5;9;oP20;a,b/a,c/a,y1,z2,z3,x4,y4,z4;-;BiGaO3;BiGaO3");
    vproto.push_back("A2B_oP24_55_2g2h_gh;2;24;55;4;15;oP24;a,b/a,c/a,x1,y1,x2,y2,x3,y3,x4,y4,x5,y5,x6,y6;-;GeAs2;GeAs2");
    vproto.push_back("A3B5_oP16_55_ch_agh;2;16;55;4;9;oP16;a,b/a,c/a,x3,y3,x4,y4,x5,y5;-;Rh5Ge3;Rh5Ge3");
    vproto.push_back("A_oP16_55_2g2h;1;16;55;3;11;oP16;a,b/a,c/a,x1,y1,x2,y2,x3,y3,x4,y4;-;C;R-carbon");
    vproto.push_back("A2B_oP6_58_g_a;2;6;58;4;5;oP6;a,b/a,c/a,x2,y2;C50;alpha-PdCl2;alpha-PdCl2");
    vproto.push_back("ABC_oP6_59_a_b_a;3;6;59;5;6;oP6;a,b/a,c/a,z1,z2,z3;E0_{5};FeOCl;FeOCl"); //DX 20180925 - added Strukturbericht
    vproto.push_back("A2B3_oP20_60_d_cd;2;20;60;4;10;oP20;a,b/a,c/a,y1,x2,y2,z2,x3,y3,z3;-;Rh2S3;Rh2S3");
    vproto.push_back("A3B_oP32_60_3d_d;2;32;60;4;15;oP32;a,b/a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4;-;WO3;WO3");
    vproto.push_back("A7B8_oP120_60_7d_8d;2;120;60;4;48;oP120;a,b/a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15;-;beta-C7H8;beta-Toluene");
    vproto.push_back("AB_oP48_61_3c_3c;2;48;61;4;21;oP48;a,b/a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6;-;Benzene;Benzene");
    vproto.push_back("A2B3_oP20_62_2c_3c;2;20;62;4;13;oP20;a,b/a,c/a,x1,z1,x2,z2,x3,z3,x4,z4,x5,z5;D5_{10};Cr3C2;Tongbaite");
    vproto.push_back("A2B4C_oP28_62_ac_2cd_c;3;28;62;5;14;oP28;a,b/a,c/a,x2,z2,x3,z3,x4,z4,x5,z5,x6,y6,z6;S1_{2};Mg2SiO4;Forsterite");
    //DX 20180619 - added info to label in part 1 vproto.push_back("A2B_oP12_62_2c_c;2;12;62;4;9;oP12;a,b/a,c/a,x1,z1,x2,z2,x3,z3;C29;SrH2;SrH2");
    vproto.push_back("A3B_oP16_62_cd_c;2;16;62;4;10;oP16;a,b/a,c/a,x1,z1,x2,z2,x3,y3,z3;D0_{20};epsilon-NiAl3;epsilon-NiAl3");
    vproto.push_back("AB2C3_oP24_62_c_d_cd;3;24;62;5;13;oP24;a,b/a,c/a,x1,z1,x2,z2,x3,y3,z3,x4,y4,z4;E9_{e};CuFe2S3;Cubanite");
    vproto.push_back("AB3_oP16_62_c_3c;2;16;62;4;11;oP16;a,b/a,c/a,x1,z1,x2,z2,x3,z3,x4,z4;D0_{8};MoO3;Molybdite");
    vproto.push_back("AB4C_oP24_62_c_2cd_c;3;24;62;5;14;oP24;a,b/a,c/a,x1,z1,x2,z2,x3,z3,x4,z4,x5,y5,z5;H0_{2};BaSO4;Barite");
    //DX 20180619 - added info to label in part 1 vproto.push_back("AB_oP8_62_c_c;2;8;62;4;7;oP8;a,b/a,c/a,x1,z1,x2,z2;B14;FeAs;Westerveldite");
    vproto.push_back("A2BC3_oC24_63_e_c_cg;3;12;63;5;8;oC24;a,b/a,c/a,y1,y2,x3,x4,y4;-;KFe2S3;Rasvumite");
    vproto.push_back("A43B5C17_oC260_63_c8fg6h_cfg_ce3f2h;3;130;63;5;59;oC260;a,b/a,c/a,y1,y2,y3,x4,y5,z5,y6,z6,y7,z7,y8,z8,y9,z9,y10,z10,y11,z11,y12,z12,y13,z13,y14,z14,y15,z15,y16,z16,x17,y17,x18,y18,x19,y19,z19,x20,y20,z20,x21,y21,z21,x22,y22,z22,x23,y23,z23,x24,y24,z24,x25,y25,z25,x26,y26,z26;-;La43Ni17Mg5;La43Ni17Mg5");
    vproto.push_back("A6B_oC28_63_efg_c;2;14;63;4;9;oC28;a,b/a,c/a,y1,x2,y3,z3,x4,y4;D2_{h};MnAl6;MnAl6");
    vproto.push_back("AB3C_oC20_63_a_cf_c;3;10;63;5;7;oC20;a,b/a,c/a,y2,y3,y4,z4;-;MgSiO3;Post-perovskite");
    vproto.push_back("AB4C_oC24_63_a_fg_c;3;12;63;5;8;oC24;a,b/a,c/a,y2,y3,z3,x4,y4;-;MgO4S;MgSO4");
    vproto.push_back("AB4C_oC24_63_c_fg_c;3;12;63;5;9;oC24;a,b/a,c/a,y1,y2,y3,z3,x4,y4;H0_{1};CaSO4;Anhydrite");
    vproto.push_back("A2B_oC24_64_2f_f;2;12;64;4;9;oC24;a,b/a,c/a,y1,z1,y2,z2,y3,z3;-;H2S;H2S");
    vproto.push_back("A2B4C_oC28_66_l_kl_a;3;14;66;5;8;oC28;a,b/a,c/a,z2,x3,y3,x4,y4;-;SrAl2Se4;SrAl2Se4");
    vproto.push_back("A3B_oC64_66_gi2lm_2l;2;32;66;4;16;oC64;a,b/a,c/a,x1,z2,x3,y3,x4,y4,x5,y5,x6,y6,x7,y7,z7;-;H3S;H3S");
    vproto.push_back("A3B_oC64_66_kl2m_bdl;2;32;66;4;14;oC64;a,b/a,c/a,z3,x4,y4,x5,y5,x6,y6,z6,x7,y7,z7;-;beta-ThI3;beta-ThI3");
    vproto.push_back("A2BC_oC16_67_ag_b_g;3;8;67;5;5;oC16;a,b/a,c/a,z3,z4;-;Al2CuIr;Al2CuIr");
    vproto.push_back("ABC2_oC16_67_b_g_ag;3;8;67;5;5;oC16;a,b/a,c/a,z3,z4;-;HoCuP2;HoCuP2");
    vproto.push_back("AB_oC8_67_a_g;2;4;67;4;4;oC8;a,b/a,c/a,z2;-;alpha-FeSe;alpha-FeSe");
    vproto.push_back("AB_oC8_67_a_g;2;4;67;4;4;oC8;a,b/a,c/a,z2;-;alpha-PbO;alpha-PbO");
    vproto.push_back("AB4_oC20_68_a_i;2;10;68;4;6;oC20;a,b/a,c/a,x2,y2,z2;-;PdSn4;PdSn4");
    vproto.push_back("AB2_oF48_70_f_fg;2;12;70;4;6;oF48;a,b/a,c/a,y1,y2,z3;D1_{f};Mn2B;Mn2B");
    vproto.push_back("A4B3_oI14_71_gh_cg;2;7;71;4;6;oI14;a,b/a,c/a,y2,y3,y4;D7_{b};Ta3B4;Ta3B4");
    vproto.push_back("ABC_oI12_71_h_j_g;3;6;71;5;6;oI12;a,b/a,c/a,y1,y2,z3;-;NbPS;NbPS");
    vproto.push_back("ABCD3_oI48_73_d_e_e_ef;4;24;73;6;10;oI48;a,b/a,c/a,y1,z2,z3,z4,x5,y5,z5;-;KAg[CO3];KAg");
    vproto.push_back("A2B_oI12_74_h_e;2;6;74;4;6;oI12;a,b/a,c/a,z1,y2,z2;-;KHg2;KHg2");
    vproto.push_back("A4B_oI20_74_beh_e;2;10;74;4;7;oI20;a,b/a,c/a,z2,z3,y4,z4;D1_{b};Al4U;Al4U");
    vproto.push_back("AB2C12D4_tP76_75_2a2b_2d_12d_4d;4;76;75;6;60;tP76;a,c/a,z1,z2,z3,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15,x16,y16,z16,x17,y17,z17,x18,y18,z18,x19,y19,z19,x20,y20,z20,x21,y21,z21,x22,y22,z22;-;BaCr2Ru4O12;BaCr2Ru4O12");
    vproto.push_back("A2BC_tP16_76_2a_a_a;3;16;76;5;14;tP16;a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4;-;LaRhC2;LaRhC2");
    vproto.push_back("A3B7_tP40_76_3a_7a;2;40;76;4;32;tP40;a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10;-;Cs3P7;Cs3P7");
    vproto.push_back("A2B6CD7_tP64_77_2d_6d_d_ab6d;4;64;77;6;49;tP64;a,c/a,z1,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15,x16,y16,z16,x17,y17,z17;-;MgB2O(OH)6;Pinnoite");
    vproto.push_back("A2B_tP48_77_8d_4d;2;48;77;4;38;tP48;a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12;-;H2S;H2S III");
    vproto.push_back("A2B7C2_tP88_78_4a_14a_4a;3;88;78;5;68;tP88;a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15,x16,y16,z16,x17,y17,z17,x18,y18,z18,x19,y19,z19,x20,y20,z20,x21,y21,z21,x22,y22,z22;-;Sr2As2O7;Sr2As2O7");
    vproto.push_back("A2BC2_tI20_79_c_2a_c;3;10;79;5;10;tI20;a,c/a,z1,z2,x3,y3,z3,x4,y4,z4;-;TlZn2Sb2;TlZn2Sb2");
    vproto.push_back("AB2_tI48_80_2b_4b;2;24;80;4;20;tI48;a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6;-;beta-NbO2;beta-NbO2");
    vproto.push_back("AB2_tP12_81_adg_2h;2;12;81;4;9;tP12;a,c/a,z3,x4,y4,z4,x5,y5,z5;-;GeSe2;GeSe2");
    vproto.push_back("A3B_tI32_82_3g_g;2;16;82;4;14;tI32;a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4;D0_{e};Ni3P;Ni3P");
    vproto.push_back("A3B2_tP10_83_adk_j;2;10;83;4;6;tP10;a,c/a,x3,y3,x4,y4;-;Ti2Ge3;Ti2Ge3");
    vproto.push_back("A2B_tP30_85_ab2g_cg;2;30;85;4;12;tP30;a,c/a,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6;-;SrBr2;SrBr2");
    vproto.push_back("AB3_tP32_86_g_3g;2;32;86;4;14;tP32;a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4;-;Ti3P;Ti3P");
    vproto.push_back("A4B_tI20_88_f_a;2;10;88;4;5;tI20;a,c/a,x2,y2,z2;-;ThCl4;ThCl4");
    vproto.push_back("AB2_tI96_88_2f_4f;2;48;88;4;20;tI96;a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6;-;alpha-NbO2;alpha-NbO2");
    vproto.push_back("A17BC4D_tP184_89_17p_p_4p_io;4;184;89;6;70;tP184;a,c/a,z1,x2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15,x16,y16,z16,x17,y17,z17,x18,y18,z18,x19,y19,z19,x20,y20,z20,x21,y21,z21,x22,y22,z22,x23,y23,z23,x24,y24,z24;-;C17FeO4Pt;C17FeO4Pt");
    vproto.push_back("A4B2C13D_tP40_90_g_d_cef2g_c;4;40;90;6;16;tP40;a,c/a,z1,z2,z3,x4,x5,x6,y6,z6,x7,y7,z7,x8,y8,z8;-;Na4Ti2Si8O22[H2O]4;Na4Ti2Si8O224");
    vproto.push_back("AB4C17D4E_tP54_90_a_g_c4g_g_c;5;54;90;7;22;tP54;a,c/a,z2,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9;-;BaCu4[VO][PO4]4;BaCu44");
    vproto.push_back("ABC_tP24_91_d_d_d;3;24;91;5;11;tP24;a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3;-;ThBC;ThBC");
    vproto.push_back("AB32CD4E8_tP184_93_i_16p_af_2p_4p;5;184;93;7;69;tP184;a,c/a,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15,x16,y16,z16,x17,y17,z17,x18,y18,z18,x19,y19,z19,x20,y20,z20,x21,y21,z21,x22,y22,z22,x23,y23,z23,x24,y24,z24,x25,y25,z25;-;AsPh4CeS8P4Me8;AsPh4CeS8P4Me8");
    vproto.push_back("A14B3C5_tP44_94_c3g_ad_bg;3;44;94;5;16;tP44;a,c/a,z3,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8;-;Na5Fe3F14;Na5Fe3F14");
    vproto.push_back("A6B2C_tP18_94_eg_c_a;3;18;94;5;7;tP18;a,c/a,z2,x3,x4,y4,z4;-;Li2MoF6;Li2MoF6");
    vproto.push_back("ABC_tP24_95_d_d_d;3;24;95;5;11;tP24;a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3;-;ThBC;ThBC");
    vproto.push_back("A2B8CD_tI24_97_d_k_a_b;4;12;97;6;5;tI24;a,c/a,x4,y4,z4;-;NaGdCu2F8;NaGdCu2F8");
    vproto.push_back("AB8C2_tI44_97_e_2k_cd;3;22;97;5;9;tI44;a,c/a,z3,x4,y4,z4,x5,y5,z5;-;Ta2Se8I;Ta2Se8I");
    vproto.push_back("A2B_tI12_98_f_a;2;6;98;4;3;tI12;a,c/a,x2;-;CdAs2;CdAs2");
    vproto.push_back("A2B8C2D_tP26_100_c_abcd_c_a;4;26;100;6;14;tP26;a,c/a,z1,z2,z3,x4,z4,x5,z5,x6,z6,x7,y7,z7;-;Ba2TiSi2O8;Fresnoite");
    vproto.push_back("A3B11C6_tP40_100_ac_bc2d_cd;3;40;100;5;19;tP40;a,c/a,z1,z2,x3,z3,x4,z4,x5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8;-;Ce3Si6N11;Ce3Si6N11");
    vproto.push_back("A7B7C2_tP32_101_bde_ade_d;3;32;101;5;16;tP32;a,c/a,z1,z2,x3,z3,x4,z4,x5,z5,x6,y6,z6,x7,y7,z7;-;gamma-MgNiSn;gamma-MgNiSn");
    vproto.push_back("A2B3_tP20_102_2c_b2c;2;20;102;4;11;tP20;a,c/a,z1,x2,z2,x3,z3,x4,z4,x5,z5;-;Gd3Al2;Gd3Al2");
    vproto.push_back("AB4_tP10_103_a_d;2;10;103;4;6;tP10;a,c/a,z1,x2,y2,z2;-;NbTe4;NbTe4");
    vproto.push_back("A5B5C4_tP28_104_ac_ac_c;3;28;104;5;13;tP28;a,c/a,z1,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5;-;Ba5In4Bi5;Ba5In4Bi5");
    vproto.push_back("AB6C4_tP22_104_a_2ac_c;3;22;104;5;11;tP22;a,c/a,z1,z2,z3,x4,y4,z4,x5,y5,z5;-;Tl4HgI6;Tl4HgI6");
    vproto.push_back("A2BC2_tP20_105_f_ac_2e;3;20;105;5;11;tP20;a,c/a,z1,z2,x3,z3,x4,z4,x5,y5,z5;-;BaGe2As2;BaGe2As2");
    vproto.push_back("A3BC3D_tP64_106_3c_c_3c_c;4;64;106;6;26;tP64;a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8;-;NaZn[OH]3;NaZn3");
    vproto.push_back("A5B7_tI24_107_ac_abd;2;12;107;4;9;tI24;a,c/a,z1,z2,z3,x4,z4,x5,z5;-;Co5Ge7;Co5Ge7");
    vproto.push_back("AB_tI4_107_a_a;2;2;107;4;4;tI4;a,c/a,z1,z2;-;GeP;GeP");
    vproto.push_back("A3B5_tI32_108_ac_a2c;2;16;108;4;10;tI32;a,c/a,z1,z2,x3,z3,x4,z4,x5,z5;-;Sr5Si3;Sr5Si3");
    vproto.push_back("ABC_tI12_109_a_a_a;3;6;109;5;5;tI12;a,c/a,z1,z2,z3;-;LaPtSi;LaPtSi");
    vproto.push_back("AB_tI8_109_a_a;2;4;109;4;4;tI8;a,c/a,z1,z2;-;NbAs;NbAs");
    vproto.push_back("A2BC8_tI176_110_2b_b_8b;3;88;110;5;35;tI176;a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11;-;Be[BH4]2;Be2");
    vproto.push_back("A2B_tP12_111_2n_adf;2;12;111;4;6;tP12;a,c/a,x4,z4,x5,z5;-;MnF2;MnF2");
    vproto.push_back("AB_tP8_111_n_n;2;8;111;4;6;tP8;a,c/a,x1,z1,x2,z2;-;VN;VN"); //DX 20180925 - prototype name should be VN not NV
    vproto.push_back("AB4C_tP12_112_b_n_e;3;12;112;5;5;tP12;a,c/a,x3,y3,z3;-;alpha-CuAlCl4;alpha-CuAlCl4");
    vproto.push_back("A2BC7D2_tP24_113_e_a_cef_e;4;24;113;6;12;tP24;a,c/a,z2,x3,z3,x4,z4,x5,z5,x6,y6,z6;S5_{3};Ca2MgSi2O7;Akermanite");
    vproto.push_back("A3B_tP32_114_3e_e;2;32;114;4;14;tP32;a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4;-;SeO3;SeO3");
    vproto.push_back("A4B_tP10_114_e_a;2;10;114;4;5;tP10;a,c/a,x2,y2,z2;-;Pd4Se;Pd4Se");
    vproto.push_back("A2B3_tP5_115_g_ag;2;5;115;4;4;tP5;a,c/a,z2,z3;-;Rh3P2;Rh3P2");
    vproto.push_back("AB2_tP12_115_j_egi;2;12;115;4;7;tP12;a,c/a,z1,z2,x3,x4,z4;-;HgI2;HgI2");
    vproto.push_back("A2B3_tP20_116_bci_fj;2;20;116;4;7;tP20;a,c/a,x3,z4,x5,y5,z5;-;Ru2Sn3;Ru2Sn3");
    vproto.push_back("A2B3_tP20_117_i_adgh;2;20;117;4;7;tP20;a,c/a,x3,x4,x5,y5,z5;D5_{12};beta-Bi2O3;beta-Bi2O3"); //DX 20180925 - added Strukturbericht designation
    vproto.push_back("A3B_tP16_118_ei_f;2;16;118;4;7;tP16;a,c/a,z1,x2,x3,y3,z3;-;RuIn3;RuIn3");
    vproto.push_back("A5B3_tP32_118_g2i_aceh;2;32;118;4;11;tP32;a,c/a,z3,x4,z5,x6,y6,z6,x7,y7,z7;-;Ir3Ga5;Ir3Ga5");
    vproto.push_back("A3B_tI24_119_b2i_af;2;12;119;4;7;tI24;a,c/a,z3,x4,z4,x5,z5;-;RbGa3;RbGa3");
    vproto.push_back("AB_tI4_119_c_a;2;2;119;4;2;tI4;a,c/a;-;GaSb;GaSb");
    vproto.push_back("A4BC2_tI28_120_i_d_e;3;14;120;5;6;tI28;a,c/a,x2,x3,y3,z3;-;KAu4Sn2;KAu4Sn2");
    vproto.push_back("A4BC4D_tP10_123_gh_a_i_d;4;10;123;6;5;tP10;a,c/a,z3,z4,z5;-;CaRbFe4As4;CaRbFe4As4"); //DX 20180925 - prototype name should be Ca not Cs
    vproto.push_back("AB4C_tP12_124_a_m_c;3;12;124;5;4;tP12;a,c/a,x3,y3;-;Nb4CoSi;Nb4CoSi");
    vproto.push_back("AB4_tP10_124_a_m;2;10;124;4;4;tP10;a,c/a,x2,y2;-;NbTe4;NbTe4");
    vproto.push_back("A4B_tP10_125_m_a;2;10;125;4;4;tP10;a,c/a,x2,z2;-;PtPb4;PtPb4");
    vproto.push_back("ABC4_tP12_125_a_b_m;3;12;125;5;4;tP12;a,c/a,x3,z3;-;KCeSe4;KCeSe4");
    vproto.push_back("A2BC4_tP28_126_cd_e_k;3;28;126;5;6;tP28;a,c/a,z3,x4,y4,z4;-;BiAl2S4;BiAl2S4");
    vproto.push_back("A4B_tP20_127_ehj_g;2;20;127;4;7;tP20;a,c/a,z1,x2,x3,x4,y4;D1_{e};ThB4;ThB4");
    vproto.push_back("A6B2C_tP18_128_eh_d_a;3;18;128;5;5;tP18;a,c/a,z3,x4,y4;-;K2SnCl6;K2SnCl6");
    vproto.push_back("A7B2C_tP40_128_egi_h_e;3;40;128;5;10;tP40;a,c/a,z1,z2,x3,x4,y4,x5,y5,z5;E9_{a};FeCu2Al7;FeCu2Al7");
    vproto.push_back("A2BC4_tP28_130_f_c_g;3;28;130;5;7;tP28;a,c/a,z1,x2,x3,y3,z3;-;CuBi2O4;CuBi2O4");
    vproto.push_back("A5B3_tP32_130_cg_cf;2;32;130;4;8;tP32;a,c/a,z1,z2,x3,x4,y4,z4;-;Ba5Si3;Ba5Si3");
    vproto.push_back("A2B2C4D_tP18_132_e_i_o_d;4;18;132;6;5;tP18;a,c/a,x3,x4,z4;-;Rb2TiCu2Se4;Rb2TiCu2S4");
    vproto.push_back("AB6C_tP16_132_d_io_a;3;16;132;5;5;tP16;a,c/a,x3,x4,z4;-;AgUF6;AgUF6");
    vproto.push_back("AB3_tP32_133_h_i2j;2;32;133;4;6;tP32;a,c/a,x1,x2,x3,x4;-;beta-V3S;beta-V3S");
    vproto.push_back("A2B_tP24_135_gh_h;2;24;135;4;7;tP24;a,c/a,x1,x2,y2,x3,y3;C47;SeO2;Downeyite");
    vproto.push_back("A4B2C_tP28_135_gh_h_d;3;28;135;5;7;tP28;a,c/a,x2,x3,y3,x4,y4;-;ZnSb2O4;ZnSb2O4");
    vproto.push_back("A2B3_tP40_137_cdf_3g;2;40;137;4;11;tP40;a,c/a,z1,z2,x3,y4,z4,y5,z5,y6,z6;D5_{9};Zn3P2;Zn3P2");
    vproto.push_back("A2B_tP6_137_d_a;2;6;137;4;3;tP6;a,c/a,z2;-;ZrO2;ZrO2");
    vproto.push_back("A4BC4_tP18_137_g_b_g;3;18;137;5;6;tP18;a,c/a,y2,z2,y3,z3;-;CeCo4B4;CeCo4B4");
    vproto.push_back("AB2_tP6_137_a_d;2;6;137;4;3;tP6;a,c/a,z2;C13;HgI2;HgI2");
    vproto.push_back("A_tP12_138_bi;1;12;138;3;4;tP12;a,c/a,x2,z2;-;C;C");
    vproto.push_back("AB_tI8_139_e_e;2;4;139;4;4;tI8;a,c/a,z1,z2;D3_{1};Hg2Cl2;Calomel");
    vproto.push_back("A3B5_tI32_140_ah_bk;2;16;140;4;5;tI32;a,c/a,x3,x4,y4;D8_{m};W5Si3;W5Si3");
    vproto.push_back("A3B5_tI32_140_ah_cl;2;16;140;4;5;tI32;a,c/a,x3,x4,z4;D8_{l};Cr5B3;Cr5B3");
    //DX 20180619 - added info to label in part 1 vproto.push_back("A2B_tI12_141_e_a;2;6;141;4;3;tI12;a,c/a,z2;C_{c};alpha-ThSi2;alpha-ThSi2");
    vproto.push_back("A_tI16_142_f;1;8;142;3;3;tI16;a,c/a,x1;-;S;S-III");
    vproto.push_back("A4B14C3_hP21_143_bd_ac4d_d;3;21;143;5;23;hP21;a,c/a,z1,z2,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9;-;Ta3Al4O13[OH];Simpsonite");
    vproto.push_back("A4B6C_hP11_143_bd_2d_a;3;11;143;5;13;hP11;a,c/a,z1,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5;-;ScRh6P4;ScRh6P4");
    vproto.push_back("AB2_hP12_143_cd_ab2d;2;12;143;4;14;hP12;a,c/a,z1,z2,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6;-;MoS2;MoS2");
    vproto.push_back("A4B_hP15_144_4a_a;2;15;144;4;17;hP15;a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5;-;IrGe4;IrGe4");
    vproto.push_back("AB_hP6_144_a_a;2;6;144;4;8;hP6;a,c/a,x1,y1,z1,x2,y2,z2;-;ZnTe;ZnTe"); //DX 20180925 - prototype name should be ZnTe not TeZn
    vproto.push_back("A2B3C3DE7_hP48_145_2a_3a_3a_a_7a;5;48;145;7;50;hP48;a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15,x16,y16,z16;-;NaCa3[CO3]2F3[H2O];Sheldrickite");
    vproto.push_back("A3BC_hR5_146_b_a_a;3;5;146;5;7;hR5;a,c/a,x1,x2,x3,y3,z3;-;gamma-Ag3SI;gamma-Ag3SI");
    vproto.push_back("ABC3_hR10_146_2a_2a_2b;3;10;146;5;12;hR10;a,c/a,x1,x2,x3,x4,x5,y5,z5,x6,y6,z6;-;FePSe3;FePSe3");
    vproto.push_back("A2B4C_hR42_148_2f_4f_f;3;42;148;5;23;hR42;a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7;S1_{3};Be2SiO4;Phenakite");
    vproto.push_back("A2B_hR18_148_2f_f;2;18;148;4;11;hR18;a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3;-;beta-PdCl2;beta-PdCl2");
    vproto.push_back("AB3_hP24_149_acgi_3l;2;24;149;4;13;hP24;a,c/a,z3,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7;-;Ti3O;Ti3O");
    vproto.push_back("A3B_hP24_153_3c_2b;2;24;153;4;13;hP24;a,c/a,x1,x2,x3,y3,z3,x4,y4,z4,x5,y5,z5;-;CrCl3;CrCl3");
    vproto.push_back("A_hP9_154_bc;1;9;154;3;6;hP9;a,c/a,x1,x2,y2,z2;-;S;S-II");
    vproto.push_back("AB2_hP9_156_b2c_3a2bc;2;9;156;4;11;hP9;a,c/a,z1,z2,z3,z4,z5,z6,z7,z8,z9;-;CdI2;CdI2");
    vproto.push_back("AB_hP12_156_2ab3c_2ab3c;2;12;156;4;14;hP12;a,c/a,z1,z2,z3,z4,z5,z6,z7,z8,z9,z10,z11,z12;-;CuI;CuI");
    vproto.push_back("AB_hP4_156_ab_ab;2;4;156;4;6;hP4;a,c/a,z1,z2,z3,z4;-;beta-CuI;beta-CuI");//DX 20180925 - shifted Wyckoff positions, orig: vproto.push_back("AB_hP4_156_ac_ac;2;4;156;4;6;hP4;a,c/a,z1,z2,z3,z4;-;beta-CuI;beta-CuI");
    vproto.push_back("A5B6C2_hP13_157_2ac_2c_b;3;13;157;5;11;hP13;a,c/a,z1,z2,z3,x4,z4,x5,z5,x6,z6;-;Ag5Pb2O6;Ag5Pb2O6");
    vproto.push_back("A3B_hP8_158_d_a;2;8;158;4;6;hP8;a,c/a,z1,x2,y2,z2;-;beta-RuCl3;beta-RuCl3");
    vproto.push_back("A2B3_hP20_159_bc_2c;2;20;159;4;12;hP20;a,c/a,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4;-;Bi2O3;Bi2O3");
    vproto.push_back("A4B3_hP28_159_ab2c_2c;2;28;159;4;16;hP28;a,c/a,z1,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6;-;alpha-Si3N4;Nierite");
    vproto.push_back("AB4C7D_hP26_159_b_ac_a2c_b;4;26;159;6;15;hP26;a,c/a,z1,z2,z3,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7;-;YbBaCo4O7;YbBaCo4O7");
    vproto.push_back("A3B_hR4_160_b_a;2;4;160;4;5;hR4;a,c/a,x1,x2,z2;-;H3S;H3S");
    vproto.push_back("A8B5_hR26_160_a3bc_a3b;2;26;160;4;19;hR26;a,c/a,x1,x2,x3,z3,x4,z4,x5,z5,x6,z6,x7,z7,x8,z8,x9,y9,z9;D8_{10};Cr5Al8;Cr5Al8"); //DX 20180925 - prototype name should be Cr5Al8 not Al8Cr5
    vproto.push_back("ABC_hR3_160_a_a_a;3;3;160;5;5;hR3;a,c/a,x1,x2,x3;F0_{2};COS;Carbonyl Sulphide");
    vproto.push_back("AB_hR10_160_5a_5a;2;10;160;4;12;hR10;a,c/a,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10;B7;SiC;Moissanite-15R");
    vproto.push_back("A2B3_hP5_164_d_ad;2;5;164;4;4;hP5;a,c/a,z2,z3;D5_{2};La2O3;La2O3");
    vproto.push_back("AB2_hP9_164_bd_c2d;2;9;164;4;6;hP9;a,c/a,z2,z3,z4,z5;-;deltaH^II-NW2;deltaH^II-NW2");
    vproto.push_back("ABC2_hP4_164_a_b_d;3;4;164;5;3;hP4;a,c/a,z3;-;CuNiSb2;CuNiSb2");
    vproto.push_back("A3B_hP24_165_bdg_f;2;24;165;4;7;hP24;a,c/a,z2,x3,x4,y4,z4;D0_{21};Cu3P;Cu3P");
    vproto.push_back("A4B3_hR7_166_2c_ac;2;7;166;4;5;hR7;a,c/a,x2,x3,x4;D7_{1};Al4C3;Al4C3");
    vproto.push_back("ABC_hR6_166_c_c_c;3;6;166;5;5;hR6;a,c/a,x1,x2,x3;-;SmSI;SmSI");
    vproto.push_back("AB3C_hR10_167_b_e_a;3;10;167;5;3;hR10;a,c/a,x3;-;PrNiO3;PrNiO3");
    vproto.push_back("ABC2_hR24_167_e_e_2e;3;24;167;5;6;hR24;a,c/a,x1,x2,x3,x4;F5_{13};KBO2;KBO2");
    vproto.push_back("A2B13C4_hP57_168_d_c6d_2d;3;57;168;5;30;hP57;a,c/a,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10;-;K2Ta4O9F4;K2Ta4O9F4");
    vproto.push_back("AB4C_hP72_168_2d_8d_2d;3;72;168;5;38;hP72;a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12;-;Al[PO4];Al");
    vproto.push_back("A2B3_hP30_169_2a_3a;2;30;169;4;17;hP30;a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5;-;alpha-Al2S3;alpha-Al2S3");
    vproto.push_back("A2B3_hP30_170_2a_3a;2;30;170;4;17;hP30;a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5;-;Al2S3;Al2S3");
    vproto.push_back("A10B2C_hP39_171_5c_c_a;3;39;171;5;21;hP39;a,c/a,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7;-;Sr[S2O6][H2O]4;Sr4");
    vproto.push_back("A10B2C_hP39_172_5c_c_a;3;39;172;5;21;hP39;a,c/a,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7;-;Sr[S2O6][H2O]4;Sr4");
    vproto.push_back("A3B_hP8_173_c_b;2;8;173;4;6;hP8;a,c/a,z1,x2,y2,z2;-;PI3;PI3");
    vproto.push_back("A4B3_hP14_173_bc_c;2;14;173;4;9;hP14;a,c/a,z1,x2,y2,z2,x3,y3,z3;-;beta-Si3N4;beta-Si3N4");
    vproto.push_back("A12B7C2_hP21_174_2j2k_ajk_cf;3;21;174;5;14;hP21;a,c/a,x4,y4,x5,y5,x6,y6,x7,y7,x8,y8,x9,y9;-;Fe12Zr2P7;Fe12Zr2P7");
    vproto.push_back("ABC_hP12_174_cj_fk_aj;3;12;174;5;8;hP12;a,c/a,x4,y4,x5,y5,x6,y6;-;GdSI;GdSI");
    vproto.push_back("A8B7C6_hP21_175_ck_aj_k;3;21;175;5;8;hP21;a,c/a,x3,y3,x4,y4,x5,y5;-;Nb7Ru6B8;Nb7Ru6B8");
    vproto.push_back("ABC_hP36_175_jk_jk_jk;3;36;175;5;14;hP36;a,c/a,x1,y1,x2,y2,x3,y3,x4,y4,x5,y5,x6,y6;-;Mg[NH];Mg");
    vproto.push_back("A3B2_hP10_176_h_bc;2;10;176;4;4;hP10;a,c/a,x3,y3;-;Er3Ru2;Er3Ru2");//DX 20180925 - shifted Wyckoff positions, orig: vproto.push_back("A3B2_hP10_176_h_bd;2;10;176;4;4;hP10;a,c/a,x3,y3;-;Er3Ru2;Er3Ru2");
    vproto.push_back("A3B3C_hP14_176_h_h_c;3;14;176;5;6;hP14;a,c/a,x2,y2,x3,y3;-;Fe3Te3Tl;Fe3Te3Tl");//DX 20180925 - shifted Wyckoff positions, orig: vproto.push_back("A3B3C_hP14_176_h_h_d;3;14;176;5;6;hP14;a,c/a,x2,y2,x3,y3;-;Fe3Te3Tl;Fe3Te3Tl");
    vproto.push_back("A3B_hP8_176_h_c;2;8;176;4;4;hP8;a,c/a,x2,y2;-;UCl3;UCl3");//DX 20180925 - shifted Wyckoff positions, orig: vproto.push_back("A3B_hP8_176_h_d;2;8;176;4;4;hP8;a,c/a,x2,y2;-;UCl3;UCl3");
    vproto.push_back("A2B_hP36_177_j2lm_n;2;36;177;4;9;hP36;a,c/a,x1,x2,x3,x4,x5,y5,z5;-;SiO2;SiO2");
    vproto.push_back("AB3_hP24_178_b_ac;2;24;178;4;7;hP24;a,c/a,x1,x2,x3,y3,z3;-;AuF3;AuF3");
    vproto.push_back("A_hP6_178_a;1;6;178;3;3;hP6;a,c/a,x1;-;Sc;Sc-V");
    vproto.push_back("AB3_hP24_179_b_ac;2;24;179;4;7;hP24;a,c/a,x1,x2,x3,y3,z3;-;AuF3;AuF3");
    vproto.push_back("A2B_hP9_181_j_c;2;9;181;4;3;hP9;a,c/a,x2;-;beta-SiO2;beta-SiO2");
    vproto.push_back("ABC_hP3_183_a_a_a;3;3;183;5;5;hP3;a,c/a,z1,z2,z3;-;AuCN;AuCN");
    vproto.push_back("AB_hP6_183_c_ab;2;6;183;4;5;hP6;a,c/a,z1,z2,z3;-;CrFe3NiSn5;CrFe3NiSn5");
    vproto.push_back("AB4C_hP72_184_d_4d_d;3;72;184;5;20;hP72;a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6;-;Al[PO4];Al");
    vproto.push_back("A3BC_hP30_185_cd_c_ab;3;30;185;5;11;hP30;a,c/a,z1,z2,x3,z3,x4,z4,x5,y5,z5;-;KNiCl3;KNiCl3");
    vproto.push_back("A3B_hP24_185_ab2c_c;2;24;185;4;10;hP24;a,c/a,z1,z2,x3,z3,x4,z4,x5,z5;-;Cu3P;Cu3P");
    vproto.push_back("A3B_hP8_185_c_a;2;8;185;4;5;hP8;a,c/a,z1,x2,z2;-;beta-RuCl3;beta-RuCl3");
    vproto.push_back("AB3_hP24_185_c_ab2c;2;24;185;4;10;hP24;a,c/a,z1,z2,x3,z3,x4,z4,x5,z5;-;Na3As;Na3As");
    vproto.push_back("A3B7_hP20_186_c_b2c;2;20;186;4;9;hP20;a,c/a,z1,x2,z2,x3,z3,x4,z4;D10_{2};Fe3Th7;Fe3Th7");
    vproto.push_back("AB3_hP4_187_e_fh;2;4;187;4;3;hP4;a,c/a,z3;-;Re3N;Re3N");
    vproto.push_back("A3BC_hP10_188_k_c_a;3;10;188;5;4;hP10;a,c/a,x3,y3;-;LiScI3;LiScI3");//DX 20180925 - shifted Wyckoff positions, orig: vproto.push_back("A3BC_hP10_188_k_a_e;3;10;188;5;4;hP10;a,c/a,x3,y3;-;LiScI3;LiScI3");
    vproto.push_back("AB9C4_hP28_188_e_kl_ak;3;28;188;5;9;hP28;a,c/a,x3,y3,x4,y4,x5,y5,z5;S3_{2};BaSi4O9;BaSi4O9"); //DX 20180925 - added Strukturbericht
    vproto.push_back("A8BC3D6_hP18_189_bfh_a_g_i;4;18;189;6;7;hP18;a,c/a,x3,x4,z5,x6,z6;E9_{b};pi-FeMg3Al8Si6;pi-FeMg3Al8Si6");
    vproto.push_back("A9BC3D5_hP18_189_fi_a_g_bh;4;18;189;6;7;hP18;a,c/a,x3,x4,z5,x6,z6;-;pi-FeMg3Al9Si5;pi-FeMg3Al9Si5");
    vproto.push_back("A2B_hP18_190_gh_bf;2;18;190;4;6;hP18;a,c/a,z2,x3,x4,y4;-;Li2Sb;Li2Sb");
    vproto.push_back("A5B3_hP16_190_bdh_g;2;16;190;4;5;hP16;a,c/a,x3,x4,y4;-;alpha-Sm3Ge5;alpha-Sm3Ge5");
    vproto.push_back("AB_hP24_190_i_afh;2;24;190;4;8;hP24;a,c/a,z2,x3,y3,x4,y4,z4;-;FeS;Troilite");
    vproto.push_back("A2B3C18D6_hP58_192_c_f_lm_l;4;58;192;6;9;hP58;a,c/a,x3,y3,x4,y4,x5,y5,z5;G3_{1};Be3Al2Si6O18;Beryl");
    vproto.push_back("AB2_hP72_192_m_j2kl;2;72;192;4;10;hP72;a,c/a,x1,x2,x3,x4,y4,x5,y5,z5;-;AlPO4;AlPO4");
    vproto.push_back("A5B3_hP16_193_dg_g;2;16;193;4;4;hP16;a,c/a,x2,x3;D8_{8};Mn5Si3;Mavlyanovite"); //DX 20180925 - added Strukturbericht
    vproto.push_back("A3B_hP16_194_gh_ac;2;16;194;4;3;hP16;a,c/a,x4;D0_{24};Ni3Ti;Ni3Ti");
    vproto.push_back("A5B2_hP28_194_ahk_ch;2;28;194;4;6;hP28;a,c/a,x3,x4,x5,z5;D8_{11};Co2Al5;Co2Al5");
    vproto.push_back("A9B3C_hP26_194_hk_h_a;3;26;194;5;6;hP26;a,c/a,x2,x3,x4,z4;E9_{c};Al9Mn3Si;Al9Mn3Si");
    vproto.push_back("A12BC4_cP34_195_2j_ab_2e;3;34;195;5;9;cP34;a,x3,x4,x5,y5,z5,x6,y6,z6;-;PrRu4P12;PrRu4P12");
    vproto.push_back("A12B2C_cF60_196_h_bc_a;3;15;196;5;4;cF60;a,x4,y4,z4;-;Cu2Fe[CN]6;Cu2Fe6");
    //DX 20180925 - wrong space group, should be 210: vproto.push_back("A12B36CD12_cF488_196_2h_6h_ac_fgh;4;122;196;6;30;cF488;a,x3,x4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13;-;MgB12H12[H2O]12;MgB12H1212");
    vproto.push_back("ABC3_cP20_198_a_a_b;3;20;198;5;6;cP20;a,x1,x2,x3,y3,z3;G3;NaClO3;Sodium Chlorate");
    vproto.push_back("A2B11_cP39_200_f_aghij;2;39;200;4;7;cP39;a,x2,x3,x4,x5,y6,z6;D8_{c};Mg2Zn11;Mg2Zn11");
    vproto.push_back("AB3C_cP60_201_be_fh_g;3;60;201;5;7;cP60;a,x2,x3,x4,x5,y5,z5;-;KSbO3;KSbO3");//DX 20180925 - shifted Wyckoff positions, orig: vproto.push_back("AB3C_cP60_201_ce_fh_g;3;60;201;5;7;cP60;a,x2,x3,x4,x5,y5,z5;-;KSbO3;KSbO3");
    vproto.push_back("A6B6C_cF104_202_h_h_c;3;26;202;5;5;cF104;a,y2,z2,y3,z3;-;KB6H6;KB6H6");
    vproto.push_back("A_cF240_202_h2i;1;60;202;3;9;cF240;a,y1,z1,x2,y2,z2,x3,y3,z3;-;C;FCC C60 Buckminsterfullerine");
    vproto.push_back("A2BCD3E6_cF208_203_e_c_d_f_g;5;52;203;7;6;cF208;a,x3,x4,x5,y5,z5;-;Na3Co(CO3)2Cl;Pyrochlore");
    vproto.push_back("A4B2C6D16E_cF232_203_e_d_f_eg_a;5;58;203;7;7;cF232;a,x3,x4,x5,x6,y6,z6;-;Na6Mg2(SO4)(CO3)4;Tychite");
    vproto.push_back("AB3C16_cF160_203_a_bc_eg;3;40;203;5;5;cF160;a,x4,x5,y5,z5;-;Rb3AsSe16;Rb3AsSe16");//DX 20180925 - shifted Wyckoff positions, orig: vproto.push_back("AB3C16_cF160_203_b_ad_eg;3;40;203;5;5;cF160;a,x4,x5,y5,z5;-;Rb3AsSe16;Rb3AsSe16");
    vproto.push_back("A2B3C6_cP264_205_2d_ab2c2d_6d;3;264;205;5;33;cP264;a,x3,x4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14;-;Ca3Al2O6;Ca3Al2O6");
    vproto.push_back("A_cP240_205_10d;1;240;205;3;31;cP240;a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10;-;C;Simple Cubic C60 Buckminsterfullerine");
    vproto.push_back("AB3C2_cI96_206_c_e_ad;3;48;206;5;6;cI96;a,x2,x3,x4,y4,z4;E9_{d};AlLi3N2;AlLi3N2");
    vproto.push_back("A17B15_cP64_207_acfk_eij;2;64;207;4;8;cP64;a,x3,x4,y5,y6,x7,y7,z7;-;Pd17Se15;Pd17Se15");
    vproto.push_back("A3B_cP16_208_j_b;2;16;208;4;2;cP16;a,x2;-;PH3;PH3");
    vproto.push_back("A6B2CD6E_cP64_208_m_ad_b_m_c;5;64;208;7;7;cP64;a,x5,y5,z5,x6,y6,z6;-;Cs2ZnFe[CN]6;Cs2ZnFe6");
    vproto.push_back("A24BC_cF104_209_j_a_b;3;26;209;5;4;cF104;a,x3,y3,z3;-;F6KP;F6KP");
    vproto.push_back("A12B36CD12_cF488_210_h_3h_a_fg;4;122;210;6;15;cF488;a,x2,y3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7;-;MgB12H12[H2O]12;MgB12H1212"); //DX 20180925 - moved this structure from SG196 to SG210
    vproto.push_back("A12B6C_cF608_210_4h_2h_e;3;152;210;5;20;cF608;a,x1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7;-;Te[OH]6;Te6");
    vproto.push_back("A2B_cI72_211_hi_i;2;36;211;4;4;cI72;a,y1,y2,y3;-;SiO2;SiO2");
    vproto.push_back("A2B_cP12_212_c_a;2;12;212;4;2;cP12;a,x2;-;SrSi2;SrSi2");
    vproto.push_back("A3B3C_cI56_214_g_h_a;3;28;214;5;3;cI56;a,y2,y3;-;Ca3PI3;Ca3PI3");
    vproto.push_back("A3BC2_cI48_214_f_a_e;3;24;214;5;3;cI48;a,x2,x3;-;Ag3AuTe2;Petzite");
    vproto.push_back("A4B9_cP52_215_ei_3efgi;2;52;215;4;11;cP52;a,x1,x2,x3,x4,x5,x6,x7,z7,x8,z8;D8_{3};gamma-Cu9Al4;gamma-brass");
    vproto.push_back("ABCD_cF16_216_c_d_b_a;4;4;216;6;1;cF16;a;-;LiMgAuSn;Quaternary Heusler"); //DX 20180925 - fixed typo in "Quaternary"
    vproto.push_back("A3B4C_cP16_218_c_e_a;3;16;218;5;2;cP16;a,x3;H2_{1};Ag3[PO4];Ag3"); //DX 20180925 - added Strukturbericht
    vproto.push_back("A7BC3D13_cF192_219_de_b_c_ah;4;48;219;6;5;cF192;a,x5,x6,y6,z6;-;Mg3B7ClO13;Boracite");
    vproto.push_back("A15B4_cI76_220_ae_c;2;38;220;4;5;cI76;a,x2,x3,y3,z3;D8_{6};Cu15Si4;Cu15Si4");
    vproto.push_back("A4B3_cI28_220_c_a;2;14;220;4;2;cI28;a,x2;D7_{3};Th3P4;Th3P4");
    vproto.push_back("A2B3C6_cP33_221_cd_ag_fh;3;33;221;5;4;cP33;a,x4,x5,x6;E9_{1};Ca3Al2O6;Ca3Al2O6");
    vproto.push_back("A5B3C16_cP96_222_ce_d_fi;3;96;222;5;6;cP96;a,x3,x4,x5,y5,z5;-;Ce5Mo3O16;Ce5Mo3O16");
    vproto.push_back("A23B6_cF116_225_bd2f_e;2;29;225;4;4;cF116;a,x3,x4,x5;D8_{a};Th6Mn23;Th6Mn23");
    vproto.push_back("A6B2C_cF36_225_e_c_a;3;9;225;5;2;cF36;a,x3;J1_{1};K2PtCl6;K2PtCl6");
    vproto.push_back("AB13_cF112_226_a_bi;2;28;226;4;3;cF112;a,y3,z3;D2_{3};NaZn13;NaZn13");
    vproto.push_back("A2B2C7_cF88_227_c_d_af;3;22;227;5;2;cF88;a,x4;-;Eu2Ir2O7;Pyrochlore Iridate");
    vproto.push_back("A3B4_cF56_227_ad_e;2;14;227;4;2;cF56;a,x3;D7_{2};Co3O4;Spinel");
    vproto.push_back("A5BCD6_cF416_228_eg_c_b_h;4;104;228;6;6;cF416;a,x3,y4,x5,y5,z5;-;CuCrCl5[NH3]6;CuCrCl56");
    vproto.push_back("A6B_cF224_228_h_c;2;56;228;4;4;cF224;a,x2,y2,z2;-;TeO6H6;TeO6H6");
    vproto.push_back("A3B10_cI52_229_e_fh;2;26;229;4;4;cI52;a,x1,x2,y3;D8_{1};gamma-Fe3Zn10;gamma-brass");
    vproto.push_back("A4B_cI10_229_c_a;2;5;229;4;1;cI10;a;-;beta-Hg4Pt;beta-Hg4Pt");
    vproto.push_back("A7B3_cI40_229_df_e;2;20;229;4;3;cI40;a,x2,x3;D8_{f};Ir3Ge7;Ir3Ge7");
    vproto.push_back("A2B3C12D3_cI160_230_a_c_h_d;4;80;230;6;4;cI160;a,x4,y4,z4;S1_{4};Co3Al2Si3O12;Garnet");

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

