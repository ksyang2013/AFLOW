// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *                                                                         *
// ***************************************************************************
// Written by Stefano Curtarolo - 2016
// Written by David Hicks - 2016

#ifndef _AFLOW_ANRL_CPP
#define _AFLOW_ANRL_CPP
#include "aflow.h"
//#include "aflow_anrl.h"

using spacegroup::SpaceGroupOptionRequired;
using std::ostream;
using std::vector;
using std::deque;
using aurostd::ExtractToStringEXPLICIT; 

// ***************************************************************************
namespace anrl { // put them in order
  // -------------------------------------------------------------------------
  // Part 1
  // -------------------------------------------------------------------------
  uint PrototypeANRL_AB2_aP12_1_4a_8a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 1
  uint PrototypeANRL_ABC2_aP16_1_4a_4a_8a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 2
  uint PrototypeANRL_A2B_aP6_2_2i_i(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 3
  uint PrototypeANRL_A_aP4_2_aci(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 4
  uint PrototypeANRL_A2B_mP12_3_bc3e_2e(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 5
  uint PrototypeANRL_A_mP4_4_2a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 6
  uint PrototypeANRL_A_mC12_5_3c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 7
  uint PrototypeANRL_A3BC_mC10_8_ab_a_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 8
  uint PrototypeANRL_A2B_mC144_9_24a_12a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 9
  uint PrototypeANRL_AB_mP4_11_e_e(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 10
  uint PrototypeANRL_ABC3_mP10_11_e_e_ef(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 11
  uint PrototypeANRL_A_mP16_11_8e(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 12
  uint PrototypeANRL_AB2_mC6_12_a_i(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 13
  uint PrototypeANRL_A_mC34_12_ah3i2j(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 14
  uint PrototypeANRL_AB3_mC16_12_g_ij(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 15
  uint PrototypeANRL_A5B2_mC14_12_a2i_i(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 16
  uint PrototypeANRL_A_mC4_12_i(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 17
  uint PrototypeANRL_ABC4_mP12_13_e_a_2g(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 18
  uint PrototypeANRL_A_mP84_13_21g(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 19
  uint PrototypeANRL_A2B_mP12_14_2e_e(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 20
  uint PrototypeANRL_A_mP32_14_8e(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 21
  uint PrototypeANRL_A_mP64_14_16e(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 22
  uint PrototypeANRL_A2B5_mC28_15_f_e2f(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 23
  uint PrototypeANRL_AB_mC8_15_c_e(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 24
  uint PrototypeANRL_A2B_mC48_15_ae3f_2f(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 25
  uint PrototypeANRL_ABC6D2_mC40_15_e_e_3f_f(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 26
  uint PrototypeANRL_ABC4_oP12_16_ag_cd_2u(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 27
  uint PrototypeANRL_AB3_oP16_18_ab_3c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 28
  uint PrototypeANRL_A2B_oP12_19_2a_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 29
  uint PrototypeANRL_A2B_oC24_20_abc_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 30
  uint PrototypeANRL_AB_oP2_25_b_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 31
  uint PrototypeANRL_AB2_oP24_28_acd_2c3d(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 32
  uint PrototypeANRL_AB3C4_oP16_31_a_ab_2ab(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 33
  uint PrototypeANRL_AB_oP8_33_a_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 34
  uint PrototypeANRL_AB3C4_oP32_33_a_3a_4a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 35
  uint PrototypeANRL_A2B_oC12_36_2a_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 36
  uint PrototypeANRL_A2BC_oC8_38_e_a_b(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 37
  uint PrototypeANRL_A2B_oC12_38_de_ab(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 38
  uint PrototypeANRL_AB4_oC20_41_a_2b(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 39
  uint PrototypeANRL_AB2_oC24_41_2a_2b(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 40
  uint PrototypeANRL_AB2_oF72_43_ab_3b(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 41
  uint PrototypeANRL_AB_oI4_44_a_b(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 42
  uint PrototypeANRL_A2B3C7D_oP13_47_t_aq_eqrs_h(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 43
  uint PrototypeANRL_AB_oP4_51_e_f(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 44
  uint PrototypeANRL_A3B2_oP20_56_ce_e(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 45
  uint PrototypeANRL_ABCD_oP16_57_d_c_d_d(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 46
  uint PrototypeANRL_AB_oP8_57_d_d(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 47
  uint PrototypeANRL_AB2_oP6_58_a_g(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 48 //49 //50
  uint PrototypeANRL_AB_oP4_59_a_b(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 51
  uint PrototypeANRL_ABC_oP6_59_a_a_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 52
  uint PrototypeANRL_A3B_oP8_59_bf_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 53
  uint PrototypeANRL_AB_oP16_61_c_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 54
  uint PrototypeANRL_A2B_oP24_61_2c_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 55
  uint PrototypeANRL_A3B2_oP20_62_3c_2c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 56
  uint PrototypeANRL_AB3C_oP20_62_c_cd_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 57
  uint PrototypeANRL_A4B_oP20_62_2cd_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 58
  uint PrototypeANRL_AB2C_oP16_62_c_2c_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 59
  uint PrototypeANRL_A2B_oP12_62_2c_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 60 //61 //62
  uint PrototypeANRL_AB_oP8_62_c_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 63 // 64 // 68 // 69
  uint PrototypeANRL_AB3_oP16_62_c_cd(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 65
  uint PrototypeANRL_A3B7_oP40_62_cd_3c2d(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 66
  uint PrototypeANRL_A_oP8_62_2c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 67
  uint PrototypeANRL_AB2C_oC16_63_c_2c_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 70
  uint PrototypeANRL_A2B_oC12_63_2c_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 71
  uint PrototypeANRL_AB_oC8_63_c_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 72
  uint PrototypeANRL_A_oC4_63_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 73
  uint PrototypeANRL_A_oC8_64_f(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 74 // 76 // 77
  uint PrototypeANRL_A2B2C_oC80_64_efg_efg_df(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 75
  uint PrototypeANRL_AB_oC8_65_j_g(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 78
  uint PrototypeANRL_A3B5_oC16_65_ah_bej(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 79
  uint PrototypeANRL_AB3_oC8_65_a_bf(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 80
  uint PrototypeANRL_AB_oF8_69_a_b(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 81
  uint PrototypeANRL_A_oF8_70_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 82
  uint PrototypeANRL_A2B_oF24_70_e_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 83
  uint PrototypeANRL_A_oF128_70_4h(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 84
  uint PrototypeANRL_AB2_oI6_71_a_i(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 85
  uint PrototypeANRL_AB2_oI6_71_a_g(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 86
  uint PrototypeANRL_A2B_oI12_72_j_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 87
  uint PrototypeANRL_AB4C_tI12_82_c_g_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 88
  uint PrototypeANRL_A2BC4_tI14_82_bc_a_g(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 89
  uint PrototypeANRL_AB_tP16_84_cej_k(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 90
  uint PrototypeANRL_A4B5_tI18_87_h_ah(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 91
  uint PrototypeANRL_AB4_tI10_87_a_h(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 92
  uint PrototypeANRL_A2B_tP12_92_b_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 93
  uint PrototypeANRL_A2B_tP36_96_3b_ab(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 94
  uint PrototypeANRL_A_tP12_96_ab(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 95
  uint PrototypeANRL_A3BC_tP5_99_bc_a_b(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 96
  uint PrototypeANRL_AB3_tP8_113_a_ce(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 97
  uint PrototypeANRL_A2BC4D_tI16_121_d_a_i_b(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 98
  uint PrototypeANRL_ABC2_tI16_122_a_b_d(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 99
  uint PrototypeANRL_AB5C_tP7_123_b_ci_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 100
  uint PrototypeANRL_AB3_tP4_123_a_ce(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 101
  uint PrototypeANRL_AB_tP2_123_a_d(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 102
  uint PrototypeANRL_ABC2_tP4_123_d_a_f(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 103
  uint PrototypeANRL_A2B3_tP10_127_g_ah(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 104
  uint PrototypeANRL_ABCD_tP8_129_c_b_a_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 105
  uint PrototypeANRL_A_tP4_129_ac(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 106
  uint PrototypeANRL_ABC_tP6_129_c_a_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 107
  uint PrototypeANRL_A2B_tP6_129_ac_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 108
  uint PrototypeANRL_AB_tP4_129_a_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 109
  uint PrototypeANRL_AB_tP4_129_c_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 110
  uint PrototypeANRL_AB_tP4_131_c_e(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 111
  uint PrototypeANRL_A_tP50_134_b2m2n(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 112
  uint PrototypeANRL_A_tP30_136_bf2ij(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 113
  uint PrototypeANRL_AB_tP8_136_g_f(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 114
  uint PrototypeANRL_A2B_tP6_136_f_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 115
  uint PrototypeANRL_sigma_tP30_136_bf2ij(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 116
  uint PrototypeANRL_A_tP4_136_f(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 117
  uint PrototypeANRL_A_tP16_138_j(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 118
  uint PrototypeANRL_A3B_tI16_139_cde_e(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 119
  uint PrototypeANRL_A_tI4_139_e(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 120
  uint PrototypeANRL_AB2C4_tI14_139_a_e_ce(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 121
  uint PrototypeANRL_A12B_tI26_139_fij_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 122
  uint PrototypeANRL_A_tI2_139_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 123 //131
  uint PrototypeANRL_A_tI8_139_h(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 124
  uint PrototypeANRL_A3B_tI8_139_bd_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 125
  uint PrototypeANRL_AB2_tI6_139_a_e(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 126
  uint PrototypeANRL_A4B5_tI18_139_i_ah(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 127
  uint PrototypeANRL_A4B_tI10_139_de_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 128
  uint PrototypeANRL_A8B_tI18_139_hi_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 129
  uint PrototypeANRL_A2B_tI6_139_d_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 130
  uint PrototypeANRL_A2B_tI12_140_h_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 132
  uint PrototypeANRL_AB3_tI16_140_b_ah(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 133
  uint PrototypeANRL_AB_tI16_140_ab_h(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 134
  uint PrototypeANRL_A4BC_tI24_141_h_b_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 135
  uint PrototypeANRL_A_tI4_141_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 136
  uint PrototypeANRL_A3B4_tI28_141_ad_h(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 137
  uint PrototypeANRL_A2B_tI12_141_e_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 138
  uint PrototypeANRL_AB_tI16_141_e_e(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 139
  uint PrototypeANRL_A2B_tI24_141_2e_e(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 140
  uint PrototypeANRL_AB_tI8_141_a_b(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 141
  uint PrototypeANRL_A2B3_tI80_141_ceh_3h(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 142
  uint PrototypeANRL_ABC4_tI96_142_e_ab_2g(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 143
  uint PrototypeANRL_A2B_hP9_147_g_ad(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 144
  uint PrototypeANRL_AB_hR16_148_cf_cf(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 145
  uint PrototypeANRL_AB3_hR8_148_c_f(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 146
  uint PrototypeANRL_AB_hR26_148_b2f_a2f(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 147
  uint PrototypeANRL_AB3C_hR10_148_c_f_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 148
  uint PrototypeANRL_A2B_hP9_150_ef_bd(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 149
  uint PrototypeANRL_A3B_hP24_151_3c_2a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 150
  uint PrototypeANRL_A2B_hP9_152_c_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 151
  uint PrototypeANRL_A_hP3_152_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 152
  uint PrototypeANRL_AB_hP6_154_a_b(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 153
  uint PrototypeANRL_AB3_hR8_155_c_de(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 154
  uint PrototypeANRL_A3B2_hR5_155_e_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 155
  uint PrototypeANRL_AB_hR6_160_b_b(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 156
  uint PrototypeANRL_AB_hR6_160_3a_3a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 157
  uint PrototypeANRL_ABC3_hR10_161_a_a_b(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 158
  uint PrototypeANRL_AB2_hP9_162_ad_k(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 159
  uint PrototypeANRL_AB2CD2_hP36_163_h_i_bf_i(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 160
  uint PrototypeANRL_A3B2_hP5_164_ad_d(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 161
  uint PrototypeANRL_AB2_hP3_164_a_d(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 162
  uint PrototypeANRL_A3B_hP24_165_adg_f(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 163
  uint PrototypeANRL_AB_hR2_166_a_b(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 164
  uint PrototypeANRL_A_hR2_166_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 165 // 172 // 175
  uint PrototypeANRL_A_hR1_166_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 166 // 170
  uint PrototypeANRL_A7B6_hR13_166_ah_3c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 167
  uint PrototypeANRL_A_hR3_166_ac(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 168
  uint PrototypeANRL_A2B3_hR5_166_c_ac(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 169
  uint PrototypeANRL_A5B2_hR7_166_a2c_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 171
  uint PrototypeANRL_A_hR12_166_2h(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 173
  uint PrototypeANRL_ABC2_hR4_166_a_b_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 174
  uint PrototypeANRL_A_hR105_166_bc9h4i(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 176
  uint PrototypeANRL_A6B_hR7_166_g_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 177
  uint PrototypeANRL_ABC3_hR10_167_a_b_e(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 178 // 179
  uint PrototypeANRL_A2B3_hR10_167_c_e(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 180
  uint PrototypeANRL_A2B_hP18_180_fi_bd(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 181
  uint PrototypeANRL_AB2_hP9_180_d_j(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 182
  uint PrototypeANRL_A2B_hP9_180_j_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 183
  uint PrototypeANRL_AB3_hP8_182_c_g(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 184
  uint PrototypeANRL_A_hP4_186_ab(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 185
  uint PrototypeANRL_AB_hP8_186_ab_ab(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 186
  uint PrototypeANRL_AB_hP4_186_b_b(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 187
  uint PrototypeANRL_AB_hP12_186_a2b_a2b(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 188
  uint PrototypeANRL_A5B3C_hP18_186_2a3b_2ab_b(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 189
  uint PrototypeANRL_AB_hP4_186_b_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 190
  uint PrototypeANRL_ABC_hP3_187_a_d_f(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 191
  uint PrototypeANRL_AB_hP2_187_d_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 192
  uint PrototypeANRL_A2B_hP9_189_fg_bc(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 193
  uint PrototypeANRL_AB4C_hP6_191_a_h_b(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 194
  uint PrototypeANRL_AB5_hP6_191_a_cg(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 195
  uint PrototypeANRL_A_hP1_191_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 196
  uint PrototypeANRL_A3B_hP4_191_bc_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 197
  uint PrototypeANRL_AB2_hP3_191_a_d(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 198
  uint PrototypeANRL_A2B_hP6_191_h_e(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 199
  uint PrototypeANRL_AB_hP6_191_f_ad(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 200
  uint PrototypeANRL_AB_hP8_194_ad_f(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 201
  uint PrototypeANRL_A_hP6_194_h(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 202
  uint PrototypeANRL_AB_hP12_194_af_bf(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 203
  uint PrototypeANRL_A_hP4_194_ac(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 204
  uint PrototypeANRL_AB3_hP8_194_c_bf(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 205
  uint PrototypeANRL_AB2_hP6_194_b_f(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 206
  uint PrototypeANRL_AB_hP4_194_c_d(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 207
  uint PrototypeANRL_ABC2_hP8_194_d_a_f(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 208
  uint PrototypeANRL_A3B_hP8_194_h_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 209
  uint PrototypeANRL_A_hP4_194_bc(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 210
  uint PrototypeANRL_AB2_hP6_194_c_f(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 211
  uint PrototypeANRL_A5B2_hP14_194_abdf_f(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 212
  uint PrototypeANRL_AB2_hP12_194_f_ah(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 213
  uint PrototypeANRL_ABC_hP6_194_c_d_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 214
  uint PrototypeANRL_A_hP4_194_f(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 215
  uint PrototypeANRL_AB2_hP6_194_c_ad(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 216
  uint PrototypeANRL_AB3C4_hP16_194_c_af_ef(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 217
  uint PrototypeANRL_A_hP2_194_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 218
  uint PrototypeANRL_AB2_hP24_194_ef_fgh(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 219
  uint PrototypeANRL_AB_hP12_194_df_ce(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 220
  uint PrototypeANRL_AB_hP4_194_c_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 221
  uint PrototypeANRL_A2B_hP12_194_cg_f(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 222
  uint PrototypeANRL_A4B_cI40_197_cde_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 223
  uint PrototypeANRL_ABC_cP12_198_a_a_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 224
  uint PrototypeANRL_A3B_cP16_198_b_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 225
  uint PrototypeANRL_A_cP8_198_2a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 226
  uint PrototypeANRL_AB_cP8_198_a_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 227 // 228
  uint PrototypeANRL_AB_cI16_199_a_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 229
  uint PrototypeANRL_AB32C48_cI162_204_a_2efg_2gh(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 230
  uint PrototypeANRL_A3B_cI32_204_g_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 231
  uint PrototypeANRL_A12B_cI26_204_g_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 232
  uint PrototypeANRL_A_cP8_205_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 233
  uint PrototypeANRL_AB_cP16_205_c_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 234
  uint PrototypeANRL_AB2_cP12_205_a_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 235
  uint PrototypeANRL_AB3C6_cI80_206_a_d_e(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 236
  uint PrototypeANRL_A_cI16_206_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 237
  uint PrototypeANRL_A_cP20_213_cd(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 238
  uint PrototypeANRL_A3B4C_cP8_215_d_e_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 239
  uint PrototypeANRL_AB4_cP5_215_a_e(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 240
  uint PrototypeANRL_AB3C4_cP8_215_a_c_e(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 241
  uint PrototypeANRL_AB5_cF24_216_a_ce(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 242
  uint PrototypeANRL_ABC_cF12_216_b_c_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 243
  uint PrototypeANRL_AB_cF8_216_c_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 244
  uint PrototypeANRL_A4B_cI10_217_c_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 245
  uint PrototypeANRL_A_cI58_217_ac2g(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 246
  uint PrototypeANRL_A5B8_cI52_217_ce_cg(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 247
  uint PrototypeANRL_A_cI16_220_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 248
  uint PrototypeANRL_A3B2_cI40_220_d_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 249
  uint PrototypeANRL_AB_cP2_221_b_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 250
  uint PrototypeANRL_AB_cP6_221_c_d(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 251
  uint PrototypeANRL_AB3C_cP5_221_a_c_b(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 252
  uint PrototypeANRL_AB27CD3_cP32_221_a_dij_b_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 253
  uint PrototypeANRL_AB3_cP4_221_a_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 254
  uint PrototypeANRL_A_cP1_221_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 255
  uint PrototypeANRL_AB11_cP36_221_c_agij(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 256
  uint PrototypeANRL_AB11CD3_cP16_221_a_dg_b_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 257
  uint PrototypeANRL_A3B_cP4_221_d_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 258
  uint PrototypeANRL_A6B_cP7_221_f_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 259
  uint PrototypeANRL_A3B_cP8_223_c_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 260
  uint PrototypeANRL_A_cP46_223_dik(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 261
  uint PrototypeANRL_A2B_cP6_224_b_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 262
  uint PrototypeANRL_A7B_cF32_225_bd_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 263
  uint PrototypeANRL_AB3_cF16_225_a_bc(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 264
  uint PrototypeANRL_A9B16C7_cF128_225_acd_2f_be(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 265
  uint PrototypeANRL_A12B_cF52_225_i_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 266
  uint PrototypeANRL_AB2_cF12_225_a_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 267
  uint PrototypeANRL_A6B23_cF116_225_e_acfh(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 268
  uint PrototypeANRL_AB2C_cF16_225_a_c_b(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 269
  uint PrototypeANRL_A_cF4_225_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 270
  uint PrototypeANRL_AB18C8_cF108_225_a_eh_f(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 271
  uint PrototypeANRL_AB_cF8_225_a_b(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 272
  uint PrototypeANRL_A2B_cF24_227_c_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 273
  uint PrototypeANRL_AB2_cF96_227_e_cf(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 274
  uint PrototypeANRL_AB_cF16_227_a_b(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 275
  uint PrototypeANRL_A_cF136_227_aeg(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 276
  uint PrototypeANRL_A2B_cF24_227_d_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 277
  uint PrototypeANRL_A_cF8_227_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 278
  uint PrototypeANRL_A2BC4_cF56_227_d_a_e(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 279
  uint PrototypeANRL_AB2_cF48_227_c_e(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 280
  uint PrototypeANRL_AB3C3_cF112_227_c_de_f(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 281
  uint PrototypeANRL_A_cI2_229_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 282
  uint PrototypeANRL_A3B_cI8_229_b_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 283
  uint PrototypeANRL_A4B3_cI14_229_c_b(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 284
  uint PrototypeANRL_A2B7_cI54_229_e_afh(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 285
  uint PrototypeANRL_AB12C3_cI32_229_a_h_b(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 286
  uint PrototypeANRL_AB4C3_cI16_229_a_c_b(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 287
  uint PrototypeANRL_A4B3_cI112_230_af_g(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 288
  // -------------------------------------------------------------------------
  // Part 2
  // -------------------------------------------------------------------------
  uint PrototypeANRL_A2B_aP6_2_aei_i(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 1
  uint PrototypeANRL_A8B5_mP13_6_a7b_3a2b(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 2
  uint PrototypeANRL_AB_mP4_6_2b_2a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 3
  uint PrototypeANRL_A2B_mP12_7_4a_2a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 4
  uint PrototypeANRL_A2B_mP18_7_6a_3a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 5
  uint PrototypeANRL_A3B_mP16_7_6a_2a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 6
  uint PrototypeANRL_A9B2_mP22_7_9a_2a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 7
  uint PrototypeANRL_A5B3_mC32_9_5a_3a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 8
  uint PrototypeANRL_AB3_mC16_9_a_3a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 9
  uint PrototypeANRL_A2B_mP6_10_mn_bg(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 10
  uint PrototypeANRL_AB3_mP16_10_mn_3m3n(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 11
  uint PrototypeANRL_ABC2_mP8_10_ac_eh_mn(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 12
  uint PrototypeANRL_AB_mP6_10_en_am(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 13
  uint PrototypeANRL_A_mP8_10_2m2n(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 14
  uint PrototypeANRL_A7B2C2_mC22_12_aij_h_i(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 15
  uint PrototypeANRL_A_mC16_12_4i(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 16
  uint PrototypeANRL_A2B_mP12_13_2g_ef(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 17
  uint PrototypeANRL_A2B_mP6_14_e_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 18
  uint PrototypeANRL_A7B8_mP120_14_14e_16e(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 19
  uint PrototypeANRL_AB3_mC16_15_e_cf(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 20
  uint PrototypeANRL_A_mC24_15_2e2f(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 21
  uint PrototypeANRL_A2B_oP12_17_abe_e(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 22
  uint PrototypeANRL_AB3_oP16_19_a_3a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 23
  uint PrototypeANRL_AB2_oC6_21_a_k(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 24
  uint PrototypeANRL_A2BC2_oF40_22_fi_ad_gh(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 25
  uint PrototypeANRL_AB_oF8_22_a_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 26
  uint PrototypeANRL_A3B_oI32_23_ij2k_k(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 27
  uint PrototypeANRL_A8B2C12D2E_oI50_23_bcfk_i_3k_j_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 28
  uint PrototypeANRL_ABC2_oI16_23_ab_i_k(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 29
  uint PrototypeANRL_ABC4_oI12_23_a_b_k(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 30
  uint PrototypeANRL_AB7CD2_oI44_24_a_b3d_c_ac(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 31
  uint PrototypeANRL_A2B_oP12_26_abc_ab(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 32 // 33
  uint PrototypeANRL_A5B_oP24_26_3a3b2c_ab(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 34
  uint PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 35
  uint PrototypeANRL_A2B_oP12_29_2a_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 36
  uint PrototypeANRL_AB2_oP12_29_a_2a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 37
  uint PrototypeANRL_ABC_oP12_29_a_a_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 38
  uint PrototypeANRL_A5B3C15_oP46_30_a2c_bc_a7c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 39
  uint PrototypeANRL_ABC3_oP20_30_2a_c_3c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 40
  uint PrototypeANRL_A13B2C2_oP34_32_a6c_c_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 41
  uint PrototypeANRL_A2B3_oP40_33_4a_6a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 42
  uint PrototypeANRL_A2B8C_oP22_34_c_4c_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 43
  uint PrototypeANRL_AB2_oP6_34_a_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 44
  uint PrototypeANRL_AB8C2_oC22_35_a_ab3e_e(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 45
  uint PrototypeANRL_AB_oC8_36_a_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 46
  uint PrototypeANRL_A2B5C2_oC36_37_d_c2d_d(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 47
  uint PrototypeANRL_A2B3_oC40_39_2d_2c2d(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 48
  uint PrototypeANRL_A9BC_oC44_39_3c3d_a_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 49
  uint PrototypeANRL_AB2C_oC16_40_a_2b_b(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 50
  uint PrototypeANRL_AB3_oC16_40_b_3b(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 51
  uint PrototypeANRL_A10B3_oF52_42_2abce_ab(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 52
  uint PrototypeANRL_AB_oF8_42_a_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 53
  uint PrototypeANRL_A2BC2_oI20_45_c_b_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 54
  uint PrototypeANRL_ABC_oI36_46_ac_bc_3b(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 55
  uint PrototypeANRL_A2B8CD_oP24_48_k_2m_d_b(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 56
  uint PrototypeANRL_A5B2_oP14_49_dehq_ab(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 57
  uint PrototypeANRL_AB2C8D_oP24_49_g_q_2qr_e(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 58
  uint PrototypeANRL_A2BC4_oP28_50_ij_ac_ijm(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 59
  uint PrototypeANRL_A3BC2_oP48_50_3m_m_2m(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 60
  uint PrototypeANRL_A2B_oP24_52_2e_cd(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 61
  uint PrototypeANRL_A3B2_oP20_52_de_cd(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 62
  uint PrototypeANRL_ABC2_oP16_53_h_e_gh(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 63
  uint PrototypeANRL_ABC3_oP20_53_e_g_hi(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 64
  uint PrototypeANRL_ABC3_oP20_54_e_d_cf(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 65
  uint PrototypeANRL_A2B_oP24_55_2g2h_gh(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 66
  uint PrototypeANRL_A3B5_oP16_55_ch_agh(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 67
  uint PrototypeANRL_A_oP16_55_2g2h(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 68
  uint PrototypeANRL_A2B_oP6_58_g_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 69
  uint PrototypeANRL_ABC_oP6_59_a_b_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 70
  uint PrototypeANRL_A2B3_oP20_60_d_cd(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 71
  uint PrototypeANRL_A3B_oP32_60_3d_d(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 72
  uint PrototypeANRL_A7B8_oP120_60_7d_8d(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 73
  uint PrototypeANRL_AB_oP48_61_3c_3c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 74
  uint PrototypeANRL_A2B3_oP20_62_2c_3c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 75
  uint PrototypeANRL_A2B4C_oP28_62_ac_2cd_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 76
  uint PrototypeANRL_A2B_oP12_62_2c_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 77
  uint PrototypeANRL_A3B_oP16_62_cd_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 78
  uint PrototypeANRL_AB2C3_oP24_62_c_d_cd(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 79
  uint PrototypeANRL_AB3_oP16_62_c_3c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 80
  uint PrototypeANRL_AB4C_oP24_62_c_2cd_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 81
  uint PrototypeANRL_AB_oP8_62_c_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 82
  uint PrototypeANRL_A2BC3_oC24_63_e_c_cg(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 83
  uint PrototypeANRL_A43B5C17_oC260_63_c8fg6h_cfg_ce3f2h(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 84
  uint PrototypeANRL_A6B_oC28_63_efg_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 85
  uint PrototypeANRL_AB3C_oC20_63_a_cf_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 86
  uint PrototypeANRL_AB4C_oC24_63_a_fg_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 87
  uint PrototypeANRL_AB4C_oC24_63_c_fg_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 88
  uint PrototypeANRL_A2B_oC24_64_2f_f(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 89
  uint PrototypeANRL_A2B4C_oC28_66_l_kl_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 90
  uint PrototypeANRL_A3B_oC64_66_gi2lm_2l(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 91
  uint PrototypeANRL_A3B_oC64_66_kl2m_bdl(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 92
  uint PrototypeANRL_A2BC_oC16_67_ag_b_g(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 93
  uint PrototypeANRL_ABC2_oC16_67_b_g_ag(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 94
  uint PrototypeANRL_AB_oC8_67_a_g(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 95 //96
  uint PrototypeANRL_AB4_oC20_68_a_i(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 97
  uint PrototypeANRL_AB2_oF48_70_f_fg(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 98
  uint PrototypeANRL_A4B3_oI14_71_gh_cg(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 99
  uint PrototypeANRL_ABC_oI12_71_h_j_g(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 100
  uint PrototypeANRL_ABCD3_oI48_73_d_e_e_ef(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 101
  uint PrototypeANRL_A2B_oI12_74_h_e(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 102
  uint PrototypeANRL_A4B_oI20_74_beh_e(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 103
  uint PrototypeANRL_AB2C12D4_tP76_75_2a2b_2d_12d_4d(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 104
  uint PrototypeANRL_A2BC_tP16_76_2a_a_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 105
  uint PrototypeANRL_A3B7_tP40_76_3a_7a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 106
  uint PrototypeANRL_A2B6CD7_tP64_77_2d_6d_d_ab6d(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 107
  uint PrototypeANRL_A2B_tP48_77_8d_4d(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 108
  uint PrototypeANRL_A2B7C2_tP88_78_4a_14a_4a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 109
  uint PrototypeANRL_A2BC2_tI20_79_c_2a_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 110
  uint PrototypeANRL_AB2_tI48_80_2b_4b(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 111
  uint PrototypeANRL_AB2_tP12_81_adg_2h(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 112
  uint PrototypeANRL_A3B_tI32_82_3g_g(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 113
  uint PrototypeANRL_A3B2_tP10_83_adk_j(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 114
  uint PrototypeANRL_A2B_tP30_85_ab2g_cg(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 115
  uint PrototypeANRL_AB3_tP32_86_g_3g(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 116
  uint PrototypeANRL_A4B_tI20_88_f_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 117
  uint PrototypeANRL_AB2_tI96_88_2f_4f(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 118
  uint PrototypeANRL_A17BC4D_tP184_89_17p_p_4p_io(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 119
  uint PrototypeANRL_A4B2C13D_tP40_90_g_d_cef2g_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 120
  uint PrototypeANRL_AB4C17D4E_tP54_90_a_g_c4g_g_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 121
  uint PrototypeANRL_ABC_tP24_91_d_d_d(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 122
  uint PrototypeANRL_AB32CD4E8_tP184_93_i_16p_af_2p_4p(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 123
  uint PrototypeANRL_A14B3C5_tP44_94_c3g_ad_bg(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 124
  uint PrototypeANRL_A6B2C_tP18_94_eg_c_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 125
  uint PrototypeANRL_ABC_tP24_95_d_d_d(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 126
  uint PrototypeANRL_A2B8CD_tI24_97_d_k_a_b(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 127
  uint PrototypeANRL_AB8C2_tI44_97_e_2k_cd(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 128
  uint PrototypeANRL_A2B_tI12_98_f_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 129
  uint PrototypeANRL_A2B8C2D_tP26_100_c_abcd_c_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 130
  uint PrototypeANRL_A3B11C6_tP40_100_ac_bc2d_cd(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 131
  uint PrototypeANRL_A7B7C2_tP32_101_bde_ade_d(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 132
  uint PrototypeANRL_A2B3_tP20_102_2c_b2c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 133
  uint PrototypeANRL_AB4_tP10_103_a_d(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 134
  uint PrototypeANRL_A5B5C4_tP28_104_ac_ac_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 135
  uint PrototypeANRL_AB6C4_tP22_104_a_2ac_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 136
  uint PrototypeANRL_A2BC2_tP20_105_f_ac_2e(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 137
  uint PrototypeANRL_A3BC3D_tP64_106_3c_c_3c_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 138
  uint PrototypeANRL_A5B7_tI24_107_ac_abd(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 139
  uint PrototypeANRL_AB_tI4_107_a_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 140
  uint PrototypeANRL_A3B5_tI32_108_ac_a2c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 141
  uint PrototypeANRL_ABC_tI12_109_a_a_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 142
  uint PrototypeANRL_AB_tI8_109_a_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 143
  uint PrototypeANRL_A2BC8_tI176_110_2b_b_8b(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 144
  uint PrototypeANRL_A2B_tP12_111_2n_adf(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 145
  uint PrototypeANRL_AB_tP8_111_n_n(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 146
  uint PrototypeANRL_AB4C_tP12_112_b_n_e(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 147
  uint PrototypeANRL_A2BC7D2_tP24_113_e_a_cef_e(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 148
  uint PrototypeANRL_A3B_tP32_114_3e_e(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 149
  uint PrototypeANRL_A4B_tP10_114_e_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 150
  uint PrototypeANRL_A2B3_tP5_115_g_ag(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 151
  uint PrototypeANRL_AB2_tP12_115_j_egi(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 152
  uint PrototypeANRL_A2B3_tP20_116_bci_fj(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 153
  uint PrototypeANRL_A2B3_tP20_117_i_adgh(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 154
  uint PrototypeANRL_A3B_tP16_118_ei_f(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 155
  uint PrototypeANRL_A5B3_tP32_118_g2i_aceh(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 156
  uint PrototypeANRL_A3B_tI24_119_b2i_af(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 157
  uint PrototypeANRL_AB_tI4_119_c_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 158
  uint PrototypeANRL_A4BC2_tI28_120_i_d_e(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 159
  uint PrototypeANRL_A4BC4D_tP10_123_gh_a_i_d(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 160
  uint PrototypeANRL_AB4C_tP12_124_a_m_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 161
  uint PrototypeANRL_AB4_tP10_124_a_m(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 162
  uint PrototypeANRL_A4B_tP10_125_m_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 163
  uint PrototypeANRL_ABC4_tP12_125_a_b_m(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 164
  uint PrototypeANRL_A2BC4_tP28_126_cd_e_k(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 165
  uint PrototypeANRL_A4B_tP20_127_ehj_g(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 166
  uint PrototypeANRL_A6B2C_tP18_128_eh_d_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 167
  uint PrototypeANRL_A7B2C_tP40_128_egi_h_e(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 168
  uint PrototypeANRL_A2BC4_tP28_130_f_c_g(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 169
  uint PrototypeANRL_A5B3_tP32_130_cg_cf(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 170
  uint PrototypeANRL_A2B2C4D_tP18_132_e_i_o_d(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 171
  uint PrototypeANRL_AB6C_tP16_132_d_io_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 172
  uint PrototypeANRL_AB3_tP32_133_h_i2j(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 173
  uint PrototypeANRL_A2B_tP24_135_gh_h(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 174
  uint PrototypeANRL_A4B2C_tP28_135_gh_h_d(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 175
  uint PrototypeANRL_A2B3_tP40_137_cdf_3g(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 176
  uint PrototypeANRL_A2B_tP6_137_d_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 177
  uint PrototypeANRL_A4BC4_tP18_137_g_b_g(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 178
  uint PrototypeANRL_AB2_tP6_137_a_d(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 179
  uint PrototypeANRL_A_tP12_138_bi(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 180
  uint PrototypeANRL_AB_tI8_139_e_e(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 181
  uint PrototypeANRL_A3B5_tI32_140_ah_bk(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 182
  uint PrototypeANRL_A3B5_tI32_140_ah_cl(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 183
  uint PrototypeANRL_A2B_tI12_141_e_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 184
  uint PrototypeANRL_A_tI16_142_f(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 185
  uint PrototypeANRL_A4B14C3_hP21_143_bd_ac4d_d(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 186
  uint PrototypeANRL_A4B6C_hP11_143_bd_2d_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 187
  uint PrototypeANRL_AB2_hP12_143_cd_ab2d(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 188
  uint PrototypeANRL_A4B_hP15_144_4a_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 189
  uint PrototypeANRL_AB_hP6_144_a_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 190
  uint PrototypeANRL_A2B3C3DE7_hP48_145_2a_3a_3a_a_7a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 191
  uint PrototypeANRL_A3BC_hR5_146_b_a_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 192
  uint PrototypeANRL_ABC3_hR10_146_2a_2a_2b(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 193
  uint PrototypeANRL_A2B4C_hR42_148_2f_4f_f(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 194
  uint PrototypeANRL_A2B_hR18_148_2f_f(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 195
  uint PrototypeANRL_AB3_hP24_149_acgi_3l(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 196
  uint PrototypeANRL_A3B_hP24_153_3c_2b(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 197
  uint PrototypeANRL_A_hP9_154_bc(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 198
  uint PrototypeANRL_AB2_hP9_156_b2c_3a2bc(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 199
  uint PrototypeANRL_AB_hP12_156_2ab3c_2ab3c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 200
  uint PrototypeANRL_AB_hP4_156_ab_ab(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 201 //DX 20180925 - changed label
  uint PrototypeANRL_A5B6C2_hP13_157_2ac_2c_b(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 202
  uint PrototypeANRL_A3B_hP8_158_d_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 203
  uint PrototypeANRL_A2B3_hP20_159_bc_2c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 204
  uint PrototypeANRL_A4B3_hP28_159_ab2c_2c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 205
  uint PrototypeANRL_AB4C7D_hP26_159_b_ac_a2c_b(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 206
  uint PrototypeANRL_A3B_hR4_160_b_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 207
  uint PrototypeANRL_A8B5_hR26_160_a3bc_a3b(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 208
  uint PrototypeANRL_ABC_hR3_160_a_a_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 209
  uint PrototypeANRL_AB_hR10_160_5a_5a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 210
  uint PrototypeANRL_A2B3_hP5_164_d_ad(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 211
  uint PrototypeANRL_AB2_hP9_164_bd_c2d(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 212
  uint PrototypeANRL_ABC2_hP4_164_a_b_d(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 213
  uint PrototypeANRL_A3B_hP24_165_bdg_f(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 214
  uint PrototypeANRL_A4B3_hR7_166_2c_ac(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 215
  uint PrototypeANRL_ABC_hR6_166_c_c_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 216
  uint PrototypeANRL_AB3C_hR10_167_b_e_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 217
  uint PrototypeANRL_ABC2_hR24_167_e_e_2e(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 218
  uint PrototypeANRL_A2B13C4_hP57_168_d_c6d_2d(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 219
  uint PrototypeANRL_AB4C_hP72_168_2d_8d_2d(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 220
  uint PrototypeANRL_A2B3_hP30_169_2a_3a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 221
  uint PrototypeANRL_A2B3_hP30_170_2a_3a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 222
  uint PrototypeANRL_A10B2C_hP39_171_5c_c_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 223
  uint PrototypeANRL_A10B2C_hP39_172_5c_c_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 224
  uint PrototypeANRL_A3B_hP8_173_c_b(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 225
  uint PrototypeANRL_A4B3_hP14_173_bc_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 226
  uint PrototypeANRL_A12B7C2_hP21_174_2j2k_ajk_cf(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 227
  uint PrototypeANRL_ABC_hP12_174_cj_fk_aj(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 228
  uint PrototypeANRL_A8B7C6_hP21_175_ck_aj_k(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 229
  uint PrototypeANRL_ABC_hP36_175_jk_jk_jk(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 230
  uint PrototypeANRL_A3B2_hP10_176_h_bc(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 231 //DX 20180925 - changed label
  uint PrototypeANRL_A3B3C_hP14_176_h_h_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 232 //DX 20180925 - changed label
  uint PrototypeANRL_A3B_hP8_176_h_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 233 //DX 20180925 - changed label
  uint PrototypeANRL_A2B_hP36_177_j2lm_n(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 234
  uint PrototypeANRL_AB3_hP24_178_b_ac(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 235
  uint PrototypeANRL_A_hP6_178_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 236
  uint PrototypeANRL_AB3_hP24_179_b_ac(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 237
  uint PrototypeANRL_A2B_hP9_181_j_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 238
  uint PrototypeANRL_ABC_hP3_183_a_a_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 239
  uint PrototypeANRL_AB_hP6_183_c_ab(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 240
  uint PrototypeANRL_AB4C_hP72_184_d_4d_d(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 241
  uint PrototypeANRL_A3BC_hP30_185_cd_c_ab(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 242
  uint PrototypeANRL_A3B_hP24_185_ab2c_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 243
  uint PrototypeANRL_A3B_hP8_185_c_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 244
  uint PrototypeANRL_AB3_hP24_185_c_ab2c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 245
  uint PrototypeANRL_A3B7_hP20_186_c_b2c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 246
  uint PrototypeANRL_AB3_hP4_187_e_fh(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 247
  uint PrototypeANRL_A3BC_hP10_188_k_c_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 248 //DX 20180925 - changed label
  uint PrototypeANRL_AB9C4_hP28_188_e_kl_ak(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 249
  uint PrototypeANRL_A8BC3D6_hP18_189_bfh_a_g_i(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 250
  uint PrototypeANRL_A9BC3D5_hP18_189_fi_a_g_bh(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 251
  uint PrototypeANRL_A2B_hP18_190_gh_bf(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 252
  uint PrototypeANRL_A5B3_hP16_190_bdh_g(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 253
  uint PrototypeANRL_AB_hP24_190_i_afh(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 254
  uint PrototypeANRL_A2B3C18D6_hP58_192_c_f_lm_l(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 255
  uint PrototypeANRL_AB2_hP72_192_m_j2kl(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 256
  uint PrototypeANRL_A5B3_hP16_193_dg_g(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 257
  uint PrototypeANRL_A3B_hP16_194_gh_ac(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 258
  uint PrototypeANRL_A5B2_hP28_194_ahk_ch(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 259
  uint PrototypeANRL_A9B3C_hP26_194_hk_h_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 260
  uint PrototypeANRL_A12BC4_cP34_195_2j_ab_2e(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 261
  uint PrototypeANRL_A12B2C_cF60_196_h_bc_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 262
  uint PrototypeANRL_ABC3_cP20_198_a_a_b(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 263
  uint PrototypeANRL_A2B11_cP39_200_f_aghij(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 264
  uint PrototypeANRL_AB3C_cP60_201_be_fh_g(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 265 //DX 20180925 - changed label
  uint PrototypeANRL_A6B6C_cF104_202_h_h_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 266
  uint PrototypeANRL_A_cF240_202_h2i(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 267
  uint PrototypeANRL_A2BCD3E6_cF208_203_e_c_d_f_g(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 268
  uint PrototypeANRL_A4B2C6D16E_cF232_203_e_d_f_eg_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 269
  uint PrototypeANRL_AB3C16_cF160_203_a_bc_eg(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 270 //DX 20180925 - changed label
  uint PrototypeANRL_A2B3C6_cP264_205_2d_ab2c2d_6d(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 271
  uint PrototypeANRL_A_cP240_205_10d(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 272
  uint PrototypeANRL_AB3C2_cI96_206_c_e_ad(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 273
  uint PrototypeANRL_A17B15_cP64_207_acfk_eij(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 274
  uint PrototypeANRL_A3B_cP16_208_j_b(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 275
  uint PrototypeANRL_A6B2CD6E_cP64_208_m_ad_b_m_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 277
  uint PrototypeANRL_A24BC_cF104_209_j_a_b(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 277
  uint PrototypeANRL_A12B36CD12_cF488_210_h_3h_a_fg(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 278
  uint PrototypeANRL_A12B6C_cF608_210_4h_2h_e(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 279
  uint PrototypeANRL_A2B_cI72_211_hi_i(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 280
  uint PrototypeANRL_A2B_cP12_212_c_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 281
  uint PrototypeANRL_A3B3C_cI56_214_g_h_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 282
  uint PrototypeANRL_A3BC2_cI48_214_f_a_e(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 283
  uint PrototypeANRL_A4B9_cP52_215_ei_3efgi(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 284
  uint PrototypeANRL_ABCD_cF16_216_c_d_b_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 285
  uint PrototypeANRL_A3B4C_cP16_218_c_e_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 286
  uint PrototypeANRL_A7BC3D13_cF192_219_de_b_c_ah(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 287
  uint PrototypeANRL_A15B4_cI76_220_ae_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 288
  uint PrototypeANRL_A4B3_cI28_220_c_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 289
  uint PrototypeANRL_A2B3C6_cP33_221_cd_ag_fh(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 290
  uint PrototypeANRL_A5B3C16_cP96_222_ce_d_fi(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 291
  uint PrototypeANRL_A23B6_cF116_225_bd2f_e(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 292
  uint PrototypeANRL_A6B2C_cF36_225_e_c_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 293
  uint PrototypeANRL_AB13_cF112_226_a_bi(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 294
  uint PrototypeANRL_A2B2C7_cF88_227_c_d_af(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 295
  uint PrototypeANRL_A3B4_cF56_227_ad_e(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 296
  uint PrototypeANRL_A5BCD6_cF416_228_eg_c_b_h(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 297
  uint PrototypeANRL_A6B_cF224_228_h_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 298
  uint PrototypeANRL_A3B10_cI52_229_e_fh(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 299
  uint PrototypeANRL_A4B_cI10_229_c_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 300
  uint PrototypeANRL_A7B3_cI40_229_df_e(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 301
  uint PrototypeANRL_A2B3C12D3_cI160_230_a_c_h_d(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG); // 302
}

// *************************************************************************** 
// the mother of all the AFLOW-NRL prototypes
namespace anrl {
  xstructure PrototypeANRL(ostream &oss,string label,string parameters,deque<string> &vatomX,deque<double> &vvolumeX,double volume_in,int mode,bool flip_option) { // COMPLETE ONE
    //   { vector<string> tokens;aurostd::string2tokens(label,tokens,"."); if(tokens.size()>1) XHOST.DEBUG=TRUE; }
    // XHOST.DEBUG=TRUE;
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) { cerr << "anrl::PrototypeANRL(ostream &oss,string label,deque<string> &vatomX,deque<double> &vvolumeX,double volume_in,int mode,bool flip_option)" << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL: label=" << label << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL: parameters=" << parameters << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL: volume_in=" << volume_in << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL: mode=" << mode << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL: flip_option=" << flip_option << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL: vatomX.size()=" << vatomX.size() << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL: vatomX ="; for(uint i=0;i<vatomX.size();i++) {cerr << " " << vatomX.at(i);} cerr << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL: vvolumeX.size()=" << vvolumeX.size() << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL: vvolumeX ="; for(uint i=0;i<vvolumeX.size();i++) {cerr << " " << vvolumeX.at(i);} cerr << endl;}
 
    deque<string> vatomX_backup(vatomX);
    xvector<double> origin(3);origin.clear();
    xstructure str("");str.lattice.clear();
    
    stringstream web;

    // ---------------------------------------------------------------------------
    // declaration
    vector<string> vproto,vlabel; 
    vector<uint>   vproto_nspecies,vproto_natoms,vproto_spacegroup,vproto_nunderscores,vproto_nparameters; 
    vector<string> vproto_Pearson_symbol,vproto_params,vproto_Strukturbericht,vproto_prototype,vproto_dialect;
        
    // ---------------------------------------------------------------------------
    // load defaults
    anrl::PrototypeANRL_LoadList(vproto,vlabel,vproto_nspecies,vproto_natoms,
				 vproto_spacegroup,vproto_nunderscores,vproto_nparameters,vproto_Pearson_symbol,vproto_params,
				 vproto_Strukturbericht,vproto_prototype,vproto_dialect);

    // ---------------------------------------------------------------------------
    // create strings
    string label_anrl="",straus;
    string label_permutations=""; deque<uint> vpermutation;
    vector<string> tokens;

    // ---------------------------------------------------------------------------
    // search for label_permutations
    aurostd::string2tokens(label,tokens,".");
    if(LDEBUG) { cerr << "anrl::PrototypeANRL: tokens.size()=" << tokens.size() << endl;}
    
    if(tokens.size()==0) { label_anrl=label; }
    if(tokens.size()==1) { label_anrl=tokens.at(0); }
    if(tokens.size()==2) { label_anrl=tokens.at(0); label_permutations=tokens.at(1); }
    
    for(uint i=0;i<label_permutations.size();i++) vpermutation.push_back(aurostd::mod(label_permutations.at(i)-65,32));

    if(LDEBUG) { cerr << "anrl::PrototypeANRL: label=" << label << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL: label.size()=" << label.size() << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL: label_anrl=" << label_anrl << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL: label_anrl.size()=" << label_anrl.size() << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL: label_permutations=" << label_permutations << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL: vpermutation.size()=" << vpermutation.size() << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL: vpermutation ="; for(uint i=0;i<vpermutation.size();i++) {cerr << " " << vpermutation.at(i);} cerr << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL: vatomX.size()=" << vatomX.size() << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL: vatomX ="; for(uint i=0;i<vatomX.size();i++) {cerr << " " << vatomX.at(i);} cerr << endl;}
    
    // ---------------------------------------------------------------------------
    // search
    bool found=FALSE;
    uint ifound=0;
    for(uint i=0;i<vlabel.size()&&!found;i++) {
      if(vlabel.at(i)==label_anrl) {  // FIX
	found=TRUE;
	ifound=i;
      }
    }
    
    // ---------------------------------------------------------------------------
    // not found
    if(!found) {
      cerr << "ERROR - anrl::PrototypeANRL: prototype not found [label=" << label << "]" << endl;
      exit(0);
    }
      
    if(LDEBUG) { cerr << "anrl::PrototypeANRL: FOUND" << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL: ifound=" << ifound << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL: vlabel.at(ifound)=" << vlabel.at(ifound) << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL: vproto_nspecies.at(ifound)=" << vproto_nspecies.at(ifound) << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL: vproto_natoms.at(ifound)=" << vproto_natoms.at(ifound) << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL: vproto_spacegroup.at(ifound)=" << vproto_spacegroup.at(ifound) << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL: vproto.at(ifound)=" << vproto.at(ifound) << endl;}
      
    // ---------------------------------------------------------------------------
    // check for vpermutation size and errors
    if(vpermutation.size()>0 && vpermutation.size()!=vproto_nspecies.at(ifound)) {
      cerr << "ERROR - anrl::PrototypeANRL: wrong number of permutation species [label=" << label << "]" << endl;
      cerr << "        anrl::PrototypeANRL: vproto_nspecies.at(ifound)=" << vproto_nspecies.at(ifound) << endl;
      cerr << "        anrl::PrototypeANRL: vpermutation.size()=" << vpermutation.size() << endl;
      cerr << "        anrl::PrototypeANRL: vpermutation ="; for(uint i=0;i<vpermutation.size();i++) {cerr << " " << vpermutation.at(i);} cerr << endl;
      exit(0);
    }
    // ---------------------------------------------------------------------------
    // check for vatomX size and errors
    if(vatomX.size()>0 && vatomX.size()!=vproto_nspecies.at(ifound)) {
      cerr << "ERROR - anrl::PrototypeANRL: wrong number of name species [label=" << label << "]" << endl;
      cerr << "        anrl::PrototypeANRL: vproto_nspecies.at(ifound)=" << vproto_nspecies.at(ifound) << endl;
      cerr << "        anrl::PrototypeANRL: vatomX.size()=" << vatomX.size() << endl;
      cerr << "        anrl::PrototypeANRL: vatomX ="; for(uint i=0;i<vatomX.size();i++) {cerr << " " << vatomX.at(i);} cerr << endl;
      exit(0);
    }
    
    for(uint i=0;i<vproto_nspecies.at(ifound);i++) { // number of species
      str.num_each_type.push_back(0);str.comp_each_type.push_back(0.0);
      str.species.push_back("");str.species_pp.push_back("");str.species_pp_type.push_back("");str.species_pp_version.push_back("");
      str.species_pp_ZVAL.push_back(0.0);
      str.species_pp_vLDAU.push_back(deque<double>());
      str.species_volume.push_back(0.0);
      str.species_mass.push_back(0.0);
    }

    // ---------------------------------------------------------------------------
    // explicit equations - DX 20180615 
    uint print_mode = 0; // no equations
    if(XHOST.vflag_pflow.flag("PROTO::EQUATIONS_ONLY")) {
      print_mode = 1; // equations only
      str.symbolic_math_representation_only=TRUE;  //DX 20180618 print symbolic math representation only
    }
    else if(XHOST.vflag_pflow.flag("PROTO::ADD_EQUATIONS")) {
      print_mode = 2; // equations + parameters
      str.constrained_symmetry_calculation=TRUE;  //DX 20180618 appends information to geometry file for calculation
    }


    // -------------------------------------------------------------------------
    // Part 1
    // -------------------------------------------------------------------------
    // ---------------------------------------------------------------------------
    // 1 // ./aflow --proto=AB2_aP12_1_4a_8a --params=5.417,1.0,1.0,90.0,90.0,90.0,0.001,0.002,0.003,0.4966,0.0001,0.5036,0.5001,0.502,0.0011,-0.0006,0.5013,0.5038,0.3857,0.3832,0.384,0.1149,0.6114,0.8846,0.8854,0.1157,0.6143,0.6153,0.8865,0.1141,0.6151,0.6132,0.6137,0.8854,0.3818,0.1149,0.1147,0.8856,0.3841,0.3857,0.1161,0.8842
    if(vlabel.at(ifound)=="AB2_aP12_1_4a_8a") {
      PrototypeANRL_AB2_aP12_1_4a_8a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 2 // ./aflow --proto=ABC2_aP16_1_4a_4a_8a --params=6.554,1.00061031431,1.92662496186,100.43475,100.46074,107.53,0.3267,0.582,0.177,0.565,-0.0132,0.4424,0.5217,0.3883,0.6767,-0.0744,0.6254,-0.0574,0.0338,0.0476,0.2599,0.0831,0.6072,0.4974,-0.0131,0.0949,0.7583,0.5449,0.1443,-0.0022,-0.0211,0.5213,0.2073,0.2907,0.5956,-0.0183,-0.0616,0.0602,0.4998,0.5068,-0.0175,0.2448,0.4596,0.0397,0.708,0.5326,0.352,0.4818,0.0,0.0,0.0,-0.078,0.569,0.7448
    if(vlabel.at(ifound)=="ABC2_aP16_1_4a_4a_8a") {
      PrototypeANRL_ABC2_aP16_1_4a_4a_8a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 3 // ./aflow --proto=A2B_aP6_2_2i_i --params=4.56,1.54824561404,1.62280701754,80.2,106.96667,98.2,0.557,0.73,0.165,0.82,0.803,0.695,0.397,0.639,0.463
    if(vlabel.at(ifound)=="A2B_aP6_2_2i_i") {
      PrototypeANRL_A2B_aP6_2_2i_i(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 4 // ./aflow --proto=A_aP4_2_aci --params=3.307,2.24130631993,0.844572119746,89.06,85.15,85.7,0.572,0.259,0.433
    if(vlabel.at(ifound)=="A_aP4_2_aci") {
      PrototypeANRL_A_aP4_2_aci(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 5 // ./aflow --proto=A2B_mP12_3_bc3e_2e --params=4.1605,0.992524936907,1.78370388174,101.3752,0.15907,0.73859,0.02399,0.752,0.18927,0.38562,0.71473,0.64074,0.48963,0.20196,0.18802,0.18244,0.0,0.69651,0.38098,0.58564,0.17797
    if(vlabel.at(ifound)=="A2B_mP12_3_bc3e_2e") {
      PrototypeANRL_A2B_mP12_3_bc3e_2e(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 6 // ./aflow --proto=A_mP4_4_2a --params=3.104,2.42042525773,1.53350515464,92.71,0.25,0.23,0.48,0.48,0.0,0.02
    if(vlabel.at(ifound)=="A_mP4_4_2a") {
      PrototypeANRL_A_mP4_4_2a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 7 // ./aflow --proto=A_mC12_5_3c --params=7.42,0.578167115903,1.90026954178,92.0,0.05,0.27,0.245,0.63,0.3,0.4,0.245,0.43,0.07
    if(vlabel.at(ifound)=="A_mC12_5_3c") {
      PrototypeANRL_A_mC12_5_3c(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 8 // ./aflow --proto=A3BC_mC10_8_ab_a_a --params=5.72204,0.9978207073,0.722908263486,90.498,0.5515,-0.0994,0.0,0.0,0.523,0.4492,0.288,0.2434,0.3729
    if(vlabel.at(ifound)=="A3BC_mC10_8_ab_a_a") {
      PrototypeANRL_A3BC_mC10_8_ab_a_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 9 // ./aflow --proto=A2B_mC144_9_24a_12a --params=18.524,0.270092852516,1.28535953358,105.82,0.5749,0.351,0.8182,0.0707,0.34,0.8476,0.7315,0.138,0.4851,0.2509,0.144,0.5152,0.4155,0.352,0.6741,-0.0873,0.352,0.6434,0.8773,0.164,-0.0787,0.416,0.168,-0.0639,0.7741,0.145,0.7538,0.2336,0.143,0.7402,0.6195,0.341,0.5847,0.0811,0.343,0.5661,-0.0034,0.011,0.6062,0.3533,0.489,0.5665,0.6498,0.005,0.6711,0.1524,0.496,0.7805,0.8636,0.499,0.7328,0.3361,0.003,0.8333,0.0052,0.493,0.7398,0.1369,0.011,-0.0732,0.4927,0.492,0.8868,0.5,0.468,0.5,0.2252,0.491,0.5898,0.2744,0.021,-0.0845,0.0507,0.041,0.5642,0.2036,0.447,0.7347,-0.0802,0.049,0.6225,0.5751,0.043,0.7955,0.4247,0.048,0.6971,0.2643,0.444,0.5386,0.8023,0.449,0.7661,0.6453,0.041,0.6027,0.8531,0.463,-0.0984,0.4493,0.466,-0.0642,0.2244,0.059,-0.0395,0.0697,0.049,0.8702
    if(vlabel.at(ifound)=="A2B_mC144_9_24a_12a") {
      PrototypeANRL_A2B_mC144_9_24a_12a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 10 // ./aflow --proto=AB_mP4_11_e_e --params=2.8837,1.42393452856,1.61854561848,82.062,0.0387,0.8252,0.5887,0.7184
    if(vlabel.at(ifound)=="AB_mP4_11_e_e") {
      PrototypeANRL_AB_mP4_11_e_e(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 11 // ./aflow --proto=ABC3_mP10_11_e_e_ef --params=4.63,1.20259179266,1.52203023758,110.21,0.121,0.1745,0.3531,0.7086,0.4009,0.1165,0.8544,0.5361,0.6943
    if(vlabel.at(ifound)=="ABC3_mP10_11_e_e_ef") {
      PrototypeANRL_ABC3_mP10_11_e_e_ef(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 12 // ./aflow --proto=A_mP16_11_8e --params=6.183,0.779880316998,1.77308749798,101.79,0.345,0.162,0.767,0.168,0.128,0.34,0.657,0.457,0.025,0.618,0.473,0.653,0.328,-0.074,0.869,0.894
    if(vlabel.at(ifound)=="A_mP16_11_8e") {
      PrototypeANRL_A_mP16_11_8e(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 13 // ./aflow --proto=AB2_mC6_12_a_i --params=7.189,0.613019891501,0.705105021561,90.04,0.6879,0.2889
    if(vlabel.at(ifound)=="AB2_mC6_12_a_i") {
      PrototypeANRL_AB2_mC6_12_a_i(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 14 // ./aflow --proto=A_mC34_12_ah3i2j --params=11.93871,0.876392843113,0.658278825769,129.00411,0.22,0.854,0.241,0.663,0.745,0.566,0.238,0.355,0.232,-0.037,0.333,0.35,0.586
    if(vlabel.at(ifound)=="A_mC34_12_ah3i2j") {
      PrototypeANRL_A_mC34_12_ah3i2j(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 15 // ./aflow --proto=AB3_mC16_12_g_ij --params=5.914,1.73047007102,1.03956712885,108.25,0.1662,0.2147,0.2263,0.2518,0.32131,0.2248
    if(vlabel.at(ifound)=="AB3_mC16_12_g_ij") {
      PrototypeANRL_AB3_mC16_12_g_ij(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 16 // ./aflow --proto=A5B2_mC14_12_a2i_i --params=9.188,0.430343926861,0.705158902917,97.56,0.14286,0.42857,0.28571,0.85714,0.42857,0.28571
    if(vlabel.at(ifound)=="A5B2_mC14_12_a2i_i") {
      PrototypeANRL_A5B2_mC14_12_a2i_i(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 17 // ./aflow --proto=A_mC4_12_i --params=5.403,0.635387747548,0.940033314825,132.32,0.106,0.173
    if(vlabel.at(ifound)=="A_mC4_12_i") {
      PrototypeANRL_A_mC4_12_i(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 18 // ./aflow --proto=ABC4_mP12_13_e_a_2g --params=8.95,0.500335195531,1.63360893855,145.35,0.5182,0.2986,0.0278,0.0003,0.2821,0.4045,0.2366
    if(vlabel.at(ifound)=="ABC4_mP12_13_e_a_2g") {
      PrototypeANRL_ABC4_mP12_13_e_a_2g(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 19 // ./aflow --proto=A_mP84_13_21g --params=9.21,0.99348534202,2.45385450597,106.1,0.30089,0.20127,0.18147,0.17387,0.03262,0.11695,0.05014,-0.05231,0.18035,-0.07589,0.78099,0.11634,0.79463,0.67872,0.1738,0.68463,0.51532,0.10402,0.56601,0.44932,0.17224,0.42424,0.27741,0.11672,0.0412,0.39067,0.07245,-0.00092,0.15881,0.04497,0.78847,0.13878,0.07346,0.7486,-0.09081,0.04464,0.53574,0.87264,0.06842,0.50833,0.63715,0.03304,0.30515,0.63715,0.06617,0.25041,0.40555,0.0442,0.146,0.38905,0.17219,0.86038,0.10055,0.17357,0.59606,0.82384,0.1694,0.41856,0.64581,0.16732,-0.05418,0.32296,0.2006
    if(vlabel.at(ifound)=="A_mP84_13_21g") {
      PrototypeANRL_A_mP84_13_21g(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 20 // ./aflow --proto=A2B_mP12_14_2e_e --params=5.1505,1.01186292593,1.03238520532,99.23,0.07,0.3317,0.3447,0.4496,0.7569,0.4792,0.2754,0.0395,0.2083
    if(vlabel.at(ifound)=="A2B_mP12_14_2e_e") {
      PrototypeANRL_A2B_mP12_14_2e_e(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 21 // ./aflow --proto=A_mP32_14_8e --params=9.31,0.866809881847,1.38023630505,93.13333,0.437,0.185,0.084,0.246,0.273,-0.023,0.24,0.102,0.828,0.05,-0.08,0.852,0.157,0.669,-0.09,0.142,0.66,0.09,0.368,0.746,0.16,0.334,0.021,0.21
    if(vlabel.at(ifound)=="A_mP32_14_8e") {
      PrototypeANRL_A_mP32_14_8e(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 22 // ./aflow --proto=A_mP64_14_16e --params=15.018,0.979691037422,0.585231056066,93.61,0.18313,0.14063,0.03451,0.22856,0.28408,0.12262,0.35548,0.31907,-0.00548,0.47826,0.28776,0.16131,0.52853,0.14438,0.09345,0.47966,0.04033,0.27102,0.36296,-0.02818,0.15123,0.22521,0.04261,0.2343,0.09552,0.48601,0.14213,0.01298,0.58883,0.27815,-0.01931,0.71476,0.12135,0.08347,0.82945,0.18553,0.19177,0.81338,0.00963,0.3102,0.73961,0.14402,0.30834,0.59137,0.04778,0.24353,0.50553,0.23353
    if(vlabel.at(ifound)=="A_mP64_14_16e") {
      PrototypeANRL_A_mP64_14_16e(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 23 // ./aflow --proto=A2B5_mC28_15_f_e2f --params=12.786,0.387533239481,0.427968090099,97.03333,0.5727,0.106,0.311,0.077,0.0958,0.0952,0.4213,0.7127,0.0726,0.3138
    if(vlabel.at(ifound)=="A2B5_mC28_15_f_e2f") {
      PrototypeANRL_A2B5_mC28_15_f_e2f(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 24 // ./aflow --proto=AB_mC8_15_c_e --params=4.6837,0.730747058949,1.0950317057,120.34,0.4184
    if(vlabel.at(ifound)=="AB_mC8_15_c_e") {
      PrototypeANRL_AB_mC8_15_c_e(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 25 // ./aflow --proto=A2B_mC48_15_ae3f_2f --params=7.1356,1.73344918437,1.00532541062,120.34,0.1163,0.266,0.1234,0.9401,0.3114,0.1038,0.3282,0.0172,0.2117,0.4782,0.14033,0.10833,0.07227,0.50682,0.15799,0.54077
    if(vlabel.at(ifound)=="A2B_mC48_15_ae3f_2f") {
      PrototypeANRL_A2B_mC48_15_ae3f_2f(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 26 // ./aflow --proto=ABC6D2_mC40_15_e_e_3f_f --params=9.79,0.901123595506,0.548518896834,105.81,0.3082,-0.0942,0.3888,0.4123,0.8659,0.1365,0.2411,0.6799,0.1468,0.4802,0.0124,0.2117,0.4057,0.7764
    if(vlabel.at(ifound)=="ABC6D2_mC40_15_e_e_3f_f") {
      PrototypeANRL_ABC6D2_mC40_15_e_e_3f_f(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 27 // ./aflow --proto=ABC4_oP12_16_ag_cd_2u --params=5.61,1.01069518717,1.61319073084,0.2,0.26,0.125,0.74,0.8,0.63
    if(vlabel.at(ifound)=="ABC4_oP12_16_ag_cd_2u") {
      PrototypeANRL_ABC4_oP12_16_ag_cd_2u(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 28 // ./aflow --proto=AB3_oP16_18_ab_3c --params=8.32,1.15865384615,0.579326923077,0.0,0.0,0.25,0.25,0.0,0.25,0.5,0.5,0.124,0.309,0.382
    if(vlabel.at(ifound)=="AB3_oP16_18_ab_3c") {
      PrototypeANRL_AB3_oP16_18_ab_3c(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 29 // ./aflow --proto=A2B_oP12_19_2a_a --params=7.764,0.909582689335,0.558088614116,0.185,0.07,0.465,0.055,0.765,-0.008,0.884,-0.011,0.391
    if(vlabel.at(ifound)=="A2B_oP12_19_2a_a") {
      PrototypeANRL_A2B_oP12_19_2a_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 30 // ./aflow --proto=A2B_oC24_20_abc_c --params=8.74,0.576659038902,0.942791762014,0.3336,0.4403,0.2453,0.1971,0.2713,0.33154,0.03589,0.81143
    if(vlabel.at(ifound)=="A2B_oC24_20_abc_c") {
      PrototypeANRL_A2B_oC24_20_abc_c(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 31 // ./aflow --proto=AB_oP2_25_b_a --params=2.8102,1.87104120703,1.0769696107,0.0,0.25
    if(vlabel.at(ifound)=="AB_oP2_25_b_a") {
      PrototypeANRL_AB_oP2_25_b_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 32 // ./aflow --proto=AB2_oP24_28_acd_2c3d --params=16.54,0.533252720677,0.269649334946,0.0,0.319,0.014,0.018,0.042,0.617,0.042,0.624,0.334,0.5,0.503,0.301,0.042,0.632,0.636,0.5,0.619,0.036,0.5
    if(vlabel.at(ifound)=="AB2_oP24_28_acd_2c3d") {
      PrototypeANRL_AB2_oP24_28_acd_2c3d(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 33 // ./aflow --proto=AB3C4_oP16_31_a_ab_2ab --params=7.43,0.869448183042,0.831763122476,0.8268,0.0,0.1514,0.4983,0.8226,0.6454,0.1436,0.1166,0.2466,0.3255,-0.0134,0.2598,0.3364,0.6184
    if(vlabel.at(ifound)=="AB3C4_oP16_31_a_ab_2ab") {
      PrototypeANRL_AB3C4_oP16_31_a_ab_2ab(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 34 // ./aflow --proto=AB_oP8_33_a_a --params=5.2857,1.11007056776,0.659950432298,0.1996,0.5867,0.2506,0.002,0.2003,0.25  
    // Change z1 to be different than z2 to get SG #33
    if(vlabel.at(ifound)=="AB_oP8_33_a_a") {
      PrototypeANRL_AB_oP8_33_a_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 35 // ./aflow --proto=AB3C4_oP32_33_a_3a_4a --params=9.11,1.01866081229,1.1613611416,0.2187,0.4807,0.2031,0.4418,0.2052,0.0015,0.4488,0.1967,0.4146,0.1422,0.9176,0.2246,0.191,0.2506,0.2228,0.3424,0.5361,0.0415,0.0069,0.5876,0.2212,0.3355,0.546,0.3761
    if(vlabel.at(ifound)=="AB3C4_oP32_33_a_3a_4a") {
      PrototypeANRL_AB3C4_oP32_33_a_3a_4a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 36 // ./aflow --proto=A2B_oC12_36_2a_a --params=4.624,1.46820934256,2.69139273356,0.333,0.0,0.061,0.134,0.395,0.366
    if(vlabel.at(ifound)=="A2B_oC12_36_2a_a") {
      PrototypeANRL_A2B_oC12_36_2a_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 37 // ./aflow --proto=A2BC_oC8_38_e_a_b --params=3.875,1.17470967742,1.59019354839,0.0,0.6144,0.155,0.2914
    if(vlabel.at(ifound)=="A2BC_oC8_38_e_a_b") {
      PrototypeANRL_A2BC_oC8_38_e_a_b(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 38 // ./aflow --proto=A2B_oC12_38_de_ab --params=4.684,1.81084543126,1.0269000854,0.06,0.5,0.17,0.56,0.17,0.0
    // Change z3 from z4 to get SG #38 (discussed in paper)
    if(vlabel.at(ifound)=="A2B_oC12_38_de_ab") {
      PrototypeANRL_A2B_oC12_38_de_ab(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 39 // ./aflow --proto=AB4_oC20_41_a_2b --params=6.388,1.00485284909,1.7778647464,0.0,0.673,0.327,0.376,0.827,0.673,0.125
    // Change z3 to get SG #41
    if(vlabel.at(ifound)=="AB4_oC20_41_a_2b") {
      PrototypeANRL_AB4_oC20_41_a_2b(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 40 // ./aflow --proto=AB2_oC24_41_2a_2b --params=6.478,1.0,1.87635072553,0.01,0.238,0.342,0.158,0.125,0.25,0.25,-0.125
    // Change z4 to get SG #41
    if(vlabel.at(ifound)=="AB2_oC24_41_2a_2b") {
      PrototypeANRL_AB2_oC24_41_2a_2b(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 41 // ./aflow --proto=AB2_oF72_43_ab_3b --params=11.66,1.91595197256,0.58833619211,0.0,0.125,0.13889,0.0,0.02222,0.08056,0.18333,0.15278,-0.01389,-0.18333,0.0625,0.125,0.27778
    if(vlabel.at(ifound)=="AB2_oF72_43_ab_3b") {
      PrototypeANRL_AB2_oF72_43_ab_3b(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 42 // ./aflow --proto=AB_oI4_44_a_b --params=4.92,0.973577235772,0.535569105691,0.0,0.425
    if(vlabel.at(ifound)=="AB_oI4_44_a_b") {
      PrototypeANRL_AB_oI4_44_a_b(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 43 // ./aflow --proto=A2B3C7D_oP13_47_t_aq_eqrs_h --params=3.8187,1.01691675177,3.05567339671,0.3554,0.1579,0.3771,0.3788,0.18445
    if(vlabel.at(ifound)=="A2B3C7D_oP13_47_t_aq_eqrs_h") {
      PrototypeANRL_A2B3C7D_oP13_47_t_aq_eqrs_h(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 44 // ./aflow --proto=AB_oP4_51_e_f --params=4.7549,0.661969757513,1.0209678437,0.8125,0.3125
    if(vlabel.at(ifound)=="AB_oP4_51_e_f") {
      PrototypeANRL_AB_oP4_51_e_f(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 45 // ./aflow --proto=A3B2_oP20_56_ce_e --params=4.911,2.53797597231,1.10201588271,0.029,0.147,0.058,0.861,0.044,0.128,0.179
    if(vlabel.at(ifound)=="A3B2_oP20_56_ce_e") {
      PrototypeANRL_A3B2_oP20_56_ce_e(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 46 // ./aflow --proto=ABCD_oP16_57_d_c_d_d --params=6.707,0.997614432682,1.13627553303,0.208,0.7704,0.2871,0.889,0.4154,0.605,0.1087
    if(vlabel.at(ifound)=="ABCD_oP16_57_d_c_d_d") {
      PrototypeANRL_ABCD_oP16_57_d_c_d_d(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 47 // ./aflow --proto=AB_oP8_57_d_d --params=6.09556,0.900425883758,0.850291031505,0.8593,0.0628,0.255,0.0096
    if(vlabel.at(ifound)=="AB_oP8_57_d_d") {
      PrototypeANRL_AB_oP8_57_d_d(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 48 // ./aflow --proto=AB2_oP6_58_a_g --params=6.24,1.03044871795,0.673076923077,0.275,0.325
    // 49 // ./aflow --proto=AB2_oP6_58_a_g --params=4.704,0.917942176871,0.601615646259,0.66667,0.25
    // 50 // ./aflow --proto=AB2_oP6_58_a_g --params=4.4446,1.22049228277,0.761913333033,0.2004,0.3787
    if(vlabel.at(ifound)=="AB2_oP6_58_a_g") {
      PrototypeANRL_AB2_oP6_58_a_g(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 51 // ./aflow --proto=AB_oP4_59_a_b --params=3.15,1.29841269841,2.20634920635,0.051,0.277
    if(vlabel.at(ifound)=="AB_oP4_59_a_b") {
      PrototypeANRL_AB_oP4_59_a_b(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 52 // ./aflow --proto=ABC_oP6_59_a_a_a --params=5.68,0.700704225352,1.01056338028,0.1499,0.4237,0.6255
    if(vlabel.at(ifound)=="ABC_oP6_59_a_a_a") {
      PrototypeANRL_ABC_oP6_59_a_a_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 53 // ./aflow --proto=A3B_oP8_59_bf_a --params=5.162,0.842115459124,0.877760557923,0.67125,0.329,0.505,0.174
    if(vlabel.at(ifound)=="A3B_oP8_59_bf_a") {
      PrototypeANRL_A3B_oP8_59_bf_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 54 // ./aflow --proto=AB_oP16_61_c_c --params=6.471,1.27538247566,1.31757070005,0.136,0.072,0.108,0.456,0.119,0.872
    if(vlabel.at(ifound)=="AB_oP16_61_c_c") {
      PrototypeANRL_AB_oP16_61_c_c(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 55 // ./aflow --proto=A2B_oP24_61_2c_c --params=9.174,0.375953782429,0.560061042075,0.0095,0.1491,0.1835,0.2314,0.111,0.5366,0.1289,0.0972,0.8628
    if(vlabel.at(ifound)=="A2B_oP24_61_2c_c") {
      PrototypeANRL_A2B_oP24_61_2c_c(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 56 // ./aflow --proto=A3B2_oP20_62_3c_2c --params=11.282,0.339443361106,0.994947704308,0.2922,0.19181,0.4504,0.877,0.6246,0.5611,-0.02937,0.17398,0.64939,-0.03603
    if(vlabel.at(ifound)=="A3B2_oP20_62_3c_2c") {
      PrototypeANRL_A3B2_oP20_62_3c_2c(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 57 // ./aflow --proto=AB3C_oP20_62_c_cd_a --params=5.4224,1.41099881971,0.996661994689,0.4877,-0.0084,0.0313,0.0586,0.288,0.537,0.213
    if(vlabel.at(ifound)=="AB3C_oP20_62_c_cd_a") {
      PrototypeANRL_AB3C_oP20_62_c_cd_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 58 // ./aflow --proto=A4B_oP20_62_2cd_c --params=5.464,0.810395314788,1.36749633968,0.22451,0.65626,0.55801,0.6466,0.05131,0.36362,0.13079,0.0579,0.06543
    if(vlabel.at(ifound)=="A4B_oP20_62_2cd_c") {
      PrototypeANRL_A4B_oP20_62_2cd_c(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 59 // ./aflow --proto=AB2C_oP16_62_c_2c_c --params=6.018,0.630741110003,2.4086075108,0.2522,0.8276,0.6221,0.095,0.8706,0.8244,0.226,0.06333
    if(vlabel.at(ifound)=="AB2C_oP16_62_c_2c_c") {
      PrototypeANRL_AB2C_oP16_62_c_2c_c(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 60 // ./aflow --proto=A2B_oP12_62_2c_c --params=4.918,0.7600650671,1.44550630338,0.038,0.282,0.674,0.562,0.202,0.611
    // 61 // ./aflow --proto=A2B_oP12_62_2c_c --params=12.735,0.468237141735,0.339615233608,0.733,0.125,0.508,0.722,0.874,0.447
    // 62 // ./aflow --proto=A2B_oP12_62_2c_c --params=7.6204,0.595008136056,1.1869718125,0.125,0.4217,0.0202,0.837,0.2377,0.0959
    // (from part 2) 77 // ./aflow --proto=A2B_oP12_62_2c_c --params=3.875,1.64232258065,1.89496774194,0.004,0.758,0.24,0.07,0.24,0.39
    if(vlabel.at(ifound)=="A2B_oP12_62_2c_c") {
      PrototypeANRL_A2B_oP12_62_2c_c(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 63 // ./aflow --proto=AB_oP8_62_c_c --params=10.42,0.349328214971,0.411708253359,0.375,0.333,0.139,0.389
    // 64 // ./aflow --proto=AB_oP8_62_c_c --params=5.24160,0.606723137973,1.12622100122,0.0056,0.1952,0.1879,0.5696
    // 68 // ./aflow --proto=AB_oP8_62_c_c --params=5.495,0.536123748863,0.737579617834,0.125,0.69,-0.18,0.125
    // 69 // ./aflow --proto=AB_oP8_62_c_c --params=11.18,0.356171735242,0.387209302326,0.3507,0.0201,0.61937,0.3806
    // (from part 2) 82 // ./aflow --proto=AB_oP8_62_c_c --params=5.454,0.609644297763,1.10542720939,0.2005,0.5741,0.0058,0.1993
    if(vlabel.at(ifound)=="AB_oP8_62_c_c") {
      PrototypeANRL_AB_oP8_62_c_c(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 65 // ./aflow --proto=AB3_oP16_62_c_cd --params=5.09,1.3257367387,0.888605108055,0.39,0.05,0.036,0.852,0.186,0.063,0.328
    if(vlabel.at(ifound)=="AB3_oP16_62_c_cd") {
      PrototypeANRL_AB3_oP16_62_c_cd(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 66 // ./aflow --proto=A3B7_oP40_62_cd_3c2d --params=4.526,1.54882898807,2.68272205038,0.4594,0.5629,0.0579,0.6261,0.2501,0.2063,0.2619,0.4165,0.0288,0.0291,0.3428,0.0565,0.0642,0.8119,0.2509,0.0657,0.0218
    if(vlabel.at(ifound)=="A3B7_oP40_62_cd_3c2d") {
      PrototypeANRL_A3B7_oP40_62_cd_3c2d(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 67 // ./aflow --proto=A_oP8_62_2c --params=6.663,0.708839861924,0.73345339937,0.464,0.292,0.181,0.658
    if(vlabel.at(ifound)=="A_oP8_62_2c") {
      PrototypeANRL_A_oP8_62_2c(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 70 // ./aflow --proto=AB2C_oC16_63_c_2c_c --params=3.577,4.56863293263,1.09538719597,0.06109,-0.0558,0.1792,0.33096
    if(vlabel.at(ifound)=="AB2C_oC16_63_c_2c_c") {
      PrototypeANRL_AB2C_oC16_63_c_2c_c(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 71 // ./aflow --proto=A2B_oC12_63_2c_c --params=3.73,3.94638069705,0.983914209115,0.061,0.75,0.396
    if(vlabel.at(ifound)=="A2B_oC12_63_2c_c") {
      PrototypeANRL_A2B_oC12_63_2c_c(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 72 // ./aflow --proto=AB_oC8_63_c_c --params=2.9782,2.64253575985,0.985360284736,0.436,0.14525
    if(vlabel.at(ifound)=="AB_oC8_63_c_c") {
      PrototypeANRL_AB_oC8_63_c_c(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 73 // ./aflow --proto=A_oC4_63_c --params=2.8444,2.06331739558,1.73379271551,0.10228
    if(vlabel.at(ifound)=="A_oC4_63_c") {
      PrototypeANRL_A_oC4_63_c(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 74 // ./aflow --proto=A_oC8_64_f --params=4.523,1.69378730931,1.0002210922,0.1549,0.081
    // 76 // ./aflow --proto=A_oC8_64_f --params=3.3136,3.16211974891,1.32070859488,0.10168,0.08056
    // 77 // ./aflow --proto=A_oC8_64_f --params=7.11906,0.654575182679,1.37596817557,0.15485,0.1175
    if(vlabel.at(ifound)=="A_oC8_64_f") {
      PrototypeANRL_A_oC8_64_f(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 75 // ./aflow --proto=A2B2C_oC80_64_efg_efg_df --params=10.922,0.866233290606,0.682933528658,0.84657,0.0946,0.9271,0.5886,0.276,-0.0792,0.2314,0.27981,-0.0113,0.1278,0.3415,0.2438,0.1245,0.175,0.2231
    if(vlabel.at(ifound)=="A2B2C_oC80_64_efg_efg_df") {
      PrototypeANRL_A2B2C_oC80_64_efg_efg_df(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 78 // ./aflow --proto=AB_oC8_65_j_g --params=5.971,1.1314687657,0.468263272484,0.28,0.22
    if(vlabel.at(ifound)=="AB_oC8_65_j_g") {
      PrototypeANRL_AB_oC8_65_j_g(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 79 // ./aflow --proto=A3B5_oC16_65_ah_bej --params=8.031,0.926410160628,0.491595069107,0.25,0.225
    if(vlabel.at(ifound)=="A3B5_oC16_65_ah_bej") {
      PrototypeANRL_A3B5_oC16_65_ah_bej(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 80 // ./aflow --proto=AB3_oC8_65_a_bf --params=5.82068,1.35259626023,0.493507631411
    if(vlabel.at(ifound)=="AB3_oC8_65_a_bf") {
      PrototypeANRL_AB3_oC8_65_a_bf(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 81 // ./aflow --proto=AB_oF8_69_a_b --params=6.08,0.903782894737,0.851973684211
    if(vlabel.at(ifound)=="AB_oF8_69_a_b") {
      PrototypeANRL_AB_oF8_69_a_b(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 82 // ./aflow --proto=A_oF8_70_a --params=3.1587,1.82613100326,3.21714629436
    if(vlabel.at(ifound)=="A_oF8_70_a") {
      PrototypeANRL_A_oF8_70_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 83 // ./aflow --proto=A2B_oF24_70_e_a --params=8.2671,0.580614725841,1.0342804611,0.4615
    if(vlabel.at(ifound)=="A2B_oF24_70_e_a") {
      PrototypeANRL_A2B_oF24_70_e_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 84 // ./aflow --proto=A_oF128_70_4h --params=10.4646,1.22947843205,2.33988876785,0.14415,0.04732,0.0486,0.29277,0.2269,0.25406,0.21598,0.28022,0.32618,0.21405,0.15761,0.37947
    if(vlabel.at(ifound)=="A_oF128_70_4h") {
      PrototypeANRL_A_oF128_70_4h(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 85 // ./aflow --proto=AB2_oI6_71_a_i --params=3.144,0.994910941476,2.44179389313,0.339
    if(vlabel.at(ifound)=="AB2_oI6_71_a_i") {
      PrototypeANRL_AB2_oI6_71_a_i(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 86 // ./aflow --proto=AB2_oI6_71_a_g --params=2.75984,2.9999963766,1.4241115427,0.35333
    if(vlabel.at(ifound)=="AB2_oI6_71_a_g") {
      PrototypeANRL_AB2_oI6_71_a_g(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 87 // ./aflow --proto=A2B_oI12_72_j_a --params=9.583,0.585829072316,0.578837524783,0.1182,0.2088
    if(vlabel.at(ifound)=="A2B_oI12_72_j_a") {
      PrototypeANRL_A2B_oI12_72_j_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 88 // ./aflow --proto=AB4C_tI12_82_c_g_a --params=4.3404,1.53216293429,0.256,0.2566,0.3722
    if(vlabel.at(ifound)=="AB4C_tI12_82_c_g_a") {
      PrototypeANRL_AB4C_tI12_82_c_g_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 89 // ./aflow --proto=A2BC4_tI14_82_bc_a_g --params=5.55,1.85585585586,0.26,0.25,0.13
    if(vlabel.at(ifound)=="A2BC4_tI14_82_bc_a_g") {
      PrototypeANRL_A2BC4_tI14_82_bc_a_g(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 90 // ./aflow --proto=AB_tP16_84_cej_k --params=6.429,1.02830922383,0.46779,0.25713,0.19361,0.30754,0.22904
    if(vlabel.at(ifound)=="AB_tP16_84_cej_k") {
      PrototypeANRL_AB_tP16_84_cej_k(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 91 // ./aflow --proto=A4B5_tI18_87_h_ah --params=10.164,0.37111373475,0.2797,-0.0589,0.3752,0.6856
    if(vlabel.at(ifound)=="A4B5_tI18_87_h_ah") {
      PrototypeANRL_A4B5_tI18_87_h_ah(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 92 // ./aflow --proto=AB4_tI10_87_a_h --params=5.72,0.623076923077,0.4,0.8
    if(vlabel.at(ifound)=="AB4_tI10_87_a_h") {
      PrototypeANRL_AB4_tI10_87_a_h(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 93 // ./aflow --proto=A2B_tP12_92_b_a --params=4.957,1.39001412144,0.3047,0.2381,0.1109,0.1826
    if(vlabel.at(ifound)=="A2B_tP12_92_b_a") {
      PrototypeANRL_A2B_tP12_92_b_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 94 // ./aflow --proto=A2B_tP36_96_3b_ab --params=7.464,1.15487674169,0.41,0.445,0.132,0.4,0.117,0.123,0.296,0.344,0.297,0.143,0.326,0.12,0.248
    if(vlabel.at(ifound)=="A2B_tP36_96_3b_ab") {
      PrototypeANRL_A2B_tP36_96_3b_ab(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 95 // ./aflow --proto=A_tP12_96_ab --params=5.51889,1.25999974633,0.0849,0.1752,0.3792,0.2742
    if(vlabel.at(ifound)=="A_tP12_96_ab") {
      PrototypeANRL_A_tP12_96_ab(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 96 // ./aflow --proto=A3BC_tP5_99_bc_a_b --params=4.046,1.02308452793,0.0,0.8973,0.4517,0.3785
    if(vlabel.at(ifound)=="A3BC_tP5_99_bc_a_b") {
      PrototypeANRL_A3BC_tP5_99_bc_a_b(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 97 // ./aflow --proto=AB3_tP8_113_a_ce --params=6.871,0.606622034638,0.206,0.1797,0.476
    if(vlabel.at(ifound)=="AB3_tP8_113_a_ce") {
      PrototypeANRL_AB3_tP8_113_a_ce(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 98 // ./aflow --proto=A2BC4D_tI16_121_d_a_i_b --params=5.46,1.96428571429,0.245,0.132
    if(vlabel.at(ifound)=="A2BC4D_tI16_121_d_a_i_b") {
      PrototypeANRL_A2BC4D_tI16_121_d_a_i_b(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 99 // ./aflow --proto=ABC2_tI16_122_a_b_d --params=5.289,1.97069389299,0.2574
    if(vlabel.at(ifound)=="ABC2_tI16_122_a_b_d") {
      PrototypeANRL_ABC2_tI16_122_a_b_d(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 100 // ./aflow --proto=AB5C_tP7_123_b_ci_a --params=4.207,1.61516520086,0.312
    if(vlabel.at(ifound)=="AB5C_tP7_123_b_ci_a") {
      PrototypeANRL_AB5C_tP7_123_b_ci_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 101 // ./aflow --proto=AB3_tP4_123_a_ce --params=4.158,0.864357864358
    if(vlabel.at(ifound)=="AB3_tP4_123_a_ce") {
      PrototypeANRL_AB3_tP4_123_a_ce(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 102 // ./aflow --proto=AB_tP2_123_a_d --params=2.8,1.31071428571
    if(vlabel.at(ifound)=="AB_tP2_123_a_d") {
      PrototypeANRL_AB_tP2_123_a_d(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 103 // ./aflow --proto=ABC2_tP4_123_d_a_f --params=3.8611,0.828649866618
    if(vlabel.at(ifound)=="ABC2_tP4_123_d_a_f") {
      PrototypeANRL_ABC2_tP4_123_d_a_f(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 104 // ./aflow --proto=A2B3_tP10_127_g_ah --params=7.3364,0.530232811733,0.3841,0.1821
    if(vlabel.at(ifound)=="A2B3_tP10_127_g_ah") {
      PrototypeANRL_A2B3_tP10_127_g_ah(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 105 // ./aflow --proto=ABCD_tP8_129_c_b_a_c --params=3.6736,2.60540069686,0.6793,0.2246
    if(vlabel.at(ifound)=="ABCD_tP8_129_c_b_a_c") {
      PrototypeANRL_ABCD_tP8_129_c_b_a_c(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 106 // ./aflow --proto=A_tP4_129_ac --params=4.897,0.69185215438,0.375
    if(vlabel.at(ifound)=="A_tP4_129_ac") {
      PrototypeANRL_A_tP4_129_ac(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 107 // ./aflow --proto=ABC_tP6_129_c_a_c --params=4.11,1.76301703163,0.6497,0.2058
    if(vlabel.at(ifound)=="ABC_tP6_129_c_a_c") {
      PrototypeANRL_ABC_tP6_129_c_a_c(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 108 // ./aflow --proto=A2B_tP6_129_ac_c --params=4.0006,1.52584612308,0.27,0.7
    if(vlabel.at(ifound)=="A2B_tP6_129_ac_c") {
      PrototypeANRL_A2B_tP6_129_ac_c(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 109 // ./aflow --proto=AB_tP4_129_a_c --params=3.9645,1.26008323874,0.2368
    if(vlabel.at(ifound)=="AB_tP4_129_a_c") {
      PrototypeANRL_AB_tP4_129_a_c(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 110 // ./aflow --proto=AB_tP4_129_c_c --params=3.107,1.90505310589,0.1,0.65
    if(vlabel.at(ifound)=="AB_tP4_129_c_c") {
      PrototypeANRL_AB_tP4_129_c_c(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 111 // ./aflow --proto=AB_tP4_131_c_e --params=4.9073,1.24500234345
    if(vlabel.at(ifound)=="AB_tP4_131_c_e") {
      PrototypeANRL_AB_tP4_131_c_e(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 112 // ./aflow --proto=A_tP50_134_b2m2n --params=8.74,0.575514874142,0.0048,0.1685,0.1305,0.628,0.1695,0.5228,0.1635,0.0753,0.3383,0.1485
    if(vlabel.at(ifound)=="A_tP50_134_b2m2n") {
      PrototypeANRL_A_tP50_134_b2m2n(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 113 // ./aflow --proto=A_tP30_136_bf2ij --params=10.59,0.532011331445,0.1033,0.3667,0.0383,0.5608,0.2354,0.3183,0.27
    if(vlabel.at(ifound)=="A_tP30_136_bf2ij") {
      PrototypeANRL_A_tP30_136_bf2ij(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 114 // ./aflow --proto=AB_tP8_136_g_f --params=4.75,0.576842105263,0.31,0.336
    if(vlabel.at(ifound)=="AB_tP8_136_g_f") {
      PrototypeANRL_AB_tP8_136_g_f(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 115 // ./aflow --proto=A2B_tP6_136_f_a --params=4.5922,0.644005052045,0.30496
    if(vlabel.at(ifound)=="A2B_tP6_136_f_a") {
      PrototypeANRL_A2B_tP6_136_f_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 116 // ./aflow --proto=sigma_tP30_136_bf2ij --params=8.7966,0.518177477662,0.39864,0.13122,0.46349,0.06609,0.73933,0.18267,0.25202
    if(vlabel.at(ifound)=="sigma_tP30_136_bf2ij") {
      PrototypeANRL_sigma_tP30_136_bf2ij(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 117 // ./aflow --proto=A_tP4_136_f --params=3.957,1.29112964367,0.098
    if(vlabel.at(ifound)=="A_tP4_136_f") {
      PrototypeANRL_A_tP4_136_f(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 118 // ./aflow --proto=A_tP16_138_j --params=8.56,0.714953271028,0.375,-0.083,0.857
    if(vlabel.at(ifound)=="A_tP16_138_j") {
      PrototypeANRL_A_tP16_138_j(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 119 // ./aflow --proto=A3B_tI16_139_cde_e --params=3.9993,4.3215062636,0.37498,0.11886
    if(vlabel.at(ifound)=="A3B_tI16_139_cde_e") {
      PrototypeANRL_A3B_tI16_139_cde_e(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 120 // ./aflow --proto=A_tI4_139_e --params=3.34916,1.94217355994,0.819
    if(vlabel.at(ifound)=="A_tI4_139_e") {
      PrototypeANRL_A_tI4_139_e(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 121 // ./aflow --proto=AB2C4_tI14_139_a_e_ce --params=3.7817,3.50337149959,0.36075,0.1824
    if(vlabel.at(ifound)=="AB2C4_tI14_139_a_e_ce") {
      PrototypeANRL_AB2C4_tI14_139_a_e_ce(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 122 // ./aflow --proto=A12B_tI26_139_fij_a --params=8.47,0.584415584416,0.361,0.278
    if(vlabel.at(ifound)=="A12B_tI26_139_fij_a") {
      PrototypeANRL_A12B_tI26_139_fij_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 123 // ./aflow --proto=A_tI2_139_a --params=4.6002,1.07523585931
    // 131 // ./aflow --proto=A_tI2_139_a --params=3.932,0.823499491353
    if(vlabel.at(ifound)=="A_tI2_139_a") {
      PrototypeANRL_A_tI2_139_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 124 // ./aflow --proto=A_tI8_139_h --params=4.33184,0.574102459925,0.17916
    if(vlabel.at(ifound)=="A_tI8_139_h") {
      PrototypeANRL_A_tI8_139_h(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 125 // ./aflow --proto=A3B_tI8_139_bd_a --params=3.8537,2.22744375535
    if(vlabel.at(ifound)=="A3B_tI8_139_bd_a") {
      PrototypeANRL_A3B_tI8_139_bd_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 126 // ./aflow --proto=AB2_tI6_139_a_e --params=3.2064,2.44754241517,0.3353
    if(vlabel.at(ifound)=="AB2_tI6_139_a_e") {
      PrototypeANRL_AB2_tI6_139_a_e(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 127 // ./aflow --proto=A4B5_tI18_139_i_ah --params=8.91,0.361391694725,0.328,0.348
    if(vlabel.at(ifound)=="A4B5_tI18_139_i_ah") {
      PrototypeANRL_A4B5_tI18_139_i_ah(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 128 // ./aflow --proto=A4B_tI10_139_de_a --params=4.53,2.45033112583,0.38
    if(vlabel.at(ifound)=="A4B_tI10_139_de_a") {
      PrototypeANRL_A4B_tI10_139_de_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 129 // ./aflow --proto=A8B_tI18_139_hi_a --params=8.312,0.468840230991,0.333,0.327
    if(vlabel.at(ifound)=="A8B_tI18_139_hi_a") {
      PrototypeANRL_A8B_tI18_139_hi_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 130 // ./aflow --proto=A2B_tI6_139_d_a --params=4.1,1.22682926829
    if(vlabel.at(ifound)=="A2B_tI6_139_d_a") {
      PrototypeANRL_A2B_tI6_139_d_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 132 // ./aflow --proto=A2B_tI12_140_h_a --params=6.04,0.804635761589,0.158
    if(vlabel.at(ifound)=="A2B_tI12_140_h_a") {
      PrototypeANRL_A2B_tI12_140_h_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 133 // ./aflow --proto=AB3_tI16_140_b_ah --params=6.017,1.44241316271,0.231
    if(vlabel.at(ifound)=="AB3_tI16_140_b_ah") {
      PrototypeANRL_AB3_tI16_140_b_ah(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 134 // ./aflow --proto=AB_tI16_140_ab_h --params=8.03,0.87297633873,0.179
    if(vlabel.at(ifound)=="AB_tI16_140_ab_h") {
      PrototypeANRL_AB_tI16_140_ab_h(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 135 // ./aflow --proto=A4BC_tI24_141_h_b_a --params=6.6042,0.905423821205,0.066,0.1951
    if(vlabel.at(ifound)=="A4BC_tI24_141_h_b_a") {
      PrototypeANRL_A4BC_tI24_141_h_b_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 136 // ./aflow --proto=A_tI4_141_a --params=5.8318,0.545611989437
    if(vlabel.at(ifound)=="A_tI4_141_a") {
      PrototypeANRL_A_tI4_141_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 137 // ./aflow --proto=A3B4_tI28_141_ad_h --params=5.765,1.63781439722,0.0278,0.2589
    if(vlabel.at(ifound)=="A3B4_tI28_141_ad_h") {
      PrototypeANRL_A3B4_tI28_141_ad_h(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 138 // ./aflow --proto=A2B_tI12_141_e_a --params=3.785,2.51360634082,0.08306
    // (from part 2) 184 // ./aflow --proto=A2B_tI12_141_e_a --params=4.126,3.47697527872,0.2915
    if(vlabel.at(ifound)=="A2B_tI12_141_e_a") {
      PrototypeANRL_A2B_tI12_141_e_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 139 // ./aflow --proto=AB_tI16_141_e_e --params=3.108,5.45045045045,0.227,0.071
    if(vlabel.at(ifound)=="AB_tI16_141_e_e") {
      PrototypeANRL_AB_tI16_141_e_e(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 140 // ./aflow --proto=A2B_tI24_141_2e_e --params=4.046,6.28917449333,0.125,0.289,-0.051
    if(vlabel.at(ifound)=="A2B_tI24_141_2e_e") {
      PrototypeANRL_A2B_tI24_141_2e_e(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 141 // ./aflow --proto=AB_tI8_141_a_b --params=3.325,3.42255639098
    if(vlabel.at(ifound)=="AB_tI8_141_a_b") {
      PrototypeANRL_AB_tI8_141_a_b(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 142 // ./aflow --proto=A2B3_tI80_141_ceh_3h --params=7.5937,4.26037373086,0.2044,0.5201,0.3324,0.516,0.2547,0.494,0.0859,0.4667,0.4164
    if(vlabel.at(ifound)=="A2B3_tI80_141_ceh_3h") {
      PrototypeANRL_A2B3_tI80_141_ceh_3h(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 143 // ./aflow --proto=ABC4_tI96_142_e_ab_2g --params=10.914,1.77396005131,0.0375,0.2482,0.3197,-0.0867,0.0923,0.1117,0.0025
    if(vlabel.at(ifound)=="ABC4_tI96_142_e_ab_2g") {
      PrototypeANRL_ABC4_tI96_142_e_ab_2g(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 144 // ./aflow --proto=A2B_hP9_147_g_ad --params=7.636,0.369264012572,0.25,0.33333,0.0,0.25
    if(vlabel.at(ifound)=="A2B_hP9_147_g_ad") {
      PrototypeANRL_A2B_hP9_147_g_ad(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 145 // RHL: ./aflow --proto=AB_hR16_148_cf_cf --params=6.29713,1.8633345667,0.11546,0.21,0.10706,0.81289,0.19519,0.1848,0.6754,0.3468
    // 145 // HEX: ./aflow --proto=AB_hR16_148_cf_cf --params=6.29713,1.8633345667,0.11546,0.21,0.10706,0.81289,0.19519,0.1848,0.6754,0.3468 --hex
    if(vlabel.at(ifound)=="AB_hR16_148_cf_cf") {
      PrototypeANRL_AB_hR16_148_cf_cf(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 146 // RHL: ./aflow --proto=AB3_hR8_148_c_f --params=7.49626,2.75900649124,0.33333,0.088,0.755,0.421
    // 146 // HEX: ./aflow --proto=AB3_hR8_148_c_f --params=7.49626,2.75900649124,0.33333,0.088,0.755,0.421 --hex
    if(vlabel.at(ifound)=="AB3_hR8_148_c_f") {
      PrototypeANRL_AB3_hR8_148_c_f(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 147 // RHL: ./aflow --proto=AB_hR26_148_b2f_a2f --params=15.659,0.335334312536,0.054,0.346,0.098,0.754,0.15699,0.6,0.555,0.84401,0.599,0.252,0.65501,0.098
    // 147 // HEX: ./aflow --proto=AB_hR26_148_b2f_a2f --params=15.659,0.335334312536,0.054,0.346,0.098,0.754,0.15699,0.6,0.555,0.84401,0.599,0.252,0.65501,0.098 --hex
    if(vlabel.at(ifound)=="AB_hR26_148_b2f_a2f") {
      PrototypeANRL_AB_hR26_148_b2f_a2f(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 148 // RHL: ./aflow --proto=AB3C_hR10_148_c_f_c --params=5.0884,2.76815894977,0.35537,0.1464,0.22174,0.56249,0.95095
    // 148 // HEX: ./aflow --proto=AB3C_hR10_148_c_f_c --params=5.0884,2.76815894977,0.35537,0.1464,0.22174,0.56249,0.95095 --hex
    if(vlabel.at(ifound)=="AB3C_hR10_148_c_f_c") {
      PrototypeANRL_AB3C_hR10_148_c_f_c(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 149 // ./aflow --proto=A2B_hP9_150_ef_bd --params=5.85,0.589743589744,0.875,0.26,0.6
    if(vlabel.at(ifound)=="A2B_hP9_150_ef_bd") {
      PrototypeANRL_A2B_hP9_150_ef_bd(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 150 // ./aflow --proto=A3B_hP24_151_3c_2a --params=6.017,2.87518697025,0.8889,0.5556,0.8889,0.1111,0.0731,0.5556,0.4444,0.0731,0.2222,0.77778,0.0731
    if(vlabel.at(ifound)=="A3B_hP24_151_3c_2a") {
      PrototypeANRL_A3B_hP24_151_3c_2a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 151 // ./aflow --proto=A2B_hP9_152_c_a --params=4.914,1.10012210012,0.4699,0.413,0.2668,0.214
    if(vlabel.at(ifound)=="A2B_hP9_152_c_a") {
      PrototypeANRL_A2B_hP9_152_c_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 152 // ./aflow --proto=A_hP3_152_a --params=4.3662,1.13453346159,0.2254
    if(vlabel.at(ifound)=="A_hP3_152_a") {
      PrototypeANRL_A_hP3_152_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 153 // ./aflow --proto=AB_hP6_154_a_b --params=4.145,2.29095295537,0.7198,0.4889
    if(vlabel.at(ifound)=="AB_hP6_154_a_b") {
      PrototypeANRL_AB_hP6_154_a_b(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 154 // RHL: ./aflow --proto=AB3_hR8_155_c_de --params=4.91608,2.53341483458,0.237,0.43,0.07
    // 154 // HEX: ./aflow --proto=AB3_hR8_155_c_de --params=4.91608,2.53341483458,0.237,0.43,0.07 --hex
    if(vlabel.at(ifound)=="AB3_hR8_155_c_de") {
      PrototypeANRL_AB3_hR8_155_c_de(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 155 // RHL: ./aflow --proto=A3B2_hR5_155_e_c --params=5.73296,1.24097324942,0.2521,0.2449
    // 155 // HEX: ./aflow --proto=A3B2_hR5_155_e_c --params=5.73296,1.24097324942,0.2521,0.2449 --hex
    if(vlabel.at(ifound)=="A3B2_hR5_155_e_c") {
      PrototypeANRL_A3B2_hR5_155_e_c(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 156 // RHL: ./aflow --proto=AB_hR6_160_b_b --params=9.619,0.327466472606,0.00019,0.26362,0.7288,0.39161
    // 156 // HEX: ./aflow --proto=AB_hR6_160_b_b --params=9.619,0.327466472606,0.00019,0.26362,0.7288,0.39161 --hex
    if(vlabel.at(ifound)=="AB_hR6_160_b_b") {
      PrototypeANRL_AB_hR6_160_b_b(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 157 // RHL: ./aflow --proto=AB_hR6_160_3a_3a --params=3.01791,7.34847294982,0.0,0.22222,0.77778,0.08333,0.30556,0.86111
    // 157 // HEX: ./aflow --proto=AB_hR6_160_3a_3a --params=3.01791,7.34847294982,0.0,0.22222,0.77778,0.08333,0.30556,0.86111 --hex
    if(vlabel.at(ifound)=="AB_hR6_160_3a_3a") {
      PrototypeANRL_AB_hR6_160_3a_3a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 158 // RHL: ./aflow --proto=ABC3_hR10_161_a_a_b --params=5.2542,2.64091583876,0.2875,0.0128,0.74643,0.14093,0.36263
    // 158 // HEX: ./aflow --proto=ABC3_hR10_161_a_a_b --params=5.2542,2.64091583876,0.2875,0.0128,0.74643,0.14093,0.36263 --hex
    if(vlabel.at(ifound)=="ABC3_hR10_161_a_a_b") {
      PrototypeANRL_ABC3_hR10_161_a_a_b(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 159 // ./aflow --proto=AB2_hP9_162_ad_k --params=4.917,0.929021761237,0.325,0.272
    if(vlabel.at(ifound)=="AB2_hP9_162_ad_k") {
      PrototypeANRL_AB2_hP9_162_ad_k(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 160 // ./aflow --proto=AB2CD2_hP36_163_h_i_bf_i --params=7.384,2.37716684724,0.01,0.833,0.33333,0.03833,0.141,0.03167,0.365,0.083
    if(vlabel.at(ifound)=="AB2CD2_hP36_163_h_i_bf_i") {
      PrototypeANRL_AB2CD2_hP36_163_h_i_bf_i(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 161 // ./aflow --proto=A3B2_hP5_164_ad_d --params=4.0282,1.21409066084,0.648,0.149
    if(vlabel.at(ifound)=="A3B2_hP5_164_ad_d") {
      PrototypeANRL_A3B2_hP5_164_ad_d(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 162 // ./aflow --proto=AB2_hP3_164_a_d --params=4.24,1.61320754717,0.252
    if(vlabel.at(ifound)=="AB2_hP3_164_a_d") {
      PrototypeANRL_AB2_hP3_164_a_d(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 163 // ./aflow --proto=A3B_hP24_165_adg_f --params=6.308,1.03994927077,0.167,0.666,0.356,0.028,0.096
    if(vlabel.at(ifound)=="A3B_hP24_165_adg_f") {
      PrototypeANRL_A3B_hP24_165_adg_f(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 164 // RHL: ./aflow --proto=AB_hR2_166_a_b --params=3.13,4.78594249201
    // 164 // HEX: ./aflow --proto=AB_hR2_166_a_b --params=3.13,4.78594249201 --hex
    if(vlabel.at(ifound)=="AB_hR2_166_a_b") {
      PrototypeANRL_AB_hR2_166_a_b(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 165 // RHL: ./aflow --proto=A_hR2_166_c --params=3.7595,2.7815666977,0.22754
    // 165 // HEX: ./aflow --proto=A_hR2_166_c --params=3.7595,2.7815666977,0.22754 --hex
    // 172 // RHL: ./aflow --proto=A_hR2_166_c --params=2.456,4.08957654723,0.16667
    // 172 // HEX: ./aflow --proto=A_hR2_166_c --params=2.456,4.08957654723,0.16667 --hex
    // 175 // RHL: ./aflow --proto=A_hR2_166_c --params=3.289,3.42991790818,0.0543 
    // 175 // HEX: ./aflow --proto=A_hR2_166_c --params=3.289,3.42991790818,0.0543 --hex
    if(vlabel.at(ifound)=="A_hR2_166_c") {
      PrototypeANRL_A_hR2_166_c(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 166 // RHL: ./aflow --proto=A_hR1_166_a --params=5.07846,0.968139947937
    // 166 // HEX: ./aflow --proto=A_hR1_166_a --params=5.07846,0.968139947937 --hex
    // 170 // RHL: ./aflow --proto=A_hR1_166_a --params=3.45741,1.92728082582
    // 170 // HEX: ./aflow --proto=A_hR1_166_a --params=3.45741,1.92728082582 --hex
    if(vlabel.at(ifound)=="A_hR1_166_a") {
      PrototypeANRL_A_hR1_166_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 167 // RHL: ./aflow --proto=A7B6_hR13_166_ah_3c --params=4.757,5.4319949548,0.167,0.346,0.448,0.09,0.59001
    // 167 // HEX: ./aflow --proto=A7B6_hR13_166_ah_3c --params=4.757,5.4319949548,0.167,0.346,0.448,0.09,0.59001 --hex
    if(vlabel.at(ifound)=="A7B6_hR13_166_ah_3c") {
      PrototypeANRL_A7B6_hR13_166_ah_3c(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 168 // RHL: ./aflow --proto=A_hR3_166_ac --params=3.62036,7.25049442597,0.22222
    // 168 // HEX: ./aflow --proto=A_hR3_166_ac --params=3.62036,7.25049442597,0.22222 --hex
    if(vlabel.at(ifound)=="A_hR3_166_ac") {
      PrototypeANRL_A_hR3_166_ac(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 169 // RHL: ./aflow --proto=A2B3_hR5_166_c_ac --params=4.36914,6.96313919902,0.399,0.208
    // 169 // HEX: ./aflow --proto=A2B3_hR5_166_c_ac --params=4.36914,6.96313919902,0.399,0.208 --hex
    if(vlabel.at(ifound)=="A2B3_hR5_166_c_ac") {
      PrototypeANRL_A2B3_hR5_166_c_ac(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 171 // RHL: ./aflow --proto=A5B2_hR7_166_a2c_c --params=3.011,6.9511790103,0.186,0.33333,0.075
    // 171 // HEX: ./aflow --proto=A5B2_hR7_166_a2c_c --params=3.011,6.9511790103,0.186,0.33333,0.075 --hex
    if(vlabel.at(ifound)=="A5B2_hR7_166_a2c_c") {
      PrototypeANRL_A5B2_hR7_166_a2c_c(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 173 // RHL: ./aflow --proto=A_hR12_166_2h --params=4.908,2.56022616137,0.0104,0.65729,0.2206,0.6323
    // 173 // HEX: ./aflow --proto=A_hR12_166_2h --params=4.908,2.56022616137,0.0104,0.65729,0.2206,0.6323 --hex
    if(vlabel.at(ifound)=="A_hR12_166_2h") {
      PrototypeANRL_A_hR12_166_2h(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 174 // RHL: ./aflow --proto=ABC2_hR4_166_a_b_c --params=3.5561,5.44557239673,0.2667
    // 174 // HEX: ./aflow --proto=ABC2_hR4_166_a_b_c --params=3.5561,5.44557239673,0.2667 --hex
    if(vlabel.at(ifound)=="ABC2_hR4_166_a_b_c") {
      PrototypeANRL_ABC2_hR4_166_a_b_c(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 176 // RHL: ./aflow --proto=A_hR105_166_bc9h4i --params=10.96,2.17974452555,0.3848,0.3843,0.21309,0.4895,0.21780,0.3873,0.56899,0.1991,0.50609,0.1983,0.68740,0.1032,0.49209,0.9933,0.66980,0.1008,0.83740,0.0025,0.16801,0.3622,0.58109,0.0976,0.3765,0.68261,0.2024,0.1673,0.55209,0.8921,0.1777,0.3473,0.0033
    // 176 // HEX: ./aflow --proto=A_hR105_166_bc9h4i --params=10.96,2.17974452555,0.3848,0.3843,0.21309,0.4895,0.21780,0.3873,0.56899,0.1991,0.50609,0.1983,0.68740,0.1032,0.49209,0.9933,0.66980,0.1008,0.83740,0.0025,0.16801,0.3622,0.58109,0.0976,0.3765,0.68261,0.2024,0.1673,0.55209,0.8921,0.1777,0.3473,0.0033 --hex
    if(vlabel.at(ifound)=="A_hR105_166_bc9h4i") {
      PrototypeANRL_A_hR105_166_bc9h4i(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 177 // RHL: ./aflow --proto=A6B_hR7_166_g_a --params=4.33304,3.13251204697,0.16667
    // 177 // HEX: ./aflow --proto=A6B_hR7_166_g_a --params=4.33304,3.13251204697,0.16667 --hex
    if(vlabel.at(ifound)=="A6B_hR7_166_g_a") {
      PrototypeANRL_A6B_hR7_166_g_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 178 // RHL: ./aflow --proto=ABC3_hR10_167_a_b_e --params=5.285,2.62039735099,0.85756666666667 
    // 178 // HEX: ./aflow --proto=ABC3_hR10_167_a_b_e --params=5.285,2.62039735099,0.85756666666667 --hex
    // 179 // RHL: ./aflow --proto=ABC3_hR10_167_a_b_e --params=4.988,3.42040898156,0.5067 
    // 179 // HEX: ./aflow --proto=ABC3_hR10_167_a_b_e --params=4.988,3.42040898156,0.5067 --hex
    if(vlabel.at(ifound)=="ABC3_hR10_167_a_b_e") {
      PrototypeANRL_ABC3_hR10_167_a_b_e(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 180 // RHL: ./aflow --proto=A2B3_hR10_167_c_e --params=4.7607,2.72957758313,0.35216,0.5561
    // 180 // HEX: ./aflow --proto=A2B3_hR10_167_c_e --params=4.7607,2.72957758313,0.35216,0.5561 --hex
    if(vlabel.at(ifound)=="A2B3_hR10_167_c_e") {
      PrototypeANRL_A2B3_hR10_167_c_e(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 181 // ./aflow --proto=A2B_hP18_180_fi_bd --params=5.198,2.54136206233,0.163,0.1141
    if(vlabel.at(ifound)=="A2B_hP18_180_fi_bd") {
      PrototypeANRL_A2B_hP18_180_fi_bd(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 182 // ./aflow --proto=AB2_hP9_180_d_j --params=4.42758,1.43826876081,0.16559
    if(vlabel.at(ifound)=="AB2_hP9_180_d_j") {
      PrototypeANRL_AB2_hP9_180_d_j(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 183 // ./aflow --proto=A2B_hP9_180_j_c --params=4.9977,1.09252256038,0.2072
    if(vlabel.at(ifound)=="A2B_hP9_180_j_c") {
      PrototypeANRL_A2B_hP9_180_j_c(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 184 // ./aflow --proto=AB3_hP8_182_c_g --params=4.8507,0.866967654153,0.3249
    if(vlabel.at(ifound)=="AB3_hP8_182_c_g") {
      PrototypeANRL_AB3_hP8_182_c_g(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 185 // ./aflow --proto=A_hP4_186_ab --params=2.47,2.75303643725,0.0,0.07143
    if(vlabel.at(ifound)=="A_hP4_186_ab") {
      PrototypeANRL_A_hP4_186_ab(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 186 // ./aflow --proto=AB_hP8_186_ab_ab --params=3.08051,3.27374363336,0.18784,0.0,0.43671,0.24982
    if(vlabel.at(ifound)=="AB_hP8_186_ab_ab") {
      PrototypeANRL_AB_hP8_186_ab_ab(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 187 // ./aflow --proto=AB_hP4_186_b_b --params=3.8227,1.63776911607,0.3748,0.0
    if(vlabel.at(ifound)=="AB_hP4_186_b_b") {
      PrototypeANRL_AB_hP4_186_b_b(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 188 // ./aflow --proto=AB_hP12_186_a2b_a2b --params=3.08129,4.90695780014,0.1254,0.0,0.29215,-0.0415,0.16675,0.8335
    if(vlabel.at(ifound)=="AB_hP12_186_a2b_a2b") {
      PrototypeANRL_AB_hP12_186_a2b_a2b(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 189 // ./aflow --proto=A5B3C_hP18_186_2a3b_2ab_b --params=3.281,6.57726302956,0.155,0.345,0.0,0.248,0.045,0.261,0.455,0.367,0.137
    if(vlabel.at(ifound)=="A5B3C_hP18_186_2a3b_2ab_b") {
      PrototypeANRL_A5B3C_hP18_186_2a3b_2ab_b(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 190 // ./aflow --proto=AB_hP4_186_b_a --params=2.51,2.66932270916,0.0,0.05
    if(vlabel.at(ifound)=="AB_hP4_186_b_a") {
      PrototypeANRL_AB_hP4_186_b_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 191 // ./aflow --proto=ABC_hP3_187_a_d_f --params=4.535,1.0769570011
    if(vlabel.at(ifound)=="ABC_hP3_187_a_d_f") {
      PrototypeANRL_ABC_hP3_187_a_d_f(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 192 // ./aflow --proto=AB_hP2_187_d_a --params=2.9065,0.975950455875
    if(vlabel.at(ifound)=="AB_hP2_187_d_a") {
      PrototypeANRL_AB_hP2_187_d_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 193 // ./aflow --proto=A2B_hP9_189_fg_bc --params=5.877,0.584822188191,0.256,0.589
    if(vlabel.at(ifound)=="A2B_hP9_189_fg_bc") {
      PrototypeANRL_A2B_hP9_189_fg_bc(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 194 // ./aflow --proto=AB4C_hP6_191_a_h_b --params=3.04436,2.20489035462,0.2413
    if(vlabel.at(ifound)=="AB4C_hP6_191_a_h_b") {
      PrototypeANRL_AB4C_hP6_191_a_h_b(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 195 // ./aflow --proto=AB5_hP6_191_a_cg --params=5.405,0.773913043478
    if(vlabel.at(ifound)=="AB5_hP6_191_a_cg") {
      PrototypeANRL_AB5_hP6_191_a_cg(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 196 // ./aflow --proto=A_hP1_191_a --params=3.2062,0.931195808122
    if(vlabel.at(ifound)=="A_hP1_191_a") {
      PrototypeANRL_A_hP1_191_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 197 // ./aflow --proto=A3B_hP4_191_bc_a --params=3.6576,1.05902777778
    if(vlabel.at(ifound)=="A3B_hP4_191_bc_a") {
      PrototypeANRL_A3B_hP4_191_bc_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 198 // ./aflow --proto=AB2_hP3_191_a_d --params=3.005,1.08276206323
    if(vlabel.at(ifound)=="AB2_hP3_191_a_d") {
      PrototypeANRL_AB2_hP3_191_a_d(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 199 // ./aflow --proto=A2B_hP6_191_h_e --params=4.237,1.71040830776,0.306,0.16
    if(vlabel.at(ifound)=="A2B_hP6_191_h_e") {
      PrototypeANRL_A2B_hP6_191_h_e(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 200 // ./aflow --proto=AB_hP6_191_f_ad --params=5.279,0.806914188293
    if(vlabel.at(ifound)=="AB_hP6_191_f_ad") {
      PrototypeANRL_AB_hP6_191_f_ad(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 201 // ./aflow --proto=AB_hP8_194_ad_f --params=3.64,3.37362637363,0.125
    if(vlabel.at(ifound)=="AB_hP8_194_ad_f") {
      PrototypeANRL_AB_hP8_194_ad_f(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 202 // ./aflow --proto=A_hP6_194_h --params=4.40445,0.568892824303,0.44799
    if(vlabel.at(ifound)=="A_hP6_194_h") {
      PrototypeANRL_A_hP6_194_h(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 203 // ./aflow --proto=AB_hP12_194_af_bf --params=3.01,4.85382059801,0.166,0.583
    if(vlabel.at(ifound)=="AB_hP12_194_af_bf") {
      PrototypeANRL_AB_hP12_194_af_bf(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 204 // ./aflow --proto=A_hP4_194_ac --params=3.77,3.2175066313
    if(vlabel.at(ifound)=="A_hP4_194_ac") {
      PrototypeANRL_A_hP4_194_ac(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 205 // ./aflow --proto=AB3_hP8_194_c_bf --params=5.088,1.76533018868,-0.083
    if(vlabel.at(ifound)=="AB3_hP8_194_c_bf") {
      PrototypeANRL_AB3_hP8_194_c_bf(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 206 // ./aflow --proto=AB2_hP6_194_b_f --params=4.895,1.58324821246,0.045
    if(vlabel.at(ifound)=="AB2_hP6_194_b_f") {
      PrototypeANRL_AB2_hP6_194_b_f(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 207 // ./aflow --proto=AB_hP4_194_c_d --params=2.50399,2.66023426611
    if(vlabel.at(ifound)=="AB_hP4_194_c_d") {
      PrototypeANRL_AB_hP4_194_c_d(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 208 // ./aflow --proto=ABC2_hP8_194_d_a_f --params=2.86,4.48251748252,0.086
    if(vlabel.at(ifound)=="ABC2_hP8_194_d_a_f") {
      PrototypeANRL_ABC2_hP8_194_d_a_f(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 209 // ./aflow --proto=A3B_hP8_194_h_c --params=5.295,0.802077431539,0.8392
    if(vlabel.at(ifound)=="A3B_hP8_194_h_c") {
      PrototypeANRL_A3B_hP8_194_h_c(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 210 // ./aflow --proto=A_hP4_194_bc --params=2.464,2.72362012987
    if(vlabel.at(ifound)=="A_hP4_194_bc") {
      PrototypeANRL_A_hP4_194_bc(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 211 // ./aflow --proto=AB2_hP6_194_c_f --params=3.161,3.8895919013,0.6275
    if(vlabel.at(ifound)=="AB2_hP6_194_c_f") {
      PrototypeANRL_AB2_hP6_194_c_f(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 212 // ./aflow --proto=A5B2_hP14_194_abdf_f --params=2.982,4.651240778,0.528,0.139
    if(vlabel.at(ifound)=="A5B2_hP14_194_abdf_f") {
      PrototypeANRL_A5B2_hP14_194_abdf_f(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 213 // ./aflow --proto=AB2_hP12_194_f_ah --params=5.223,1.64005360904,0.06286,0.83048
    if(vlabel.at(ifound)=="AB2_hP12_194_f_ah") {
      PrototypeANRL_AB2_hP12_194_f_ah(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 214 // ./aflow --proto=ABC_hP6_194_c_d_a --params=2.752,2.56468023256
    if(vlabel.at(ifound)=="ABC_hP6_194_c_d_a") {
      PrototypeANRL_ABC_hP6_194_c_d_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 215 // ./aflow --proto=A_hP4_194_f --params=2.508,1.66786283892,0.05995
    if(vlabel.at(ifound)=="A_hP4_194_f") {
      PrototypeANRL_A_hP4_194_f(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 216 // ./aflow --proto=AB2_hP6_194_c_ad --params=4.186,1.22527472527
    if(vlabel.at(ifound)=="AB2_hP6_194_c_ad") {
      PrototypeANRL_AB2_hP6_194_c_ad(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 217 // ./aflow --proto=AB3C4_hP16_194_c_af_ef --params=2.988,7.82195448461,0.1543,0.605,0.0539
    if(vlabel.at(ifound)=="AB3C4_hP16_194_c_af_ef") {
      PrototypeANRL_AB3C4_hP16_194_c_af_ef(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 218 // ./aflow --proto=A_hP2_194_c --params=3.2093,1.62359393014
    if(vlabel.at(ifound)=="A_hP2_194_c") {
      PrototypeANRL_A_hP2_194_c(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 219 // ./aflow --proto=AB2_hP24_194_ef_fgh --params=4.824,3.28067993367,0.04598,0.84417,0.12514,0.16429
    if(vlabel.at(ifound)=="AB2_hP24_194_ef_fgh") {
      PrototypeANRL_AB2_hP24_194_ef_fgh(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 220 // ./aflow --proto=AB_hP12_194_df_ce --params=3.976,4.12022132797,0.0637,0.10724
    if(vlabel.at(ifound)=="AB_hP12_194_df_ce") {
      PrototypeANRL_AB_hP12_194_df_ce(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 221 // ./aflow --proto=AB_hP4_194_c_a --params=3.619,1.39375518099
    if(vlabel.at(ifound)=="AB_hP4_194_c_a") {
      PrototypeANRL_AB_hP4_194_c_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 222 // ./aflow --proto=A2B_hP12_194_cg_f --params=5.052,1.63697545527,0.062
    if(vlabel.at(ifound)=="A2B_hP12_194_cg_f") {
      PrototypeANRL_A2B_hP12_194_cg_f(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 223 // ./aflow --proto=A4B_cI40_197_cde_c --params=8.4295,0.1668,0.3345,0.6476,0.7484
    if(vlabel.at(ifound)=="A4B_cI40_197_cde_c") {
      PrototypeANRL_A4B_cI40_197_cde_c(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 224 // ./aflow --proto=ABC_cP12_198_a_a_a --params=5.881,-0.024,0.39,0.875
    if(vlabel.at(ifound)=="ABC_cP12_198_a_a_a") {
      PrototypeANRL_ABC_cP12_198_a_a_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 225 // ./aflow --proto=A3B_cP16_198_b_a --params=5.1305,0.2107,0.3689,0.2671,0.1159
    if(vlabel.at(ifound)=="A3B_cP16_198_b_a") {
      PrototypeANRL_A3B_cP16_198_b_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 226 // ./aflow --proto=A_cP8_198_2a --params=5.65,0.0699,-0.0378
    if(vlabel.at(ifound)=="A_cP8_198_2a") {
      PrototypeANRL_A_cP8_198_2a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 227 // ./aflow --proto=AB_cP8_198_a_a --params=5.63,-0.042,0.067
    // 228 // ./aflow --proto=AB_cP8_198_a_a --params=4.48688,0.13652,0.8424
    if(vlabel.at(ifound)=="AB_cP8_198_a_a") {
      PrototypeANRL_AB_cP8_198_a_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 229 // ./aflow --proto=AB_cI16_199_a_a --params=6.3557,0.294,0.0347
    if(vlabel.at(ifound)=="AB_cI16_199_a_a") {
      PrototypeANRL_AB_cI16_199_a_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 230 // ./aflow --proto=AB32C48_cI162_204_a_2efg_2gh --params=14.16,0.8203,0.5998,0.1836,0.2942,0.8806,0.0908,0.8499,0.1748,0.6993,0.686,0.0969,0.332
    if(vlabel.at(ifound)=="AB32C48_cI162_204_a_2efg_2gh") {
      PrototypeANRL_AB32C48_cI162_204_a_2efg_2gh(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 231 // ./aflow --proto=A3B_cI32_204_g_c --params=7.58,0.3431,0.8497
    if(vlabel.at(ifound)=="A3B_cI32_204_g_c") {
      PrototypeANRL_A3B_cI32_204_g_c(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 232 // ./aflow --proto=A12B_cI26_204_g_a --params=7.58,0.184,0.691
    if(vlabel.at(ifound)=="A12B_cI26_204_g_a") {
      PrototypeANRL_A12B_cI26_204_g_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 233 // ./aflow --proto=A_cP8_205_c --params=5.65,0.05569
    if(vlabel.at(ifound)=="A_cP8_205_c") {
      PrototypeANRL_A_cP8_205_c(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 234 // ./aflow --proto=AB_cP16_205_c_c --params=6.4162,0.1527,0.6297
    if(vlabel.at(ifound)=="AB_cP16_205_c_c") {
      PrototypeANRL_AB_cP16_205_c_c(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 235 // ./aflow --proto=AB2_cP12_205_a_c --params=5.417,0.3851
    if(vlabel.at(ifound)=="AB2_cP12_205_a_c") {
      PrototypeANRL_AB2_cP12_205_a_c(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 236 // ./aflow --proto=AB3C6_cI80_206_a_d_e --params=9.4,-0.0344,0.338,0.1,0.125
    if(vlabel.at(ifound)=="AB3C6_cI80_206_a_d_e") {
      PrototypeANRL_AB3C6_cI80_206_a_d_e(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 237 // ./aflow --proto=A_cI16_206_c --params=4.11971,0.1001
    if(vlabel.at(ifound)=="A_cI16_206_c") {
      PrototypeANRL_A_cI16_206_c(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 238 // ./aflow --proto=A_cP20_213_cd --params=6.315,0.06361,0.20224
    if(vlabel.at(ifound)=="A_cP20_213_cd") {
      PrototypeANRL_A_cP20_213_cd(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 239 // ./aflow --proto=A3B4C_cP8_215_d_e_a --params=5.3912,0.2372
    if(vlabel.at(ifound)=="A3B4C_cP8_215_d_e_a") {
      PrototypeANRL_A3B4C_cP8_215_d_e_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 240 // ./aflow --proto=AB4_cP5_215_a_e --params=3.878,0.265
    if(vlabel.at(ifound)=="AB4_cP5_215_a_e") {
      PrototypeANRL_AB4_cP5_215_a_e(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 241 // ./aflow --proto=AB3C4_cP8_215_a_c_e --params=5.28,0.25
    if(vlabel.at(ifound)=="AB3C4_cP8_215_a_c_e") {
      PrototypeANRL_AB3C4_cP8_215_a_c_e(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 242 // ./aflow --proto=AB5_cF24_216_a_ce --params=6.1,0.625
    if(vlabel.at(ifound)=="AB5_cF24_216_a_ce") {
      PrototypeANRL_AB5_cF24_216_a_ce(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 243 // ./aflow --proto=ABC_cF12_216_b_c_a --params=6.24
    if(vlabel.at(ifound)=="ABC_cF12_216_b_c_a") {
      PrototypeANRL_ABC_cF12_216_b_c_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 244 // ./aflow --proto=AB_cF8_216_c_a --params=5.4093
    if(vlabel.at(ifound)=="AB_cF8_216_c_a") {
      PrototypeANRL_AB_cF8_216_c_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 245 // ./aflow --proto=A4B_cI10_217_c_a --params=5.45858,0.165
    if(vlabel.at(ifound)=="A4B_cI10_217_c_a") {
      PrototypeANRL_A4B_cI10_217_c_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 246 // ./aflow --proto=A_cI58_217_ac2g --params=8.911,0.31787,-0.08958,0.28194,0.64294,0.03457
    if(vlabel.at(ifound)=="A_cI58_217_ac2g") {
      PrototypeANRL_A_cI58_217_ac2g(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 247 // ./aflow --proto=A5B8_cI52_217_ce_cg --params=8.8664,0.32774,0.10781,0.64421,0.68844,0.03674
    if(vlabel.at(ifound)=="A5B8_cI52_217_ce_cg") {
      PrototypeANRL_A5B8_cI52_217_ce_cg(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 248 // ./aflow --proto=A_cI16_220_c --params=5.2716,0.049
    if(vlabel.at(ifound)=="A_cI16_220_c") {
      PrototypeANRL_A_cI16_220_c(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 249 // ./aflow --proto=A3B2_cI40_220_d_c --params=8.135,0.0492,0.2896
    if(vlabel.at(ifound)=="A3B2_cI40_220_d_c") {
      PrototypeANRL_A3B2_cI40_220_d_c(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 250 // ./aflow --proto=AB_cP2_221_b_a --params=4.07925
    if(vlabel.at(ifound)=="AB_cP2_221_b_a") {
      PrototypeANRL_AB_cP2_221_b_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 251 // ./aflow --proto=AB_cP6_221_c_d --params=4.2101
    if(vlabel.at(ifound)=="AB_cP6_221_c_d") {
      PrototypeANRL_AB_cP6_221_c_d(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 252 // ./aflow --proto=AB3C_cP5_221_a_c_b --params=3.795
    if(vlabel.at(ifound)=="AB3C_cP5_221_a_c_b") {
      PrototypeANRL_AB3C_cP5_221_a_c_b(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 253 // ./aflow --proto=AB27CD3_cP32_221_a_dij_b_c --params=7.04,0.245,0.26
    if(vlabel.at(ifound)=="AB27CD3_cP32_221_a_dij_b_c") {
      PrototypeANRL_AB27CD3_cP32_221_a_dij_b_c(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 254 // ./aflow --proto=AB3_cP4_221_a_c --params=3.7402
    if(vlabel.at(ifound)=="AB3_cP4_221_a_c") {
      PrototypeANRL_AB3_cP4_221_a_c(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 255 // ./aflow --proto=A_cP1_221_a --params=3.34
    if(vlabel.at(ifound)=="A_cP1_221_a") {
      PrototypeANRL_A_cP1_221_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 256 // ./aflow --proto=AB11_cP36_221_c_agij --params=9.6,0.345,0.225,0.115
    if(vlabel.at(ifound)=="AB11_cP36_221_c_agij") {
      PrototypeANRL_AB11_cP36_221_c_agij(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 257 // ./aflow --proto=AB11CD3_cP16_221_a_dg_b_c --params=5.74,0.245
    if(vlabel.at(ifound)=="AB11CD3_cP16_221_a_dg_b_c") {
      PrototypeANRL_AB11CD3_cP16_221_a_dg_b_c(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 258 // ./aflow --proto=A3B_cP4_221_d_a --params=3.734
    if(vlabel.at(ifound)=="A3B_cP4_221_d_a") {
      PrototypeANRL_A3B_cP4_221_d_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 259 // ./aflow --proto=A6B_cP7_221_f_a --params=4.145,0.2117
    if(vlabel.at(ifound)=="A6B_cP7_221_f_a") {
      PrototypeANRL_A6B_cP7_221_f_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 260 // ./aflow --proto=A3B_cP8_223_c_a --params=4.556
    if(vlabel.at(ifound)=="A3B_cP8_223_c_a") {
      PrototypeANRL_A3B_cP8_223_c_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 261 // ./aflow --proto=A_cP46_223_dik --params=10.355,0.1837,0.1172,0.3077
    if(vlabel.at(ifound)=="A_cP46_223_dik") {
      PrototypeANRL_A_cP46_223_dik(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 262 // ./aflow --proto=A2B_cP6_224_b_a --params=4.267
    if(vlabel.at(ifound)=="A2B_cP6_224_b_a") {
      PrototypeANRL_A2B_cP6_224_b_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 263 // ./aflow --proto=A7B_cF32_225_bd_a --params=9.45
    if(vlabel.at(ifound)=="A7B_cF32_225_bd_a") {
      PrototypeANRL_A7B_cF32_225_bd_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 264 // ./aflow --proto=AB3_cF16_225_a_bc --params=5.853
    if(vlabel.at(ifound)=="AB3_cF16_225_a_bc") {
      PrototypeANRL_AB3_cF16_225_a_bc(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 265 // ./aflow --proto=A9B16C7_cF128_225_acd_2f_be --params=11.48,0.25,0.875,0.625
    if(vlabel.at(ifound)=="A9B16C7_cF128_225_acd_2f_be") {
      PrototypeANRL_A9B16C7_cF128_225_acd_2f_be(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 266 // ./aflow --proto=A12B_cF52_225_i_a --params=7.468,0.666
    if(vlabel.at(ifound)=="A12B_cF52_225_i_a") {
      PrototypeANRL_A12B_cF52_225_i_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 267 // ./aflow --proto=AB2_cF12_225_a_c --params=5.4631
    if(vlabel.at(ifound)=="AB2_cF12_225_a_c") {
      PrototypeANRL_AB2_cF12_225_a_c(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 268 // ./aflow --proto=A6B23_cF116_225_e_acfh --params=10.65,0.2765,0.6191,0.6699
    if(vlabel.at(ifound)=="A6B23_cF116_225_e_acfh") {
      PrototypeANRL_A6B23_cF116_225_e_acfh(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 269 // ./aflow --proto=AB2C_cF16_225_a_c_b --params=5.95
    if(vlabel.at(ifound)=="AB2C_cF16_225_a_c_b") {
      PrototypeANRL_AB2C_cF16_225_a_c_b(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 270 // ./aflow --proto=A_cF4_225_a --params=3.61491
    if(vlabel.at(ifound)=="A_cF4_225_a") {
      PrototypeANRL_A_cF4_225_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 271 // ./aflow --proto=AB18C8_cF108_225_a_eh_f --params=10.56,0.325,0.65833,0.66
    if(vlabel.at(ifound)=="AB18C8_cF108_225_a_eh_f") {
      PrototypeANRL_AB18C8_cF108_225_a_eh_f(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 272 // ./aflow --proto=AB_cF8_225_a_b --params=5.63931
    if(vlabel.at(ifound)=="AB_cF8_225_a_b") {
      PrototypeANRL_AB_cF8_225_a_b(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 273 // ./aflow --proto=A2B_cF24_227_c_a --params=7.166
    if(vlabel.at(ifound)=="A2B_cF24_227_c_a") {
      PrototypeANRL_A2B_cF24_227_c_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 274 // ./aflow --proto=AB2_cF96_227_e_cf --params=11.278,0.215,0.44
    if(vlabel.at(ifound)=="AB2_cF96_227_e_cf") {
      PrototypeANRL_AB2_cF96_227_e_cf(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 275 // ./aflow --proto=AB_cF16_227_a_b --params=7.483
    if(vlabel.at(ifound)=="AB_cF16_227_a_b") {
      PrototypeANRL_AB_cF16_227_a_b(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 276 // ./aflow --proto=A_cF136_227_aeg --params=14.864,0.2624,0.1824,0.3701
    if(vlabel.at(ifound)=="A_cF136_227_aeg") {
      PrototypeANRL_A_cF136_227_aeg(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 277 // ./aflow --proto=A2B_cF24_227_d_a --params=7.02
    if(vlabel.at(ifound)=="A2B_cF24_227_d_a") {
      PrototypeANRL_A2B_cF24_227_d_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 278 // ./aflow --proto=A_cF8_227_a --params=3.55
    if(vlabel.at(ifound)=="A_cF8_227_a") {
      PrototypeANRL_A_cF8_227_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 279 // ./aflow --proto=A2BC4_cF56_227_d_a_e --params=8.0832,0.7376
    if(vlabel.at(ifound)=="A2BC4_cF56_227_d_a_e") {
      PrototypeANRL_A2BC4_cF56_227_d_a_e(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 280 // ./aflow --proto=AB2_cF48_227_c_e --params=8.6,0.245
    if(vlabel.at(ifound)=="AB2_cF48_227_c_e") {
      PrototypeANRL_AB2_cF48_227_c_e(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 281 // ./aflow --proto=AB3C3_cF112_227_c_de_f --params=11.087,0.7047,0.323
    if(vlabel.at(ifound)=="AB3C3_cF112_227_c_de_f") {
      PrototypeANRL_AB3C3_cF112_227_c_de_f(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 282 // ./aflow --proto=A_cI2_229_a --params=3.155
    if(vlabel.at(ifound)=="A_cI2_229_a") {
      PrototypeANRL_A_cI2_229_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 283 // ./aflow --proto=A3B_cI8_229_b_a --params=2.984
    if(vlabel.at(ifound)=="A3B_cI8_229_b_a") {
      PrototypeANRL_A3B_cI8_229_b_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 284 // ./aflow --proto=A4B3_cI14_229_c_b --params=6.226
    if(vlabel.at(ifound)=="A4B3_cI14_229_c_b") {
      PrototypeANRL_A4B3_cI14_229_c_b(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 285 // ./aflow --proto=A2B7_cI54_229_e_afh --params=11.618,0.6862,0.1704,0.6503
    if(vlabel.at(ifound)=="A2B7_cI54_229_e_afh") {
      PrototypeANRL_A2B7_cI54_229_e_afh(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 286 // ./aflow --proto=AB12C3_cI32_229_a_h_b --params=7.04,0.7625
    if(vlabel.at(ifound)=="AB12C3_cI32_229_a_h_b") {
      PrototypeANRL_AB12C3_cI32_229_a_h_b(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 287 // ./aflow --proto=AB4C3_cI16_229_a_c_b --params=5.74
    if(vlabel.at(ifound)=="AB4C3_cI16_229_a_c_b") {
      PrototypeANRL_AB4C3_cI16_229_a_c_b(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 288 // ./aflow --proto=A4B3_cI112_230_af_g --params=11.411,0.0,0.625
    if(vlabel.at(ifound)=="A4B3_cI112_230_af_g") {
      PrototypeANRL_A4B3_cI112_230_af_g(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    }
    // -------------------------------------------------------------------------
    // Part 2
    // -------------------------------------------------------------------------
    // 1 // ./aflow --proto=A2B_aP6_2_aei_i --params=2.7804,1.00438785786,1.53100273342,78.28,76.53,70.42,0.259,0.415,0.663,0.192,0.183,0.255
    if(vlabel.at(ifound)=="A2B_aP6_2_aei_i")
      PrototypeANRL_A2B_aP6_2_aei_i(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 2 // ./aflow --proto=A8B5_mP13_6_a7b_3a2b --params=6.5369054321,0.490874879533,1.43786810261,109.592,0.5119,0.7172,0.4546,0.4285,0.044,0.5911,0.37,0.0053,0.3887,0.2195,0.6277,0.0001,0.1818,0.4663,0.7456,0.5256,0.1826,0.7901,0.7824,0.2724,0.0,0.0,0.8188,0.8026,0.0123,0.2051
    if(vlabel.at(ifound)=="A8B5_mP13_6_a7b_3a2b")
      PrototypeANRL_A8B5_mP13_6_a7b_3a2b(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 3 // ./aflow --proto=AB_mP4_6_2b_2a --params=3.580975428,1.00027925162,1.00167550965,90.04,0.0,0.0,0.518,0.507,0.026,0.501,0.529,0.027
    if(vlabel.at(ifound)=="AB_mP4_6_2b_2a")
      PrototypeANRL_AB_mP4_6_2b_2a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 4 // ./aflow --proto=A2B_mP12_7_4a_2a --params=5.0942,0.627360527659,1.04603274312,90.38,0.498,-0.062,0.472,-0.023,0.574,0.143,0.777,-0.052,0.799,0.271,0.151,0.261,-0.001,0.808,0.36,0.494,0.649,0.658
    if(vlabel.at(ifound)=="A2B_mP12_7_4a_2a")
      PrototypeANRL_A2B_mP12_7_4a_2a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 5 // ./aflow --proto=A2B_mP18_7_6a_3a --params=6.5499839723,1.91328244277,1.22717557252,127.75,0.6753,0.4828,0.3393,0.4711,0.6498,0.3067,0.4174,0.3433,0.0893,0.2883,0.1551,0.3432,0.1472,0.1112,0.5582,0.0,0.0746,0.0,0.5693,0.07852,0.0933,0.1184,0.41567,0.2938,0.8473,0.27155,0.6656
    if(vlabel.at(ifound)=="A2B_mP18_7_6a_3a")
      PrototypeANRL_A2B_mP18_7_6a_3a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 6 // ./aflow --proto=A3B_mP16_7_6a_2a --params=5.2778002048,0.97690325515,1.45210125431,91.762,0.5044,0.292,0.01,0.5764,0.215,0.586,0.0,0.209,0.0,0.0864,0.29,0.58,0.2874,0.0717,0.287,0.7924,0.4201,0.301,0.2874,0.014,0.0012,0.7994,0.528,0.078
    if(vlabel.at(ifound)=="A3B_mP16_7_6a_2a")
      PrototypeANRL_A3B_mP16_7_6a_2a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 7 // ./aflow --proto=A9B2_mP22_7_9a_2a --params=6.4165085102,0.999298672153,1.36910105356,93.39,0.4944,0.2369,0.0135,0.2965,0.2623,0.5568,0.4955,0.4453,0.2912,0.2817,0.0557,0.2805,0.8067,0.5491,0.0506,0.0,0.0309,0.0,0.6774,0.1539,0.7442,0.096,0.6343,0.3258,0.8638,0.2389,0.2526,0.6359,0.1217,0.4513,0.156,0.3761,0.11377
    if(vlabel.at(ifound)=="A9B2_mP22_7_9a_2a")
      PrototypeANRL_A9B2_mP22_7_9a_2a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 8 // ./aflow --proto=A5B3_mC32_9_5a_3a --params=8.12077,0.718445418353,1.12797801194,115.809,0.009,-0.003,0.269,0.129,0.341,0.45,0.37,0.119,0.066,0.142,0.351,0.147,0.356,0.135,0.348,0.0,0.5182,0.0,0.136,0.2,0.309,0.365,0.2924,0.196
    if(vlabel.at(ifound)=="A5B3_mC32_9_5a_3a")
      PrototypeANRL_A5B3_mC32_9_5a_3a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 9 // ./aflow --proto=AB3_mC16_9_a_3a --params=3.5367,2.67777872028,0.986823875364,93.018,0.50691,0.35551,0.49997,0.28659,0.26502,0.27913,0.57076,0.06256,0.41996,0.42173,0.07037,0.56036
    if(vlabel.at(ifound)=="AB3_mC16_9_a_3a")
      PrototypeANRL_AB3_mC16_9_a_3a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 10 // ./aflow --proto=A2B_mP6_10_mn_bg --params=4.012,0.819541375872,2.93668993021,97.03,0.843,0.126,0.558,0.644
    if(vlabel.at(ifound)=="A2B_mP6_10_mn_bg")
      PrototypeANRL_A2B_mP6_10_mn_bg(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 11 // ./aflow --proto=AB3_mP16_10_mn_3m3n --params=3.678,0.71098423056,1.23817292007,91.0,0.263,0.339,0.08,0.059,0.795,0.341,0.438,-0.069,0.763,0.161,0.705,0.841,0.58,0.441,0.062,0.431
    if(vlabel.at(ifound)=="AB3_mP16_10_mn_3m3n")
      PrototypeANRL_AB3_mP16_10_mn_3m3n(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 12 // ./aflow --proto=ABC2_mP8_10_ac_eh_mn --params=5.1177434753,0.86107854631,1.45123094959,90.021,0.6089,0.24179,0.1277,0.24913
    if(vlabel.at(ifound)=="ABC2_mP8_10_ac_eh_mn")
      PrototypeANRL_ABC2_mP8_10_ac_eh_mn(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 13 // ./aflow --proto=AB_mP6_10_en_am --params=5.1700416367,0.61508704062,1.49709864605,104.5,0.234,0.66,0.263,0.336
    if(vlabel.at(ifound)=="AB_mP6_10_en_am")
      PrototypeANRL_AB_mP6_10_en_am(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 14 // ./aflow --proto=A_mP8_10_2m2n --params=4.7302,0.527461840937,0.86332501797,106.1,0.1175,0.6746,0.5344,0.3333,0.1131,0.8977,0.4209,0.1319
    if(vlabel.at(ifound)=="A_mP8_10_2m2n")
      PrototypeANRL_A_mP8_10_2m2n(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 15 // ./aflow --proto=A7B2C2_mC22_12_aij_h_i --params=6.65,1.29563909774,0.704661654135,102.2,0.30503,0.38654,0.22171,0.22108,-0.08762,0.23655,0.15499,0.71826
    if(vlabel.at(ifound)=="A7B2C2_mC22_12_aij_h_i")
      PrototypeANRL_A7B2C2_mC22_12_aij_h_i(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 16 // ./aflow --proto=A_mC16_12_4i --params=9.089,0.274617669711,0.451534822313,96.96,-0.0572,0.1206,0.4419,0.3467,0.7858,-0.0594,0.2715,0.4149
    if(vlabel.at(ifound)=="A_mC16_12_4i")
      PrototypeANRL_A_mC16_12_4i(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 17 // ./aflow --proto=A2B_mP12_13_2g_ef --params=5.6255,0.61198115723,1.23780997245,127.44,0.1808,-0.004,0.155,0.346,0.225,0.345,0.273,0.573
    if(vlabel.at(ifound)=="A2B_mP12_13_2g_ef")
      PrototypeANRL_A2B_mP12_13_2g_ef(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 18 // ./aflow --proto=A2B_mP6_14_e_a --params=5.5496,0.695689779444,1.15512829753,107.151,0.255,0.2573,0.3141
    if(vlabel.at(ifound)=="A2B_mP6_14_e_a")
      PrototypeANRL_A2B_mP6_14_e_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 19 // ./aflow --proto=A7B8_mP120_14_14e_16e --params=7.5889,0.766725085322,3.55545599494,106.136,0.25252,0.5751,0.3402,0.76924,0.888,0.47034,0.13489,0.3866,0.33222,0.8218,0.8146,0.42747,0.12506,0.2343,0.29191,0.78423,0.9464,0.38285,0.23276,0.2673,0.25885,0.69189,0.1533,0.38027,0.34844,0.4547,0.26591,0.63874,0.2281,0.42248,0.35819,0.6069,0.3061,0.67752,0.0969,0.46714,0.26777,0.7386,0.3843,0.8119,0.7461,0.51889,0.0602,0.3617,0.3547,0.8843,0.6723,0.4287,0.0438,0.1068,0.2871,0.822,0.8945,0.354,0.2273,0.1621,0.2315,0.6653,0.243,0.3497,0.4219,0.4796,0.2431,0.5754,0.3699,0.4209,0.4385,0.7352,0.3104,0.6409,0.1506,0.496,0.9409,0.7672,0.5381,0.7891,0.5835,0.5099,0.7334,0.7952,0.5403,0.2081,0.8842,0.3711,0.3975,0.7664,0.4019,0.2077,0.6717,0.4087
    if(vlabel.at(ifound)=="A7B8_mP120_14_14e_16e")
      PrototypeANRL_A7B8_mP120_14_14e_16e(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 20 // ./aflow --proto=AB3_mC16_15_e_cf --params=3.282,2.64137720902,0.970749542962,91.9,0.86,0.577,0.068,0.169
    if(vlabel.at(ifound)=="AB3_mC16_15_e_cf")
      PrototypeANRL_AB3_mC16_15_e_cf(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 21 // ./aflow --proto=A_mC24_15_2e2f --params=4.939,0.569143551326,0.838023891476,142.47,0.1012,0.3684,0.226,0.0672,0.2464,0.3443,0.1958,0.2227
    if(vlabel.at(ifound)=="A_mC24_15_2e2f")
      PrototypeANRL_A_mC24_15_2e2f(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 22 // ./aflow --proto=A2B_oP12_17_abe_e --params=7.0499703347,1.11347517732,0.614184397158,0.893,0.878,0.379,0.225,0.522,0.202,0.275,0.022
    if(vlabel.at(ifound)=="A2B_oP12_17_abe_e")
      PrototypeANRL_A2B_oP12_17_abe_e(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 23 // ./aflow --proto=AB3_oP16_19_a_3a --params=5.668,0.758997882851,0.528757939308,0.584,0.123,0.027,0.31,0.159,0.417,0.257,0.073,0.603,0.983,0.124,0.227
    if(vlabel.at(ifound)=="AB3_oP16_19_a_3a")
      PrototypeANRL_AB3_oP16_19_a_3a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 24 // ./aflow --proto=AB2_oC6_21_a_k --params=3.3982513,1.3943496174,1.40170688639,0.268
    if(vlabel.at(ifound)=="AB2_oC6_21_a_k")
      PrototypeANRL_AB2_oC6_21_a_k(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 25 // ./aflow --proto=A2BC2_oF40_22_fi_ad_gh --params=6.4793068924,1.3977127159,1.54737394473,0.3134,0.3627,0.1144,0.0719
    if(vlabel.at(ifound)=="A2BC2_oF40_22_fi_ad_gh")
      PrototypeANRL_A2BC2_oF40_22_fi_ad_gh(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 26 // ./aflow --proto=AB_oF8_22_a_c --params=5.5400291632,0.990433212996,0.937725631773
    if(vlabel.at(ifound)=="AB_oF8_22_a_c")
      PrototypeANRL_AB_oF8_22_a_c(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 27 // ./aflow --proto=A3B_oI32_23_ij2k_k --params=5.82463,1.24369101557,1.32254065924,0.04851,0.45153,0.75475,0.49255,0.20405,0.4421,0.23223,0.28616,0.76005,0.8216,0.36488
    if(vlabel.at(ifound)=="A3B_oI32_23_ij2k_k")
      PrototypeANRL_A3B_oI32_23_ij2k_k(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 28 // ./aflow --proto=A8B2C12D2E_oI50_23_bcfk_i_3k_j_a --params=10.767,0.502554100492,1.4969815176,0.2511,0.3298,0.1693,0.2465,0.0107,0.1695,0.1308,0.2443,0.0826,0.3792,0.7558,0.0801,0.1294,0.7488,0.2546
    if(vlabel.at(ifound)=="A8B2C12D2E_oI50_23_bcfk_i_3k_j_a")
      PrototypeANRL_A8B2C12D2E_oI50_23_bcfk_i_3k_j_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 29 // ./aflow --proto=ABC2_oI16_23_ab_i_k --params=5.3999384419,1.15740740741,2.00555555555,0.28,0.25,0.2,0.115
    if(vlabel.at(ifound)=="ABC2_oI16_23_ab_i_k")
      PrototypeANRL_ABC2_oI16_23_ab_i_k(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 30 // ./aflow --proto=ABC4_oI12_23_a_b_k --params=5.2501580231,1.06666666665,1.7219047619,0.21,0.2,0.115
    if(vlabel.at(ifound)=="ABC4_oI12_23_a_b_k")
      PrototypeANRL_ABC4_oI12_23_a_b_k(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 31 // ./aflow --proto=AB7CD2_oI44_24_a_b3d_c_ac --params=7.0501914381,1.035035461,1.41546099291,0.7511,0.2496,0.1361,-0.0003,0.5002,0.7501,0.2213,0.3356,0.5662,0.0681,0.1362,-0.0648,0.0709,0.1385
    if(vlabel.at(ifound)=="AB7CD2_oI44_24_a_b3d_c_ac")
      PrototypeANRL_AB7CD2_oI44_24_a_b3d_c_ac(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 32 // ./aflow --proto=A2B_oP12_26_abc_ab --params=4.6806,0.627034995513,1.05710806307,0.455,0.858,0.179,0.623,0.048,0.545,0.375,0.355,0.751,0.119,0.213
    // 33 // ./aflow --proto=A2B_oP12_26_abc_ab --params=5.0725824349,0.881353258932,1.48474034935,0.12,0.0,0.2484,0.461,0.246,0.289,0.3781,0.0862,0.253,0.652,0.12
    if(vlabel.at(ifound)=="A2B_oP12_26_abc_ab")
      PrototypeANRL_A2B_oP12_26_abc_ab(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 34 // ./aflow --proto=A5B_oP24_26_3a3b2c_ab --params=6.465013742,1.07113689096,1.87440061873,0.1628,0.2907,0.112,0.0,0.1452,0.8256,0.4859,0.1227,0.1146,0.0096,0.3286,0.1358,0.1565,0.2905,0.7214,0.2818,0.2489,0.2854,0.3947,0.249,0.0963,0.5436
    if(vlabel.at(ifound)=="A5B_oP24_26_3a3b2c_ab")
      PrototypeANRL_A5B_oP24_26_3a3b2c_ab(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 35 // ./aflow --proto=A6B4C16D_oP108_27_abcd4e_4e_16e_e --params=13.0409854622,0.999877309088,0.703110981603,0.0,0.005,0.03,0.05,0.373,0.379,0.296,0.122,0.125,0.249,0.369,0.129,0.28,0.122,0.376,0.26,0.058,0.751,0.449,0.237,0.574,0.086,0.257,0.012,0.554,0.483,0.247,0.028,0.729,0.183,0.412,0.162,0.684,0.26,0.23,0.164,0.679,0.651,0.299,0.6217,0.4,0.254,0.238,0.433,0.091,0.44,0.445,0.395,0.455,0.245,0.09,0.304,0.404,0.057,0.126,0.093,0.054,0.102,0.243,0.402,0.336,0.099,0.466,0.124,0.392,0.467,0.154,0.103,0.253,0.191,0.047,0.385,0.425,0.052,0.111,0.41,0.74634,0.2513,0.28218
    if(vlabel.at(ifound)=="A6B4C16D_oP108_27_abcd4e_4e_16e_e")
      PrototypeANRL_A6B4C16D_oP108_27_abcd4e_4e_16e_e(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 36 // ./aflow --proto=A2B_oP12_29_2a_a --params=5.2594682584,0.963498098863,0.965209125484,0.639,0.068,0.0,0.771,0.537,0.106,0.53,0.267,0.356
    if(vlabel.at(ifound)=="A2B_oP12_29_2a_a")
      PrototypeANRL_A2B_oP12_29_2a_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 37 // ./aflow --proto=AB2_oP12_29_a_2a --params=5.4179557747,1.0,1.0,0.5049,0.2419,0.0,0.615,0.135,0.3834,0.615,0.635,0.1134
    if(vlabel.at(ifound)=="AB2_oP12_29_a_2a")
      PrototypeANRL_AB2_oP12_29_a_2a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 38 // ./aflow --proto=ABC_oP12_29_a_a_a --params=5.2594682584,0.963498098863,0.965209125484,0.61885,0.63065,0.11668,0.50496,0.24091,0.0,0.61734,0.13129,0.37996
    if(vlabel.at(ifound)=="ABC_oP12_29_a_a_a")
      PrototypeANRL_ABC_oP12_29_a_a_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 39 // ./aflow --proto=A5B3C15_oP46_30_a2c_bc_a7c --params=5.4630061801,1.00183049606,3.84605528098,0.5,0.107,0.523,0.2442,0.513,0.019,0.3675,0.026,0.016,0.1074,0.003,0.02,0.1949,0.05,0.074,0.097,0.687,0.22,0.0844,0.254,0.226,0.413,0.52,0.105,0.505,0.718,0.31,0.305,0.763,0.271,0.695,0.752,0.289
    if(vlabel.at(ifound)=="A5B3C15_oP46_30_a2c_bc_a7c")
      PrototypeANRL_A5B3C15_oP46_30_a2c_bc_a7c(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 40 // ./aflow --proto=ABC3_oP20_30_2a_c_3c --params=4.480982758,1.71211783083,3.19772372237,0.8543,0.5,0.3011,0.721,0.4096,0.3607,0.7305,0.1733,0.5847,0.8629,0.0496,0.6302,0.8523,0.2991
    if(vlabel.at(ifound)=="ABC3_oP20_30_2a_c_3c")
      PrototypeANRL_ABC3_oP20_30_2a_c_3c(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 41 // ./aflow --proto=A13B2C2_oP34_32_a6c_c_c --params=8.4261784988,1.1450866366,0.701103726564,0.0,0.808,0.55,0.8946,0.614,0.7595,-0.0255,0.8197,0.8044,0.642,0.5367,0.6481,0.3894,0.7689,-0.0703,0.2917,0.8392,0.688,0.2793,0.67666,0.60773,0.0844,0.8643,0.8235,0.4147
    if(vlabel.at(ifound)=="A13B2C2_oP34_32_a6c_c_c")
      PrototypeANRL_A13B2C2_oP34_32_a6c_c_c(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 42 // ./aflow --proto=A2B3_oP40_33_4a_6a --params=4.8437,1.71975968784,1.84873134174,0.6787,0.8416,0.0,0.1846,0.3432,0.7868,0.8115,0.6489,0.6972,0.6677,0.4696,0.9993,0.329,0.8313,0.8927,0.0248,0.4908,0.6292,0.4717,0.6647,0.6381,0.5145,0.6728,0.1212,0.8608,0.3301,0.8662,0.336,0.4992,0.9
    if(vlabel.at(ifound)=="A2B3_oP40_33_4a_6a")
      PrototypeANRL_A2B3_oP40_33_4a_6a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 43 // ./aflow --proto=A2B8C_oP22_34_c_4c_a --params=6.2740008991,2.04032515143,1.38284985656,0.5,0.605,0.8135,0.499,0.7349,0.6796,0.009,0.743,-0.0925,0.3064,0.2388,0.5935,0.2196,0.7131,0.6459,0.513
    if(vlabel.at(ifound)=="A2B8C_oP22_34_c_4c_a")
      PrototypeANRL_A2B8C_oP22_34_c_4c_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 44 // ./aflow --proto=AB2_oP6_34_a_c --params=5.8327827022,1.12083390481,0.548158688784,0.5,0.6881,0.8565,0.0097
    if(vlabel.at(ifound)=="AB2_oP6_34_a_c")
      PrototypeANRL_AB2_oP6_34_a_c(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 45 // ./aflow --proto=AB8C2_oC22_35_a_ab3e_e --params=4.534758023,1.1408839779,2.72580110498,0.0,0.5838,-0.042,0.189,0.5461,0.0961,0.122,0.2982,0.1399,0.1866,0.0038
    if(vlabel.at(ifound)=="AB8C2_oC22_35_a_ab3e_e")
      PrototypeANRL_AB8C2_oC22_35_a_ab3e_e(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 46 // ./aflow --proto=AB_oC8_36_a_a --params=5.825,0.945115879828,0.922403433476,0.25,0.0,0.081,0.83
    if(vlabel.at(ifound)=="AB_oC8_36_a_a")
      PrototypeANRL_AB_oC8_36_a_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 47 // ./aflow --proto=A2B5C2_oC36_37_d_c2d_d --params=5.8071602545,2.51110728429,0.821939039086,0.0,0.654,0.0584,0.0469,0.6705,0.0718,0.471,0.0932,0.1377,0.4004,0.1552,0.14836,0.0571
    if(vlabel.at(ifound)=="A2B5C2_oC36_37_d_c2d_d")
      PrototypeANRL_A2B5C2_oC36_37_d_c2d_d(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 48 // ./aflow --proto=A2B3_oC40_39_2d_2c2d --params=7.47831,2.30292940517,0.749515599113,0.3449,0.5,0.0058,0.7654,0.173,0.5398,0.3889,0.5843,0.6154,0.3917,0.28104,0.16619,0.0206,0.10455,0.6073,0.0144
    if(vlabel.at(ifound)=="A2B3_oC40_39_2d_2c2d")
      PrototypeANRL_A2B3_oC40_39_2d_2c2d(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 49 // ./aflow --proto=A9BC_oC44_39_3c3d_a_c --params=6.2123847492,1.90215711527,2.62509658728,0.5,0.2673,0.7493,0.6176,0.8879,0.6246,0.6059,0.6191,0.7476,0.2266,0.0857,0.2459,0.1789,0.5932,0.0649,0.18,0.5948,0.4281
    if(vlabel.at(ifound)=="A9BC_oC44_39_3c3d_a_c")
      PrototypeANRL_A9BC_oC44_39_3c3d_a_c(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 50 // ./aflow --proto=AB2C_oC16_40_a_2b_b --params=6.4580033647,1.33199132859,1.69634561784,0.0,0.3467,0.1773,0.8661,0.3301,0.7281,0.0093
    if(vlabel.at(ifound)=="AB2C_oC16_40_a_2b_b")
      PrototypeANRL_AB2C_oC16_40_a_2b_b(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 51 // ./aflow --proto=AB3_oC16_40_b_3b --params=4.3850022361,5.92335058951,0.997331752147,0.83109,0.0,0.70417,0.0021,0.56971,0.4966,-0.07002,0.4978
    if(vlabel.at(ifound)=="AB3_oC16_40_b_3b")
      PrototypeANRL_AB3_oC16_40_b_3b(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 52 // ./aflow --proto=A10B3_oF52_42_2abce_ab --params=7.4494846573,1.70036689767,1.04688136975,0.76,0.27,0.0,0.31,0.06,0.79,0.6,0.17,0.11,0.07
    if(vlabel.at(ifound)=="A10B3_oF52_42_2abce_ab")
      PrototypeANRL_A10B3_oF52_42_2abce_ab(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 53 // ./aflow --proto=AB_oF8_42_a_a --params=2.5000573158,1.33999999997,1.73599999999,0.0,0.333
    if(vlabel.at(ifound)=="AB_oF8_42_a_a")
      PrototypeANRL_AB_oF8_42_a_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 54 // ./aflow --proto=A2BC2_oI20_45_c_b_c --params=5.9678340183,1.97721179624,0.981568364609,0.25,0.2729,0.629,0.493,0.2296,0.8629,0.474
    if(vlabel.at(ifound)=="A2BC2_oI20_45_c_b_c")
      PrototypeANRL_A2BC2_oI20_45_c_b_c(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 55 // ./aflow --proto=ABC_oI36_46_ac_bc_3b --params=6.9969353511,1.54780620265,0.898527940538,0.25,0.5253,0.0054,0.7207,0.7706,-0.0021,-0.0823,0.2996,0.7963,0.5295,0.6236,0.1199,0.006,0.3325,0.4952
    if(vlabel.at(ifound)=="ABC_oI36_46_ac_bc_3b")
      PrototypeANRL_ABC_oI36_46_ac_bc_3b(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 56 // ./aflow --proto=A2B8CD_oP24_48_k_2m_d_b --params=6.3302356554,1.0,1.50710900473,0.4978,0.56,0.629,0.114,0.642,0.558,0.601
    if(vlabel.at(ifound)=="A2B8CD_oP24_48_k_2m_d_b")
      PrototypeANRL_A2B8CD_oP24_48_k_2m_d_b(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 57 // ./aflow --proto=A5B2_oP14_49_dehq_ab --params=3.6705354001,1.69539132804,2.12544314151,0.002,0.681
    // lowering tolerance value for aflowSYM to 0.001 resolves SG #49
    if(vlabel.at(ifound)=="A5B2_oP14_49_dehq_ab")
      PrototypeANRL_A5B2_oP14_49_dehq_ab(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 58 // ./aflow --proto=AB2C8D_oP24_49_g_q_2qr_e --params=6.3302356554,1.0,1.50710900473,0.5184,0.7996,0.258,-0.0697,0.3652,0.6307,0.2617,0.1943,0.8286
    if(vlabel.at(ifound)=="AB2C8D_oP24_49_g_q_2qr_e")
      PrototypeANRL_AB2C8D_oP24_49_g_q_2qr_e(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 59 // ./aflow --proto=A2BC4_oP28_50_ij_ac_ijm --params=5.5348069961,2.26684733515,0.987895212291,0.8726,0.445,0.3867,-0.076,0.5,0.743,0.226
    if(vlabel.at(ifound)=="A2BC4_oP28_50_ij_ac_ijm")
      PrototypeANRL_A2BC4_oP28_50_ij_ac_ijm(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 60 // ./aflow --proto=A3BC2_oP48_50_3m_m_2m --params=11.0940438033,1.50045069408,0.472480620154,0.0515,0.8152,0.238,0.602,0.5567,0.315,0.604,0.8447,0.116,0.61286,0.66266,0.2344,0.11964,0.66839,0.2479,0.63183,0.00412,0.2439
    if(vlabel.at(ifound)=="A3BC2_oP48_50_3m_m_2m")
      PrototypeANRL_A3BC2_oP48_50_3m_m_2m(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 61 // ./aflow --proto=A2B_oP24_52_2e_cd --params=7.2235008877,1.34576036547,1.32098013428,0.3175,0.6759,0.8271,0.1762,0.5576,0.5093,0.0419,0.8142
    if(vlabel.at(ifound)=="A2B_oP24_52_2e_cd")
      PrototypeANRL_A2B_oP24_52_2e_cd(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 62 // ./aflow --proto=A3B2_oP20_52_de_cd --params=15.5832022616,0.434585119311,0.426069989123,0.443,0.4294,0.0001,0.6539,0.064,0.0788
    if(vlabel.at(ifound)=="A3B2_oP20_52_de_cd")
      PrototypeANRL_A3B2_oP20_52_de_cd(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 63 // ./aflow --proto=ABC2_oP16_53_h_e_gh --params=11.6904618594,0.282132455361,0.78890717994,0.79514,0.3193,0.1198,0.3549,0.224,0.7493
    if(vlabel.at(ifound)=="ABC2_oP16_53_h_e_gh")
      PrototypeANRL_ABC2_oP16_53_h_e_gh(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 64 // ./aflow --proto=ABC3_oP20_53_e_g_hi --params=14.3630679002,0.312469539787,0.535821207264,0.6826,0.2856,0.3458,0.2708,0.6247,0.6057,0.3575
    if(vlabel.at(ifound)=="ABC3_oP20_53_e_g_hi")
      PrototypeANRL_ABC3_oP20_53_e_g_hi(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 65 // ./aflow --proto=ABC3_oP20_54_e_d_cf --params=5.3467489374,0.956627452436,1.81710213776,0.8667,0.6417,0.8902,0.4055,0.2686,0.5503
    if(vlabel.at(ifound)=="ABC3_oP20_54_e_d_cf")
      PrototypeANRL_ABC3_oP20_54_e_d_cf(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 66 // ./aflow --proto=A2B_oP24_55_2g2h_gh --params=10.1600191136,1.45275590551,0.366929133855,0.0628,0.4022,0.1014,0.1118,0.2024,0.2667,0.226,0.0384,0.3532,0.2953,0.4192,0.1378
    if(vlabel.at(ifound)=="A2B_oP24_55_2g2h_gh")
      PrototypeANRL_A2B_oP24_55_2g2h_gh(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 67 // ./aflow --proto=A3B5_oP16_55_ch_agh --params=5.4199981729,1.90405904061,0.730627306278,0.348,0.22,0.112,0.152,0.17,0.393
    if(vlabel.at(ifound)=="A3B5_oP16_55_ch_agh")
      PrototypeANRL_A3B5_oP16_55_ch_agh(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 68 // ./aflow --proto=A_oP16_55_2g2h --params=7.7886,0.613101199189,0.320442698303,0.6731,-0.037,0.8435,0.8087,-0.0454,0.8613,0.5704,0.8926
    if(vlabel.at(ifound)=="A_oP16_55_2g2h")
      PrototypeANRL_A_oP16_55_2g2h(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 69 // ./aflow --proto=A2B_oP6_58_g_a --params=3.7572,2.89952624295,0.89063664431,0.3326,0.63309
    if(vlabel.at(ifound)=="A2B_oP6_58_g_a")
      PrototypeANRL_A2B_oP6_58_g_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 70 // ./aflow --proto=ABC_oP6_59_a_b_a --params=3.301,1.14298697364,2.39612238716,0.32961,-0.04795,0.89243
    if(vlabel.at(ifound)=="ABC_oP6_59_a_b_a")
      PrototypeANRL_ABC_oP6_59_a_b_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 71 // ./aflow --proto=A2B3_oP20_60_d_cd --params=8.4598035458,0.706855791961,0.724586288418,0.547,0.394,0.75,0.033,0.348,0.611,0.396
    if(vlabel.at(ifound)=="A2B3_oP20_60_d_cd")
      PrototypeANRL_A2B3_oP20_60_d_cd(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 72 // ./aflow --proto=A3B_oP32_60_3d_d --params=7.3397836195,1.05524748967,1.03197678378,0.5016,0.7205,0.0322,0.2167,0.7591,0.2582,0.2197,0.5016,0.013,0.248,0.783,0.0291
    if(vlabel.at(ifound)=="A3B_oP32_60_3d_d")
      PrototypeANRL_A3B_oP32_60_3d_d(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 73 // ./aflow --proto=A7B8_oP120_60_7d_8d --params=13.85,0.824548736462,0.53357400722,0.382,0.3932,0.5783,0.2764,0.38,0.5411,0.2264,0.4631,0.4411,0.1289,0.4507,0.4069,0.0794,0.3554,0.4717,0.1276,0.2723,0.571,0.2251,0.2844,0.6054,0.2616,0.5329,0.393,0.094,0.5116,0.3344,0.0088,0.3466,0.4468,0.0917,0.2026,0.6185,0.2593,0.2232,0.6779,0.3972,0.4784,0.595,0.4192,0.3624,0.4724,0.4002,0.3492,0.6895
    if(vlabel.at(ifound)=="A7B8_oP120_60_7d_8d")
      PrototypeANRL_A7B8_oP120_60_7d_8d(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 74 // ./aflow --proto=AB_oP48_61_3c_3c --params=6.914,1.08128435059,1.38313566676,0.1297,0.5762,0.40803,0.1235,0.6328,0.54518,0.0057,0.4432,0.36289,0.2172,0.6275,0.346,0.2068,0.7225,0.5756,0.0095,0.4051,0.2704
    if(vlabel.at(ifound)=="AB_oP48_61_3c_3c")
      PrototypeANRL_AB_oP48_61_3c_3c(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 75 // ./aflow --proto=A2B3_oP20_62_2c_3c --params=5.5329,0.511305102207,2.07339731425,0.1008,0.2055,0.2432,-0.0464,0.0157,0.4015,0.1808,0.7737,0.8691,-0.0688
    if(vlabel.at(ifound)=="A2B3_oP20_62_2c_3c")
      PrototypeANRL_A2B3_oP20_62_2c_3c(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 76 // ./aflow --proto=A2B4C_oP28_62_ac_2cd_c --params=10.193,0.586382811734,0.466202295693,0.2774,-0.0085,0.0913,0.7657,0.4474,0.2215,0.094,0.4262,0.1628,0.0331,0.2777
    if(vlabel.at(ifound)=="A2B4C_oP28_62_ac_2cd_c")
      PrototypeANRL_A2B4C_oP28_62_ac_2cd_c(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    //DX IN PART 1 SECTION // 77 // ./aflow --proto=A2B_oP12_62_2c_c --params=3.875,1.64232258065,1.89496774194,0.004,0.758,0.24,0.07,0.24,0.39
    //DX IN PART 1 SECTION if(vlabel.at(ifound)=="A2B_oP12_62_2c_c")
    //DX IN PART 1 SECTION  PrototypeANRL_A2B_oP12_62_2c_c(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 78 // ./aflow --proto=A3B_oP16_62_cd_c --params=6.5982,1.11416750023,0.727789397108,0.011,0.415,0.369,0.555,0.174,0.053,0.856
    if(vlabel.at(ifound)=="A3B_oP16_62_cd_c")
      PrototypeANRL_A3B_oP16_62_cd_c(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 79 // ./aflow --proto=AB2C3_oP24_62_c_d_cd --params=6.231,1.78414379714,1.03787514043,0.123,0.0823,0.2579,0.4127,0.1366,0.087,0.5853,0.267,0.0846,-0.088
    if(vlabel.at(ifound)=="AB2C3_oP24_62_c_d_cd")
      PrototypeANRL_AB2C3_oP24_62_c_d_cd(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 80 // ./aflow --proto=AB3_oP16_62_c_3c --params=13.855,0.266791771923,0.28601948755,0.398,0.425,0.074,-0.026,0.414,-0.066,0.276,0.49
    if(vlabel.at(ifound)=="AB3_oP16_62_c_3c")
      PrototypeANRL_AB3_oP16_62_c_3c(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 81 // ./aflow --proto=AB4C_oP24_62_c_2cd_c --params=8.884,0.614362899595,0.805155335434,0.8154,0.3419,0.5878,0.6062,0.3192,0.5515,0.437,0.6914,0.4186,0.4702,0.819
    if(vlabel.at(ifound)=="AB4C_oP24_62_c_2cd_c")
      PrototypeANRL_AB4C_oP24_62_c_2cd_c(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    //DX IN PART 1 SECTION // 82 // ./aflow --proto=AB_oP8_62_c_c --params=5.454,0.609644297763,1.10542720939,0.2005,0.5741,0.0058,0.1993
    //DX IN PART 1 SECTION if(vlabel.at(ifound)=="AB_oP8_62_c_c")
    //DX IN PART 1 SECTION  PrototypeANRL_AB_oP8_62_c_c(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 83 // ./aflow --proto=A2BC3_oC24_63_e_c_cg --params=9.049,1.21770361366,0.600176815118,0.6699,0.1191,0.8502,0.2174,0.3859
    if(vlabel.at(ifound)=="A2BC3_oC24_63_e_c_cg")
      PrototypeANRL_A2BC3_oC24_63_e_c_cg(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 84 // ./aflow --proto=A43B5C17_oC260_63_c8fg6h_cfg_ce3f2h --params=10.11895,1.73974572461,4.16742843872,0.38998,0.0432,0.5987,0.1275,0.13992,0.17412,0.7278,0.20648,0.19723,0.62632,0.1388,0.02884,0.39065,0.0988,0.06788,0.55569,0.59876,0.1328,0.28079,0.54761,0.0678,0.6912,0.07797,0.09224,0.45337,0.1675,0.43146,0.0091,0.18607,0.20422,0.3338,0.3774,0.18828,0.32898,0.174,0.19344,0.2014,0.10123,0.30731,0.03499,0.20655,0.30495,0.49596,0.12681,0.19578,0.33018,0.026,0.31697,0.03483,0.04522,0.2812,0.17213,0.16828,0.2807,0.35936,0.09067
    if(vlabel.at(ifound)=="A43B5C17_oC260_63_c8fg6h_cfg_ce3f2h")
      PrototypeANRL_A43B5C17_oC260_63_c8fg6h_cfg_ce3f2h(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 85 // ./aflow --proto=A6B_oC28_63_efg_c --params=7.5551,0.860266574896,1.17435904224,0.45686,0.32602,0.13917,0.10039,0.31768,0.28622
    if(vlabel.at(ifound)=="A6B_oC28_63_efg_c")
      PrototypeANRL_A6B_oC28_63_efg_c(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 86 // ./aflow --proto=AB3C_oC20_63_a_cf_c --params=2.456,3.27442996743,2.48086319218,0.077,0.747,0.631,-0.064
    if(vlabel.at(ifound)=="AB3C_oC20_63_a_cf_c")
      PrototypeANRL_AB3C_oC20_63_a_cf_c(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 87 // ./aflow --proto=AB4C_oC24_63_a_fg_c --params=5.182,1.52315708221,1.25549980702,0.37,0.25,0.06,0.25,0.47
    if(vlabel.at(ifound)=="AB4C_oC24_63_a_fg_c")
      PrototypeANRL_AB4C_oC24_63_a_fg_c(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 88 // ./aflow --proto=AB4C_oC24_63_c_fg_c --params=6.995,0.892780557541,0.999714081487,0.6524,0.15556,0.7025,-0.0819,0.1699,0.0162
    if(vlabel.at(ifound)=="AB4C_oC24_63_c_fg_c")
      PrototypeANRL_AB4C_oC24_63_c_fg_c(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 89 // ./aflow --proto=A2B_oC24_64_2f_f --params=2.9986,1.44314013206,2.534682852,0.372,0.27,0.6,0.092,0.371,0.615
    if(vlabel.at(ifound)=="A2B_oC24_64_2f_f")
      PrototypeANRL_A2B_oC24_64_2f_f(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 90 // ./aflow --proto=A2B4C_oC28_66_l_kl_a --params=6.2700590867,1.72567783094,1.73046251993,0.833,0.005,0.268,0.737,0.42
    if(vlabel.at(ifound)=="A2B4C_oC28_66_l_kl_a")
      PrototypeANRL_A2B4C_oC28_66_l_kl_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 91 // ./aflow --proto=A3B_oC64_66_gi2lm_2l --params=8.157,1.00294225818,0.59212945936,0.54552,0.3287,0.39296,0.16145,0.33837,0.89039,0.24136,0.07877,0.42337,0.74171,0.33285,0.66798,0.24948
    if(vlabel.at(ifound)=="A3B_oC64_66_gi2lm_2l")
      PrototypeANRL_A3B_oC64_66_gi2lm_2l(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 92 // ./aflow --proto=A3B_oC64_66_kl2m_bdl --params=8.735,2.32364052662,1.67842014883,0.1826,-0.0318,0.1994,0.327,0.1716,0.2894,0.451,0.1302,0.1133,0.3773,0.3708
    if(vlabel.at(ifound)=="A3B_oC64_66_kl2m_bdl")
      PrototypeANRL_A3B_oC64_66_kl2m_bdl(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 93 // ./aflow --proto=A2BC_oC16_67_ag_b_g --params=5.0644067238,1.60320657112,1.02379259962,0.3305,0.8198
    if(vlabel.at(ifound)=="A2BC_oC16_67_ag_b_g")
      PrototypeANRL_A2BC_oC16_67_ag_b_g(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 94 // ./aflow --proto=ABC2_oC16_67_b_g_ag --params=5.2729860526,1.00606865161,1.82912952778,0.765,0.3403
    if(vlabel.at(ifound)=="ABC2_oC16_67_b_g_ag")
      PrototypeANRL_ABC2_oC16_67_b_g_ag(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 95 // ./aflow --proto=AB_oC8_67_a_g --params=5.32495,0.997010300566,1.02896740814,0.26686
    // lowering tolerance value for aflowSYM to 0.001 resolves SG #67
    // 96 // ./aflow --proto=AB_oC8_67_a_g --params=5.6124272786,0.999376380873,0.8895303257,0.7642
    if(vlabel.at(ifound)=="AB_oC8_67_a_g")
      PrototypeANRL_AB_oC8_67_a_g(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 97 // ./aflow --proto=AB4_oC20_68_a_i --params=6.44222,1.77662048176,0.991771470083,0.3313,0.1233,0.0801
    if(vlabel.at(ifound)=="AB4_oC20_68_a_i")
      PrototypeANRL_AB4_oC20_68_a_i(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 98 // ./aflow --proto=AB2_oF48_70_f_fg --params=4.2082,3.45504015969,1.7326647973,0.2495,0.54337,0.29445
    if(vlabel.at(ifound)=="AB2_oF48_70_f_fg")
      PrototypeANRL_AB2_oF48_70_f_fg(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 99 // ./aflow --proto=A4B3_oI14_71_gh_cg --params=3.29,4.25531914894,0.951367781155,0.375,0.18,0.444
    if(vlabel.at(ifound)=="A4B3_oI14_71_gh_cg")
      PrototypeANRL_A4B3_oI14_71_gh_cg(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 100 // ./aflow --proto=ABC_oI12_71_h_j_g --params=3.438,3.4554973822,1.37434554974,0.212,0.1232,0.235
    if(vlabel.at(ifound)=="ABC_oI12_71_h_j_g")
      PrototypeANRL_ABC_oI12_71_h_j_g(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 101 // ./aflow --proto=ABCD3_oI48_73_d_e_e_ef --params=5.7750411261,1.02926406926,3.53056277057,0.63427,0.3734,0.18032,0.311,0.1349,0.6124,0.0967
    if(vlabel.at(ifound)=="ABCD3_oI48_73_d_e_e_ef")
      PrototypeANRL_ABCD3_oI48_73_d_e_e_ef(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 102 // ./aflow --proto=A2B_oI12_74_h_e --params=5.158858099,1.56976744188,1.70096899227,-0.047,0.56,0.663
    if(vlabel.at(ifound)=="A2B_oI12_74_h_e")
      PrototypeANRL_A2B_oI12_74_h_e(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 103 // ./aflow --proto=A4B_oI20_74_beh_e --params=4.39,1.42369020501,3.12528473804,-0.111,0.111,-0.033,0.314
    if(vlabel.at(ifound)=="A4B_oI20_74_beh_e")
      PrototypeANRL_A4B_oI20_74_beh_e(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 104 // ./aflow --proto=AB2C12D4_tP76_75_2a2b_2d_12d_4d --params=9.8880614494,0.922431229769,0.3333,0.0,0.5003,0.1666,0.167,0.348,0.1666,0.333,0.152,0.0,0.208,0.152,0.8333,0.458,0.152,0.8333,0.292,0.348,0.6666,0.042,0.348,0.6666,0.208,0.152,0.5,0.458,0.152,0.5,0.292,0.348,0.3333,0.042,0.348,0.3333,0.208,0.152,0.1666,0.458,0.152,0.1666,0.292,0.348,0.0,0.042,0.348,0.0,0.167,0.348,0.8333,0.333,0.152,0.6666,0.167,0.348,0.5,0.333,0.152,0.3333
    if(vlabel.at(ifound)=="AB2C12D4_tP76_75_2a2b_2d_12d_4d")
      PrototypeANRL_AB2C12D4_tP76_75_2a2b_2d_12d_4d(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 105 // ./aflow --proto=A2BC_tP16_76_2a_a_a --params=3.9810644174,3.85606631503,0.143,0.173,0.1347,0.353,0.141,0.2027,0.3461,0.3475,0.5937,0.1519,0.1578,0.0
    if(vlabel.at(ifound)=="A2BC_tP16_76_2a_a_a")
      PrototypeANRL_A2BC_tP16_76_2a_a_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 106 // ./aflow --proto=A3B7_tP40_76_3a_7a --params=9.046,1.84766747734,0.7435,0.3852,0.0,0.4169,0.733,0.8359,0.026,0.8404,-0.0086,0.79,0.6,0.8105,0.443,0.095,-0.053,0.106,0.473,0.8928,0.357,0.024,0.061,0.629,0.794,0.017,-0.002,0.341,0.7049,0.011,0.29,0.84
    if(vlabel.at(ifound)=="A3B7_tP40_76_3a_7a")
      PrototypeANRL_A3B7_tP40_76_3a_7a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 107 // ./aflow --proto=A2B6CD7_tP64_77_2d_6d_d_ab6d --params=7.6159522234,1.07480314961,0.0,0.035,0.1167,0.1203,0.418,0.392,0.3817,0.132,0.362,0.099,0.309,0.136,0.393,0.161,0.192,0.355,0.442,0.312,0.133,0.03,0.028,0.14,0.14,0.47,0.35,0.32,0.2491,0.7483,0.273,0.2602,0.0188,0.013,0.2316,0.4834,0.184,0.1644,0.2577,0.529,0.3406,0.2365,0.012,0.0182,0.2119,0.275,0.4808,0.3008,0.268
    if(vlabel.at(ifound)=="A2B6CD7_tP64_77_2d_6d_d_ab6d")
      PrototypeANRL_A2B6CD7_tP64_77_2d_6d_d_ab6d(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 108 // ./aflow --proto=A2B_tP48_77_8d_4d --params=13.47,0.304380103935,0.193,0.14,0.163,0.091,0.059,0.837,0.693,0.14,0.163,0.591,0.059,0.837,0.193,0.64,0.163,0.193,0.64,0.837,0.693,0.64,0.837,0.591,0.559,0.163,0.11,0.14,0.0,0.61,0.14,0.0,0.11,0.64,0.0,0.61,0.64,0.0
    if(vlabel.at(ifound)=="A2B_tP48_77_8d_4d")
      PrototypeANRL_A2B_tP48_77_8d_4d(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 109 // ./aflow --proto=A2B7C2_tP88_78_4a_14a_4a --params=7.1088964505,3.60337042297,0.14561,0.23414,0.69696,0.22057,0.48276,0.79885,0.18673,0.10525,0.29124,0.24906,0.35443,0.19001,0.4586,0.3299,0.02876,0.0936,0.3899,0.1424,0.4187,0.2147,0.16739,0.1085,0.2261,0.23449,0.2999,0.2626,0.32782,0.0777,0.3681,0.7514,0.2725,0.3751,0.65898,0.0556,0.3427,0.02051,0.3251,0.3037,0.52047,0.4811,0.0546,0.59334,0.0026,0.0009,0.31632,0.3813,0.3334,0.82115,0.2817,0.0576,0.71776,0.168,0.0503,-0.08181,0.26434,0.22687,0.42638,0.02701,0.34822,0.57318,0.15408,0.39584,-0.07135,0.37469,0.25769,0.06464
    if(vlabel.at(ifound)=="A2B7C2_tP88_78_4a_14a_4a")
      PrototypeANRL_A2B7C2_tP88_78_4a_14a_4a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 110 // ./aflow --proto=A2BC2_tI20_79_c_2a_c --params=8.6489709127,0.842525147414,0.0,0.4896,0.337,0.164,0.2196,0.1519,0.1578,0.0
    if(vlabel.at(ifound)=="A2BC2_tI20_79_c_2a_c")
      PrototypeANRL_A2BC2_tI20_79_c_2a_c(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 111 // ./aflow --proto=AB2_tI48_80_2b_4b --params=9.693,0.617455896007,0.2621,0.5076,0.0299,0.2455,0.4909,0.4804,0.3974,0.1497,0.0077,0.1102,0.3642,-0.0098,0.6086,0.3609,0.5064,0.65,0.1038,0.2484
    if(vlabel.at(ifound)=="AB2_tI48_80_2b_4b")
      PrototypeANRL_AB2_tI48_80_2b_4b(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 112 // ./aflow --proto=AB2_tP12_81_adg_2h --params=5.3391020438,1.87980670176,0.25,0.2739,0.234,0.128,0.2289,0.23,0.6373
    if(vlabel.at(ifound)=="AB2_tP12_81_adg_2h")
      PrototypeANRL_AB2_tP12_81_adg_2h(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 113 // ./aflow --proto=A3B_tI32_82_3g_g --params=8.954,0.495711413893,0.0775,0.1117,0.2391,0.3649,0.0321,0.9765,0.1689,0.22,0.7524,0.2862,0.0487,0.4807
    if(vlabel.at(ifound)=="A3B_tI32_82_3g_g")
      PrototypeANRL_A3B_tI32_82_3g_g(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 114 // ./aflow --proto=A3B2_tP10_83_adk_j --params=6.2840064744,0.638128580522,0.375,0.191,0.109,0.314
    if(vlabel.at(ifound)=="A3B2_tP10_83_adk_j")
      PrototypeANRL_A3B2_tP10_83_adk_j(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 115 // ./aflow --proto=A2B_tP30_85_ab2g_cg --params=11.6179902668,0.614391461525,0.6517,0.5428,0.6612,0.5963,0.6531,0.541,0.1258,0.5856,0.1045,0.2524
    if(vlabel.at(ifound)=="A2B_tP30_85_ab2g_cg")
      PrototypeANRL_A2B_tP30_85_ab2g_cg(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 116 // ./aflow --proto=AB3_tP32_86_g_3g --params=9.997999333,0.499899979999,0.0439,0.20812,0.5354,0.11009,0.22151,0.0295,0.14275,0.66613,0.7153,0.53342,0.06957,0.7593
    if(vlabel.at(ifound)=="AB3_tP32_86_g_3g")
      PrototypeANRL_AB3_tP32_86_g_3g(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 117 // ./aflow --proto=A4B_tI20_88_f_a --params=6.407994829,2.01685393259,0.147,0.017,0.298
    if(vlabel.at(ifound)=="A4B_tI20_88_f_a")
      PrototypeANRL_A4B_tI20_88_f_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 118 // ./aflow --proto=AB2_tI96_88_2f_4f --params=13.66,0.436603221083,0.1155,0.1249,0.4746,0.1356,0.125,0.0267,-0.0134,0.1262,-0.0046,-0.0251,0.1252,0.5,0.2739,0.1245,-0.0002,0.2631,0.1241,0.5043
    if(vlabel.at(ifound)=="AB2_tI96_88_2f_4f")
      PrototypeANRL_AB2_tI96_88_2f_4f(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 119 // ./aflow --proto=A17BC4D_tP184_89_17p_p_4p_io --params=16.2503,0.81957871547,0.86772,0.10962,0.1709,0.485,0.7893,0.2334,0.4627,0.7103,0.2146,0.4383,0.6125,0.2885,0.416,0.5631,0.354,0.4242,0.6323,0.3206,0.4579,0.7209,0.234,0.226,0.661,0.3111,0.2275,0.696,0.3148,0.2616,0.794,0.236,0.287,0.818,0.179,0.265,0.733,-0.0696,0.3435,-0.084,-0.078,0.2512,-0.097,-0.004,0.378,0.549,-0.01,0.309,0.593,-0.003,0.238,0.558,-0.006,0.165,0.6,0.2676,0.3463,0.6943,0.1965,0.4869,0.8803,0.096,0.4889,0.7614,0.1143,0.375,0.0177,-0.0146,0.3768,0.8623
    if(vlabel.at(ifound)=="A17BC4D_tP184_89_17p_p_4p_io")
      PrototypeANRL_A17BC4D_tP184_89_17p_p_4p_io(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 120 // ./aflow --proto=A4B2C13D_tP40_90_g_d_cef2g_c --params=7.3739946979,1.45226471387,-0.0654,0.7769,0.83072,0.6476,0.7846,0.648,0.751,0.001,0.7502,0.5354,0.32757,0.7289,0.8811,0.2635
    if(vlabel.at(ifound)=="A4B2C13D_tP40_90_g_d_cef2g_c")
      PrototypeANRL_A4B2C13D_tP40_90_g_d_cef2g_c(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 121 // ./aflow --proto=AB4C17D4E_tP54_90_a_g_c4g_g_c --params=9.5601062765,0.748953974895,0.2035,-0.023,0.7639,0.5152,0.407,0.6525,0.8664,0.4465,0.8786,0.838,0.2708,0.6666,0.8964,0.0943,0.6838,0.656,0.2468,0.7207,0.8114,0.2609
    if(vlabel.at(ifound)=="AB4C17D4E_tP54_90_a_g_c4g_g_c")
      PrototypeANRL_AB4C17D4E_tP54_90_a_g_c4g_g_c(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 122 // ./aflow --proto=ABC_tP24_91_d_d_d --params=3.7620082462,6.71079213192,0.303,0.202,0.019,0.296,0.189,0.08,0.2975,0.1983,0.1795
    if(vlabel.at(ifound)=="ABC_tP24_91_d_d_d")
      PrototypeANRL_ABC_tP24_91_d_d_d(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 123 // ./aflow --proto=AB32CD4E8_tP184_93_i_16p_af_2p_4p --params=11.461,2.92112381119,0.62279,0.369,0.14,0.5403,0.251,0.129,0.6173,0.138,0.647,0.7882,0.215,0.58,0.8645,0.138,0.517,0.591,0.132,0.58,0.5541,0.244,0.605,0.536,0.345,0.554,0.5495,0.343,0.485,0.5842,0.237,0.463,0.6031,0.008,0.366,0.6554,0.077,0.371,0.6892,0.082,0.273,0.7151,0.023,0.17,0.7031,-0.04,0.165,0.6678,-0.05,0.265,0.6408,0.2272,0.1019,0.5636,0.2601,0.5802,0.812,0.1021,0.205,0.5436,0.1936,-0.0676,0.5558,0.404,0.6758,0.8049,0.2809,0.418,0.7912
    if(vlabel.at(ifound)=="AB32CD4E8_tP184_93_i_16p_af_2p_4p")
      PrototypeANRL_AB32CD4E8_tP184_93_i_16p_af_2p_4p(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 124 // ./aflow --proto=A14B3C5_tP44_94_c3g_ad_bg --params=7.3451315192,1.41592920353,0.8225,-0.0016,0.7086,-0.0024,0.6281,-0.0287,0.6575,0.6242,0.54,0.75,0.5448,0.7312,0.7348,0.7403
    if(vlabel.at(ifound)=="A14B3C5_tP44_94_c3g_ad_bg")
      PrototypeANRL_A14B3C5_tP44_94_c3g_ad_bg(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 125 // ./aflow --proto=A6B2C_tP18_94_eg_c_a --params=4.6863209101,1.96124874637,0.6623,0.7093,0.684,0.707,0.6579
    if(vlabel.at(ifound)=="A6B2C_tP18_94_eg_c_a")
      PrototypeANRL_A6B2C_tP18_94_eg_c_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 126 // ./aflow --proto=ABC_tP24_95_d_d_d --params=3.7620082462,6.71079213192,0.303,0.202,-0.019,0.296,0.189,-0.08,0.2975,0.1983,0.8205
    if(vlabel.at(ifound)=="ABC_tP24_95_d_d_d")
      PrototypeANRL_ABC_tP24_95_d_d_d(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 127 // ./aflow --proto=A2B8CD_tI24_97_d_k_a_b --params=5.4068544677,1.92010356944,0.1697,0.3128,0.1237
    if(vlabel.at(ifound)=="A2B8CD_tI24_97_d_k_a_b")
      PrototypeANRL_A2B8CD_tI24_97_d_k_a_b(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 128 // ./aflow --proto=AB8C2_tI44_97_e_2k_cd --params=9.5317026755,1.33879580768,0.1553,0.0449,0.284,0.3693,0.312,0.1212,0.1191
    if(vlabel.at(ifound)=="AB8C2_tI44_97_e_2k_cd")
      PrototypeANRL_AB8C2_tI44_97_e_2k_cd(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 129 // ./aflow --proto=A2B_tI12_98_f_a --params=7.953376649,0.587954231111,0.44
    if(vlabel.at(ifound)=="A2B_tI12_98_f_a")
      PrototypeANRL_A2B_tI12_98_f_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 130 // ./aflow --proto=A2B8C2D_tP26_100_c_abcd_c_a --params=8.527,0.611047261639,0.7904,0.4646,0.3707,0.32701,0.0,0.1259,0.7949,0.1282,0.4871,0.2924,0.5772,0.3571
    if(vlabel.at(ifound)=="A2B8C2D_tP26_100_c_abcd_c_a")
      PrototypeANRL_A2B8C2D_tP26_100_c_abcd_c_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 131 // ./aflow --proto=A3B11C6_tP40_100_ac_bc2d_cd --params=10.1311200179,0.47843253381,0.0,0.0672,0.18111,0.0154,0.6536,0.6912,0.6174,0.0432,0.5814,0.678,0.6402,0.7288,0.574,0.1742,0.7098,0.5782,0.5334
    if(vlabel.at(ifound)=="A3B11C6_tP40_100_ac_bc2d_cd")
      PrototypeANRL_A3B11C6_tP40_100_ac_bc2d_cd(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 132 // ./aflow --proto=A7B7C2_tP32_101_bde_ade_d --params=9.8510697809,0.697188102729,0.0,0.4692,0.17136,0.73945,0.2254,0.3402,0.30926,0.0086,0.2384,0.5244,0.2259,0.0352,0.3449,0.0281
    if(vlabel.at(ifound)=="A7B7C2_tP32_101_bde_ade_d")
      PrototypeANRL_A7B7C2_tP32_101_bde_ade_d(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 133 // ./aflow --proto=A2B3_tP20_102_2c_b2c --params=8.3289849893,0.909833113223,0.5,0.604,0.439,0.623,0.031,0.795,0.725,0.848,0.251
    if(vlabel.at(ifound)=="A2B3_tP20_102_2c_b2c")
      PrototypeANRL_A2B3_tP20_102_2c_b2c(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 134 // ./aflow --proto=AB4_tP10_103_a_d --params=6.5509768136,1.04518394138,0.0,0.144,0.3276,0.242
    if(vlabel.at(ifound)=="AB4_tP10_103_a_d")
      PrototypeANRL_AB4_tP10_103_a_d(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 135 // ./aflow --proto=A5B5C4_tP28_104_ac_ac_c --params=10.6225961282,0.84830508475,0.5,0.8821,0.8116,0.6057,0.3261,0.60942,0.80921,0.00978,0.8116,-0.072,0.1681
    if(vlabel.at(ifound)=="A5B5C4_tP28_104_ac_ac_c")
      PrototypeANRL_A5B5C4_tP28_104_ac_ac_c(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 136 // ./aflow --proto=AB6C4_tP22_104_a_2ac_c --params=9.3940153509,0.981690440703,0.786,0.5,0.0649,0.8297,0.6458,0.286,0.6491,0.8588,0.036
    if(vlabel.at(ifound)=="AB6C4_tP22_104_a_2ac_c")
      PrototypeANRL_AB6C4_tP22_104_a_2ac_c(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 137 // ./aflow --proto=A2BC2_tP20_105_f_ac_2e --params=7.7858653925,1.11276650398,0.0,0.0241,0.3351,0.3664,0.1632,0.6068,0.3453,0.229,0.2558
    if(vlabel.at(ifound)=="A2BC2_tP20_105_f_ac_2e")
      PrototypeANRL_A2BC2_tP20_105_f_ac_2e(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 138 // ./aflow --proto=A3BC3D_tP64_106_3c_c_3c_c --params=10.8312004139,0.492105992062,0.836,0.614,0.2,0.531,0.845,0.11,0.858,0.825,0.08,0.5929,0.62036,0.0343,0.8146,0.6194,0.0688,0.5681,0.8298,0.1322,0.8669,-0.0994,0.0636,0.69285,-0.05677,0.0
    if(vlabel.at(ifound)=="A3BC3D_tP64_106_3c_c_3c_c")
      PrototypeANRL_A3BC3D_tP64_106_3c_c_3c_c(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 139 // ./aflow --proto=A5B7_tI24_107_ac_abd --params=7.6400197048,0.760471204184,0.0,0.056,0.04,0.22,0.0,0.243,0.29
    if(vlabel.at(ifound)=="A5B7_tI24_107_ac_abd")
      PrototypeANRL_A5B7_tI24_107_ac_abd(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 140 // ./aflow --proto=AB_tI4_107_a_a --params=3.5440505103,1.57477426639,0.0,0.427
    if(vlabel.at(ifound)=="AB_tI4_107_a_a")
      PrototypeANRL_AB_tI4_107_a_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 141 // ./aflow --proto=A3B5_tI32_108_ac_a2c --params=8.0549870847,1.94761018001,0.007,0.75,0.109,0.257,0.676,0.114,0.676,0.4
    if(vlabel.at(ifound)=="A3B5_tI32_108_ac_a2c")
      PrototypeANRL_A3B5_tI32_108_ac_a2c(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 142 // ./aflow --proto=ABC_tI12_109_a_a_a --params=4.2490694941,3.42174629325,0.081,0.666,0.5
    if(vlabel.at(ifound)=="ABC_tI12_109_a_a_a")
      PrototypeANRL_ABC_tI12_109_a_a_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 143 // ./aflow --proto=AB_tI8_109_a_a --params=3.4517145504,3.38383984705,0.5416,0.5
    if(vlabel.at(ifound)=="AB_tI8_109_a_a")
      PrototypeANRL_AB_tI8_109_a_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 144 // ./aflow --proto=A2BC8_tI176_110_2b_b_8b --params=13.62003,0.668135092213,0.5491,0.8362,0.259,0.8065,0.6386,0.3704,0.1998,0.0869,0.0,0.5687,-0.0955,0.3098,0.6012,0.7807,0.2082,0.5033,0.7931,0.3482,0.6413,0.5083,0.4159,0.8436,0.6054,0.2752,0.8374,0.7122,0.3904,0.7286,0.645,0.3515,0.815,0.5899,0.4761
    if(vlabel.at(ifound)=="A2BC8_tI176_110_2b_b_8b")
      PrototypeANRL_A2BC8_tI176_110_2b_b_8b(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 145 // ./aflow --proto=A2B_tP12_111_2n_adf --params=5.1219931862,1.02616165562,0.205,0.28,0.301,0.622
    if(vlabel.at(ifound)=="A2B_tP12_111_2n_adf")
      PrototypeANRL_A2B_tP12_111_2n_adf(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 146 // ./aflow --proto=AB_tP8_111_n_n --params=4.1305540686,0.998959140196,0.2522,0.7473,0.24415,0.24404
    // lowering tolerance value for aflowSYM to 0.001 resolves SG #111
    if(vlabel.at(ifound)=="AB_tP8_111_n_n")
      PrototypeANRL_AB_tP8_111_n_n(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 147 // ./aflow --proto=AB4C_tP12_112_b_n_e --params=5.4410982776,1.85862633019,0.2334,0.2761,0.6295
    if(vlabel.at(ifound)=="AB4C_tP12_112_b_n_e")
      PrototypeANRL_AB4C_tP12_112_b_n_e(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 148 // ./aflow --proto=A2BC7D2_tP24_113_e_a_cef_e --params=7.8338,0.639306594501,0.8201,0.8324,0.4935,0.6407,0.7471,0.6396,0.0642,0.0798,0.1862,0.7856
    if(vlabel.at(ifound)=="A2BC7D2_tP24_113_e_a_cef_e")
      PrototypeANRL_A2BC7D2_tP24_113_e_a_cef_e(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 149 // ./aflow --proto=A3B_tP32_114_3e_e --params=9.6362543341,0.547945205485,0.6,0.743,0.809,0.836,0.555,0.604,0.881,0.899,0.844,0.5125,0.7273,0.563
    if(vlabel.at(ifound)=="A3B_tP32_114_3e_e")
      PrototypeANRL_A3B_tP32_114_3e_e(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 150 // ./aflow --proto=A4B_tP10_114_e_a --params=5.2323591487,1.07923706139,0.626,0.768,0.846
    if(vlabel.at(ifound)=="A4B_tP10_114_e_a")
      PrototypeANRL_A4B_tP10_114_e_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 151 // ./aflow --proto=A2B3_tP5_115_g_ag --params=3.3269188443,1.84881274424,0.253,0.6308
    if(vlabel.at(ifound)=="A2B3_tP5_115_g_ag")
      PrototypeANRL_A2B3_tP5_115_g_ag(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 152 // ./aflow --proto=AB2_tP12_115_j_egi --params=8.7863400452,0.701854022742,0.0288,0.0315,0.26398,0.24918,0.24945
    if(vlabel.at(ifound)=="AB2_tP12_115_j_egi")
      PrototypeANRL_AB2_tP12_115_j_egi(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 153 // ./aflow --proto=A2B3_tP20_116_bci_fj --params=6.1720115185,1.606448477,0.177,0.625,0.655,0.216,0.582
    if(vlabel.at(ifound)=="A2B3_tP20_116_bci_fj")
      PrototypeANRL_A2B3_tP20_116_bci_fj(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 154 // ./aflow --proto=A2B3_tP20_117_i_adgh --params=7.7289660931,0.727131582349,0.73,0.73,0.75,0.52,0.25
    if(vlabel.at(ifound)=="A2B3_tP20_117_i_adgh")
      PrototypeANRL_A2B3_tP20_117_i_adgh(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 155 // ./aflow --proto=A3B_tP16_118_ei_f --params=6.9983025398,1.03510852635,0.237,0.15,0.343,0.149,0.509
    if(vlabel.at(ifound)=="A3B_tP16_118_ei_f")
      PrototypeANRL_A3B_tP16_118_ei_f(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 156 // ./aflow --proto=A5B3_tP32_118_g2i_aceh --params=5.8229835854,2.43860552978,0.6709,0.675,0.5861,0.73,0.85,0.5515,0.84,0.8,0.15
    if(vlabel.at(ifound)=="A5B3_tP32_118_g2i_aceh")
      PrototypeANRL_A5B3_tP32_118_g2i_aceh(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 157 // ./aflow --proto=A3B_tI24_119_b2i_af --params=6.315,2.37529691211,0.372,0.2068,0.2229,0.3067,0.3917
    if(vlabel.at(ifound)=="A3B_tI24_119_b2i_af")
      PrototypeANRL_A3B_tI24_119_b2i_af(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 158 // ./aflow --proto=AB_tI4_119_c_a --params=5.4790101504,0.558496075934
    if(vlabel.at(ifound)=="AB_tI4_119_c_a")
      PrototypeANRL_AB_tI4_119_c_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 159 // ./aflow --proto=A4BC2_tI28_120_i_d_e --params=8.8470588481,0.924381146154,0.856,0.6452,0.6575,0.0851
    if(vlabel.at(ifound)=="A4BC2_tI28_120_i_d_e")
      PrototypeANRL_A4BC2_tI28_120_i_d_e(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 160 // ./aflow --proto=A4BC4D_tP10_123_gh_a_i_d --params=3.8757,3.38106664603,0.3336,0.1193,0.2246
    if(vlabel.at(ifound)=="A4BC4D_tP10_123_gh_a_i_d")
      PrototypeANRL_A4BC4D_tP10_123_gh_a_i_d(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 161 // ./aflow --proto=AB4C_tP12_124_a_m_c --params=6.1884808485,0.816448537725,0.162,0.662
    if(vlabel.at(ifound)=="AB4C_tP12_124_a_m_c")
      PrototypeANRL_AB4C_tP12_124_a_m_c(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 162 // ./aflow --proto=AB4_tP10_124_a_m --params=6.4989671731,1.05200800123,0.1425,0.3361
    if(vlabel.at(ifound)=="AB4_tP10_124_a_m")
      PrototypeANRL_AB4_tP10_124_a_m(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 163 // ./aflow --proto=A4B_tP10_125_m_a --params=6.6398746049,0.899096385539,0.425,0.255
    if(vlabel.at(ifound)=="A4B_tP10_125_m_a")
      PrototypeANRL_A4B_tP10_125_m_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 164 // ./aflow --proto=ABC4_tP12_125_a_b_m --params=6.3759876428,1.30630489336,0.3822,0.2163
    if(vlabel.at(ifound)=="ABC4_tP12_125_a_b_m")
      PrototypeANRL_ABC4_tP12_125_a_b_m(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 165 // ./aflow --proto=A2BC4_tP28_126_cd_e_k --params=7.4920241479,1.58609183128,0.11808,0.0865,0.5924,0.1254
    if(vlabel.at(ifound)=="A2BC4_tP28_126_cd_e_k")
      PrototypeANRL_A2BC4_tP28_126_cd_e_k(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 166 // ./aflow --proto=A4B_tP20_127_ehj_g --params=7.256,0.56684123484,0.2,0.31,0.1,0.2,0.04
    if(vlabel.at(ifound)=="A4B_tP20_127_ehj_g")
      PrototypeANRL_A4B_tP20_127_ehj_g(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 167 // ./aflow --proto=A6B2C_tP18_128_eh_d_a --params=7.057532571,1.41383170154,0.2523,0.2217,0.2511
    if(vlabel.at(ifound)=="A6B2C_tP18_128_eh_d_a")
      PrototypeANRL_A6B2C_tP18_128_eh_d_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 168 // ./aflow --proto=A7B2C_tP40_128_egi_h_e --params=6.336,2.34690656566,0.366,0.2008,0.165,0.278,0.088,0.198,0.42,0.1
    if(vlabel.at(ifound)=="A7B2C_tP40_128_egi_h_e")
      PrototypeANRL_A7B2C_tP40_128_egi_h_e(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 169 // ./aflow --proto=A2BC4_tP28_130_f_c_g --params=8.5103337343,0.683196239716,0.58,0.5815,0.045,0.136,0.597
    if(vlabel.at(ifound)=="A2BC4_tP28_130_f_c_g")
      PrototypeANRL_A2BC4_tP28_130_f_c_g(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 170 // ./aflow --proto=A5B3_tP32_130_cg_cf --params=8.465,1.94329592439,0.2271,0.0095,0.1482,0.57997,0.07997,0.10688
    if(vlabel.at(ifound)=="A5B3_tP32_130_cg_cf")
      PrototypeANRL_A5B3_tP32_130_cg_cf(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 171 // ./aflow --proto=A2B2C4D_tP18_132_e_i_o_d --params=5.6046,2.34700067801,0.2369,0.26316,0.34803
    if(vlabel.at(ifound)=="A2B2C4D_tP18_132_e_i_o_d")
      PrototypeANRL_A2B2C4D_tP18_132_e_i_o_d(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 172 // ./aflow --proto=AB6C_tP16_132_d_io_a --params=5.4229923801,1.46597824081,0.3,0.2,0.333
    if(vlabel.at(ifound)=="AB6C_tP16_132_d_io_a")
      PrototypeANRL_AB6C_tP16_132_d_io_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 173 // ./aflow --proto=AB3_tP32_133_h_i2j --params=9.3810096033,0.497068542799,-0.0329,0.8986,0.658,0.0472
    if(vlabel.at(ifound)=="AB3_tP32_133_h_i2j")
      PrototypeANRL_AB3_tP32_133_h_i2j(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 174 // ./aflow --proto=A2B_tP24_135_gh_h --params=8.3218,0.607332548247,0.36248,-0.05789,0.17358,0.13396,0.20929
    if(vlabel.at(ifound)=="A2B_tP24_135_gh_h")
      PrototypeANRL_A2B_tP24_135_gh_h(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 175 // ./aflow --proto=A4B2C_tP28_135_gh_h_d --params=8.4909023894,0.697208809336,0.169,0.114,0.386,0.167,0.175
    if(vlabel.at(ifound)=="A4B2C_tP28_135_gh_h_d")
      PrototypeANRL_A4B2C_tP28_135_gh_h_d(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 176 // ./aflow --proto=A2B3_tP40_137_cdf_3g --params=8.097,1.41410398913,0.0,0.011,0.511,0.533,0.147,0.467,0.864,0.5,0.603
    if(vlabel.at(ifound)=="A2B3_tP40_137_cdf_3g")
      PrototypeANRL_A2B3_tP40_137_cdf_3g(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 177 // ./aflow --proto=A2B_tP6_137_d_a --params=3.64008007,1.4478021978,0.565
    if(vlabel.at(ifound)=="A2B_tP6_137_d_a")
      PrototypeANRL_A2B_tP6_137_d_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 178 // ./aflow --proto=A4BC4_tP18_137_g_b_g --params=5.0593101041,1.39612571654,0.08,0.1,0.503,0.384
    if(vlabel.at(ifound)=="A4BC4_tP18_137_g_b_g")
      PrototypeANRL_A4BC4_tP18_137_g_b_g(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 179 // ./aflow --proto=AB2_tP6_137_a_d --params=4.3675,2.8551803091,0.389
    if(vlabel.at(ifound)=="AB2_tP6_137_a_d")
      PrototypeANRL_AB2_tP6_137_a_d(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 180 // ./aflow --proto=A_tP12_138_bi --params=3.388,1.77420306966,0.086,0.107
    if(vlabel.at(ifound)=="A_tP12_138_bi")
      PrototypeANRL_A_tP12_138_bi(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 181 // ./aflow --proto=AB_tI8_139_e_e --params=4.4795,2.43451278044,0.3356,0.119
    if(vlabel.at(ifound)=="AB_tI8_139_e_e")
      PrototypeANRL_AB_tI8_139_e_e(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 182 // ./aflow --proto=A3B5_tI32_140_ah_bk --params=9.64,0.515560165975,0.17,0.074,0.223
    if(vlabel.at(ifound)=="A3B5_tI32_140_ah_bk")
      PrototypeANRL_A3B5_tI32_140_ah_bk(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 183 // ./aflow --proto=A3B5_tI32_140_ah_cl --params=5.46,1.91575091575,0.625,0.166,0.15
    if(vlabel.at(ifound)=="A3B5_tI32_140_ah_cl")
      PrototypeANRL_A3B5_tI32_140_ah_cl(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    //DX IN PART 1 // 184 // ./aflow --proto=A2B_tI12_141_e_a --params=4.126,3.47697527872,0.2915
    //DX IN PART 1 if(vlabel.at(ifound)=="A2B_tI12_141_e_a")
    //DX IN PART 1  PrototypeANRL_A2B_tI12_141_e_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 185 // ./aflow --proto=A_tI16_142_f --params=8.5939,0.420984651904,0.1405
    if(vlabel.at(ifound)=="A_tI16_142_f")
      PrototypeANRL_A_tI16_142_f(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 186 // ./aflow --proto=A4B14C3_hP21_143_bd_ac4d_d --params=7.3813879247,0.611841213925,0.0,0.305,0.497,0.453,0.09,0.302,0.087,0.605,0.536,0.253,0.059,0.583,0.141,0.429,0.08,0.444,0.296,0.071,0.2215,0.2757,0.806
    if(vlabel.at(ifound)=="A4B14C3_hP21_143_bd_ac4d_d")
      PrototypeANRL_A4B14C3_hP21_143_bd_ac4d_d(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 187 // ./aflow --proto=A4B6C_hP11_143_bd_2d_a --params=6.9672687469,0.526119402977,0.0004,-0.0007,0.8181,0.1915,0.4998,0.5475,0.48,0.0003,0.1859,0.799,0.5002
    // lowering tolerance value for aflowSYM to 0.001 resolves SG #143
    if(vlabel.at(ifound)=="A4B6C_hP11_143_bd_2d_a")
      PrototypeANRL_A4B6C_hP11_143_bd_2d_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 188 // ./aflow --proto=AB2_hP12_143_cd_ab2d --params=6.5003859369,0.944615384626,0.0,0.5,0.25,0.0548,0.2679,0.25,0.33333,0.16667,0.5,0.0,0.5,0.0
    if(vlabel.at(ifound)=="AB2_hP12_143_cd_ab2d")
      PrototypeANRL_AB2_hP12_143_cd_ab2d(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 189 // ./aflow --proto=A4B_hP15_144_4a_a --params=6.2151204722,1.25245374096,0.4904,0.2194,0.2268,0.2226,0.4873,0.1142,0.0775,0.0012,0.0,0.6097,0.0014,0.0018,0.3178,0.0008,0.5062
    if(vlabel.at(ifound)=="A4B_hP15_144_4a_a")
      PrototypeANRL_A4B_hP15_144_4a_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 190 // ./aflow --proto=AB_hP6_144_a_a --params=4.0452497415,2.30951792334,0.33,0.16,0.197,0.13,0.32,0.0
    if(vlabel.at(ifound)=="AB_hP6_144_a_a")
      PrototypeANRL_AB_hP6_144_a_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 191 // ./aflow --proto=A2B3C3DE7_hP48_145_2a_3a_3a_a_7a --params=6.7260126433,2.23669342849,0.567,0.317,0.168,0.432,0.752,0.169,0.6161,-0.0009,0.0,-0.0011,0.6284,0.0045,0.3706,0.371,-0.0043,-0.003,0.266,0.01,0.271,0.002,-0.005,0.73,0.734,-0.02,0.001,0.297,0.165,0.664,0.348,0.09,0.668,0.333,0.238,0.353,0.273,0.167,0.321,0.654,-0.092,0.336,0.666,0.761,0.646,-0.079,0.165,0.203,0.202,0.831
    if(vlabel.at(ifound)=="A2B3C3DE7_hP48_145_2a_3a_3a_a_7a")
      PrototypeANRL_A2B3C3DE7_hP48_145_2a_3a_3a_a_7a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 192 // RHL: ./aflow --proto=A3BC_hR5_146_b_a_a --params=6.8858547087,1.22475238897,0.47,0.0,0.49,-0.002,-0.14399
    // 192 // HEX: ./aflow --proto=A3BC_hR5_146_b_a_a --params=6.8858547087,1.22475238897,0.47,0.0,0.49,-0.002,-0.14399 --hex
    if(vlabel.at(ifound)=="A3BC_hR5_146_b_a_a")
      PrototypeANRL_A3BC_hR5_146_b_a_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 193 // RHL: ./aflow --proto=ABC3_hR10_146_2a_2a_2b --params=6.2648898445,3.16041500397,0.3911,0.7256,0.1132,0.0,0.1454,-0.1886,0.474,0.6151,-0.0128,0.3247
    // 193 // HEX: ./aflow --proto=ABC3_hR10_146_2a_2a_2b --params=6.2648898445,3.16041500397,0.3911,0.7256,0.1132,0.0,0.1454,-0.1886,0.474,0.6151,-0.0128,0.3247 --hex
    if(vlabel.at(ifound)=="ABC3_hR10_146_2a_2a_2b")
      PrototypeANRL_ABC3_hR10_146_2a_2a_2b(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 194 // RHL: ./aflow --proto=A2B4C_hR42_148_2f_4f_f --params=12.4401,0.661634552777,0.20546,-0.56835,0.60961,0.87265,0.10211,-0.72078,0.66221,-0.37042,-0.03997,0.41696,-0.2507,0.0834,0.70548,0.00275,0.03736,0.37614,-0.32772,0.70785,0.53815,-0.23404,-0.05466
    // 194 // HEX: ./aflow --proto=A2B4C_hR42_148_2f_4f_f --params=12.4401,0.661634552777,0.20546,-0.56835,0.60961,0.87265,0.10211,-0.72078,0.66221,-0.37042,-0.03997,0.41696,-0.2507,0.0834,0.70548,0.00275,0.03736,0.37614,-0.32772,0.70785,0.53815,-0.23404,-0.05466 --hex
    if(vlabel.at(ifound)=="A2B4C_hR42_148_2f_4f_f")
      PrototypeANRL_A2B4C_hR42_148_2f_4f_f(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 195 // RHL: ./aflow --proto=A2B_hR18_148_2f_f --params=13.0471,0.659280606418,0.7407,-0.76315,0.0238,1.14609,0.38403,-0.59217,0.08714,0.06386,0.3263
    // 195 // HEX: ./aflow --proto=A2B_hR18_148_2f_f --params=13.0471,0.659280606418,0.7407,-0.76315,0.0238,1.14609,0.38403,-0.59217,0.08714,0.06386,0.3263 --hex
    if(vlabel.at(ifound)=="A2B_hR18_148_2f_f")
      PrototypeANRL_A2B_hR18_148_2f_f(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 196 // ./aflow --proto=AB3_hP24_149_acgi_3l --params=5.1418156677,2.78268310709,0.33333,0.33333,0.0,0.33333,0.421,0.33333,0.33333,0.246,0.0,0.33333,0.088
    if(vlabel.at(ifound)=="AB3_hP24_149_acgi_3l")
      PrototypeANRL_AB3_hP24_149_acgi_3l(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 197 // ./aflow --proto=A3B_hP24_153_3c_2b --params=6.017,2.87518697025,0.1111,0.4444,0.1111,0.2222,0.09357,0.4444,0.8888,0.09357,0.77778,0.55558,0.09357
    if(vlabel.at(ifound)=="A3B_hP24_153_3c_2b")
      PrototypeANRL_A3B_hP24_153_3c_2b(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 198 // ./aflow --proto=A_hP9_154_bc --params=6.9082,0.616557134999,0.876,0.23,0.534,0.051
    if(vlabel.at(ifound)=="A_hP9_154_bc")
      PrototypeANRL_A_hP9_154_bc(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 199 // ./aflow --proto=AB2_hP9_156_b2c_3a2bc --params=4.239807232,4.83608490569,0.0,0.66667,0.33333,0.75,0.16667,0.5,0.08333,0.41667,0.83333
    if(vlabel.at(ifound)=="AB2_hP9_156_b2c_3a2bc")
      PrototypeANRL_AB2_hP9_156_b2c_3a2bc(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 200 // ./aflow --proto=AB_hP12_156_2ab3c_2ab3c --params=4.2499813346,4.90823529409,0.375,0.70833,0.5,0.83333,0.04167,0.16667,0.45833,0.79167,0.125,0.33333,0.66667,0.0
    if(vlabel.at(ifound)=="AB_hP12_156_2ab3c_2ab3c")
      PrototypeANRL_AB_hP12_156_2ab3c_2ab3c(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // --------------------------------------------------------------------------- //DX 20180925 - changed label/params - START
    // 201 // ./aflow --proto=AB_hP4_156_ab_ab --params=4.2794836776,1.67515774714,0.636,0.0,0.894,0.5 
    if(found && vlabel.at(ifound)=="AB_hP4_156_ab_ab")
      PrototypeANRL_AB_hP4_156_ab_ab(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // --------------------------------------------------------------------------- //DX 20180925 - changed label/params - END
    // 202 // ./aflow --proto=A5B6C2_hP13_157_2ac_2c_b --params=5.939,1.08233709379,0.264,0.736,0.022,0.5,0.522,0.603,0.215,0.397,0.829
    if(vlabel.at(ifound)=="A5B6C2_hP13_157_2ac_2c_b")
      PrototypeANRL_A5B6C2_hP13_157_2ac_2c_b(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 203 // ./aflow --proto=A3B_hP8_158_d_a --params=6.12,0.924509803922,0.0,0.318,0.027,0.237
    if(vlabel.at(ifound)=="A3B_hP8_158_d_a")
      PrototypeANRL_A3B_hP8_158_d_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 204 // ./aflow --proto=A2B3_hP20_159_bc_2c --params=7.7488577892,0.813266227893,0.0,0.337,0.154,0.021,0.451,0.061,0.287,0.149,0.282,0.083
    if(vlabel.at(ifound)=="A2B3_hP20_159_bc_2c")
      PrototypeANRL_A2B3_hP20_159_bc_2c(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 205 // ./aflow --proto=A4B3_hP28_159_ab2c_2c --params=7.7479899913,0.724961280333,0.0,0.3649,0.0424,0.3891,0.0408,0.3169,0.3198,0.2712,0.0821,0.5089,0.3172,0.1712,0.2563,0.0274
    if(vlabel.at(ifound)=="A4B3_hP28_159_ab2c_2c")
      PrototypeANRL_A4B3_hP28_159_ab2c_2c(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 206 // ./aflow --proto=AB4C7D_hP26_159_b_ac_a2c_b --params=6.2653109789,1.6324735851,0.0,0.3057,0.0613,0.4379,0.3425,0.1575,0.2472,0.0015,0.5149,0.3102,0.3347,0.1157,0.0618
    if(vlabel.at(ifound)=="AB4C7D_hP26_159_b_ac_a2c_b")
      PrototypeANRL_AB4C7D_hP26_159_b_ac_a2c_b(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 207 // RHL: ./aflow --proto=A3B_hR4_160_b_a --params=4.405,0.610442678774,0.0,0.52073,-0.02217
    // 207 // HEX: ./aflow --proto=A3B_hR4_160_b_a --params=4.405,0.610442678774,0.0,0.52073,-0.02217 --hex
    if(vlabel.at(ifound)=="A3B_hR4_160_b_a")
      PrototypeANRL_A3B_hR4_160_b_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 208 // RHL: ./aflow --proto=A8B5_hR26_160_a3bc_a3b --params=12.71838,0.624030733474,0.672,0.194,0.654,0.01201,0.349,0.58199,0.722,0.35601,1.003,-0.20599,0.998,-0.66,0.355,1.00601,1.033,-0.339,0.28799
    // 208 // HEX: ./aflow --proto=A8B5_hR26_160_a3bc_a3b --params=12.71838,0.624030733474,0.672,0.194,0.654,0.01201,0.349,0.58199,0.722,0.35601,1.003,-0.20599,0.998,-0.66,0.355,1.00601,1.033,-0.339,0.28799 --hex
    if(vlabel.at(ifound)=="A8B5_hR26_160_a3bc_a3b")
      PrototypeANRL_A8B5_hR26_160_a3bc_a3b(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 209 // RHL: ./aflow --proto=ABC_hR3_160_a_a_a --params=6.1703,0.949908432329,0.0,0.79356,0.25763
    // 209 // HEX: ./aflow --proto=ABC_hR3_160_a_a_a --params=6.1703,0.949908432329,0.0,0.79356,0.25763 --hex
    if(vlabel.at(ifound)=="ABC_hR3_160_a_a_a")
      PrototypeANRL_ABC_hR3_160_a_a_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 210 // RHL: ./aflow --proto=AB_hR10_160_5a_5a --params=3.09,12.2653721683,0.0,0.13333,0.4,0.6,0.86667,0.05,0.18333,0.45,0.65,-0.08333
    // 210 // HEX: ./aflow --proto=AB_hR10_160_5a_5a --params=3.09,12.2653721683,0.0,0.13333,0.4,0.6,0.86667,0.05,0.18333,0.45,0.65,-0.08333 --hex
    if(vlabel.at(ifound)=="AB_hR10_160_5a_5a")
      PrototypeANRL_AB_hR10_160_5a_5a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 211 // ./aflow --proto=A2B3_hP5_164_d_ad --params=3.9381,1.55813717275,0.2467,0.647
    if(vlabel.at(ifound)=="A2B3_hP5_164_d_ad")
      PrototypeANRL_A2B3_hP5_164_d_ad(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 212 // ./aflow --proto=AB2_hP9_164_bd_c2d --params=2.89,7.90657439446,0.0607,0.154,0.27263,0.39403
    if(vlabel.at(ifound)=="AB2_hP9_164_bd_c2d")
      PrototypeANRL_AB2_hP9_164_bd_c2d(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 213 // ./aflow --proto=ABC2_hP4_164_a_b_d --params=4.0482,1.26777333136,0.271
    if(vlabel.at(ifound)=="ABC2_hP4_164_a_b_d")
      PrototypeANRL_ABC2_hP4_164_a_b_d(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 214 // ./aflow --proto=A3B_hP24_165_bdg_f --params=7.07,1.00919377652,0.17,0.38,0.69,0.07,0.08
    if(vlabel.at(ifound)=="A3B_hP24_165_bdg_f")
      PrototypeANRL_A3B_hP24_165_bdg_f(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 215 // RHL: ./aflow --proto=A4B3_hR7_166_2c_ac --params=3.335,7.48635682159,0.29422,0.12967,0.2168
    // 215 // HEX: ./aflow --proto=A4B3_hR7_166_2c_ac --params=3.335,7.48635682159,0.29422,0.12967,0.2168 --hex
    if(vlabel.at(ifound)=="A4B3_hR7_166_2c_ac")
      PrototypeANRL_A4B3_hR7_166_2c_ac(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 216 // RHL: ./aflow --proto=ABC_hR6_166_c_c_c --params=3.8548,7.95397945419,0.1159,0.3017,0.3815
    // 216 // HEX: ./aflow --proto=ABC_hR6_166_c_c_c --params=3.8548,7.95397945419,0.1159,0.3017,0.3815 --hex
    if(vlabel.at(ifound)=="ABC_hR6_166_c_c_c")
      PrototypeANRL_ABC_hR6_166_c_c_c(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 217 // RHL: ./aflow --proto=AB3C_hR10_167_b_e_a --params=5.4577,2.40134122433,0.6946
    // 217 // HEX: ./aflow --proto=AB3C_hR10_167_b_e_a --params=5.4577,2.40134122433,0.6946 --hex
    if(vlabel.at(ifound)=="AB3C_hR10_167_b_e_a")
      PrototypeANRL_AB3C_hR10_167_b_e_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 218 // RHL: ./aflow --proto=ABC2_hR24_167_e_e_2e --params=12.76,0.575235109718,1.1389,0.8113,1.0343,0.3584
    // 218 // HEX: ./aflow --proto=ABC2_hR24_167_e_e_2e --params=12.76,0.575235109718,1.1389,0.8113,1.0343,0.3584 --hex
    if(vlabel.at(ifound)=="ABC2_hR24_167_e_e_2e")
      PrototypeANRL_ABC2_hR24_167_e_e_2e(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 219 // ./aflow --proto=A2B13C4_hP57_168_d_c6d_2d --params=15.9361426085,0.244226907628,0.0,0.4965,0.179,0.555,0.035,0.187,0.007,0.349,0.034,0.018,0.168,0.388,0.01,0.195,0.569,0.02,0.093,0.454,0.539,0.273,0.095,0.555,0.2634,0.09,0.086,0.0789,0.4368,0.071
    if(vlabel.at(ifound)=="A2B13C4_hP57_168_d_c6d_2d")
      PrototypeANRL_A2B13C4_hP57_168_d_c6d_2d(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 220 // ./aflow --proto=AB4C_hP72_168_2d_8d_2d --params=13.759939558,0.609738372094,0.4498,0.113,0.6253,0.1278,0.466,0.1246,0.4154,0.2053,0.0456,0.2066,0.4189,0.5576,0.4218,0.0894,0.8248,0.1514,0.4907,0.3249,0.3746,0.0106,0.0948,0.0083,0.3667,0.4879,0.5693,0.1616,0.0357,0.1505,0.5625,0.5943,0.4451,0.1173,0.0,0.1291,0.4589,0.4993
    if(vlabel.at(ifound)=="AB4C_hP72_168_2d_8d_2d")
      PrototypeANRL_AB4C_hP72_168_2d_8d_2d(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 221 // ./aflow --proto=A2B3_hP30_169_2a_3a --params=6.4300116544,2.78071539658,0.013,0.3579,0.12736,0.3339,0.3226,0.29886,0.3347,0.0,0.004,0.0119,0.3343,0.0,0.338,0.0064,0.33823
    if(vlabel.at(ifound)=="A2B3_hP30_169_2a_3a")
      PrototypeANRL_A2B3_hP30_169_2a_3a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 222 // ./aflow --proto=A2B3_hP30_170_2a_3a --params=6.4300116544,2.78071539658,0.013,0.3579,0.87264,0.3339,0.3226,0.70114,0.3347,0.0,-0.004,0.0119,0.3343,0.0,0.338,0.0064,0.66177
    if(vlabel.at(ifound)=="A2B3_hP30_170_2a_3a")
      PrototypeANRL_A2B3_hP30_170_2a_3a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 223 // ./aflow --proto=A10B2C_hP39_171_5c_c_a --params=6.3199906634,3.05221518986,0.0,-0.004,0.253,0.23633,0.015,0.248,0.11034,0.681,0.192,0.16733,0.448,0.188,0.29533,0.449,0.25,0.048,0.373,0.065,0.50467
    if(vlabel.at(ifound)=="A10B2C_hP39_171_5c_c_a")
      PrototypeANRL_A10B2C_hP39_171_5c_c_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 224 // ./aflow --proto=A10B2C_hP39_172_5c_c_a --params=6.3199906634,3.05221518986,0.0,-0.004,0.253,0.76367,0.015,0.248,0.88966,0.681,0.192,0.83267,0.448,0.188,0.70467,0.449,0.25,-0.048,0.373,0.065,0.49533
    if(vlabel.at(ifound)=="A10B2C_hP39_172_5c_c_a")
      PrototypeANRL_A10B2C_hP39_172_5c_c_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 225 // ./aflow --proto=A3B_hP8_173_c_b --params=7.1329719992,1.03939436422,0.0,0.0337,0.3475,0.146
    if(vlabel.at(ifound)=="A3B_hP8_173_c_b")
      PrototypeANRL_A3B_hP8_173_c_b(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 226 // ./aflow --proto=A4B3_hP14_173_bc_c --params=7.603038022,0.382612126795,0.0,0.3284,0.0313,0.05,0.2314,0.4063,0.013
    if(vlabel.at(ifound)=="A4B3_hP14_173_bc_c")
      PrototypeANRL_A4B3_hP14_173_bc_c(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 227 // ./aflow --proto=A12B7C2_hP21_174_2j2k_ajk_cf --params=9.0004021308,0.399102242177,0.4309,0.3719,0.1189,0.2772,0.4163,0.1204,0.0495,0.4359,0.2232,0.124,0.2889,0.4096
    if(vlabel.at(ifound)=="A12B7C2_hP21_174_2j2k_ajk_cf")
      PrototypeANRL_A12B7C2_hP21_174_2j2k_ajk_cf(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 228 // ./aflow --proto=ABC_hP12_174_cj_fk_aj --params=10.7303215747,0.395153774459,0.30167,0.15433,0.03467,0.51733,0.14967,0.31433
    if(vlabel.at(ifound)=="ABC_hP12_174_cj_fk_aj")
      PrototypeANRL_ABC_hP12_174_cj_fk_aj(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 229 // ./aflow --proto=A8B7C6_hP21_175_ck_aj_k --params=9.5057625379,0.329093950194,0.36335,0.08627,0.0661,0.221,0.15405,0.51373
    if(vlabel.at(ifound)=="A8B7C6_hP21_175_ck_aj_k")
      PrototypeANRL_A8B7C6_hP21_175_ck_aj_k(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 230 // ./aflow --proto=ABC_hP36_175_jk_jk_jk --params=11.5799622371,0.317895264089,0.255,0.0598,0.1323,0.416,0.3483,0.0762,0.2334,0.5727,0.4279,0.0715,0.1367,0.5046
    if(vlabel.at(ifound)=="ABC_hP36_175_jk_jk_jk")
      PrototypeANRL_ABC_hP36_175_jk_jk_jk(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // --------------------------------------------------------------------------- //DX 20180925 - changed label/params - START
    // 231 // ./aflow --proto=A3B2_hP10_176_h_bc --params=7.8700439404,0.50025412961,0.3847,0.0915
    if(found && vlabel.at(ifound)=="A3B2_hP10_176_h_bc")
      PrototypeANRL_A3B2_hP10_176_h_bc(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 232 // ./aflow --proto=A3B3C_hP14_176_h_h_c --params=9.3500327107,0.451336898402,0.1493,0.1701,0.357,0.0462
    if(found && vlabel.at(ifound)=="A3B3C_hP14_176_h_h_c")
      PrototypeANRL_A3B3C_hP14_176_h_h_c(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 233 // ./aflow --proto=A3B_hP8_176_h_c --params=7.4429335392,0.580545478976,0.375,0.083 
    if(found && vlabel.at(ifound)=="A3B_hP8_176_h_c")
      PrototypeANRL_A3B_hP8_176_h_c(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // --------------------------------------------------------------------------- //DX 20180925 - changed label/params - END
    // 234 // ./aflow --proto=A2B_hP36_177_j2lm_n --params=12.7835,0.291064262526,0.61855,0.39242,0.79257,0.44445,0.52169,0.86952,0.16458
    if(vlabel.at(ifound)=="A2B_hP36_177_j2lm_n")
      PrototypeANRL_A2B_hP36_177_j2lm_n(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 235 // ./aflow --proto=AB3_hP24_178_b_ac --params=5.14898393,3.15789473684,0.8361,0.7601,0.5338,0.3099,-0.0053
    if(vlabel.at(ifound)=="AB3_hP24_178_b_ac")
      PrototypeANRL_AB3_hP24_178_b_ac(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 236 // ./aflow --proto=A_hP6_178_a --params=2.355,4.43566878981,0.461
    if(vlabel.at(ifound)=="A_hP6_178_a")
      PrototypeANRL_A_hP6_178_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 237 // ./aflow --proto=AB3_hP24_179_b_ac --params=5.14898393,3.15789473684,0.8361,0.7601,0.5338,0.3099,0.0053
    if(vlabel.at(ifound)=="AB3_hP24_179_b_ac")
      PrototypeANRL_AB3_hP24_179_b_ac(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 238 // ./aflow --proto=A2B_hP9_181_j_c --params=4.9977,1.09252256038,0.2072
    if(vlabel.at(ifound)=="A2B_hP9_181_j_c")
      PrototypeANRL_A2B_hP9_181_j_c(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 239 // ./aflow --proto=ABC_hP3_183_a_a_a --params=3.3908401495,1.49568037743,0.608,0.0,0.226
    if(vlabel.at(ifound)=="ABC_hP3_183_a_a_a")
      PrototypeANRL_ABC_hP3_183_a_a_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 240 // ./aflow --proto=AB_hP6_183_c_ab --params=5.3175214551,0.83247442072,0.0,0.513,0.01
    if(vlabel.at(ifound)=="AB_hP6_183_c_ab")
      PrototypeANRL_AB_hP6_183_c_ab(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 241 // ./aflow --proto=AB4C_hP72_184_d_4d_d --params=13.7178848276,0.616168537688,0.45652,0.12053,0.0,0.42,0.2069,0.071,0.1224,0.4519,0.2982,0.3629,0.0019,0.0649,0.1538,0.5737,0.0602,0.12298,0.4525,0.12746
    if(vlabel.at(ifound)=="AB4C_hP72_184_d_4d_d")
      PrototypeANRL_AB4C_hP72_184_d_4d_d(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 242 // ./aflow --proto=A3BC_hP30_185_cd_c_ab --params=11.795090915,0.502416278086,0.0,0.377,0.1598,0.2396,0.6647,0.1706,0.1732,0.5056,0.1148
    if(vlabel.at(ifound)=="A3BC_hP30_185_cd_c_ab")
      PrototypeANRL_A3BC_hP30_185_cd_c_ab(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 243 // ./aflow --proto=A3B_hP24_185_ab2c_c --params=6.9593,1.02639633296,0.3213,0.1998,0.2806,0.0765,0.3761,0.4246,0.3322,0.75
    if(vlabel.at(ifound)=="A3B_hP24_185_ab2c_c")
      PrototypeANRL_A3B_hP24_185_ab2c_c(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 244 // ./aflow --proto=A3B_hP8_185_c_a --params=6.1197165237,0.924509803925,0.0,0.305,0.245
    if(vlabel.at(ifound)=="A3B_hP8_185_c_a")
      PrototypeANRL_A3B_hP8_185_c_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 245 // ./aflow --proto=AB3_hP24_185_c_ab2c --params=8.7838,1.02449964708,0.2684,0.2311,0.3321,0.25,0.3153,0.5863,0.3518,-0.0769
    if(vlabel.at(ifound)=="AB3_hP24_185_c_ab2c")
      PrototypeANRL_AB3_hP24_185_c_ab2c(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 246 // ./aflow --proto=A3B7_hP20_186_c_b2c --params=9.85,0.624365482234,0.06,0.815,0.31,0.126,0.25,0.544,0.31
    if(vlabel.at(ifound)=="A3B7_hP20_186_c_b2c")
      PrototypeANRL_A3B7_hP20_186_c_b2c(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 247 // ./aflow --proto=AB3_hP4_187_e_fh --params=2.8065,2.53411722786,0.198
    if(vlabel.at(ifound)=="AB3_hP4_187_e_fh")
      PrototypeANRL_AB3_hP4_187_e_fh(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // --------------------------------------------------------------------------- //DX 20180925 - changed label/params - START
    // 248 // ./aflow --proto=A3BC_hP10_188_k_c_a --params=7.2864263258,0.928891999832,0.0041,0.32685
    if(found && vlabel.at(ifound)=="A3BC_hP10_188_k_c_a")
      PrototypeANRL_A3BC_hP10_188_k_c_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // --------------------------------------------------------------------------- //DX 20180925 - changed label/params - END
    // 249 // ./aflow --proto=AB9C4_hP28_188_e_kl_ak --params=6.4953629976,1.43896355826,0.07103,0.48306,0.12023,0.75436,0.22923,0.00127,0.6032
    if(vlabel.at(ifound)=="AB9C4_hP28_188_e_kl_ak")
      PrototypeANRL_AB9C4_hP28_188_e_kl_ak(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 250 // ./aflow --proto=A8BC3D6_hP18_189_bfh_a_g_i --params=6.62,1.20241691843,0.403,0.444,0.231,0.75,0.222
    if(vlabel.at(ifound)=="A8BC3D6_hP18_189_bfh_a_g_i")
      PrototypeANRL_A8BC3D6_hP18_189_bfh_a_g_i(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 251 // ./aflow --proto=A9BC3D5_hP18_189_fi_a_g_bh --params=6.6,1.19696969697,0.378,0.43,0.266,0.755,0.236
    if(vlabel.at(ifound)=="A9BC3D5_hP18_189_fi_a_g_bh")
      PrototypeANRL_A9BC3D5_hP18_189_fi_a_g_bh(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 252 // ./aflow --proto=A2B_hP18_190_gh_bf --params=7.946,0.789076264787,0.0225,0.294,0.612,0.01
    if(vlabel.at(ifound)=="A2B_hP18_190_gh_bf")
      PrototypeANRL_A2B_hP18_190_gh_bf(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 253 // ./aflow --proto=A5B3_hP16_190_bdh_g --params=6.9236970109,1.22634969237,0.3313,0.0628,0.6682
    if(vlabel.at(ifound)=="A5B3_hP16_190_bdh_g")
      PrototypeANRL_A5B3_hP16_190_bdh_g(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 254 // ./aflow --proto=AB_hP24_190_i_afh --params=5.9699820408,1.96984924623,0.52,0.6683,0.6653,0.3786,0.3233,0.623
    if(vlabel.at(ifound)=="AB_hP24_190_i_afh")
      PrototypeANRL_AB_hP24_190_i_afh(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 255 // ./aflow --proto=A2B3C18D6_hP58_192_c_f_lm_l --params=9.214,0.997829390059,0.3103,0.2369,0.3876,0.1159,0.4985,0.1456,0.1453
    if(vlabel.at(ifound)=="A2B3C18D6_hP58_192_c_f_lm_l")
      PrototypeANRL_A2B3C18D6_hP58_192_c_f_lm_l(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 256 // ./aflow --proto=AB2_hP72_192_m_j2kl --params=13.7705812079,0.608458538779,0.6373,0.7895,0.4221,0.6688,0.5445,0.1221,0.4551,0.686
    if(vlabel.at(ifound)=="AB2_hP72_192_m_j2kl")
      PrototypeANRL_AB2_hP72_192_m_j2kl(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 257 // ./aflow --proto=A5B3_hP16_193_dg_g --params=6.9104160691,0.696671490596,0.2358,0.5992
    if(vlabel.at(ifound)=="A5B3_hP16_193_dg_g")
      PrototypeANRL_A5B3_hP16_193_dg_g(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 258 // ./aflow --proto=A3B_hP16_194_gh_ac --params=5.096,1.6295133438,-0.16667
    if(vlabel.at(ifound)=="A3B_hP16_194_gh_ac")
      PrototypeANRL_A3B_hP16_194_gh_ac(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 259 // ./aflow --proto=A5B2_hP28_194_ahk_ch --params=7.656,0.991771159875,0.533,0.872,0.196,0.439
    if(vlabel.at(ifound)=="A5B2_hP28_194_ahk_ch")
      PrototypeANRL_A5B2_hP28_194_ahk_ch(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 260 // ./aflow --proto=A9B3C_hP26_194_hk_h_a --params=7.513,1.03087980833,0.458,0.12,0.201,-0.067
    if(vlabel.at(ifound)=="A9B3C_hP26_194_hk_h_a")
      PrototypeANRL_A9B3C_hP26_194_hk_h_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 261 // ./aflow --proto=A12BC4_cP34_195_2j_ab_2e --params=8.0357772599,0.2493,0.7507,0.14225,0.5,0.35579,0.0002,0.14286,0.35831
    // lowering tolerance value for aflowSYM to 0.001 resolves SG #195
    if(vlabel.at(ifound)=="A12BC4_cP34_195_2j_ab_2e")
      PrototypeANRL_A12BC4_cP34_195_2j_ab_2e(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 262 // ./aflow --proto=A12B2C_cF60_196_h_bc_a --params=9.9799933334,0.0,0.0625,0.25
    if(vlabel.at(ifound)=="A12B2C_cF60_196_h_bc_a")
      PrototypeANRL_A12B2C_cF60_196_h_bc_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    //// 263 // ./aflow --proto=A12B36CD12_cF488_196_2h_6h_ac_fgh --params=16.4321054599,0.12536,0.37536,0.55423,0.08845,0.00025,0.83845,0.30423,0.24975,0.0008,0.09,0.3557,0.84,0.2508,0.3943,0.0123,0.0411,0.1708,0.2911,0.2623,0.0792,0.0347,0.2125,0.0888,0.4625,0.2847,0.1612,0.125,0.01382,0.23628
    //// lowering tolerance value for aflowSYM to 0.001 resolves SG #196
    //if(vlabel.at(ifound)=="A12B36CD12_cF488_196_2h_6h_ac_fgh")
    //  PrototypeANRL_A12B36CD12_cF488_196_2h_6h_ac_fgh(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    //// ---------------------------------------------------------------------------
    // 263 // ./aflow --proto=ABC3_cP20_198_a_a_b --params=6.57,0.417,0.064,0.303,0.592,0.5
    if(vlabel.at(ifound)=="ABC3_cP20_198_a_a_b")
      PrototypeANRL_ABC3_cP20_198_a_a_b(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 264 // ./aflow --proto=A2B11_cP39_200_f_aghij --params=8.5520223662,0.18,0.34,0.265,0.278,0.157,0.257
    if(vlabel.at(ifound)=="A2B11_cP39_200_f_aghij")
      PrototypeANRL_A2B11_cP39_200_f_aghij(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // --------------------------------------------------------------------------- //DX 20180925 - changed label/params - START
    // 265 // ./aflow --proto=AB3C_cP60_201_be_fh_g --params=9.5599167841,0.3889,0.6111,0.0972,0.0389,0.0972,0.75
    if(found && vlabel.at(ifound)=="AB3C_cP60_201_be_fh_g")
      PrototypeANRL_AB3C_cP60_201_be_fh_g(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // --------------------------------------------------------------------------- //DX 20180925 - changed label/params - START
    // 266 // ./aflow --proto=A6B6C_cF104_202_h_h_c --params=10.6100296668,0.5827,0.6359,0.638,0.72
    if(vlabel.at(ifound)=="A6B6C_cF104_202_h_h_c")
      PrototypeANRL_A6B6C_cF104_202_h_h_c(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 267 // ./aflow --proto=A_cF240_202_h2i --params=14.26,0.249,0.052,0.105,0.085,0.22,0.185,0.052,0.165
    if(vlabel.at(ifound)=="A_cF240_202_h2i")
      PrototypeANRL_A_cF240_202_h2i(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 268 // ./aflow --proto=A2BCD3E6_cF208_203_e_c_d_f_g --params=13.9898,0.2826,-0.099,0.2257,0.2665,0.3531
    if(vlabel.at(ifound)=="A2BCD3E6_cF208_203_e_c_d_f_g")
      PrototypeANRL_A2BCD3E6_cF208_203_e_c_d_f_g(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 269 // ./aflow --proto=A4B2C6D16E_cF232_203_e_d_f_eg_a --params=13.9038,0.28207,0.06362,0.34379,0.26626,0.22529,0.35333
    if(vlabel.at(ifound)=="A4B2C6D16E_cF232_203_e_d_f_eg_a")
      PrototypeANRL_A4B2C6D16E_cF232_203_e_d_f_eg_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // --------------------------------------------------------------------------- //DX 20180925 - changed label/params - START
    // 270 // ./aflow --proto=AB3C16_cF160_203_a_bc_eg --params=16.6600284675,0.20516,0.01201,0.111,0.42978
    if(found && vlabel.at(ifound)=="AB3C16_cF160_203_a_bc_eg")
      PrototypeANRL_AB3C16_cF160_203_a_bc_eg(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // --------------------------------------------------------------------------- //DX 20180925 - changed label/params - END
    // 271 // ./aflow --proto=A2B3C6_cP264_205_2d_ab2c2d_6d --params=15.263,0.2561,0.375,0.2526,0.0133,0.0197,0.2444,0.2335,0.0046,0.1386,0.3763,0.1272,0.38,0.3838,0.1209,0.2777,0.1241,0.0103,0.4835,0.1315,0.2536,0.2664,0.2841,0.1049,0.235,0.4047,0.2921,0.3491,-0.0385,-0.1074,0.1509,-0.0104,-0.0242
    if(vlabel.at(ifound)=="A2B3C6_cP264_205_2d_ab2c2d_6d")
      PrototypeANRL_A2B3C6_cP264_205_2d_ab2c2d_6d(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 272 // ./aflow --proto=A_cP240_205_10d --params=14.04078,0.2294,-0.0325,0.101,0.2467,-0.054,0.0061,0.2081,0.0646,0.1289,0.2066,0.8599,-0.036,0.171,-0.0963,0.159,0.2236,0.1122,-0.0371,0.2439,0.0192,-0.0636,0.2053,0.1349,0.0616,0.1503,0.7983,0.0202,0.1323,0.8207,0.1186
    if(vlabel.at(ifound)=="A_cP240_205_10d")
      PrototypeANRL_A_cP240_205_10d(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 273 // ./aflow --proto=AB3C2_cI96_206_c_e_ad --params=9.46,0.115,0.205,0.16,0.382,0.11
    if(vlabel.at(ifound)=="AB3C2_cI96_206_c_e_ad")
      PrototypeANRL_AB3C2_cI96_206_c_e_ad(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 274 // ./aflow --proto=A17B15_cP64_207_acfk_eij --params=10.6058825779,0.2422,0.2622,0.3319,0.2701,0.142,0.1539,0.3498
    if(vlabel.at(ifound)=="A17B15_cP64_207_acfk_eij")
      PrototypeANRL_A17B15_cP64_207_acfk_eij(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 275 // ./aflow --proto=A3B_cP16_208_j_b --params=6.31,0.184
    if(vlabel.at(ifound)=="A3B_cP16_208_j_b")
      PrototypeANRL_A3B_cP16_208_j_b(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 276 // ./aflow --proto=A6B2CD6E_cP64_208_m_ad_b_m_c --params=10.3701312618,0.057,0.25,0.23,0.235,0.25,0.458
    if(vlabel.at(ifound)=="A6B2CD6E_cP64_208_m_ad_b_m_c")
      PrototypeANRL_A6B2CD6E_cP64_208_m_ad_b_m_c(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 277 // ./aflow --proto=A24BC_cF104_209_j_a_b --params=7.7099775082,0.043,0.109,0.165
    if(vlabel.at(ifound)=="A24BC_cF104_209_j_a_b")
      PrototypeANRL_A24BC_cF104_209_j_a_b(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // --------------------------------------------------------------------------- //DX 20180925 - new proto (used to be SG #196) - START
    // 278 // ./aflow --proto=A12B36CD12_cF488_210_h_3h_a_fg --params=16.4321054599,0.12536,0.51382,0.55423,0.58845,0.50025,0.5008,0.59,0.3557,0.5123,0.5411,0.1708,0.5347,0.2125,0.5888
    if(found && vlabel.at(ifound)=="A12B36CD12_cF488_210_h_3h_a_fg")
      PrototypeANRL_A12B36CD12_cF488_210_h_3h_a_fg(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // --------------------------------------------------------------------------- //DX 20180925 - new proto (used to be SG #196) - END
    // 279 // ./aflow --proto=A12B6C_cF608_210_4h_2h_e --params=16.4321054599,0.3771,0.0157,0.2009,0.1224,0.5287,0.1425,0.0068,0.0949,0.2596,0.2225,0.6117,0.1785,0.0224,0.0928,0.2607,0.1616,0.0194,0.0793,0.3503
    if(vlabel.at(ifound)=="A12B6C_cF608_210_4h_2h_e")
      PrototypeANRL_A12B6C_cF608_210_4h_2h_e(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 280 // ./aflow --proto=A2B_cI72_211_hi_i --params=9.68882,0.37338,0.15866,0.38235
    if(vlabel.at(ifound)=="A2B_cI72_211_hi_i")
      PrototypeANRL_A2B_cI72_211_hi_i(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 281 // ./aflow --proto=A2B_cP12_212_c_a --params=6.54,0.428
    if(vlabel.at(ifound)=="A2B_cP12_212_c_a")
      PrototypeANRL_A2B_cP12_212_c_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 282 // ./aflow --proto=A3B3C_cI56_214_g_h_a --params=12.31504,0.108,0.384
    if(vlabel.at(ifound)=="A3B3C_cI56_214_g_h_a")
      PrototypeANRL_A3B3C_cI56_214_g_h_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 283 // ./aflow --proto=A3BC2_cI48_214_f_a_e --params=10.38,0.266,0.365
    if(vlabel.at(ifound)=="A3BC2_cI48_214_f_a_e")
      PrototypeANRL_A3BC2_cI48_214_f_a_e(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 284 // ./aflow --proto=A4B9_cP52_215_ei_3efgi --params=8.7068,0.1157,0.8296,0.3253,0.6066,0.3534,0.8549,0.8113,0.5332,0.3153,0.0322
    if(vlabel.at(ifound)=="A4B9_cP52_215_ei_3efgi")
      PrototypeANRL_A4B9_cP52_215_ei_3efgi(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 285 // ./aflow --proto=ABCD_cF16_216_c_d_b_a --params=6.465
    if(vlabel.at(ifound)=="ABCD_cF16_216_c_d_b_a")
      PrototypeANRL_ABCD_cF16_216_c_d_b_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 286 // ./aflow --proto=A3B4C_cP16_218_c_e_a --params=6.0258147002,0.6486
    if(vlabel.at(ifound)=="A3B4C_cP16_218_c_e_a")
      PrototypeANRL_A3B4C_cP16_218_c_e_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 287 // ./aflow --proto=A7BC3D13_cF192_219_de_b_c_ah --params=12.0986,0.0808,0.0987,0.0214,0.1821
    if(vlabel.at(ifound)=="A7BC3D13_cF192_219_de_b_c_ah")
      PrototypeANRL_A7BC3D13_cF192_219_de_b_c_ah(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 288 // ./aflow --proto=A15B4_cI76_220_ae_c --params=9.718,-0.042,0.12,0.16,-0.04
    if(vlabel.at(ifound)=="A15B4_cI76_220_ae_c")
      PrototypeANRL_A15B4_cI76_220_ae_c(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 289 // ./aflow --proto=A4B3_cI28_220_c_a --params=8.6,0.08333
    if(vlabel.at(ifound)=="A4B3_cI28_220_c_a")
      PrototypeANRL_A4B3_cI28_220_c_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 290 // ./aflow --proto=A2B3C6_cP33_221_cd_ag_fh --params=7.624,0.3,0.2,0.2
    if(vlabel.at(ifound)=="A2B3C6_cP33_221_cd_ag_fh")
      PrototypeANRL_A2B3C6_cP33_221_cd_ag_fh(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 291 // ./aflow --proto=A5B3C16_cP96_222_ce_d_fi --params=11.1199004528,0.0,0.125,0.084,0.166,0.625
    if(vlabel.at(ifound)=="A5B3C16_cP96_222_ce_d_fi")
      PrototypeANRL_A5B3C16_cP96_222_ce_d_fi(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 292 // ./aflow --proto=A23B6_cF116_225_bd2f_e --params=12.523,0.203,0.178,0.378
    if(vlabel.at(ifound)=="A23B6_cF116_225_bd2f_e")
      PrototypeANRL_A23B6_cF116_225_bd2f_e(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 293 // ./aflow --proto=A6B2C_cF36_225_e_c_a --params=9.725,0.24
    if(vlabel.at(ifound)=="A6B2C_cF36_225_e_c_a")
      PrototypeANRL_A6B2C_cF36_225_e_c_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 294 // ./aflow --proto=AB13_cF112_226_a_bi --params=12.2836,0.1806,0.1192
    if(vlabel.at(ifound)=="AB13_cF112_226_a_bi")
      PrototypeANRL_AB13_cF112_226_a_bi(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 295 // ./aflow --proto=A2B2C7_cF88_227_c_d_af --params=10.2663,0.4157
    if(vlabel.at(ifound)=="A2B2C7_cF88_227_c_d_af")
      PrototypeANRL_A2B2C7_cF88_227_c_d_af(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 296 // ./aflow --proto=A3B4_cF56_227_ad_e --params=8.0835,0.2642
    if(vlabel.at(ifound)=="A3B4_cF56_227_ad_e")
      PrototypeANRL_A3B4_cF56_227_ad_e(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 297 // ./aflow --proto=A5BCD6_cF416_228_eg_c_b_h --params=22.2649775014,0.19,0.425,0.05,0.18,0.28
    if(vlabel.at(ifound)=="A5BCD6_cF416_228_eg_c_b_h")
      PrototypeANRL_A5BCD6_cF416_228_eg_c_b_h(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 298 // ./aflow --proto=A6B_cF224_228_h_c --params=15.4800297791,0.043,0.138,0.278
    if(vlabel.at(ifound)=="A6B_cF224_228_h_c")
      PrototypeANRL_A6B_cF224_228_h_c(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 299 // ./aflow --proto=A3B10_cI52_229_e_fh --params=8.9822,0.3538,0.13505,0.3045
    if(vlabel.at(ifound)=="A3B10_cI52_229_e_fh")
      PrototypeANRL_A3B10_cI52_229_e_fh(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 300 // ./aflow --proto=A4B_cI10_229_c_a --params=6.186
    if(vlabel.at(ifound)=="A4B_cI10_229_c_a")
      PrototypeANRL_A4B_cI10_229_c_a(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 301 // ./aflow --proto=A7B3_cI40_229_df_e --params=8.735,0.342,0.156
    if(vlabel.at(ifound)=="A7B3_cI40_229_df_e")
      PrototypeANRL_A7B3_cI40_229_df_e(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    // ---------------------------------------------------------------------------
    // 302 // ./aflow --proto=A2B3C12D3_cI160_230_a_c_h_d --params=11.4597,0.3471,0.4664,0.0512
    if(vlabel.at(ifound)=="A2B3C12D3_cI160_230_a_c_h_d")
      PrototypeANRL_A2B3C12D3_cI160_230_a_c_h_d(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    

    // ---------------------------------------------------------------------------
    

    // DAVID PUT THEM HERE WITH THE ORDER OF THE PAPER

    //if(vlabel.at(ifound)=="YYY") {
    //  PrototypeANRL_YYY(oss,web,str,parameters,vproto.at(ifound),print_mode,LDEBUG);
    //}

    // after you generate them you should try them with aflowSG  and if they are small with the pearson

    // DX NOT NEEDED: Parameters are always in RHL if(XHOST.vflag_pflow.flag("PROTO::RHL")) cout << "DAVID WE GOT RHL"<< endl;
    if(XHOST.vflag_pflow.flag("PROTO::HEX") && vproto_Pearson_symbol[ifound][1] == 'R') {
      vector<double> vparameters;
      aurostd::string2tokens(parameters,vparameters,",");
      uint i=0;
      double a=vparameters.at(i++);
      double covera=vparameters.at(i++);
      double c=covera*a;
      str=rhl2hex(str,a,c);
    }

    for(uint iat=0;iat<str.atoms.size();iat++) {
      str.atoms.at(iat).name_is_given=TRUE;
      str.atoms.at(iat).number=iat;//iat;    // reference position for convasp
      str.atoms.at(iat).basis=iat;//iat;     // position in the basis
      if(print_mode!=1){ //equations only //DX 20180618
      str.atoms.at(iat).cpos=F2C(str.lattice,str.atoms.at(iat).fpos);
      }
      str.num_each_type.at(str.atoms.at(iat).type)++;
      //     str.comp_each_type.at(str.atoms.at(iat).type)+=1.0; inside code
      str.species.at(str.atoms.at(iat).type)=str.atoms.at(iat).name;	
    }
    
    // ---------------------------------------------------------------------------
    // DONE
    if(print_mode!=1){ //equations only //DX 20180618
    xvector<double> data(6);
    data=Getabc_angles(str.lattice,DEGREES);
    str.a=data[1];str.b=data[2];str.c=data[3];str.alpha=data[4];str.beta=data[5];str.gamma=data[6];
    clear(str.origin);
    }
    //  if(vpflow.flag("STDPRIMCELL")) {cout << "EUREKA"<< endl;} //cout << GetStandardPrimitive(xstructure(cin,IOAFLOW_AUTO));

    // ---------------------------------------------------------------------------
    // NOW PLAY WITH PERMUTATIONS and ATOMX
    if(vpermutation.size()>0 || vatomX.size()>0) {
      if(LDEBUG) { cerr << "anrl::PrototypeANRL: PERMUTATIONS" << endl;}
      if(LDEBUG) { cerr << "anrl::PrototypeANRL: vpermutation.size()=" << vpermutation.size() << endl;}
      if(LDEBUG) { cerr << "anrl::PrototypeANRL: vpermutation ="; for(uint i=0;i<vpermutation.size();i++) {cerr << " " << vpermutation.at(i);} cerr << endl;}     
      if(LDEBUG) { cerr << "anrl::PrototypeANRL: ATOMX" << endl;}
      if(LDEBUG) { cerr << "anrl::PrototypeANRL: vatomX.size()=" << vatomX.size() << endl;}
      if(LDEBUG) { cerr << "anrl::PrototypeANRL: vatomX ="; for(uint i=0;i<vatomX.size();i++) {cerr << " " << vatomX.at(i);} cerr << endl;}
      if(print_mode!=1){ //equations only //DX 20180618
        std::deque<_atom> atoms;
        atoms=str.atoms;
        // STRIP ALL ATOMS
        while(str.atoms.size()>0) { str.RemoveAtom(0); }
        // ADD MODIFIED ATOMS
        for(uint i=0;i<atoms.size();i++) {
	  uint type=atoms.at(i).type;
	  if(vpermutation.size()>0)  { atoms.at(i).type=vpermutation.at(type); }  // PERMUTATIONS 
	  if(vpermutation.size()>0 || vatomX.size()>0) { atoms.at(i).name=vatomX.at(atoms.at(i).type); }  // PERMUTATIONS AND ATOMX
	  //	atoms.at(i).name=aurostd::mod(label_permutations.at(type)-65,32)+65;
	  str.AddAtom(atoms.at(i));
        }
      }
      str.SpeciesPutAlphabetic();
    }
    
    if(XHOST.vflag_control.flag("WWW")) {
      cout << web.str() << endl;
      exit(0);
    }

    return str;
  }
}


// ***************************************************************************
namespace anrl {
  bool PrototypeANRL_Consistency(ostream &oss,
				 uint vparameters_size,
				 uint nparameters,string prototype,string label,string Strukturbericht,
				 string Pearson_symbol,uint spacegroup,string params,uint print_mode) { //DX 20180710 - added print_mode info
    if(vparameters_size!=nparameters && print_mode!=1) {                                                //DX 20180710 - if equations only (print_mode==1), we do not need the parameters
      oss << "anrl::PrototypeANRL" << endl;
      oss << " Prototype                   : " << prototype << endl;
      oss << " AFLOW prototype label       : " << label << endl;
      oss << " Strukturbericht Designation : " << Strukturbericht << endl;
      oss << " Pearson Symbol              : " << Pearson_symbol << endl;
      oss << " Space group number          : " << GetSpaceGroupName(spacegroup) << endl;
      oss << " Space group symbol          : " << spacegroup << endl;
      oss << " AFLOW prototype command     : aflow --proto=" << label << endl;
      oss << "                                     --params=" << params << endl;
      return FALSE;
    }
    return TRUE;
  }
}

// *************************************************************************** 
// the mother of all the AFLOW-NRL prototypes
namespace anrl {
  bool vproto2tokens(string proto_line,
		     string& label,uint& nspecies,uint& natoms,uint& spacegroup,uint& nunderscores,uint& nparameters,
		     string& Pearson_symbol,string& params,string& Strukturbericht,string& prototype,string& dialect) { 

    vector<string> tokens;
    uint j=0;
    
    aurostd::string2tokens(proto_line,tokens,";");
    label=tokens.at(j++);
    nspecies=aurostd::string2utype<uint>(tokens.at(j++));
    natoms=aurostd::string2utype<uint>(tokens.at(j++));
    spacegroup=aurostd::string2utype<uint>(tokens.at(j++));
    nunderscores=aurostd::string2utype<uint>(tokens.at(j++));
    nparameters=aurostd::string2utype<uint>(tokens.at(j++));
    Pearson_symbol=tokens.at(j++);
    params=tokens.at(j++);
    Strukturbericht=tokens.at(j++);
    prototype=tokens.at(j++);
    dialect=tokens.at(j++);

    return TRUE;
  }
}

// *************************************************************************** 
// RHL to HEX transformation
namespace anrl {
  xstructure rhl2hex(xstructure& str, double& a, double& c) { 

    xstructure hex_str; //make new xstructure object
    hex_str=str;

    hex_str.atoms.clear();

    xvector<double> xn(3);   xn(1)=1.0;xn(2)=0.0;xn(3)=0.0;
    xvector<double> yn(3);   yn(1)=0.0;yn(2)=1.0;yn(3)=0.0;
    xvector<double> zn(3);   zn(1)=0.0;zn(2)=0.0;zn(3)=1.0;
    xvector<double> a1(3),a2(3),a3(3);

    xmatrix<double> rhl_lattice, hex_lattice, rtransf, htransf;
     
    //HEX lattice
    a1=(1.0/2.0)*a*xn-(sqrt(3.0)/2.0)*a*yn;
    a2=(1.0/2.0)*a*xn+(sqrt(3.0)/2.0)*a*yn;
    a3=c*zn;
    hex_str.lattice(1,1)=a1(1);hex_str.lattice(1,2)=a1(2);hex_str.lattice(1,3)=a1(3);
    hex_str.lattice(2,1)=a2(1);hex_str.lattice(2,2)=a2(2);hex_str.lattice(2,3)=a2(3);
    hex_str.lattice(3,1)=a3(1);hex_str.lattice(3,2)=a3(2);hex_str.lattice(3,3)=a3(3);

    hex_str.FixLattices(); // Reciprocal/f2c/c2f

    //RHL Transformation matrix
    rtransf(1,1)=(1.0/2.0);rtransf(1,2)=-(1.0/(2.0*sqrt(3.0)));rtransf(1,3)=(1.0/3.0);
    rtransf(2,1)=0.0;rtransf(2,2)=(1.0/sqrt(3.0));rtransf(2,3)=(1.0/3.0);
    rtransf(3,1)=-(1.0/2.0);rtransf(3,2)=-(1.0/(2.0*sqrt(3.0)));rtransf(3,3)=(1.0/3.0);

    //HEX Transformtion matrix
    htransf(1,1)=(1.0/2.0);htransf(1,2)=-(sqrt(3.0)/2.0);htransf(1,3)=0.0;
    htransf(2,1)=(1.0/2.0);htransf(2,2)=(sqrt(3.0)/2.0);htransf(2,3)=0.0;
    htransf(3,1)=0.0;htransf(3,2)=0.0;htransf(3,3)=1.0;

    xvector<double> c1(3); c1(1)=(2.0/3.0); c1(2)=(1.0/3.0); c1(3)=(1.0/3.0); //centering translation
    xvector<double> c2(3); c2(1)=(1.0/3.0); c2(2)=(2.0/3.0); c2(3)=(2.0/3.0); //centering translation

    for(uint a=0;a<str.atoms.size();a++) {
      _atom tmp;
      tmp.name=str.atoms[a].name;
      tmp.type=str.atoms[a].type;
      tmp.basis=str.atoms[a].basis;
      xvector<double> center_pos;
      center_pos=trasp(inverse(htransf))*(trasp(rtransf)*str.atoms[a].fpos); // Method for transforming RHL to HEX
      tmp.fpos=center_pos;
      hex_str.comp_each_type.at(tmp.type)+=1.0;                        
      hex_str.atoms.push_back(tmp);
      //add centering c1
      tmp.fpos=center_pos+c1; 
      hex_str.comp_each_type.at(tmp.type)+=1.0;                        
      hex_str.atoms.push_back(tmp);
      //add centering c2
      tmp.fpos=center_pos+c2; 
      hex_str.comp_each_type.at(tmp.type)+=1.0;                        
      hex_str.atoms.push_back(tmp);
    }
    return hex_str;
  }
}

#endif // _AFLOW_ANRL_CPP

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *                                                                         *
// ***************************************************************************
// Written by Stefano Curtarolo - 2016
