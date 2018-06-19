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
  uint PrototypeANRL_AB2_aP12_1_4a_8a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 1
  uint PrototypeANRL_ABC2_aP16_1_4a_4a_8a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 2
  uint PrototypeANRL_A2B_aP6_2_2i_i(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 3
  uint PrototypeANRL_A_aP4_2_aci(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 4
  uint PrototypeANRL_A2B_mP12_3_bc3e_2e(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 5
  uint PrototypeANRL_A_mP4_4_2a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 6
  uint PrototypeANRL_A_mC12_5_3c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 7
  uint PrototypeANRL_A3BC_mC10_8_ab_a_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 8
  uint PrototypeANRL_A2B_mC144_9_24a_12a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 9
  uint PrototypeANRL_AB_mP4_11_e_e(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 10
  uint PrototypeANRL_ABC3_mP10_11_e_e_ef(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 11
  uint PrototypeANRL_A_mP16_11_8e(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 12
  uint PrototypeANRL_AB2_mC6_12_a_i(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 13
  uint PrototypeANRL_A_mC34_12_ah3i2j(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 14
  uint PrototypeANRL_AB3_mC16_12_g_ij(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 15
  uint PrototypeANRL_A5B2_mC14_12_a2i_i(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 16
  uint PrototypeANRL_A_mC4_12_i(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 17
  uint PrototypeANRL_ABC4_mP12_13_e_a_2g(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 18
  uint PrototypeANRL_A_mP84_13_21g(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 19
  uint PrototypeANRL_A2B_mP12_14_2e_e(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 20
  uint PrototypeANRL_A_mP32_14_8e(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 21
  uint PrototypeANRL_A_mP64_14_16e(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 22
  uint PrototypeANRL_A2B5_mC28_15_f_e2f(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 23
  uint PrototypeANRL_AB_mC8_15_c_e(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 24
  uint PrototypeANRL_A2B_mC48_15_ae3f_2f(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 25
  uint PrototypeANRL_ABC6D2_mC40_15_e_e_3f_f(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 26
  uint PrototypeANRL_ABC4_oP12_16_ag_cd_2u(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 27
  uint PrototypeANRL_AB3_oP16_18_ab_3c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 28
  uint PrototypeANRL_A2B_oP12_19_2a_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 29
  uint PrototypeANRL_A2B_oC24_20_abc_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 30
  uint PrototypeANRL_AB_oP2_25_b_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 31
  uint PrototypeANRL_AB2_oP24_28_acd_2c3d(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 32
  uint PrototypeANRL_AB3C4_oP16_31_a_ab_2ab(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 33
  uint PrototypeANRL_AB_oP8_33_a_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 34
  uint PrototypeANRL_AB3C4_oP32_33_a_3a_4a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 35
  uint PrototypeANRL_A2B_oC12_36_2a_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 36
  uint PrototypeANRL_A2BC_oC8_38_e_a_b(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 37
  uint PrototypeANRL_A2B_oC12_38_de_ab(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 38
  uint PrototypeANRL_AB4_oC20_41_a_2b(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 39
  uint PrototypeANRL_AB2_oC24_41_2a_2b(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 40
  uint PrototypeANRL_AB2_oF72_43_ab_3b(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 41
  uint PrototypeANRL_AB_oI4_44_a_b(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 42
  uint PrototypeANRL_A2B3C7D_oP13_47_t_aq_eqrs_h(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 43
  uint PrototypeANRL_AB_oP4_51_e_f(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 44
  uint PrototypeANRL_A3B2_oP20_56_ce_e(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 45
  uint PrototypeANRL_ABCD_oP16_57_d_c_d_d(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 46
  uint PrototypeANRL_AB_oP8_57_d_d(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 47
  uint PrototypeANRL_AB2_oP6_58_a_g(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 48 //49 //50
  uint PrototypeANRL_AB_oP4_59_a_b(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 51
  uint PrototypeANRL_ABC_oP6_59_a_a_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 52
  uint PrototypeANRL_A3B_oP8_59_bf_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 53
  uint PrototypeANRL_AB_oP16_61_c_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 54
  uint PrototypeANRL_A2B_oP24_61_2c_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 55
  uint PrototypeANRL_A3B2_oP20_62_3c_2c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 56
  uint PrototypeANRL_AB3C_oP20_62_c_cd_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 57
  uint PrototypeANRL_A4B_oP20_62_2cd_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 58
  uint PrototypeANRL_AB2C_oP16_62_c_2c_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 59
  uint PrototypeANRL_A2B_oP12_62_2c_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 60 //61 //62
  uint PrototypeANRL_AB_oP8_62_c_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 63 // 64 // 68 // 69
  uint PrototypeANRL_AB3_oP16_62_c_cd(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 65
  uint PrototypeANRL_A3B7_oP40_62_cd_3c2d(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 66
  uint PrototypeANRL_A_oP8_62_2c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 67
  uint PrototypeANRL_AB2C_oC16_63_c_2c_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 70
  uint PrototypeANRL_A2B_oC12_63_2c_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 71
  uint PrototypeANRL_AB_oC8_63_c_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 72
  uint PrototypeANRL_A_oC4_63_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 73
  uint PrototypeANRL_A_oC8_64_f(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 74 // 76 // 77
  uint PrototypeANRL_A2B2C_oC80_64_efg_efg_df(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 75
  uint PrototypeANRL_AB_oC8_65_j_g(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 78
  uint PrototypeANRL_A3B5_oC16_65_ah_bej(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 79
  uint PrototypeANRL_AB3_oC8_65_a_bf(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 80
  uint PrototypeANRL_AB_oF8_69_a_b(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 81
  uint PrototypeANRL_A_oF8_70_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 82
  uint PrototypeANRL_A2B_oF24_70_e_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 83
  uint PrototypeANRL_A_oF128_70_4h(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 84
  uint PrototypeANRL_AB2_oI6_71_a_i(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 85
  uint PrototypeANRL_AB2_oI6_71_a_g(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 86
  uint PrototypeANRL_A2B_oI12_72_j_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 87
  uint PrototypeANRL_AB4C_tI12_82_c_g_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 88
  uint PrototypeANRL_A2BC4_tI14_82_bc_a_g(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 89
  uint PrototypeANRL_AB_tP16_84_cej_k(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 90
  uint PrototypeANRL_A4B5_tI18_87_h_ah(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 91
  uint PrototypeANRL_AB4_tI10_87_a_h(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 92
  uint PrototypeANRL_A2B_tP12_92_b_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 93
  uint PrototypeANRL_A2B_tP36_96_3b_ab(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 94
  uint PrototypeANRL_A_tP12_96_ab(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 95
  uint PrototypeANRL_A3BC_tP5_99_bc_a_b(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 96
  uint PrototypeANRL_AB3_tP8_113_a_ce(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 97
  uint PrototypeANRL_A2BC4D_tI16_121_d_a_i_b(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 98
  uint PrototypeANRL_ABC2_tI16_122_a_b_d(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 99
  uint PrototypeANRL_AB5C_tP7_123_b_ci_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 100
  uint PrototypeANRL_AB3_tP4_123_a_ce(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 101
  uint PrototypeANRL_AB_tP2_123_a_d(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 102
  uint PrototypeANRL_ABC2_tP4_123_d_a_f(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 103
  uint PrototypeANRL_A2B3_tP10_127_g_ah(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 104
  uint PrototypeANRL_ABCD_tP8_129_c_b_a_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 105
  uint PrototypeANRL_A_tP4_129_ac(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 106
  uint PrototypeANRL_ABC_tP6_129_c_a_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 107
  uint PrototypeANRL_A2B_tP6_129_ac_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 108
  uint PrototypeANRL_AB_tP4_129_a_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 109
  uint PrototypeANRL_AB_tP4_129_c_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 110
  uint PrototypeANRL_AB_tP4_131_c_e(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 111
  uint PrototypeANRL_A_tP50_134_b2m2n(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 112
  uint PrototypeANRL_A_tP30_136_bf2ij(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 113
  uint PrototypeANRL_AB_tP8_136_g_f(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 114
  uint PrototypeANRL_A2B_tP6_136_f_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 115
  uint PrototypeANRL_sigma_tP30_136_bf2ij(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 116
  uint PrototypeANRL_A_tP4_136_f(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 117
  uint PrototypeANRL_A_tP16_138_j(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 118
  uint PrototypeANRL_A3B_tI16_139_cde_e(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 119
  uint PrototypeANRL_A_tI4_139_e(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 120
  uint PrototypeANRL_AB2C4_tI14_139_a_e_ce(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 121
  uint PrototypeANRL_A12B_tI26_139_fij_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 122
  uint PrototypeANRL_A_tI2_139_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 123 //131
  uint PrototypeANRL_A_tI8_139_h(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 124
  uint PrototypeANRL_A3B_tI8_139_bd_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 125
  uint PrototypeANRL_AB2_tI6_139_a_e(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 126
  uint PrototypeANRL_A4B5_tI18_139_i_ah(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 127
  uint PrototypeANRL_A4B_tI10_139_de_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 128
  uint PrototypeANRL_A8B_tI18_139_hi_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 129
  uint PrototypeANRL_A2B_tI6_139_d_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 130
  uint PrototypeANRL_A2B_tI12_140_h_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 132
  uint PrototypeANRL_AB3_tI16_140_b_ah(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 133
  uint PrototypeANRL_AB_tI16_140_ab_h(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 134
  uint PrototypeANRL_A4BC_tI24_141_h_b_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 135
  uint PrototypeANRL_A_tI4_141_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 136
  uint PrototypeANRL_A3B4_tI28_141_ad_h(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 137
  uint PrototypeANRL_A2B_tI12_141_e_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 138
  uint PrototypeANRL_AB_tI16_141_e_e(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 139
  uint PrototypeANRL_A2B_tI24_141_2e_e(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 140
  uint PrototypeANRL_AB_tI8_141_a_b(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 141
  uint PrototypeANRL_A2B3_tI80_141_ceh_3h(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 142
  uint PrototypeANRL_ABC4_tI96_142_e_ab_2g(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 143
  uint PrototypeANRL_A2B_hP9_147_g_ad(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 144
  uint PrototypeANRL_AB_hR16_148_cf_cf(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 145
  uint PrototypeANRL_AB3_hR8_148_c_f(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 146
  uint PrototypeANRL_AB_hR26_148_b2f_a2f(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 147
  uint PrototypeANRL_AB3C_hR10_148_c_f_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 148
  uint PrototypeANRL_A2B_hP9_150_ef_bd(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 149
  uint PrototypeANRL_A3B_hP24_151_3c_2a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 150
  uint PrototypeANRL_A2B_hP9_152_c_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 151
  uint PrototypeANRL_A_hP3_152_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 152
  uint PrototypeANRL_AB_hP6_154_a_b(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 153
  uint PrototypeANRL_AB3_hR8_155_c_de(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 154
  uint PrototypeANRL_A3B2_hR5_155_e_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 155
  uint PrototypeANRL_AB_hR6_160_b_b(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 156
  uint PrototypeANRL_AB_hR6_160_3a_3a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 157
  uint PrototypeANRL_ABC3_hR10_161_a_a_b(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 158
  uint PrototypeANRL_AB2_hP9_162_ad_k(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 159
  uint PrototypeANRL_AB2CD2_hP36_163_h_i_bf_i(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 160
  uint PrototypeANRL_A3B2_hP5_164_ad_d(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 161
  uint PrototypeANRL_AB2_hP3_164_a_d(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 162
  uint PrototypeANRL_A3B_hP24_165_adg_f(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 163
  uint PrototypeANRL_AB_hR2_166_a_b(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 164
  uint PrototypeANRL_A_hR2_166_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 165 // 172 // 175
  uint PrototypeANRL_A_hR1_166_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 166 // 170
  uint PrototypeANRL_A7B6_hR13_166_ah_3c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 167
  uint PrototypeANRL_A_hR3_166_ac(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 168
  uint PrototypeANRL_A2B3_hR5_166_c_ac(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 169
  uint PrototypeANRL_A5B2_hR7_166_a2c_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 171
  uint PrototypeANRL_A_hR12_166_2h(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 173
  uint PrototypeANRL_ABC2_hR4_166_a_b_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 174
  uint PrototypeANRL_A_hR105_166_bc9h4i(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 176
  uint PrototypeANRL_A6B_hR7_166_g_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 177
  uint PrototypeANRL_ABC3_hR10_167_a_b_e(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 178 // 179
  uint PrototypeANRL_A2B3_hR10_167_c_e(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 180
  uint PrototypeANRL_A2B_hP18_180_fi_bd(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 181
  uint PrototypeANRL_AB2_hP9_180_d_j(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 182
  uint PrototypeANRL_A2B_hP9_180_j_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 183
  uint PrototypeANRL_AB3_hP8_182_c_g(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 184
  uint PrototypeANRL_A_hP4_186_ab(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 185
  uint PrototypeANRL_AB_hP8_186_ab_ab(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 186
  uint PrototypeANRL_AB_hP4_186_b_b(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 187
  uint PrototypeANRL_AB_hP12_186_a2b_a2b(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 188
  uint PrototypeANRL_A5B3C_hP18_186_2a3b_2ab_b(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 189
  uint PrototypeANRL_AB_hP4_186_b_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 190
  uint PrototypeANRL_ABC_hP3_187_a_d_f(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 191
  uint PrototypeANRL_AB_hP2_187_d_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 192
  uint PrototypeANRL_A2B_hP9_189_fg_bc(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 193
  uint PrototypeANRL_AB4C_hP6_191_a_h_b(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 194
  uint PrototypeANRL_AB5_hP6_191_a_cg(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 195
  uint PrototypeANRL_A_hP1_191_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 196
  uint PrototypeANRL_A3B_hP4_191_bc_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 197
  uint PrototypeANRL_AB2_hP3_191_a_d(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 198
  uint PrototypeANRL_A2B_hP6_191_h_e(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 199
  uint PrototypeANRL_AB_hP6_191_f_ad(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 200
  uint PrototypeANRL_AB_hP8_194_ad_f(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 201
  uint PrototypeANRL_A_hP6_194_h(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 202
  uint PrototypeANRL_AB_hP12_194_af_bf(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 203
  uint PrototypeANRL_A_hP4_194_ac(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 204
  uint PrototypeANRL_AB3_hP8_194_c_bf(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 205
  uint PrototypeANRL_AB2_hP6_194_b_f(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 206
  uint PrototypeANRL_AB_hP4_194_c_d(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 207
  uint PrototypeANRL_ABC2_hP8_194_d_a_f(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 208
  uint PrototypeANRL_A3B_hP8_194_h_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 209
  uint PrototypeANRL_A_hP4_194_bc(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 210
  uint PrototypeANRL_AB2_hP6_194_c_f(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 211
  uint PrototypeANRL_A5B2_hP14_194_abdf_f(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 212
  uint PrototypeANRL_AB2_hP12_194_f_ah(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 213
  uint PrototypeANRL_ABC_hP6_194_c_d_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 214
  uint PrototypeANRL_A_hP4_194_f(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 215
  uint PrototypeANRL_AB2_hP6_194_c_ad(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 216
  uint PrototypeANRL_AB3C4_hP16_194_c_af_ef(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 217
  uint PrototypeANRL_A_hP2_194_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 218
  uint PrototypeANRL_AB2_hP24_194_ef_fgh(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 219
  uint PrototypeANRL_AB_hP12_194_df_ce(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 220
  uint PrototypeANRL_AB_hP4_194_c_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 221
  uint PrototypeANRL_A2B_hP12_194_cg_f(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 222
  uint PrototypeANRL_A4B_cI40_197_cde_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 223
  uint PrototypeANRL_ABC_cP12_198_a_a_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 224
  uint PrototypeANRL_A3B_cP16_198_b_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 225
  uint PrototypeANRL_A_cP8_198_2a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 226
  uint PrototypeANRL_AB_cP8_198_a_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 227 // 228
  uint PrototypeANRL_AB_cI16_199_a_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 229
  uint PrototypeANRL_AB32C48_cI162_204_a_2efg_2gh(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 230
  uint PrototypeANRL_A3B_cI32_204_g_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 231
  uint PrototypeANRL_A12B_cI26_204_g_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 232
  uint PrototypeANRL_A_cP8_205_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 233
  uint PrototypeANRL_AB_cP16_205_c_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 234
  uint PrototypeANRL_AB2_cP12_205_a_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 235
  uint PrototypeANRL_AB3C6_cI80_206_a_d_e(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 236
  uint PrototypeANRL_A_cI16_206_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 237
  uint PrototypeANRL_A_cP20_213_cd(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 238
  uint PrototypeANRL_A3B4C_cP8_215_d_e_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 239
  uint PrototypeANRL_AB4_cP5_215_a_e(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 240
  uint PrototypeANRL_AB3C4_cP8_215_a_c_e(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 241
  uint PrototypeANRL_AB5_cF24_216_a_ce(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 242
  uint PrototypeANRL_ABC_cF12_216_b_c_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 243
  uint PrototypeANRL_AB_cF8_216_c_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 244
  uint PrototypeANRL_A4B_cI10_217_c_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 245
  uint PrototypeANRL_A_cI58_217_ac2g(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 246
  uint PrototypeANRL_A5B8_cI52_217_ce_cg(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 247
  uint PrototypeANRL_A_cI16_220_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 248
  uint PrototypeANRL_A3B2_cI40_220_d_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 249
  uint PrototypeANRL_AB_cP2_221_b_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 250
  uint PrototypeANRL_AB_cP6_221_c_d(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 251
  uint PrototypeANRL_AB3C_cP5_221_a_c_b(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 252
  uint PrototypeANRL_AB27CD3_cP32_221_a_dij_b_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 253
  uint PrototypeANRL_AB3_cP4_221_a_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 254
  uint PrototypeANRL_A_cP1_221_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 255
  uint PrototypeANRL_AB11_cP36_221_c_agij(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 256
  uint PrototypeANRL_AB11CD3_cP16_221_a_dg_b_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 257
  uint PrototypeANRL_A3B_cP4_221_d_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 258
  uint PrototypeANRL_A6B_cP7_221_f_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 259
  uint PrototypeANRL_A3B_cP8_223_c_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 260
  uint PrototypeANRL_A_cP46_223_dik(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 261
  uint PrototypeANRL_A2B_cP6_224_b_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 262
  uint PrototypeANRL_A7B_cF32_225_bd_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 263
  uint PrototypeANRL_AB3_cF16_225_a_bc(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 264
  uint PrototypeANRL_A9B16C7_cF128_225_acd_2f_be(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 265
  uint PrototypeANRL_A12B_cF52_225_i_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 266
  uint PrototypeANRL_AB2_cF12_225_a_c(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 267
  uint PrototypeANRL_A6B23_cF116_225_e_acfh(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 268
  uint PrototypeANRL_AB2C_cF16_225_a_c_b(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 269
  uint PrototypeANRL_A_cF4_225_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 270
  uint PrototypeANRL_AB18C8_cF108_225_a_eh_f(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 271
  uint PrototypeANRL_AB_cF8_225_a_b(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 272
  uint PrototypeANRL_A2B_cF24_227_c_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 273
  uint PrototypeANRL_AB2_cF96_227_e_cf(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 274
  uint PrototypeANRL_AB_cF16_227_a_b(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 275
  uint PrototypeANRL_A_cF136_227_aeg(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 276
  uint PrototypeANRL_A2B_cF24_227_d_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 277
  uint PrototypeANRL_A_cF8_227_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 278
  uint PrototypeANRL_A2BC4_cF56_227_d_a_e(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 279
  uint PrototypeANRL_AB2_cF48_227_c_e(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 280
  uint PrototypeANRL_AB3C3_cF112_227_c_de_f(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 281
  uint PrototypeANRL_A_cI2_229_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 282
  uint PrototypeANRL_A3B_cI8_229_b_a(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 283
  uint PrototypeANRL_A4B3_cI14_229_c_b(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 284
  uint PrototypeANRL_A2B7_cI54_229_e_afh(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 285
  uint PrototypeANRL_AB12C3_cI32_229_a_h_b(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 286
  uint PrototypeANRL_AB4C3_cI16_229_a_c_b(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 287
  uint PrototypeANRL_A4B3_cI112_230_af_g(ostream &oss,stringstream &web,xstructure& str,string parameters,string proto_line,bool LDEBUG); // 288
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
    // 1 // ./aflow --proto=AB2_aP12_1_4a_8a --params=5.417,1.0,1.0,90.0,90.0,90.0,0.001,0.002,0.003,0.4966,0.0001,0.5036,0.5001,0.502,0.0011,-0.0006,0.5013,0.5038,0.3857,0.3832,0.384,0.1149,0.6114,0.8846,0.8854,0.1157,0.6143,0.6153,0.8865,0.1141,0.6151,0.6132,0.6137,0.8854,0.3818,0.1149,0.1147,0.8856,0.3841,0.3857,0.1161,0.8842
    if(vlabel.at(ifound)=="AB2_aP12_1_4a_8a") {
      PrototypeANRL_AB2_aP12_1_4a_8a(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 2 // ./aflow --proto=ABC2_aP16_1_4a_4a_8a --params=6.554,1.00061031431,1.92662496186,100.43475,100.46074,107.53,0.3267,0.582,0.177,0.565,-0.0132,0.4424,0.5217,0.3883,0.6767,-0.0744,0.6254,-0.0574,0.0338,0.0476,0.2599,0.0831,0.6072,0.4974,-0.0131,0.0949,0.7583,0.5449,0.1443,-0.0022,-0.0211,0.5213,0.2073,0.2907,0.5956,-0.0183,-0.0616,0.0602,0.4998,0.5068,-0.0175,0.2448,0.4596,0.0397,0.708,0.5326,0.352,0.4818,0.0,0.0,0.0,-0.078,0.569,0.7448
    if(vlabel.at(ifound)=="ABC2_aP16_1_4a_4a_8a") {
      PrototypeANRL_ABC2_aP16_1_4a_4a_8a(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 3 // ./aflow --proto=A2B_aP6_2_2i_i --params=4.56,1.54824561404,1.62280701754,80.2,106.96667,98.2,0.557,0.73,0.165,0.82,0.803,0.695,0.397,0.639,0.463
    if(vlabel.at(ifound)=="A2B_aP6_2_2i_i") {
      PrototypeANRL_A2B_aP6_2_2i_i(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 4 // ./aflow --proto=A_aP4_2_aci --params=3.307,2.24130631993,0.844572119746,89.06,85.15,85.7,0.572,0.259,0.433
    if(vlabel.at(ifound)=="A_aP4_2_aci") {
      PrototypeANRL_A_aP4_2_aci(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 5 // ./aflow --proto=A2B_mP12_3_bc3e_2e --params=4.1605,0.992524936907,1.78370388174,101.3752,0.15907,0.73859,0.02399,0.752,0.18927,0.38562,0.71473,0.64074,0.48963,0.20196,0.18802,0.18244,0.0,0.69651,0.38098,0.58564,0.17797
    if(vlabel.at(ifound)=="A2B_mP12_3_bc3e_2e") {
      PrototypeANRL_A2B_mP12_3_bc3e_2e(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 6 // ./aflow --proto=A_mP4_4_2a --params=3.104,2.42042525773,1.53350515464,92.71,0.25,0.23,0.48,0.48,0.0,0.02
    if(vlabel.at(ifound)=="A_mP4_4_2a") {
      PrototypeANRL_A_mP4_4_2a(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 7 // ./aflow --proto=A_mC12_5_3c --params=7.42,0.578167115903,1.90026954178,92.0,0.05,0.27,0.245,0.63,0.3,0.4,0.245,0.43,0.07
    if(vlabel.at(ifound)=="A_mC12_5_3c") {
      PrototypeANRL_A_mC12_5_3c(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 8 // ./aflow --proto=A3BC_mC10_8_ab_a_a --params=5.72204,0.9978207073,0.722908263486,90.498,0.5515,-0.0994,0.0,0.0,0.523,0.4492,0.288,0.2434,0.3729
    if(vlabel.at(ifound)=="A3BC_mC10_8_ab_a_a") {
      PrototypeANRL_A3BC_mC10_8_ab_a_a(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 9 // ./aflow --proto=A2B_mC144_9_24a_12a --params=18.524,0.270092852516,1.28535953358,105.82,0.5749,0.351,0.8182,0.0707,0.34,0.8476,0.7315,0.138,0.4851,0.2509,0.144,0.5152,0.4155,0.352,0.6741,-0.0873,0.352,0.6434,0.8773,0.164,-0.0787,0.416,0.168,-0.0639,0.7741,0.145,0.7538,0.2336,0.143,0.7402,0.6195,0.341,0.5847,0.0811,0.343,0.5661,-0.0034,0.011,0.6062,0.3533,0.489,0.5665,0.6498,0.005,0.6711,0.1524,0.496,0.7805,0.8636,0.499,0.7328,0.3361,0.003,0.8333,0.0052,0.493,0.7398,0.1369,0.011,-0.0732,0.4927,0.492,0.8868,0.5,0.468,0.5,0.2252,0.491,0.5898,0.2744,0.021,-0.0845,0.0507,0.041,0.5642,0.2036,0.447,0.7347,-0.0802,0.049,0.6225,0.5751,0.043,0.7955,0.4247,0.048,0.6971,0.2643,0.444,0.5386,0.8023,0.449,0.7661,0.6453,0.041,0.6027,0.8531,0.463,-0.0984,0.4493,0.466,-0.0642,0.2244,0.059,-0.0395,0.0697,0.049,0.8702
    if(vlabel.at(ifound)=="A2B_mC144_9_24a_12a") {
      PrototypeANRL_A2B_mC144_9_24a_12a(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 10 // ./aflow --proto=AB_mP4_11_e_e --params=2.8837,1.42393452856,1.61854561848,82.062,0.0387,0.8252,0.5887,0.7184
    if(vlabel.at(ifound)=="AB_mP4_11_e_e") {
      PrototypeANRL_AB_mP4_11_e_e(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 11 // ./aflow --proto=ABC3_mP10_11_e_e_ef --params=4.63,1.20259179266,1.52203023758,110.21,0.121,0.1745,0.3531,0.7086,0.4009,0.1165,0.8544,0.5361,0.6943
    if(vlabel.at(ifound)=="ABC3_mP10_11_e_e_ef") {
      PrototypeANRL_ABC3_mP10_11_e_e_ef(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 12 // ./aflow --proto=A_mP16_11_8e --params=6.183,0.779880316998,1.77308749798,101.79,0.345,0.162,0.767,0.168,0.128,0.34,0.657,0.457,0.025,0.618,0.473,0.653,0.328,-0.074,0.869,0.894
    if(vlabel.at(ifound)=="A_mP16_11_8e") {
      PrototypeANRL_A_mP16_11_8e(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 13 // ./aflow --proto=AB2_mC6_12_a_i --params=7.189,0.613019891501,0.705105021561,90.04,0.6879,0.2889
    if(vlabel.at(ifound)=="AB2_mC6_12_a_i") {
      PrototypeANRL_AB2_mC6_12_a_i(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 14 // ./aflow --proto=A_mC34_12_ah3i2j --params=11.93871,0.876392843113,0.658278825769,129.00411,0.22,0.854,0.241,0.663,0.745,0.566,0.238,0.355,0.232,-0.037,0.333,0.35,0.586
    if(vlabel.at(ifound)=="A_mC34_12_ah3i2j") {
      PrototypeANRL_A_mC34_12_ah3i2j(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 15 // ./aflow --proto=AB3_mC16_12_g_ij --params=5.914,1.73047007102,1.03956712885,108.25,0.1662,0.2147,0.2263,0.2518,0.32131,0.2248
    if(vlabel.at(ifound)=="AB3_mC16_12_g_ij") {
      PrototypeANRL_AB3_mC16_12_g_ij(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 16 // ./aflow --proto=A5B2_mC14_12_a2i_i --params=9.188,0.430343926861,0.705158902917,97.56,0.14286,0.42857,0.28571,0.85714,0.42857,0.28571
    if(vlabel.at(ifound)=="A5B2_mC14_12_a2i_i") {
      PrototypeANRL_A5B2_mC14_12_a2i_i(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 17 // ./aflow --proto=A_mC4_12_i --params=5.403,0.635387747548,0.940033314825,132.32,0.106,0.173
    if(vlabel.at(ifound)=="A_mC4_12_i") {
      PrototypeANRL_A_mC4_12_i(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 18 // ./aflow --proto=ABC4_mP12_13_e_a_2g --params=8.95,0.500335195531,1.63360893855,145.35,0.5182,0.2986,0.0278,0.0003,0.2821,0.4045,0.2366
    if(vlabel.at(ifound)=="ABC4_mP12_13_e_a_2g") {
      PrototypeANRL_ABC4_mP12_13_e_a_2g(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 19 // ./aflow --proto=A_mP84_13_21g --params=9.21,0.99348534202,2.45385450597,106.1,0.30089,0.20127,0.18147,0.17387,0.03262,0.11695,0.05014,-0.05231,0.18035,-0.07589,0.78099,0.11634,0.79463,0.67872,0.1738,0.68463,0.51532,0.10402,0.56601,0.44932,0.17224,0.42424,0.27741,0.11672,0.0412,0.39067,0.07245,-0.00092,0.15881,0.04497,0.78847,0.13878,0.07346,0.7486,-0.09081,0.04464,0.53574,0.87264,0.06842,0.50833,0.63715,0.03304,0.30515,0.63715,0.06617,0.25041,0.40555,0.0442,0.146,0.38905,0.17219,0.86038,0.10055,0.17357,0.59606,0.82384,0.1694,0.41856,0.64581,0.16732,-0.05418,0.32296,0.2006
    if(vlabel.at(ifound)=="A_mP84_13_21g") {
      PrototypeANRL_A_mP84_13_21g(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 20 // ./aflow --proto=A2B_mP12_14_2e_e --params=5.1505,1.01186292593,1.03238520532,99.23,0.07,0.3317,0.3447,0.4496,0.7569,0.4792,0.2754,0.0395,0.2083
    if(vlabel.at(ifound)=="A2B_mP12_14_2e_e") {
      PrototypeANRL_A2B_mP12_14_2e_e(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 21 // ./aflow --proto=A_mP32_14_8e --params=9.31,0.866809881847,1.38023630505,93.13333,0.437,0.185,0.084,0.246,0.273,-0.023,0.24,0.102,0.828,0.05,-0.08,0.852,0.157,0.669,-0.09,0.142,0.66,0.09,0.368,0.746,0.16,0.334,0.021,0.21
    if(vlabel.at(ifound)=="A_mP32_14_8e") {
      PrototypeANRL_A_mP32_14_8e(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 22 // ./aflow --proto=A_mP64_14_16e --params=15.018,0.979691037422,0.585231056066,93.61,0.18313,0.14063,0.03451,0.22856,0.28408,0.12262,0.35548,0.31907,-0.00548,0.47826,0.28776,0.16131,0.52853,0.14438,0.09345,0.47966,0.04033,0.27102,0.36296,-0.02818,0.15123,0.22521,0.04261,0.2343,0.09552,0.48601,0.14213,0.01298,0.58883,0.27815,-0.01931,0.71476,0.12135,0.08347,0.82945,0.18553,0.19177,0.81338,0.00963,0.3102,0.73961,0.14402,0.30834,0.59137,0.04778,0.24353,0.50553,0.23353
    if(vlabel.at(ifound)=="A_mP64_14_16e") {
      PrototypeANRL_A_mP64_14_16e(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 23 // ./aflow --proto=A2B5_mC28_15_f_e2f --params=12.786,0.387533239481,0.427968090099,97.03333,0.5727,0.106,0.311,0.077,0.0958,0.0952,0.4213,0.7127,0.0726,0.3138
    if(vlabel.at(ifound)=="A2B5_mC28_15_f_e2f") {
      PrototypeANRL_A2B5_mC28_15_f_e2f(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 24 // ./aflow --proto=AB_mC8_15_c_e --params=4.6837,0.730747058949,1.0950317057,120.34,0.4184
    if(vlabel.at(ifound)=="AB_mC8_15_c_e") {
      PrototypeANRL_AB_mC8_15_c_e(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 25 // ./aflow --proto=A2B_mC48_15_ae3f_2f --params=7.1356,1.73344918437,1.00532541062,120.34,0.1163,0.266,0.1234,0.9401,0.3114,0.1038,0.3282,0.0172,0.2117,0.4782,0.14033,0.10833,0.07227,0.50682,0.15799,0.54077
    if(vlabel.at(ifound)=="A2B_mC48_15_ae3f_2f") {
      PrototypeANRL_A2B_mC48_15_ae3f_2f(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 26 // ./aflow --proto=ABC6D2_mC40_15_e_e_3f_f --params=9.79,0.901123595506,0.548518896834,105.81,0.3082,-0.0942,0.3888,0.4123,0.8659,0.1365,0.2411,0.6799,0.1468,0.4802,0.0124,0.2117,0.4057,0.7764
    if(vlabel.at(ifound)=="ABC6D2_mC40_15_e_e_3f_f") {
      PrototypeANRL_ABC6D2_mC40_15_e_e_3f_f(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 27 // ./aflow --proto=ABC4_oP12_16_ag_cd_2u --params=5.61,1.01069518717,1.61319073084,0.2,0.26,0.125,0.74,0.8,0.63
    if(vlabel.at(ifound)=="ABC4_oP12_16_ag_cd_2u") {
      PrototypeANRL_ABC4_oP12_16_ag_cd_2u(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 28 // ./aflow --proto=AB3_oP16_18_ab_3c --params=8.32,1.15865384615,0.579326923077,0.0,0.0,0.25,0.25,0.0,0.25,0.5,0.5,0.124,0.309,0.382
    if(vlabel.at(ifound)=="AB3_oP16_18_ab_3c") {
      PrototypeANRL_AB3_oP16_18_ab_3c(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 29 // ./aflow --proto=A2B_oP12_19_2a_a --params=7.764,0.909582689335,0.558088614116,0.185,0.07,0.465,0.055,0.765,-0.008,0.884,-0.011,0.391
    if(vlabel.at(ifound)=="A2B_oP12_19_2a_a") {
      PrototypeANRL_A2B_oP12_19_2a_a(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 30 // ./aflow --proto=A2B_oC24_20_abc_c --params=8.74,0.576659038902,0.942791762014,0.3336,0.4403,0.2453,0.1971,0.2713,0.33154,0.03589,0.81143
    if(vlabel.at(ifound)=="A2B_oC24_20_abc_c") {
      PrototypeANRL_A2B_oC24_20_abc_c(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 31 // ./aflow --proto=AB_oP2_25_b_a --params=2.8102,1.87104120703,1.0769696107,0.0,0.25
    if(vlabel.at(ifound)=="AB_oP2_25_b_a") {
      PrototypeANRL_AB_oP2_25_b_a(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 32 // ./aflow --proto=AB2_oP24_28_acd_2c3d --params=16.54,0.533252720677,0.269649334946,0.0,0.319,0.014,0.018,0.042,0.617,0.042,0.624,0.334,0.5,0.503,0.301,0.042,0.632,0.636,0.5,0.619,0.036,0.5
    if(vlabel.at(ifound)=="AB2_oP24_28_acd_2c3d") {
      PrototypeANRL_AB2_oP24_28_acd_2c3d(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 33 // ./aflow --proto=AB3C4_oP16_31_a_ab_2ab --params=7.43,0.869448183042,0.831763122476,0.8268,0.0,0.1514,0.4983,0.8226,0.6454,0.1436,0.1166,0.2466,0.3255,-0.0134,0.2598,0.3364,0.6184
    if(vlabel.at(ifound)=="AB3C4_oP16_31_a_ab_2ab") {
      PrototypeANRL_AB3C4_oP16_31_a_ab_2ab(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 34 // ./aflow --proto=AB_oP8_33_a_a --params=5.2857,1.11007056776,0.659950432298,0.1996,0.5867,0.2506,0.002,0.2003,0.25  
    // Change z1 to be different than z2 to get SG #33
    if(vlabel.at(ifound)=="AB_oP8_33_a_a") {
      PrototypeANRL_AB_oP8_33_a_a(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 35 // ./aflow --proto=AB3C4_oP32_33_a_3a_4a --params=9.11,1.01866081229,1.1613611416,0.2187,0.4807,0.2031,0.4418,0.2052,0.0015,0.4488,0.1967,0.4146,0.1422,0.9176,0.2246,0.191,0.2506,0.2228,0.3424,0.5361,0.0415,0.0069,0.5876,0.2212,0.3355,0.546,0.3761
    if(vlabel.at(ifound)=="AB3C4_oP32_33_a_3a_4a") {
      PrototypeANRL_AB3C4_oP32_33_a_3a_4a(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 36 // ./aflow --proto=A2B_oC12_36_2a_a --params=4.624,1.46820934256,2.69139273356,0.333,0.0,0.061,0.134,0.395,0.366
    if(vlabel.at(ifound)=="A2B_oC12_36_2a_a") {
      PrototypeANRL_A2B_oC12_36_2a_a(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 37 // ./aflow --proto=A2BC_oC8_38_e_a_b --params=3.875,1.17470967742,1.59019354839,0.0,0.6144,0.155,0.2914
    if(vlabel.at(ifound)=="A2BC_oC8_38_e_a_b") {
      PrototypeANRL_A2BC_oC8_38_e_a_b(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 38 // ./aflow --proto=A2B_oC12_38_de_ab --params=4.684,1.81084543126,1.0269000854,0.06,0.5,0.17,0.56,0.17,0.0
    // Change z3 from z4 to get SG #38 (discussed in paper)
    if(vlabel.at(ifound)=="A2B_oC12_38_de_ab") {
      PrototypeANRL_A2B_oC12_38_de_ab(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 39 // ./aflow --proto=AB4_oC20_41_a_2b --params=6.388,1.00485284909,1.7778647464,0.0,0.673,0.327,0.376,0.827,0.673,0.125
    // Change z3 to get SG #41
    if(vlabel.at(ifound)=="AB4_oC20_41_a_2b") {
      PrototypeANRL_AB4_oC20_41_a_2b(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 40 // ./aflow --proto=AB2_oC24_41_2a_2b --params=6.478,1.0,1.87635072553,0.01,0.238,0.342,0.158,0.125,0.25,0.25,-0.125
    // Change z4 to get SG #41
    if(vlabel.at(ifound)=="AB2_oC24_41_2a_2b") {
      PrototypeANRL_AB2_oC24_41_2a_2b(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 41 // ./aflow --proto=AB2_oF72_43_ab_3b --params=11.66,1.91595197256,0.58833619211,0.0,0.125,0.13889,0.0,0.02222,0.08056,0.18333,0.15278,-0.01389,-0.18333,0.0625,0.125,0.27778
    if(vlabel.at(ifound)=="AB2_oF72_43_ab_3b") {
      PrototypeANRL_AB2_oF72_43_ab_3b(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 42 // ./aflow --proto=AB_oI4_44_a_b --params=4.92,0.973577235772,0.535569105691,0.0,0.425
    if(vlabel.at(ifound)=="AB_oI4_44_a_b") {
      PrototypeANRL_AB_oI4_44_a_b(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 43 // ./aflow --proto=A2B3C7D_oP13_47_t_aq_eqrs_h --params=3.8187,1.01691675177,3.05567339671,0.3554,0.1579,0.3771,0.3788,0.18445
    if(vlabel.at(ifound)=="A2B3C7D_oP13_47_t_aq_eqrs_h") {
      PrototypeANRL_A2B3C7D_oP13_47_t_aq_eqrs_h(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 44 // ./aflow --proto=AB_oP4_51_e_f --params=4.7549,0.661969757513,1.0209678437,0.8125,0.3125
    if(vlabel.at(ifound)=="AB_oP4_51_e_f") {
      PrototypeANRL_AB_oP4_51_e_f(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 45 // ./aflow --proto=A3B2_oP20_56_ce_e --params=4.911,2.53797597231,1.10201588271,0.029,0.147,0.058,0.861,0.044,0.128,0.179
    if(vlabel.at(ifound)=="A3B2_oP20_56_ce_e") {
      PrototypeANRL_A3B2_oP20_56_ce_e(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 46 // ./aflow --proto=ABCD_oP16_57_d_c_d_d --params=6.707,0.997614432682,1.13627553303,0.208,0.7704,0.2871,0.889,0.4154,0.605,0.1087
    if(vlabel.at(ifound)=="ABCD_oP16_57_d_c_d_d") {
      PrototypeANRL_ABCD_oP16_57_d_c_d_d(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 47 // ./aflow --proto=AB_oP8_57_d_d --params=6.09556,0.900425883758,0.850291031505,0.8593,0.0628,0.255,0.0096
    if(vlabel.at(ifound)=="AB_oP8_57_d_d") {
      PrototypeANRL_AB_oP8_57_d_d(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 48 // ./aflow --proto=AB2_oP6_58_a_g --params=6.24,1.03044871795,0.673076923077,0.275,0.325
    // 49 // ./aflow --proto=AB2_oP6_58_a_g --params=4.704,0.917942176871,0.601615646259,0.66667,0.25
    // 50 // ./aflow --proto=AB2_oP6_58_a_g --params=4.4446,1.22049228277,0.761913333033,0.2004,0.3787
    if(vlabel.at(ifound)=="AB2_oP6_58_a_g") {
      PrototypeANRL_AB2_oP6_58_a_g(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 51 // ./aflow --proto=AB_oP4_59_a_b --params=3.15,1.29841269841,2.20634920635,0.051,0.277
    if(vlabel.at(ifound)=="AB_oP4_59_a_b") {
      PrototypeANRL_AB_oP4_59_a_b(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 52 // ./aflow --proto=ABC_oP6_59_a_a_a --params=5.68,0.700704225352,1.01056338028,0.1499,0.4237,0.6255
    if(vlabel.at(ifound)=="ABC_oP6_59_a_a_a") {
      PrototypeANRL_ABC_oP6_59_a_a_a(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 53 // ./aflow --proto=A3B_oP8_59_bf_a --params=5.162,0.842115459124,0.877760557923,0.67125,0.329,0.505,0.174
    if(vlabel.at(ifound)=="A3B_oP8_59_bf_a") {
      PrototypeANRL_A3B_oP8_59_bf_a(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 54 // ./aflow --proto=AB_oP16_61_c_c --params=6.471,1.27538247566,1.31757070005,0.136,0.072,0.108,0.456,0.119,0.872
    if(vlabel.at(ifound)=="AB_oP16_61_c_c") {
      PrototypeANRL_AB_oP16_61_c_c(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 55 // ./aflow --proto=A2B_oP24_61_2c_c --params=9.174,0.375953782429,0.560061042075,0.0095,0.1491,0.1835,0.2314,0.111,0.5366,0.1289,0.0972,0.8628
    if(vlabel.at(ifound)=="A2B_oP24_61_2c_c") {
      PrototypeANRL_A2B_oP24_61_2c_c(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 56 // ./aflow --proto=A3B2_oP20_62_3c_2c --params=11.282,0.339443361106,0.994947704308,0.2922,0.19181,0.4504,0.877,0.6246,0.5611,-0.02937,0.17398,0.64939,-0.03603
    if(vlabel.at(ifound)=="A3B2_oP20_62_3c_2c") {
      PrototypeANRL_A3B2_oP20_62_3c_2c(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 57 // ./aflow --proto=AB3C_oP20_62_c_cd_a --params=5.4224,1.41099881971,0.996661994689,0.4877,-0.0084,0.0313,0.0586,0.288,0.537,0.213
    if(vlabel.at(ifound)=="AB3C_oP20_62_c_cd_a") {
      PrototypeANRL_AB3C_oP20_62_c_cd_a(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 58 // ./aflow --proto=A4B_oP20_62_2cd_c --params=5.464,0.810395314788,1.36749633968,0.22451,0.65626,0.55801,0.6466,0.05131,0.36362,0.13079,0.0579,0.06543
    if(vlabel.at(ifound)=="A4B_oP20_62_2cd_c") {
      PrototypeANRL_A4B_oP20_62_2cd_c(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 59 // ./aflow --proto=AB2C_oP16_62_c_2c_c --params=6.018,0.630741110003,2.4086075108,0.2522,0.8276,0.6221,0.095,0.8706,0.8244,0.226,0.06333
    if(vlabel.at(ifound)=="AB2C_oP16_62_c_2c_c") {
      PrototypeANRL_AB2C_oP16_62_c_2c_c(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 60 // ./aflow --proto=A2B_oP12_62_2c_c --params=4.918,0.7600650671,1.44550630338,0.038,0.282,0.674,0.562,0.202,0.611
    // 61 // ./aflow --proto=A2B_oP12_62_2c_c --params=12.735,0.468237141735,0.339615233608,0.733,0.125,0.508,0.722,0.874,0.447
    // 62 // ./aflow --proto=A2B_oP12_62_2c_c --params=7.6204,0.595008136056,1.1869718125,0.125,0.4217,0.0202,0.837,0.2377,0.0959
    if(vlabel.at(ifound)=="A2B_oP12_62_2c_c") {
      PrototypeANRL_A2B_oP12_62_2c_c(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 63 // ./aflow --proto=AB_oP8_62_c_c --params=10.42,0.349328214971,0.411708253359,0.375,0.333,0.139,0.389
    // 64 // ./aflow --proto=AB_oP8_62_c_c --params=5.24160,0.606723137973,1.12622100122,0.0056,0.1952,0.1879,0.5696
    // 68 // ./aflow --proto=AB_oP8_62_c_c --params=5.495,0.536123748863,0.737579617834,0.125,0.69,-0.18,0.125
    // 69 // ./aflow --proto=AB_oP8_62_c_c --params=11.18,0.356171735242,0.387209302326,0.3507,0.0201,0.61937,0.3806
    if(vlabel.at(ifound)=="AB_oP8_62_c_c") {
      PrototypeANRL_AB_oP8_62_c_c(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 65 // ./aflow --proto=AB3_oP16_62_c_cd --params=5.09,1.3257367387,0.888605108055,0.39,0.05,0.036,0.852,0.186,0.063,0.328
    if(vlabel.at(ifound)=="AB3_oP16_62_c_cd") {
      PrototypeANRL_AB3_oP16_62_c_cd(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 66 // ./aflow --proto=A3B7_oP40_62_cd_3c2d --params=4.526,1.54882898807,2.68272205038,0.4594,0.5629,0.0579,0.6261,0.2501,0.2063,0.2619,0.4165,0.0288,0.0291,0.3428,0.0565,0.0642,0.8119,0.2509,0.0657,0.0218
    if(vlabel.at(ifound)=="A3B7_oP40_62_cd_3c2d") {
      PrototypeANRL_A3B7_oP40_62_cd_3c2d(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 67 // ./aflow --proto=A_oP8_62_2c --params=6.663,0.708839861924,0.73345339937,0.464,0.292,0.181,0.658
    if(vlabel.at(ifound)=="A_oP8_62_2c") {
      PrototypeANRL_A_oP8_62_2c(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 70 // ./aflow --proto=AB2C_oC16_63_c_2c_c --params=3.577,4.56863293263,1.09538719597,0.06109,-0.0558,0.1792,0.33096
    if(vlabel.at(ifound)=="AB2C_oC16_63_c_2c_c") {
      PrototypeANRL_AB2C_oC16_63_c_2c_c(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 71 // ./aflow --proto=A2B_oC12_63_2c_c --params=3.73,3.94638069705,0.983914209115,0.061,0.75,0.396
    if(vlabel.at(ifound)=="A2B_oC12_63_2c_c") {
      PrototypeANRL_A2B_oC12_63_2c_c(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 72 // ./aflow --proto=AB_oC8_63_c_c --params=2.9782,2.64253575985,0.985360284736,0.436,0.14525
    if(vlabel.at(ifound)=="AB_oC8_63_c_c") {
      PrototypeANRL_AB_oC8_63_c_c(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 73 // ./aflow --proto=A_oC4_63_c --params=2.8444,2.06331739558,1.73379271551,0.10228
    if(vlabel.at(ifound)=="A_oC4_63_c") {
      PrototypeANRL_A_oC4_63_c(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 74 // ./aflow --proto=A_oC8_64_f --params=4.523,1.69378730931,1.0002210922,0.1549,0.081
    // 76 // ./aflow --proto=A_oC8_64_f --params=3.3136,3.16211974891,1.32070859488,0.10168,0.08056
    // 77 // ./aflow --proto=A_oC8_64_f --params=7.11906,0.654575182679,1.37596817557,0.15485,0.1175
    if(vlabel.at(ifound)=="A_oC8_64_f") {
      PrototypeANRL_A_oC8_64_f(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 75 // ./aflow --proto=A2B2C_oC80_64_efg_efg_df --params=10.922,0.866233290606,0.682933528658,0.84657,0.0946,0.9271,0.5886,0.276,-0.0792,0.2314,0.27981,-0.0113,0.1278,0.3415,0.2438,0.1245,0.175,0.2231
    if(vlabel.at(ifound)=="A2B2C_oC80_64_efg_efg_df") {
      PrototypeANRL_A2B2C_oC80_64_efg_efg_df(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 78 // ./aflow --proto=AB_oC8_65_j_g --params=5.971,1.1314687657,0.468263272484,0.28,0.22
    if(vlabel.at(ifound)=="AB_oC8_65_j_g") {
      PrototypeANRL_AB_oC8_65_j_g(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 79 // ./aflow --proto=A3B5_oC16_65_ah_bej --params=8.031,0.926410160628,0.491595069107,0.25,0.225
    if(vlabel.at(ifound)=="A3B5_oC16_65_ah_bej") {
      PrototypeANRL_A3B5_oC16_65_ah_bej(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 80 // ./aflow --proto=AB3_oC8_65_a_bf --params=5.82068,1.35259626023,0.493507631411
    if(vlabel.at(ifound)=="AB3_oC8_65_a_bf") {
      PrototypeANRL_AB3_oC8_65_a_bf(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 81 // ./aflow --proto=AB_oF8_69_a_b --params=6.08,0.903782894737,0.851973684211
    if(vlabel.at(ifound)=="AB_oF8_69_a_b") {
      PrototypeANRL_AB_oF8_69_a_b(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 82 // ./aflow --proto=A_oF8_70_a --params=3.1587,1.82613100326,3.21714629436
    if(vlabel.at(ifound)=="A_oF8_70_a") {
      PrototypeANRL_A_oF8_70_a(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 83 // ./aflow --proto=A2B_oF24_70_e_a --params=8.2671,0.580614725841,1.0342804611,0.4615
    if(vlabel.at(ifound)=="A2B_oF24_70_e_a") {
      PrototypeANRL_A2B_oF24_70_e_a(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 84 // ./aflow --proto=A_oF128_70_4h --params=10.4646,1.22947843205,2.33988876785,0.14415,0.04732,0.0486,0.29277,0.2269,0.25406,0.21598,0.28022,0.32618,0.21405,0.15761,0.37947
    if(vlabel.at(ifound)=="A_oF128_70_4h") {
      PrototypeANRL_A_oF128_70_4h(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 85 // ./aflow --proto=AB2_oI6_71_a_i --params=3.144,0.994910941476,2.44179389313,0.339
    if(vlabel.at(ifound)=="AB2_oI6_71_a_i") {
      PrototypeANRL_AB2_oI6_71_a_i(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 86 // ./aflow --proto=AB2_oI6_71_a_g --params=2.75984,2.9999963766,1.4241115427,0.35333
    if(vlabel.at(ifound)=="AB2_oI6_71_a_g") {
      PrototypeANRL_AB2_oI6_71_a_g(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 87 // ./aflow --proto=A2B_oI12_72_j_a --params=9.583,0.585829072316,0.578837524783,0.1182,0.2088
    if(vlabel.at(ifound)=="A2B_oI12_72_j_a") {
      PrototypeANRL_A2B_oI12_72_j_a(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 88 // ./aflow --proto=AB4C_tI12_82_c_g_a --params=4.3404,1.53216293429,0.256,0.2566,0.3722
    if(vlabel.at(ifound)=="AB4C_tI12_82_c_g_a") {
      PrototypeANRL_AB4C_tI12_82_c_g_a(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 89 // ./aflow --proto=A2BC4_tI14_82_bc_a_g --params=5.55,1.85585585586,0.26,0.25,0.13
    if(vlabel.at(ifound)=="A2BC4_tI14_82_bc_a_g") {
      PrototypeANRL_A2BC4_tI14_82_bc_a_g(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 90 // ./aflow --proto=AB_tP16_84_cej_k --params=6.429,1.02830922383,0.46779,0.25713,0.19361,0.30754,0.22904
    if(vlabel.at(ifound)=="AB_tP16_84_cej_k") {
      PrototypeANRL_AB_tP16_84_cej_k(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 91 // ./aflow --proto=A4B5_tI18_87_h_ah --params=10.164,0.37111373475,0.2797,-0.0589,0.3752,0.6856
    if(vlabel.at(ifound)=="A4B5_tI18_87_h_ah") {
      PrototypeANRL_A4B5_tI18_87_h_ah(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 92 // ./aflow --proto=AB4_tI10_87_a_h --params=5.72,0.623076923077,0.4,0.8
    if(vlabel.at(ifound)=="AB4_tI10_87_a_h") {
      PrototypeANRL_AB4_tI10_87_a_h(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 93 // ./aflow --proto=A2B_tP12_92_b_a --params=4.957,1.39001412144,0.3047,0.2381,0.1109,0.1826
    if(vlabel.at(ifound)=="A2B_tP12_92_b_a") {
      PrototypeANRL_A2B_tP12_92_b_a(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 94 // ./aflow --proto=A2B_tP36_96_3b_ab --params=7.464,1.15487674169,0.41,0.445,0.132,0.4,0.117,0.123,0.296,0.344,0.297,0.143,0.326,0.12,0.248
    if(vlabel.at(ifound)=="A2B_tP36_96_3b_ab") {
      PrototypeANRL_A2B_tP36_96_3b_ab(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 95 // ./aflow --proto=A_tP12_96_ab --params=5.51889,1.25999974633,0.0849,0.1752,0.3792,0.2742
    if(vlabel.at(ifound)=="A_tP12_96_ab") {
      PrototypeANRL_A_tP12_96_ab(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 96 // ./aflow --proto=A3BC_tP5_99_bc_a_b --params=4.046,1.02308452793,0.0,0.8973,0.4517,0.3785
    if(vlabel.at(ifound)=="A3BC_tP5_99_bc_a_b") {
      PrototypeANRL_A3BC_tP5_99_bc_a_b(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 97 // ./aflow --proto=AB3_tP8_113_a_ce --params=6.871,0.606622034638,0.206,0.1797,0.476
    if(vlabel.at(ifound)=="AB3_tP8_113_a_ce") {
      PrototypeANRL_AB3_tP8_113_a_ce(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 98 // ./aflow --proto=A2BC4D_tI16_121_d_a_i_b --params=5.46,1.96428571429,0.245,0.132
    if(vlabel.at(ifound)=="A2BC4D_tI16_121_d_a_i_b") {
      PrototypeANRL_A2BC4D_tI16_121_d_a_i_b(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 99 // ./aflow --proto=ABC2_tI16_122_a_b_d --params=5.289,1.97069389299,0.2574
    if(vlabel.at(ifound)=="ABC2_tI16_122_a_b_d") {
      PrototypeANRL_ABC2_tI16_122_a_b_d(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 100 // ./aflow --proto=AB5C_tP7_123_b_ci_a --params=4.207,1.61516520086,0.312
    if(vlabel.at(ifound)=="AB5C_tP7_123_b_ci_a") {
      PrototypeANRL_AB5C_tP7_123_b_ci_a(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 101 // ./aflow --proto=AB3_tP4_123_a_ce --params=4.158,0.864357864358
    if(vlabel.at(ifound)=="AB3_tP4_123_a_ce") {
      PrototypeANRL_AB3_tP4_123_a_ce(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 102 // ./aflow --proto=AB_tP2_123_a_d --params=2.8,1.31071428571
    if(vlabel.at(ifound)=="AB_tP2_123_a_d") {
      PrototypeANRL_AB_tP2_123_a_d(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 103 // ./aflow --proto=ABC2_tP4_123_d_a_f --params=3.8611,0.828649866618
    if(vlabel.at(ifound)=="ABC2_tP4_123_d_a_f") {
      PrototypeANRL_ABC2_tP4_123_d_a_f(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 104 // ./aflow --proto=A2B3_tP10_127_g_ah --params=7.3364,0.530232811733,0.3841,0.1821
    if(vlabel.at(ifound)=="A2B3_tP10_127_g_ah") {
      PrototypeANRL_A2B3_tP10_127_g_ah(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 105 // ./aflow --proto=ABCD_tP8_129_c_b_a_c --params=3.6736,2.60540069686,0.6793,0.2246
    if(vlabel.at(ifound)=="ABCD_tP8_129_c_b_a_c") {
      PrototypeANRL_ABCD_tP8_129_c_b_a_c(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 106 // ./aflow --proto=A_tP4_129_ac --params=4.897,0.69185215438,0.375
    if(vlabel.at(ifound)=="A_tP4_129_ac") {
      PrototypeANRL_A_tP4_129_ac(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 107 // ./aflow --proto=ABC_tP6_129_c_a_c --params=4.11,1.76301703163,0.6497,0.2058
    if(vlabel.at(ifound)=="ABC_tP6_129_c_a_c") {
      PrototypeANRL_ABC_tP6_129_c_a_c(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 108 // ./aflow --proto=A2B_tP6_129_ac_c --params=4.0006,1.52584612308,0.27,0.7
    if(vlabel.at(ifound)=="A2B_tP6_129_ac_c") {
      PrototypeANRL_A2B_tP6_129_ac_c(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 109 // ./aflow --proto=AB_tP4_129_a_c --params=3.9645,1.26008323874,0.2368
    if(vlabel.at(ifound)=="AB_tP4_129_a_c") {
      PrototypeANRL_AB_tP4_129_a_c(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 110 // ./aflow --proto=AB_tP4_129_c_c --params=3.107,1.90505310589,0.1,0.65
    if(vlabel.at(ifound)=="AB_tP4_129_c_c") {
      PrototypeANRL_AB_tP4_129_c_c(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 111 // ./aflow --proto=AB_tP4_131_c_e --params=4.9073,1.24500234345
    if(vlabel.at(ifound)=="AB_tP4_131_c_e") {
      PrototypeANRL_AB_tP4_131_c_e(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 112 // ./aflow --proto=A_tP50_134_b2m2n --params=8.74,0.575514874142,0.0048,0.1685,0.1305,0.628,0.1695,0.5228,0.1635,0.0753,0.3383,0.1485
    if(vlabel.at(ifound)=="A_tP50_134_b2m2n") {
      PrototypeANRL_A_tP50_134_b2m2n(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 113 // ./aflow --proto=A_tP30_136_bf2ij --params=10.59,0.532011331445,0.1033,0.3667,0.0383,0.5608,0.2354,0.3183,0.27
    if(vlabel.at(ifound)=="A_tP30_136_bf2ij") {
      PrototypeANRL_A_tP30_136_bf2ij(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 114 // ./aflow --proto=AB_tP8_136_g_f --params=4.75,0.576842105263,0.31,0.336
    if(vlabel.at(ifound)=="AB_tP8_136_g_f") {
      PrototypeANRL_AB_tP8_136_g_f(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 115 // ./aflow --proto=A2B_tP6_136_f_a --params=4.5922,0.644005052045,0.30496
    if(vlabel.at(ifound)=="A2B_tP6_136_f_a") {
      PrototypeANRL_A2B_tP6_136_f_a(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 116 // ./aflow --proto=sigma_tP30_136_bf2ij --params=8.7966,0.518177477662,0.39864,0.13122,0.46349,0.06609,0.73933,0.18267,0.25202
    if(vlabel.at(ifound)=="sigma_tP30_136_bf2ij") {
      PrototypeANRL_sigma_tP30_136_bf2ij(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 117 // ./aflow --proto=A_tP4_136_f --params=3.957,1.29112964367,0.098
    if(vlabel.at(ifound)=="A_tP4_136_f") {
      PrototypeANRL_A_tP4_136_f(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 118 // ./aflow --proto=A_tP16_138_j --params=8.56,0.714953271028,0.375,-0.083,0.857
    if(vlabel.at(ifound)=="A_tP16_138_j") {
      PrototypeANRL_A_tP16_138_j(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 119 // ./aflow --proto=A3B_tI16_139_cde_e --params=3.9993,4.3215062636,0.37498,0.11886
    if(vlabel.at(ifound)=="A3B_tI16_139_cde_e") {
      PrototypeANRL_A3B_tI16_139_cde_e(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 120 // ./aflow --proto=A_tI4_139_e --params=3.34916,1.94217355994,0.819
    if(vlabel.at(ifound)=="A_tI4_139_e") {
      PrototypeANRL_A_tI4_139_e(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 121 // ./aflow --proto=AB2C4_tI14_139_a_e_ce --params=3.7817,3.50337149959,0.36075,0.1824
    if(vlabel.at(ifound)=="AB2C4_tI14_139_a_e_ce") {
      PrototypeANRL_AB2C4_tI14_139_a_e_ce(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 122 // ./aflow --proto=A12B_tI26_139_fij_a --params=8.47,0.584415584416,0.361,0.278
    if(vlabel.at(ifound)=="A12B_tI26_139_fij_a") {
      PrototypeANRL_A12B_tI26_139_fij_a(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 123 // ./aflow --proto=A_tI2_139_a --params=4.6002,1.07523585931
    // 131 // ./aflow --proto=A_tI2_139_a --params=3.932,0.823499491353
    if(vlabel.at(ifound)=="A_tI2_139_a") {
      PrototypeANRL_A_tI2_139_a(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 124 // ./aflow --proto=A_tI8_139_h --params=4.33184,0.574102459925,0.17916
    if(vlabel.at(ifound)=="A_tI8_139_h") {
      PrototypeANRL_A_tI8_139_h(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 125 // ./aflow --proto=A3B_tI8_139_bd_a --params=3.8537,2.22744375535
    if(vlabel.at(ifound)=="A3B_tI8_139_bd_a") {
      PrototypeANRL_A3B_tI8_139_bd_a(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 126 // ./aflow --proto=AB2_tI6_139_a_e --params=3.2064,2.44754241517,0.3353
    if(vlabel.at(ifound)=="AB2_tI6_139_a_e") {
      PrototypeANRL_AB2_tI6_139_a_e(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 127 // ./aflow --proto=A4B5_tI18_139_i_ah --params=8.91,0.361391694725,0.328,0.348
    if(vlabel.at(ifound)=="A4B5_tI18_139_i_ah") {
      PrototypeANRL_A4B5_tI18_139_i_ah(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 128 // ./aflow --proto=A4B_tI10_139_de_a --params=4.53,2.45033112583,0.38
    if(vlabel.at(ifound)=="A4B_tI10_139_de_a") {
      PrototypeANRL_A4B_tI10_139_de_a(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 129 // ./aflow --proto=A8B_tI18_139_hi_a --params=8.312,0.468840230991,0.333,0.327
    if(vlabel.at(ifound)=="A8B_tI18_139_hi_a") {
      PrototypeANRL_A8B_tI18_139_hi_a(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 130 // ./aflow --proto=A2B_tI6_139_d_a --params=4.1,1.22682926829
    if(vlabel.at(ifound)=="A2B_tI6_139_d_a") {
      PrototypeANRL_A2B_tI6_139_d_a(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 132 // ./aflow --proto=A2B_tI12_140_h_a --params=6.04,0.804635761589,0.158
    if(vlabel.at(ifound)=="A2B_tI12_140_h_a") {
      PrototypeANRL_A2B_tI12_140_h_a(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 133 // ./aflow --proto=AB3_tI16_140_b_ah --params=6.017,1.44241316271,0.231
    if(vlabel.at(ifound)=="AB3_tI16_140_b_ah") {
      PrototypeANRL_AB3_tI16_140_b_ah(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 134 // ./aflow --proto=AB_tI16_140_ab_h --params=8.03,0.87297633873,0.179
    if(vlabel.at(ifound)=="AB_tI16_140_ab_h") {
      PrototypeANRL_AB_tI16_140_ab_h(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 135 // ./aflow --proto=A4BC_tI24_141_h_b_a --params=6.6042,0.905423821205,0.066,0.1951
    if(vlabel.at(ifound)=="A4BC_tI24_141_h_b_a") {
      PrototypeANRL_A4BC_tI24_141_h_b_a(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 136 // ./aflow --proto=A_tI4_141_a --params=5.8318,0.545611989437
    if(vlabel.at(ifound)=="A_tI4_141_a") {
      PrototypeANRL_A_tI4_141_a(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 137 // ./aflow --proto=A3B4_tI28_141_ad_h --params=5.765,1.63781439722,0.0278,0.2589
    if(vlabel.at(ifound)=="A3B4_tI28_141_ad_h") {
      PrototypeANRL_A3B4_tI28_141_ad_h(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 138 // ./aflow --proto=A2B_tI12_141_e_a --params=3.785,2.51360634082,0.20806
    if(vlabel.at(ifound)=="A2B_tI12_141_e_a") {
      PrototypeANRL_A2B_tI12_141_e_a(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 139 // ./aflow --proto=AB_tI16_141_e_e --params=3.108,5.45045045045,0.227,0.071
    if(vlabel.at(ifound)=="AB_tI16_141_e_e") {
      PrototypeANRL_AB_tI16_141_e_e(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 140 // ./aflow --proto=A2B_tI24_141_2e_e --params=4.046,6.28917449333,0.125,0.289,-0.051
    if(vlabel.at(ifound)=="A2B_tI24_141_2e_e") {
      PrototypeANRL_A2B_tI24_141_2e_e(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 141 // ./aflow --proto=AB_tI8_141_a_b --params=3.325,3.42255639098
    if(vlabel.at(ifound)=="AB_tI8_141_a_b") {
      PrototypeANRL_AB_tI8_141_a_b(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 142 // ./aflow --proto=A2B3_tI80_141_ceh_3h --params=7.5937,4.26037373086,0.2044,0.5201,0.3324,0.516,0.2547,0.494,0.0859,0.4667,0.4164
    if(vlabel.at(ifound)=="A2B3_tI80_141_ceh_3h") {
      PrototypeANRL_A2B3_tI80_141_ceh_3h(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 143 // ./aflow --proto=ABC4_tI96_142_e_ab_2g --params=10.914,1.77396005131,0.0375,0.2482,0.3197,-0.0867,0.0923,0.1117,0.0025
    if(vlabel.at(ifound)=="ABC4_tI96_142_e_ab_2g") {
      PrototypeANRL_ABC4_tI96_142_e_ab_2g(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 144 // ./aflow --proto=A2B_hP9_147_g_ad --params=7.636,0.369264012572,0.25,0.33333,0.0,0.25
    if(vlabel.at(ifound)=="A2B_hP9_147_g_ad") {
      PrototypeANRL_A2B_hP9_147_g_ad(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 145 // RHL: ./aflow --proto=AB_hR16_148_cf_cf --params=6.29713,1.8633345667,0.11546,0.21,0.10706,0.81289,0.19519,0.1848,0.6754,0.3468
    // 145 // HEX: ./aflow --proto=AB_hR16_148_cf_cf --params=6.29713,1.8633345667,0.11546,0.21,0.10706,0.81289,0.19519,0.1848,0.6754,0.3468 --hex
    if(vlabel.at(ifound)=="AB_hR16_148_cf_cf") {
      PrototypeANRL_AB_hR16_148_cf_cf(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 146 // RHL: ./aflow --proto=AB3_hR8_148_c_f --params=7.49626,2.75900649124,0.33333,0.088,0.755,0.421
    // 146 // HEX: ./aflow --proto=AB3_hR8_148_c_f --params=7.49626,2.75900649124,0.33333,0.088,0.755,0.421 --hex
    if(vlabel.at(ifound)=="AB3_hR8_148_c_f") {
      PrototypeANRL_AB3_hR8_148_c_f(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 147 // RHL: ./aflow --proto=AB_hR26_148_b2f_a2f --params=15.659,0.335334312536,0.054,0.346,0.098,0.754,0.15699,0.6,0.555,0.84401,0.599,0.252,0.65501,0.098
    // 147 // HEX: ./aflow --proto=AB_hR26_148_b2f_a2f --params=15.659,0.335334312536,0.054,0.346,0.098,0.754,0.15699,0.6,0.555,0.84401,0.599,0.252,0.65501,0.098 --hex
    if(vlabel.at(ifound)=="AB_hR26_148_b2f_a2f") {
      PrototypeANRL_AB_hR26_148_b2f_a2f(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 148 // RHL: ./aflow --proto=AB3C_hR10_148_c_f_c --params=5.0884,2.76815894977,0.35537,0.1464,0.22174,0.56249,0.95095
    // 148 // HEX: ./aflow --proto=AB3C_hR10_148_c_f_c --params=5.0884,2.76815894977,0.35537,0.1464,0.22174,0.56249,0.95095 --hex
    if(vlabel.at(ifound)=="AB3C_hR10_148_c_f_c") {
      PrototypeANRL_AB3C_hR10_148_c_f_c(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 149 // ./aflow --proto=A2B_hP9_150_ef_bd --params=5.85,0.589743589744,0.875,0.26,0.6
    if(vlabel.at(ifound)=="A2B_hP9_150_ef_bd") {
      PrototypeANRL_A2B_hP9_150_ef_bd(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 150 // ./aflow --proto=A3B_hP24_151_3c_2a --params=6.017,2.87518697025,0.8889,0.5556,0.8889,0.1111,0.0731,0.5556,0.4444,0.0731,0.2222,0.77778,0.0731
    if(vlabel.at(ifound)=="A3B_hP24_151_3c_2a") {
      PrototypeANRL_A3B_hP24_151_3c_2a(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 151 // ./aflow --proto=A2B_hP9_152_c_a --params=4.914,1.10012210012,0.4699,0.413,0.2668,0.214
    if(vlabel.at(ifound)=="A2B_hP9_152_c_a") {
      PrototypeANRL_A2B_hP9_152_c_a(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 152 // ./aflow --proto=A_hP3_152_a --params=4.3662,1.13453346159,0.2254
    if(vlabel.at(ifound)=="A_hP3_152_a") {
      PrototypeANRL_A_hP3_152_a(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 153 // ./aflow --proto=AB_hP6_154_a_b --params=4.145,2.29095295537,0.7198,0.4889
    if(vlabel.at(ifound)=="AB_hP6_154_a_b") {
      PrototypeANRL_AB_hP6_154_a_b(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 154 // RHL: ./aflow --proto=AB3_hR8_155_c_de --params=4.91608,2.53341483458,0.237,0.43,0.07
    // 154 // HEX: ./aflow --proto=AB3_hR8_155_c_de --params=4.91608,2.53341483458,0.237,0.43,0.07 --hex
    if(vlabel.at(ifound)=="AB3_hR8_155_c_de") {
      PrototypeANRL_AB3_hR8_155_c_de(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 155 // RHL: ./aflow --proto=A3B2_hR5_155_e_c --params=5.73296,1.24097324942,0.2521,0.2449
    // 155 // HEX: ./aflow --proto=A3B2_hR5_155_e_c --params=5.73296,1.24097324942,0.2521,0.2449 --hex
    if(vlabel.at(ifound)=="A3B2_hR5_155_e_c") {
      PrototypeANRL_A3B2_hR5_155_e_c(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 156 // RHL: ./aflow --proto=AB_hR6_160_b_b --params=9.619,0.327466472606,0.00019,0.26362,0.7288,0.39161
    // 156 // HEX: ./aflow --proto=AB_hR6_160_b_b --params=9.619,0.327466472606,0.00019,0.26362,0.7288,0.39161 --hex
    if(vlabel.at(ifound)=="AB_hR6_160_b_b") {
      PrototypeANRL_AB_hR6_160_b_b(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 157 // RHL: ./aflow --proto=AB_hR6_160_3a_3a --params=3.01791,7.34847294982,0.0,0.22222,0.77778,0.08333,0.30556,0.86111
    // 157 // HEX: ./aflow --proto=AB_hR6_160_3a_3a --params=3.01791,7.34847294982,0.0,0.22222,0.77778,0.08333,0.30556,0.86111 --hex
    if(vlabel.at(ifound)=="AB_hR6_160_3a_3a") {
      PrototypeANRL_AB_hR6_160_3a_3a(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 158 // RHL: ./aflow --proto=ABC3_hR10_161_a_a_b --params=5.2542,2.64091583876,0.2875,0.0128,0.74643,0.14093,0.36263
    // 158 // HEX: ./aflow --proto=ABC3_hR10_161_a_a_b --params=5.2542,2.64091583876,0.2875,0.0128,0.74643,0.14093,0.36263 --hex
    if(vlabel.at(ifound)=="ABC3_hR10_161_a_a_b") {
      PrototypeANRL_ABC3_hR10_161_a_a_b(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 159 // ./aflow --proto=AB2_hP9_162_ad_k --params=4.917,0.929021761237,0.325,0.272
    if(vlabel.at(ifound)=="AB2_hP9_162_ad_k") {
      PrototypeANRL_AB2_hP9_162_ad_k(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 160 // ./aflow --proto=AB2CD2_hP36_163_h_i_bf_i --params=7.384,2.37716684724,0.01,0.833,0.33333,0.03833,0.141,0.03167,0.365,0.083
    if(vlabel.at(ifound)=="AB2CD2_hP36_163_h_i_bf_i") {
      PrototypeANRL_AB2CD2_hP36_163_h_i_bf_i(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 161 // ./aflow --proto=A3B2_hP5_164_ad_d --params=4.0282,1.21409066084,0.648,0.149
    if(vlabel.at(ifound)=="A3B2_hP5_164_ad_d") {
      PrototypeANRL_A3B2_hP5_164_ad_d(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 162 // ./aflow --proto=AB2_hP3_164_a_d --params=4.24,1.61320754717,0.252
    if(vlabel.at(ifound)=="AB2_hP3_164_a_d") {
      PrototypeANRL_AB2_hP3_164_a_d(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 163 // ./aflow --proto=A3B_hP24_165_adg_f --params=6.308,1.03994927077,0.167,0.666,0.356,0.028,0.096
    if(vlabel.at(ifound)=="A3B_hP24_165_adg_f") {
      PrototypeANRL_A3B_hP24_165_adg_f(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 164 // RHL: ./aflow --proto=AB_hR2_166_a_b --params=3.13,4.78594249201
    // 164 // HEX: ./aflow --proto=AB_hR2_166_a_b --params=3.13,4.78594249201 --hex
    if(vlabel.at(ifound)=="AB_hR2_166_a_b") {
      PrototypeANRL_AB_hR2_166_a_b(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 165 // RHL: ./aflow --proto=A_hR2_166_c --params=3.7595,2.7815666977,0.22754
    // 165 // HEX: ./aflow --proto=A_hR2_166_c --params=3.7595,2.7815666977,0.22754 --hex
    // 172 // RHL: ./aflow --proto=A_hR2_166_c --params=2.456,4.08957654723,0.16667
    // 172 // HEX: ./aflow --proto=A_hR2_166_c --params=2.456,4.08957654723,0.16667 --hex
    // 175 // RHL: ./aflow --proto=A_hR2_166_c --params=3.289,3.42991790818,0.0543 
    // 175 // HEX: ./aflow --proto=A_hR2_166_c --params=3.289,3.42991790818,0.0543 --hex
    if(vlabel.at(ifound)=="A_hR2_166_c") {
      PrototypeANRL_A_hR2_166_c(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 166 // RHL: ./aflow --proto=A_hR1_166_a --params=5.07846,0.968139947937
    // 166 // HEX: ./aflow --proto=A_hR1_166_a --params=5.07846,0.968139947937 --hex
    // 170 // RHL: ./aflow --proto=A_hR1_166_a --params=3.45741,1.92728082582
    // 170 // HEX: ./aflow --proto=A_hR1_166_a --params=3.45741,1.92728082582 --hex
    if(vlabel.at(ifound)=="A_hR1_166_a") {
      PrototypeANRL_A_hR1_166_a(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 167 // RHL: ./aflow --proto=A7B6_hR13_166_ah_3c --params=4.757,5.4319949548,0.167,0.346,0.448,0.09,0.59001
    // 167 // HEX: ./aflow --proto=A7B6_hR13_166_ah_3c --params=4.757,5.4319949548,0.167,0.346,0.448,0.09,0.59001 --hex
    if(vlabel.at(ifound)=="A7B6_hR13_166_ah_3c") {
      PrototypeANRL_A7B6_hR13_166_ah_3c(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 168 // RHL: ./aflow --proto=A_hR3_166_ac --params=3.62036,7.25049442597,0.22222
    // 168 // HEX: ./aflow --proto=A_hR3_166_ac --params=3.62036,7.25049442597,0.22222 --hex
    if(vlabel.at(ifound)=="A_hR3_166_ac") {
      PrototypeANRL_A_hR3_166_ac(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 169 // RHL: ./aflow --proto=A2B3_hR5_166_c_ac --params=4.36914,6.96313919902,0.399,0.208
    // 169 // HEX: ./aflow --proto=A2B3_hR5_166_c_ac --params=4.36914,6.96313919902,0.399,0.208 --hex
    if(vlabel.at(ifound)=="A2B3_hR5_166_c_ac") {
      PrototypeANRL_A2B3_hR5_166_c_ac(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 171 // RHL: ./aflow --proto=A5B2_hR7_166_a2c_c --params=3.011,6.9511790103,0.186,0.33333,0.075
    // 171 // HEX: ./aflow --proto=A5B2_hR7_166_a2c_c --params=3.011,6.9511790103,0.186,0.33333,0.075 --hex
    if(vlabel.at(ifound)=="A5B2_hR7_166_a2c_c") {
      PrototypeANRL_A5B2_hR7_166_a2c_c(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 173 // RHL: ./aflow --proto=A_hR12_166_2h --params=4.908,2.56022616137,0.0104,0.65729,0.2206,0.6323
    // 173 // HEX: ./aflow --proto=A_hR12_166_2h --params=4.908,2.56022616137,0.0104,0.65729,0.2206,0.6323 --hex
    if(vlabel.at(ifound)=="A_hR12_166_2h") {
      PrototypeANRL_A_hR12_166_2h(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 174 // RHL: ./aflow --proto=ABC2_hR4_166_a_b_c --params=3.5561,5.44557239673,0.2667
    // 174 // HEX: ./aflow --proto=ABC2_hR4_166_a_b_c --params=3.5561,5.44557239673,0.2667 --hex
    if(vlabel.at(ifound)=="ABC2_hR4_166_a_b_c") {
      PrototypeANRL_ABC2_hR4_166_a_b_c(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 176 // RHL: ./aflow --proto=A_hR105_166_bc9h4i --params=10.96,2.17974452555,0.3848,0.3843,0.21309,0.4895,0.21780,0.3873,0.56899,0.1991,0.50609,0.1983,0.68740,0.1032,0.49209,0.9933,0.66980,0.1008,0.83740,0.0025,0.16801,0.3622,0.58109,0.0976,0.3765,0.68261,0.2024,0.1673,0.55209,0.8921,0.1777,0.3473,0.0033
    // 176 // HEX: ./aflow --proto=A_hR105_166_bc9h4i --params=10.96,2.17974452555,0.3848,0.3843,0.21309,0.4895,0.21780,0.3873,0.56899,0.1991,0.50609,0.1983,0.68740,0.1032,0.49209,0.9933,0.66980,0.1008,0.83740,0.0025,0.16801,0.3622,0.58109,0.0976,0.3765,0.68261,0.2024,0.1673,0.55209,0.8921,0.1777,0.3473,0.0033 --hex
    if(vlabel.at(ifound)=="A_hR105_166_bc9h4i") {
      PrototypeANRL_A_hR105_166_bc9h4i(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 177 // RHL: ./aflow --proto=A6B_hR7_166_g_a --params=4.33304,3.13251204697,0.16667
    // 177 // HEX: ./aflow --proto=A6B_hR7_166_g_a --params=4.33304,3.13251204697,0.16667 --hex
    if(vlabel.at(ifound)=="A6B_hR7_166_g_a") {
      PrototypeANRL_A6B_hR7_166_g_a(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 178 // RHL: ./aflow --proto=ABC3_hR10_167_a_b_e --params=5.285,2.62039735099,0.85756666666667 
    // 178 // HEX: ./aflow --proto=ABC3_hR10_167_a_b_e --params=5.285,2.62039735099,0.85756666666667 --hex
    // 179 // RHL: ./aflow --proto=ABC3_hR10_167_a_b_e --params=4.988,3.42040898156,0.5067 
    // 179 // HEX: ./aflow --proto=ABC3_hR10_167_a_b_e --params=4.988,3.42040898156,0.5067 --hex
    if(vlabel.at(ifound)=="ABC3_hR10_167_a_b_e") {
      PrototypeANRL_ABC3_hR10_167_a_b_e(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 180 // RHL: ./aflow --proto=A2B3_hR10_167_c_e --params=4.7607,2.72957758313,0.35216,0.5561
    // 180 // HEX: ./aflow --proto=A2B3_hR10_167_c_e --params=4.7607,2.72957758313,0.35216,0.5561 --hex
    if(vlabel.at(ifound)=="A2B3_hR10_167_c_e") {
      PrototypeANRL_A2B3_hR10_167_c_e(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 181 // ./aflow --proto=A2B_hP18_180_fi_bd --params=5.198,2.54136206233,0.163,0.1141
    if(vlabel.at(ifound)=="A2B_hP18_180_fi_bd") {
      PrototypeANRL_A2B_hP18_180_fi_bd(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 182 // ./aflow --proto=AB2_hP9_180_d_j --params=4.42758,1.43826876081,0.16559
    if(vlabel.at(ifound)=="AB2_hP9_180_d_j") {
      PrototypeANRL_AB2_hP9_180_d_j(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 183 // ./aflow --proto=A2B_hP9_180_j_c --params=4.9977,1.09252256038,0.2072
    if(vlabel.at(ifound)=="A2B_hP9_180_j_c") {
      PrototypeANRL_A2B_hP9_180_j_c(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 184 // ./aflow --proto=AB3_hP8_182_c_g --params=4.8507,0.866967654153,0.3249
    if(vlabel.at(ifound)=="AB3_hP8_182_c_g") {
      PrototypeANRL_AB3_hP8_182_c_g(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 185 // ./aflow --proto=A_hP4_186_ab --params=2.47,2.75303643725,0.0,0.07143
    if(vlabel.at(ifound)=="A_hP4_186_ab") {
      PrototypeANRL_A_hP4_186_ab(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 186 // ./aflow --proto=AB_hP8_186_ab_ab --params=3.08051,3.27374363336,0.18784,0.0,0.43671,0.24982
    if(vlabel.at(ifound)=="AB_hP8_186_ab_ab") {
      PrototypeANRL_AB_hP8_186_ab_ab(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 187 // ./aflow --proto=AB_hP4_186_b_b --params=3.8227,1.63776911607,0.3748,0.0
    if(vlabel.at(ifound)=="AB_hP4_186_b_b") {
      PrototypeANRL_AB_hP4_186_b_b(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 188 // ./aflow --proto=AB_hP12_186_a2b_a2b --params=3.08129,4.90695780014,0.1254,0.0,0.29215,-0.0415,0.16675,0.8335
    if(vlabel.at(ifound)=="AB_hP12_186_a2b_a2b") {
      PrototypeANRL_AB_hP12_186_a2b_a2b(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 189 // ./aflow --proto=A5B3C_hP18_186_2a3b_2ab_b --params=3.281,6.57726302956,0.155,0.345,0.0,0.248,0.045,0.261,0.455,0.367,0.137
    if(vlabel.at(ifound)=="A5B3C_hP18_186_2a3b_2ab_b") {
      PrototypeANRL_A5B3C_hP18_186_2a3b_2ab_b(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 190 // ./aflow --proto=AB_hP4_186_b_a --params=2.51,2.66932270916,0.0,0.05
    if(vlabel.at(ifound)=="AB_hP4_186_b_a") {
      PrototypeANRL_AB_hP4_186_b_a(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 191 // ./aflow --proto=ABC_hP3_187_a_d_f --params=4.535,1.0769570011
    if(vlabel.at(ifound)=="ABC_hP3_187_a_d_f") {
      PrototypeANRL_ABC_hP3_187_a_d_f(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 192 // ./aflow --proto=AB_hP2_187_d_a --params=2.9065,0.975950455875
    if(vlabel.at(ifound)=="AB_hP2_187_d_a") {
      PrototypeANRL_AB_hP2_187_d_a(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 193 // ./aflow --proto=A2B_hP9_189_fg_bc --params=5.877,0.584822188191,0.256,0.589
    if(vlabel.at(ifound)=="A2B_hP9_189_fg_bc") {
      PrototypeANRL_A2B_hP9_189_fg_bc(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 194 // ./aflow --proto=AB4C_hP6_191_a_h_b --params=3.04436,2.20489035462,0.2413
    if(vlabel.at(ifound)=="AB4C_hP6_191_a_h_b") {
      PrototypeANRL_AB4C_hP6_191_a_h_b(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 195 // ./aflow --proto=AB5_hP6_191_a_cg --params=5.405,0.773913043478
    if(vlabel.at(ifound)=="AB5_hP6_191_a_cg") {
      PrototypeANRL_AB5_hP6_191_a_cg(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 196 // ./aflow --proto=A_hP1_191_a --params=3.2062,0.931195808122
    if(vlabel.at(ifound)=="A_hP1_191_a") {
      PrototypeANRL_A_hP1_191_a(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 197 // ./aflow --proto=A3B_hP4_191_bc_a --params=3.6576,1.05902777778
    if(vlabel.at(ifound)=="A3B_hP4_191_bc_a") {
      PrototypeANRL_A3B_hP4_191_bc_a(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 198 // ./aflow --proto=AB2_hP3_191_a_d --params=3.005,1.08276206323
    if(vlabel.at(ifound)=="AB2_hP3_191_a_d") {
      PrototypeANRL_AB2_hP3_191_a_d(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 199 // ./aflow --proto=A2B_hP6_191_h_e --params=4.237,1.71040830776,0.306,0.16
    if(vlabel.at(ifound)=="A2B_hP6_191_h_e") {
      PrototypeANRL_A2B_hP6_191_h_e(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 200 // ./aflow --proto=AB_hP6_191_f_ad --params=5.279,0.806914188293
    if(vlabel.at(ifound)=="AB_hP6_191_f_ad") {
      PrototypeANRL_AB_hP6_191_f_ad(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 201 // ./aflow --proto=AB_hP8_194_ad_f --params=3.64,3.37362637363,0.125
    if(vlabel.at(ifound)=="AB_hP8_194_ad_f") {
      PrototypeANRL_AB_hP8_194_ad_f(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 202 // ./aflow --proto=A_hP6_194_h --params=4.40445,0.568892824303,0.44799
    if(vlabel.at(ifound)=="A_hP6_194_h") {
      PrototypeANRL_A_hP6_194_h(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 203 // ./aflow --proto=AB_hP12_194_af_bf --params=3.01,4.85382059801,0.166,0.583
    if(vlabel.at(ifound)=="AB_hP12_194_af_bf") {
      PrototypeANRL_AB_hP12_194_af_bf(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 204 // ./aflow --proto=A_hP4_194_ac --params=3.77,3.2175066313
    if(vlabel.at(ifound)=="A_hP4_194_ac") {
      PrototypeANRL_A_hP4_194_ac(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 205 // ./aflow --proto=AB3_hP8_194_c_bf --params=5.088,1.76533018868,-0.083
    if(vlabel.at(ifound)=="AB3_hP8_194_c_bf") {
      PrototypeANRL_AB3_hP8_194_c_bf(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 206 // ./aflow --proto=AB2_hP6_194_b_f --params=4.895,1.58324821246,0.045
    if(vlabel.at(ifound)=="AB2_hP6_194_b_f") {
      PrototypeANRL_AB2_hP6_194_b_f(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 207 // ./aflow --proto=AB_hP4_194_c_d --params=2.50399,2.66023426611
    if(vlabel.at(ifound)=="AB_hP4_194_c_d") {
      PrototypeANRL_AB_hP4_194_c_d(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 208 // ./aflow --proto=ABC2_hP8_194_d_a_f --params=2.86,4.48251748252,0.086
    if(vlabel.at(ifound)=="ABC2_hP8_194_d_a_f") {
      PrototypeANRL_ABC2_hP8_194_d_a_f(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 209 // ./aflow --proto=A3B_hP8_194_h_c --params=5.295,0.802077431539,0.8392
    if(vlabel.at(ifound)=="A3B_hP8_194_h_c") {
      PrototypeANRL_A3B_hP8_194_h_c(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 210 // ./aflow --proto=A_hP4_194_bc --params=2.464,2.72362012987
    if(vlabel.at(ifound)=="A_hP4_194_bc") {
      PrototypeANRL_A_hP4_194_bc(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 211 // ./aflow --proto=AB2_hP6_194_c_f --params=3.161,3.8895919013,0.6275
    if(vlabel.at(ifound)=="AB2_hP6_194_c_f") {
      PrototypeANRL_AB2_hP6_194_c_f(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 212 // ./aflow --proto=A5B2_hP14_194_abdf_f --params=2.982,4.651240778,0.528,0.139
    if(vlabel.at(ifound)=="A5B2_hP14_194_abdf_f") {
      PrototypeANRL_A5B2_hP14_194_abdf_f(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 213 // ./aflow --proto=AB2_hP12_194_f_ah --params=5.223,1.64005360904,0.06286,0.83048
    if(vlabel.at(ifound)=="AB2_hP12_194_f_ah") {
      PrototypeANRL_AB2_hP12_194_f_ah(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 214 // ./aflow --proto=ABC_hP6_194_c_d_a --params=2.752,2.56468023256
    if(vlabel.at(ifound)=="ABC_hP6_194_c_d_a") {
      PrototypeANRL_ABC_hP6_194_c_d_a(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 215 // ./aflow --proto=A_hP4_194_f --params=2.508,1.66786283892,0.05995
    if(vlabel.at(ifound)=="A_hP4_194_f") {
      PrototypeANRL_A_hP4_194_f(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 216 // ./aflow --proto=AB2_hP6_194_c_ad --params=4.186,1.22527472527
    if(vlabel.at(ifound)=="AB2_hP6_194_c_ad") {
      PrototypeANRL_AB2_hP6_194_c_ad(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 217 // ./aflow --proto=AB3C4_hP16_194_c_af_ef --params=2.988,7.82195448461,0.1543,0.605,0.0539
    if(vlabel.at(ifound)=="AB3C4_hP16_194_c_af_ef") {
      PrototypeANRL_AB3C4_hP16_194_c_af_ef(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 218 // ./aflow --proto=A_hP2_194_c --params=3.2093,1.62359393014
    if(vlabel.at(ifound)=="A_hP2_194_c") {
      PrototypeANRL_A_hP2_194_c(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 219 // ./aflow --proto=AB2_hP24_194_ef_fgh --params=4.824,3.28067993367,0.04598,0.84417,0.12514,0.16429
    if(vlabel.at(ifound)=="AB2_hP24_194_ef_fgh") {
      PrototypeANRL_AB2_hP24_194_ef_fgh(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 220 // ./aflow --proto=AB_hP12_194_df_ce --params=3.976,4.12022132797,0.0637,0.10724
    if(vlabel.at(ifound)=="AB_hP12_194_df_ce") {
      PrototypeANRL_AB_hP12_194_df_ce(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 221 // ./aflow --proto=AB_hP4_194_c_a --params=3.619,1.39375518099
    if(vlabel.at(ifound)=="AB_hP4_194_c_a") {
      PrototypeANRL_AB_hP4_194_c_a(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 222 // ./aflow --proto=A2B_hP12_194_cg_f --params=5.052,1.63697545527,0.062
    if(vlabel.at(ifound)=="A2B_hP12_194_cg_f") {
      PrototypeANRL_A2B_hP12_194_cg_f(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 223 // ./aflow --proto=A4B_cI40_197_cde_c --params=8.4295,0.1668,0.3345,0.6476,0.7484
    if(vlabel.at(ifound)=="A4B_cI40_197_cde_c") {
      PrototypeANRL_A4B_cI40_197_cde_c(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 224 // ./aflow --proto=ABC_cP12_198_a_a_a --params=5.881,-0.024,0.39,0.875
    if(vlabel.at(ifound)=="ABC_cP12_198_a_a_a") {
      PrototypeANRL_ABC_cP12_198_a_a_a(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 225 // ./aflow --proto=A3B_cP16_198_b_a --params=5.1305,0.2107,0.3689,0.2671,0.1159
    if(vlabel.at(ifound)=="A3B_cP16_198_b_a") {
      PrototypeANRL_A3B_cP16_198_b_a(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 226 // ./aflow --proto=A_cP8_198_2a --params=5.65,0.0699,-0.0378
    if(vlabel.at(ifound)=="A_cP8_198_2a") {
      PrototypeANRL_A_cP8_198_2a(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 227 // ./aflow --proto=AB_cP8_198_a_a --params=5.63,-0.042,0.067
    // 228 // ./aflow --proto=AB_cP8_198_a_a --params=4.48688,0.13652,0.8424
    if(vlabel.at(ifound)=="AB_cP8_198_a_a") {
      PrototypeANRL_AB_cP8_198_a_a(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 229 // ./aflow --proto=AB_cI16_199_a_a --params=6.3557,0.294,0.0347
    if(vlabel.at(ifound)=="AB_cI16_199_a_a") {
      PrototypeANRL_AB_cI16_199_a_a(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 230 // ./aflow --proto=AB32C48_cI162_204_a_2efg_2gh --params=14.16,0.8203,0.5998,0.1836,0.2942,0.8806,0.0908,0.8499,0.1748,0.6993,0.686,0.0969,0.332
    if(vlabel.at(ifound)=="AB32C48_cI162_204_a_2efg_2gh") {
      PrototypeANRL_AB32C48_cI162_204_a_2efg_2gh(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 231 // ./aflow --proto=A3B_cI32_204_g_c --params=7.58,0.3431,0.8497
    if(vlabel.at(ifound)=="A3B_cI32_204_g_c") {
      PrototypeANRL_A3B_cI32_204_g_c(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 232 // ./aflow --proto=A12B_cI26_204_g_a --params=7.58,0.184,0.691
    if(vlabel.at(ifound)=="A12B_cI26_204_g_a") {
      PrototypeANRL_A12B_cI26_204_g_a(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 233 // ./aflow --proto=A_cP8_205_c --params=5.65,0.05569
    if(vlabel.at(ifound)=="A_cP8_205_c") {
      PrototypeANRL_A_cP8_205_c(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 234 // ./aflow --proto=AB_cP16_205_c_c --params=6.4162,0.1527,0.6297
    if(vlabel.at(ifound)=="AB_cP16_205_c_c") {
      PrototypeANRL_AB_cP16_205_c_c(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 235 // ./aflow --proto=AB2_cP12_205_a_c --params=5.417,0.3851
    if(vlabel.at(ifound)=="AB2_cP12_205_a_c") {
      PrototypeANRL_AB2_cP12_205_a_c(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 236 // ./aflow --proto=AB3C6_cI80_206_a_d_e --params=9.4,-0.0344,0.338,0.1,0.125
    if(vlabel.at(ifound)=="AB3C6_cI80_206_a_d_e") {
      PrototypeANRL_AB3C6_cI80_206_a_d_e(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 237 // ./aflow --proto=A_cI16_206_c --params=4.11971,0.1001
    if(vlabel.at(ifound)=="A_cI16_206_c") {
      PrototypeANRL_A_cI16_206_c(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 238 // ./aflow --proto=A_cP20_213_cd --params=6.315,0.06361,0.20224
    if(vlabel.at(ifound)=="A_cP20_213_cd") {
      PrototypeANRL_A_cP20_213_cd(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 239 // ./aflow --proto=A3B4C_cP8_215_d_e_a --params=5.3912,0.2372
    if(vlabel.at(ifound)=="A3B4C_cP8_215_d_e_a") {
      PrototypeANRL_A3B4C_cP8_215_d_e_a(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 240 // ./aflow --proto=AB4_cP5_215_a_e --params=3.878,0.265
    if(vlabel.at(ifound)=="AB4_cP5_215_a_e") {
      PrototypeANRL_AB4_cP5_215_a_e(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 241 // ./aflow --proto=AB3C4_cP8_215_a_c_e --params=5.28,0.25
    if(vlabel.at(ifound)=="AB3C4_cP8_215_a_c_e") {
      PrototypeANRL_AB3C4_cP8_215_a_c_e(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 242 // ./aflow --proto=AB5_cF24_216_a_ce --params=6.1,0.625
    if(vlabel.at(ifound)=="AB5_cF24_216_a_ce") {
      PrototypeANRL_AB5_cF24_216_a_ce(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 243 // ./aflow --proto=ABC_cF12_216_b_c_a --params=6.24
    if(vlabel.at(ifound)=="ABC_cF12_216_b_c_a") {
      PrototypeANRL_ABC_cF12_216_b_c_a(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 244 // ./aflow --proto=AB_cF8_216_c_a --params=5.4093
    if(vlabel.at(ifound)=="AB_cF8_216_c_a") {
      PrototypeANRL_AB_cF8_216_c_a(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 245 // ./aflow --proto=A4B_cI10_217_c_a --params=5.45858,0.165
    if(vlabel.at(ifound)=="A4B_cI10_217_c_a") {
      PrototypeANRL_A4B_cI10_217_c_a(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 246 // ./aflow --proto=A_cI58_217_ac2g --params=8.911,0.31787,-0.08958,0.28194,0.64294,0.03457
    if(vlabel.at(ifound)=="A_cI58_217_ac2g") {
      PrototypeANRL_A_cI58_217_ac2g(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 247 // ./aflow --proto=A5B8_cI52_217_ce_cg --params=8.8664,0.32774,0.10781,0.64421,0.68844,0.03674
    if(vlabel.at(ifound)=="A5B8_cI52_217_ce_cg") {
      PrototypeANRL_A5B8_cI52_217_ce_cg(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 248 // ./aflow --proto=A_cI16_220_c --params=5.2716,0.049
    if(vlabel.at(ifound)=="A_cI16_220_c") {
      PrototypeANRL_A_cI16_220_c(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 249 // ./aflow --proto=A3B2_cI40_220_d_c --params=8.135,0.0492,0.2896
    if(vlabel.at(ifound)=="A3B2_cI40_220_d_c") {
      PrototypeANRL_A3B2_cI40_220_d_c(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 250 // ./aflow --proto=AB_cP2_221_b_a --params=4.07925
    if(vlabel.at(ifound)=="AB_cP2_221_b_a") {
      PrototypeANRL_AB_cP2_221_b_a(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 251 // ./aflow --proto=AB_cP6_221_c_d --params=4.2101
    if(vlabel.at(ifound)=="AB_cP6_221_c_d") {
      PrototypeANRL_AB_cP6_221_c_d(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 252 // ./aflow --proto=AB3C_cP5_221_a_c_b --params=3.795
    if(vlabel.at(ifound)=="AB3C_cP5_221_a_c_b") {
      PrototypeANRL_AB3C_cP5_221_a_c_b(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 253 // ./aflow --proto=AB27CD3_cP32_221_a_dij_b_c --params=7.04,0.245,0.26
    if(vlabel.at(ifound)=="AB27CD3_cP32_221_a_dij_b_c") {
      PrototypeANRL_AB27CD3_cP32_221_a_dij_b_c(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 254 // ./aflow --proto=AB3_cP4_221_a_c --params=3.7402
    if(vlabel.at(ifound)=="AB3_cP4_221_a_c") {
      PrototypeANRL_AB3_cP4_221_a_c(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 255 // ./aflow --proto=A_cP1_221_a --params=3.34
    if(vlabel.at(ifound)=="A_cP1_221_a") {
      PrototypeANRL_A_cP1_221_a(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 256 // ./aflow --proto=AB11_cP36_221_c_agij --params=9.6,0.345,0.225,0.115
    if(vlabel.at(ifound)=="AB11_cP36_221_c_agij") {
      PrototypeANRL_AB11_cP36_221_c_agij(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 257 // ./aflow --proto=AB11CD3_cP16_221_a_dg_b_c --params=5.74,0.245
    if(vlabel.at(ifound)=="AB11CD3_cP16_221_a_dg_b_c") {
      PrototypeANRL_AB11CD3_cP16_221_a_dg_b_c(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 258 // ./aflow --proto=A3B_cP4_221_d_a --params=3.734
    if(vlabel.at(ifound)=="A3B_cP4_221_d_a") {
      PrototypeANRL_A3B_cP4_221_d_a(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 259 // ./aflow --proto=A6B_cP7_221_f_a --params=4.145,0.2117
    if(vlabel.at(ifound)=="A6B_cP7_221_f_a") {
      PrototypeANRL_A6B_cP7_221_f_a(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 260 // ./aflow --proto=A3B_cP8_223_c_a --params=4.556
    if(vlabel.at(ifound)=="A3B_cP8_223_c_a") {
      PrototypeANRL_A3B_cP8_223_c_a(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 261 // ./aflow --proto=A_cP46_223_dik --params=10.355,0.1837,0.1172,0.3077
    if(vlabel.at(ifound)=="A_cP46_223_dik") {
      PrototypeANRL_A_cP46_223_dik(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 262 // ./aflow --proto=A2B_cP6_224_b_a --params=4.267
    if(vlabel.at(ifound)=="A2B_cP6_224_b_a") {
      PrototypeANRL_A2B_cP6_224_b_a(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 263 // ./aflow --proto=A7B_cF32_225_bd_a --params=9.45
    if(vlabel.at(ifound)=="A7B_cF32_225_bd_a") {
      PrototypeANRL_A7B_cF32_225_bd_a(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 264 // ./aflow --proto=AB3_cF16_225_a_bc --params=5.853
    if(vlabel.at(ifound)=="AB3_cF16_225_a_bc") {
      PrototypeANRL_AB3_cF16_225_a_bc(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 265 // ./aflow --proto=A9B16C7_cF128_225_acd_2f_be --params=11.48,0.25,0.875,0.625
    if(vlabel.at(ifound)=="A9B16C7_cF128_225_acd_2f_be") {
      PrototypeANRL_A9B16C7_cF128_225_acd_2f_be(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 266 // ./aflow --proto=A12B_cF52_225_i_a --params=7.468,0.666
    if(vlabel.at(ifound)=="A12B_cF52_225_i_a") {
      PrototypeANRL_A12B_cF52_225_i_a(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 267 // ./aflow --proto=AB2_cF12_225_a_c --params=5.4631
    if(vlabel.at(ifound)=="AB2_cF12_225_a_c") {
      PrototypeANRL_AB2_cF12_225_a_c(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 268 // ./aflow --proto=A6B23_cF116_225_e_acfh --params=10.65,0.2765,0.6191,0.6699
    if(vlabel.at(ifound)=="A6B23_cF116_225_e_acfh") {
      PrototypeANRL_A6B23_cF116_225_e_acfh(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 269 // ./aflow --proto=AB2C_cF16_225_a_c_b --params=5.95
    if(vlabel.at(ifound)=="AB2C_cF16_225_a_c_b") {
      PrototypeANRL_AB2C_cF16_225_a_c_b(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 270 // ./aflow --proto=A_cF4_225_a --params=3.61491
    if(vlabel.at(ifound)=="A_cF4_225_a") {
      PrototypeANRL_A_cF4_225_a(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 271 // ./aflow --proto=AB18C8_cF108_225_a_eh_f --params=10.56,0.325,0.65833,0.66
    if(vlabel.at(ifound)=="AB18C8_cF108_225_a_eh_f") {
      PrototypeANRL_AB18C8_cF108_225_a_eh_f(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 272 // ./aflow --proto=AB_cF8_225_a_b --params=5.63931
    if(vlabel.at(ifound)=="AB_cF8_225_a_b") {
      PrototypeANRL_AB_cF8_225_a_b(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 273 // ./aflow --proto=A2B_cF24_227_c_a --params=7.166
    if(vlabel.at(ifound)=="A2B_cF24_227_c_a") {
      PrototypeANRL_A2B_cF24_227_c_a(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 274 // ./aflow --proto=AB2_cF96_227_e_cf --params=11.278,0.215,0.44
    if(vlabel.at(ifound)=="AB2_cF96_227_e_cf") {
      PrototypeANRL_AB2_cF96_227_e_cf(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 275 // ./aflow --proto=AB_cF16_227_a_b --params=7.483
    if(vlabel.at(ifound)=="AB_cF16_227_a_b") {
      PrototypeANRL_AB_cF16_227_a_b(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 276 // ./aflow --proto=A_cF136_227_aeg --params=14.864,0.2624,0.1824,0.3701
    if(vlabel.at(ifound)=="A_cF136_227_aeg") {
      PrototypeANRL_A_cF136_227_aeg(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 277 // ./aflow --proto=A2B_cF24_227_d_a --params=7.02
    if(vlabel.at(ifound)=="A2B_cF24_227_d_a") {
      PrototypeANRL_A2B_cF24_227_d_a(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 278 // ./aflow --proto=A_cF8_227_a --params=3.55
    if(vlabel.at(ifound)=="A_cF8_227_a") {
      PrototypeANRL_A_cF8_227_a(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 279 // ./aflow --proto=A2BC4_cF56_227_d_a_e --params=8.0832,0.7376
    if(vlabel.at(ifound)=="A2BC4_cF56_227_d_a_e") {
      PrototypeANRL_A2BC4_cF56_227_d_a_e(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 280 // ./aflow --proto=AB2_cF48_227_c_e --params=8.6,0.245
    if(vlabel.at(ifound)=="AB2_cF48_227_c_e") {
      PrototypeANRL_AB2_cF48_227_c_e(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 281 // ./aflow --proto=AB3C3_cF112_227_c_de_f --params=11.087,0.7047,0.323
    if(vlabel.at(ifound)=="AB3C3_cF112_227_c_de_f") {
      PrototypeANRL_AB3C3_cF112_227_c_de_f(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 282 // ./aflow --proto=A_cI2_229_a --params=3.155
    if(vlabel.at(ifound)=="A_cI2_229_a") {
      PrototypeANRL_A_cI2_229_a(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 283 // ./aflow --proto=A3B_cI8_229_b_a --params=2.984
    if(vlabel.at(ifound)=="A3B_cI8_229_b_a") {
      PrototypeANRL_A3B_cI8_229_b_a(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 284 // ./aflow --proto=A4B3_cI14_229_c_b --params=6.226
    if(vlabel.at(ifound)=="A4B3_cI14_229_c_b") {
      PrototypeANRL_A4B3_cI14_229_c_b(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 285 // ./aflow --proto=A2B7_cI54_229_e_afh --params=11.618,0.6862,0.1704,0.6503
    if(vlabel.at(ifound)=="A2B7_cI54_229_e_afh") {
      PrototypeANRL_A2B7_cI54_229_e_afh(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 286 // ./aflow --proto=AB12C3_cI32_229_a_h_b --params=7.04,0.7625
    if(vlabel.at(ifound)=="AB12C3_cI32_229_a_h_b") {
      PrototypeANRL_AB12C3_cI32_229_a_h_b(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 287 // ./aflow --proto=AB4C3_cI16_229_a_c_b --params=5.74
    if(vlabel.at(ifound)=="AB4C3_cI16_229_a_c_b") {
      PrototypeANRL_AB4C3_cI16_229_a_c_b(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }
    // ---------------------------------------------------------------------------
    // 288 // ./aflow --proto=A4B3_cI112_230_af_g --params=11.411,0.0,0.625
    if(vlabel.at(ifound)=="A4B3_cI112_230_af_g") {
      PrototypeANRL_A4B3_cI112_230_af_g(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
    }

    // ---------------------------------------------------------------------------
    

    // DAVID PUT THEM HERE WITH THE ORDER OF THE PAPER

    //if(vlabel.at(ifound)=="YYY") {
    //  PrototypeANRL_YYY(oss,web,str,parameters,vproto.at(ifound),LDEBUG);
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
      str.atoms.at(iat).cpos=F2C(str.lattice,str.atoms.at(iat).fpos);
      str.num_each_type.at(str.atoms.at(iat).type)++;
      //     str.comp_each_type.at(str.atoms.at(iat).type)+=1.0; inside code
      str.species.at(str.atoms.at(iat).type)=str.atoms.at(iat).name;	
    }
    
    // ---------------------------------------------------------------------------
    // DONE
    xvector<double> data(6);
    data=Getabc_angles(str.lattice,DEGREES);
    str.a=data[1];str.b=data[2];str.c=data[3];str.alpha=data[4];str.beta=data[5];str.gamma=data[6];
    clear(str.origin);

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
				 string Pearson_symbol,uint spacegroup,string params) {
    if(vparameters_size!=nparameters) {
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
