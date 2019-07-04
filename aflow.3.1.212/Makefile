# ***************************************************************************
# *                                                                         *
# *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
# *                                                                         *
# ***************************************************************************
#
#  Copyright 2003-2018 - Stefano Curtarolo - AFLOW.ORG consortium
#
#  This file is part of AFLOW software.
#
#  AFLOW is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
# 
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# ***************************************************************************

VERSNUMBER=3.1.212

# MAKEFILE FOR AFLOW3
# FOR FROZSL -lmkl_lapack95_ilp64 -lmkl_core -lmkl_intel_lp64 -lmkl_sequential

# for aflow
CPP=g++
#-O2 -frepo
# for sqllite
CC=gcc

#CPP=g++-fsf-7 
#CPP=g++-fsf-4.9 -std=gnu++11
# -std=gnu++11
# -O3
#CPP=g++ -std=c++0x-4.7 -std=c++0x
# -std=gnu++0x
#CCFLAGS= -DCYGWIN -Wall -W -fno-inline
CCFLAGS=-Wall -W -fno-inline
CCFLAGS_MT=$(CCFLAGS) -O2
#CCFLAGS_GSA=$(CCFLAGS) -O1	#CO OBSOLETE 170601
INCLUDE = -I ./
#INCLUDE = -I ./  -I /opt/intel/mkl100/include

# -fno-inline -fno-implicit-inline-templates-frepo
OPTS=
#-O1
#-O3
# -march=opteron
TODAY=-DTODAY=\"$(TIME)\"
VERS=-DAFLOW_VERSION=\"$(VERSNUMBER)\" $(TODAY)

# -DHOSTNAME=\"`hostname`\" -DYYY=\"2003-2012\"

#BIN_MARYLOU=/fslhome/fslcollab8/group/bin
BIN_MARYLOU=/fslhome/fslcollab8/bin
#BIN_RANGER=/share/home/00470/tg457283/bin
BIN=~/bin
BIN_KRAKEN=/nics/a/proj/aflow/bin
BIN_OHAD=/home/aflow/bin

EXTRA_FREE=../AFLOW3_FREE/EXTRA
ICSD_LIB=$(EXTRA_FREE)/ICSD
NIST_LIB=/common/NIST/2014
AFLOW_LIB=/common/AFLOW/LIBS
#NIST_LIB=/common/NIST/2010
PROJECT_DIRECTORY_GNDSTATE=/common/LIB2

# -static
OPTSA=-O3
# -static
DAY=$(shell date +%y)$(shell date +%m)$(shell date +%d)
#TIME="20$(shell date +%y)-$(shell date +%m)-$(shell date +%d) $(shell date +%H):$(shell date +%M):$(shell date +%S)"
TIME="20$(shell date +%y)-$(shell date +%m)-$(shell date +%d)"
LIBS=-lpthread
# -static
#CO START 170614 easy integration
LIBS_AFLOW=-ldl $(LIBS)
LIBS_STATIC=-static
LIBS_AFLOW_STATIC=$(LIBS_AFLOW) $(LIBS_STATIC)
#CO END 170614 easy integration
# -lglut -lGLU -lGL -D_HOST_=\"`hostname`\"
#AFLOW_VERSION := $(shell cat aflow.cpp | /home/auro/bin/getvalue AFLOW_VERSION)-$(TIME)
#AFLOW_VERSIONM := $(shell cat Makefile | /home/auro/bin/getvalue MAKEFILE)-$(TIME)
#DEBUG=-DXSTR_DEBUG
DEBUG=
BACKUP=/common/AFLOW/SRC/

# Library MKL dynamic link
#LIB_MKL= -L /opt/intel/mkl100/lib/32 -lmkl_lapack -lmkl_core -lmkl_intel -lmkl_sequential
# Library MKL static link
#LIB_MKL= -L /opt/intel/mkl100/lib/32/ /opt/intel/mkl100/lib/32/libmkl_lapack.a \
# /opt/intel/mkl100/lib/32/libmkl_core.a \
# /opt/intel/mkl100/lib/32/libguide.a

AUROSTD_CPP = ./AUROSTD/aurostd.cpp ./AUROSTD/aurostd_main.cpp ./AUROSTD/aurostd_argv.cpp ./AUROSTD/aurostd_xoption.cpp ./AUROSTD/aurostd_boot.cpp ./AUROSTD/aurostd_crc64.cpp ./AUROSTD/aurostd_xscalar.cpp ./AUROSTD/aurostd_xcomplex.cpp ./AUROSTD/aurostd_xvector.cpp ./AUROSTD/aurostd_xmatrix.cpp ./AUROSTD/aurostd_xtensor.cpp ./AUROSTD/aurostd_xcombos.cpp ./AUROSTD/aurostd_xerror.cpp
AUROSTD_HPP = ./AUROSTD/aurostd.h ./AUROSTD/aurostd_argv.h ./AUROSTD/aurostd_xcomplex.h ./AUROSTD/aurostd_xmatrix.h ./AUROSTD/aurostd_xoption.h ./AUROSTD/aurostd_xrandom.h ./AUROSTD/aurostd_xscalar.h ./AUROSTD/aurostd_xtensor.h ./AUROSTD/aurostd_xvector.h ./AUROSTD/aurostd_xcombos.h AUROSTD/aurostd_xerror.h

#CO 170614
SQLITE_C = ./SQLITE/sqlite3.c

#AUROSTD_CPP = aurostd.cpp
#AUROSTD_HPP = aurostd.h

#DX and CO - START
AFLOW_HPP = aflow.h aflow_pflow.h aflow_agl_debye.h aflow_apennsy.h aflow_chull.h aflow_pocc.h aflow_pocc_old.h aflow_xgndstate.h aflow_xphases.h aflow_xvaspin.h APL/aflow_apl.h APL/aflow_qha_operations.h
#DX and CO - END
AFLOWLIB_HPP = aflowlib.h


AFLOW_CLASS_SRC = aflow_xclasses.cpp 
AFLOW_FUNCS_SRC = aflow_aflowrc.cpp aflow_init.cpp aflow_xatom.cpp aflow_kbin.cpp aflow_ialien.cpp aflow_kalien.cpp aflow_kvasp.cpp aflow_ivasp.cpp aflow_avasp.cpp aflow_ovasp.cpp aflow_kaims.cpp aflow_iaims.cpp aflow_oaims.cpp aflow_lattice.cpp aflow_spacegroup.cpp aflow_mix.cpp aflow_nomix.cpp aflow_qsub.cpp aflow_pthreads.cpp aflow_symmetry.cpp aflow_neighbours.cpp aflow_surface.cpp aflow_defects.cpp aflow_wyckoff.cpp aflow_xproto.cpp aflow_xproto_gus.cpp aflow_matlab.cpp aflow_matlab_funcs.cpp aflow_gnuplot_funcs.cpp aflow_ifrozsl.cpp
# aflow_ifrozsl.cpp
PFLOW_FUNCS_SRC = aflow_aconvasp_main.cpp aflow_pflow_funcs.cpp aflow_pflow_print.cpp aflow_aconvasp.cpp 
ESTRUCTURE_FUNCS_SRC = aflow_estructure.cpp 
#[OBSOLETE] [KESONG] POCCUPATION_FUNCS_SRC = aflow_poccupation_params.cpp aflow_poccupation_edos.cpp aflow_poccupation_forcefield.cpp 
POCCUPATION_FUNCS_SRC = aflow_pocc_old.cpp
APENNSY_FUNCS_SRC = aflow_apennsy_gndstate.cpp aflow_apennsy_structures.cpp aflow_apennsy.cpp aflow_apennsy_phases.cpp aflow_apennsy_main.cpp aflow_apennsy_vaspin.cpp
LIBRARIES_FUNCS_SRC = aflowlib_libraries.cpp 
CALCULATED_SRC = aflow_data_calculated.cpp
AFLOWLIB_SRC = aflowlib_web_interface.cpp aflowlib_web_outreach.cpp 
#DX and CO - START
APL_SRC = APL/aflow_apl_kphonons.cpp APL/aflow_qha_operations.cpp APL/aflow_apl_function_fitting.cpp APL/aflow_apl_dirphoncalc.cpp APL/aflow_apl_doscalc.cpp APL/aflow_apl_gsa.cpp APL/aflow_apl_hroutines.cpp APL/aflow_apl_logger.cpp APL/aflow_apl_lrphoncalc.cpp APL/aflow_apl_ltetdos.cpp APL/aflow_apl_mpmesh.cpp APL/aflow_apl_pathbuilder.cpp APL/aflow_apl_pdisc.cpp APL/aflow_apl_phoncalc.cpp APL/aflow_apl_rsmdos.cpp APL/aflow_apl_shellhandle.cpp APL/aflow_apl_supercell.cpp APL/aflow_apl_thermalpc.cpp APL/aflow_aapl_cluster.cpp APL/aflow_aapl_ifcs.cpp APL/aflow_aapl_setup.cpp APL/aflow_aapl_tcond.cpp APL/aflow_qha_aflowin_creator.cpp APL/aflow_qha_dm_pdos_save.cpp APL/aflow_apl_hsqpoints.cpp APL/aflow_apl_uniform_mesh.cpp APL/aflow_qha_gruneisen.cpp APL/aflow_apl_group_velocity.cpp APL/aflow_apl_atomic_disp.cpp APL/aflow_scqha_gruneisen.cpp APL/aflow_scqha_eos.cpp APL/aflow_qha_energies.cpp APL/aflow_qha3phonons_eos.cpp APL/aflow_qha_eos.cpp APL/aflow_scqha_T_freqs.cpp
#DX and CO - END

AGL_SRC = aflow_agl_debye.cpp aflow_agl_eqn_state.cpp aflow_agl_get_ev.cpp aflow_agl_plot_write.cpp aflow_agl_polynomial.cpp aflow_agl_rungibbs.cpp aflow_agl_thermal.cpp aflow_agl_hugoniot.cpp
AEL_SRC = aflow_ael_elasticity.cpp aflow_ael_get_stress.cpp aflow_ael_elastic_fit.cpp 

ANRL_SRC = aflow_anrl.cpp ANRL/aflow_anrl_list.cpp ANRL/aflow_anrl_AB2_aP12_1_4a_8a.cpp ANRL/aflow_anrl_ABC2_aP16_1_4a_4a_8a.cpp ANRL/aflow_anrl_A2B_aP6_2_2i_i.cpp ANRL/aflow_anrl_A_aP4_2_aci.cpp ANRL/aflow_anrl_A2B_mP12_3_bc3e_2e.cpp ANRL/aflow_anrl_A_mP4_4_2a.cpp ANRL/aflow_anrl_A_mC12_5_3c.cpp ANRL/aflow_anrl_A3BC_mC10_8_ab_a_a.cpp ANRL/aflow_anrl_A2B_mC144_9_24a_12a.cpp ANRL/aflow_anrl_AB_mP4_11_e_e.cpp ANRL/aflow_anrl_ABC3_mP10_11_e_e_ef.cpp ANRL/aflow_anrl_A_mP16_11_8e.cpp ANRL/aflow_anrl_AB2_mC6_12_a_i.cpp ANRL/aflow_anrl_A_mC34_12_ah3i2j.cpp ANRL/aflow_anrl_AB3_mC16_12_g_ij.cpp ANRL/aflow_anrl_A5B2_mC14_12_a2i_i.cpp ANRL/aflow_anrl_A_mC4_12_i.cpp ANRL/aflow_anrl_ABC4_mP12_13_e_a_2g.cpp ANRL/aflow_anrl_A_mP84_13_21g.cpp ANRL/aflow_anrl_A2B_mP12_14_2e_e.cpp ANRL/aflow_anrl_A_mP32_14_8e.cpp ANRL/aflow_anrl_A_mP64_14_16e.cpp ANRL/aflow_anrl_A2B5_mC28_15_f_e2f.cpp ANRL/aflow_anrl_AB_mC8_15_c_e.cpp ANRL/aflow_anrl_A2B_mC48_15_ae3f_2f.cpp ANRL/aflow_anrl_ABC6D2_mC40_15_e_e_3f_f.cpp ANRL/aflow_anrl_ABC4_oP12_16_ag_cd_2u.cpp ANRL/aflow_anrl_AB3_oP16_18_ab_3c.cpp ANRL/aflow_anrl_A2B_oP12_19_2a_a.cpp ANRL/aflow_anrl_A2B_oC24_20_abc_c.cpp ANRL/aflow_anrl_AB_oP2_25_b_a.cpp ANRL/aflow_anrl_AB2_oP24_28_acd_2c3d.cpp ANRL/aflow_anrl_AB3C4_oP16_31_a_ab_2ab.cpp ANRL/aflow_anrl_AB_oP8_33_a_a.cpp ANRL/aflow_anrl_AB3C4_oP32_33_a_3a_4a.cpp ANRL/aflow_anrl_A2B_oC12_36_2a_a.cpp ANRL/aflow_anrl_A2BC_oC8_38_e_a_b.cpp ANRL/aflow_anrl_A2B_oC12_38_de_ab.cpp ANRL/aflow_anrl_AB4_oC20_41_a_2b.cpp ANRL/aflow_anrl_AB2_oC24_41_2a_2b.cpp ANRL/aflow_anrl_AB2_oF72_43_ab_3b.cpp ANRL/aflow_anrl_AB_oI4_44_a_b.cpp ANRL/aflow_anrl_A2B3C7D_oP13_47_t_aq_eqrs_h.cpp ANRL/aflow_anrl_AB_oP4_51_e_f.cpp ANRL/aflow_anrl_A3B2_oP20_56_ce_e.cpp ANRL/aflow_anrl_ABCD_oP16_57_d_c_d_d.cpp ANRL/aflow_anrl_AB_oP8_57_d_d.cpp ANRL/aflow_anrl_AB2_oP6_58_a_g.cpp ANRL/aflow_anrl_AB2_oP6_58_a_g.cpp ANRL/aflow_anrl_AB2_oP6_58_a_g.cpp ANRL/aflow_anrl_AB_oP4_59_a_b.cpp ANRL/aflow_anrl_ABC_oP6_59_a_a_a.cpp ANRL/aflow_anrl_A3B_oP8_59_bf_a.cpp ANRL/aflow_anrl_AB_oP16_61_c_c.cpp ANRL/aflow_anrl_A2B_oP24_61_2c_c.cpp ANRL/aflow_anrl_A3B2_oP20_62_3c_2c.cpp ANRL/aflow_anrl_AB3C_oP20_62_c_cd_a.cpp ANRL/aflow_anrl_A4B_oP20_62_2cd_c.cpp ANRL/aflow_anrl_AB2C_oP16_62_c_2c_c.cpp ANRL/aflow_anrl_A2B_oP12_62_2c_c.cpp ANRL/aflow_anrl_A2B_oP12_62_2c_c.cpp ANRL/aflow_anrl_A2B_oP12_62_2c_c.cpp ANRL/aflow_anrl_AB_oP8_62_c_c.cpp ANRL/aflow_anrl_AB_oP8_62_c_c.cpp ANRL/aflow_anrl_AB3_oP16_62_c_cd.cpp ANRL/aflow_anrl_A3B7_oP40_62_cd_3c2d.cpp ANRL/aflow_anrl_A_oP8_62_2c.cpp ANRL/aflow_anrl_AB_oP8_62_c_c.cpp ANRL/aflow_anrl_AB_oP8_62_c_c.cpp ANRL/aflow_anrl_AB2C_oC16_63_c_2c_c.cpp ANRL/aflow_anrl_A2B_oC12_63_2c_c.cpp ANRL/aflow_anrl_AB_oC8_63_c_c.cpp ANRL/aflow_anrl_A_oC4_63_c.cpp ANRL/aflow_anrl_A_oC8_64_f.cpp ANRL/aflow_anrl_A2B2C_oC80_64_efg_efg_df.cpp ANRL/aflow_anrl_A_oC8_64_f.cpp ANRL/aflow_anrl_A_oC8_64_f.cpp ANRL/aflow_anrl_AB_oC8_65_j_g.cpp ANRL/aflow_anrl_A3B5_oC16_65_ah_bej.cpp ANRL/aflow_anrl_AB3_oC8_65_a_bf.cpp ANRL/aflow_anrl_AB_oF8_69_a_b.cpp ANRL/aflow_anrl_A_oF8_70_a.cpp ANRL/aflow_anrl_A2B_oF24_70_e_a.cpp ANRL/aflow_anrl_A_oF128_70_4h.cpp ANRL/aflow_anrl_AB2_oI6_71_a_i.cpp ANRL/aflow_anrl_AB2_oI6_71_a_g.cpp ANRL/aflow_anrl_A2B_oI12_72_j_a.cpp ANRL/aflow_anrl_AB4C_tI12_82_c_g_a.cpp ANRL/aflow_anrl_A2BC4_tI14_82_bc_a_g.cpp ANRL/aflow_anrl_AB_tP16_84_cej_k.cpp ANRL/aflow_anrl_A4B5_tI18_87_h_ah.cpp ANRL/aflow_anrl_AB4_tI10_87_a_h.cpp ANRL/aflow_anrl_A2B_tP12_92_b_a.cpp ANRL/aflow_anrl_A2B_tP36_96_3b_ab.cpp ANRL/aflow_anrl_A_tP12_96_ab.cpp ANRL/aflow_anrl_A3BC_tP5_99_bc_a_b.cpp ANRL/aflow_anrl_AB3_tP8_113_a_ce.cpp ANRL/aflow_anrl_A2BC4D_tI16_121_d_a_i_b.cpp ANRL/aflow_anrl_ABC2_tI16_122_a_b_d.cpp ANRL/aflow_anrl_AB5C_tP7_123_b_ci_a.cpp ANRL/aflow_anrl_AB3_tP4_123_a_ce.cpp ANRL/aflow_anrl_AB_tP2_123_a_d.cpp ANRL/aflow_anrl_ABC2_tP4_123_d_a_f.cpp ANRL/aflow_anrl_A2B3_tP10_127_g_ah.cpp ANRL/aflow_anrl_ABCD_tP8_129_c_b_a_c.cpp ANRL/aflow_anrl_A_tP4_129_ac.cpp ANRL/aflow_anrl_ABC_tP6_129_c_a_c.cpp ANRL/aflow_anrl_A2B_tP6_129_ac_c.cpp ANRL/aflow_anrl_AB_tP4_129_a_c.cpp ANRL/aflow_anrl_AB_tP4_129_c_c.cpp ANRL/aflow_anrl_AB_tP4_131_c_e.cpp ANRL/aflow_anrl_A_tP50_134_b2m2n.cpp ANRL/aflow_anrl_A_tP30_136_bf2ij.cpp ANRL/aflow_anrl_AB_tP8_136_g_f.cpp ANRL/aflow_anrl_A2B_tP6_136_f_a.cpp ANRL/aflow_anrl_sigma_tP30_136_bf2ij.cpp ANRL/aflow_anrl_A_tP4_136_f.cpp ANRL/aflow_anrl_A_tP16_138_j.cpp ANRL/aflow_anrl_A3B_tI16_139_cde_e.cpp ANRL/aflow_anrl_A_tI4_139_e.cpp ANRL/aflow_anrl_AB2C4_tI14_139_a_e_ce.cpp ANRL/aflow_anrl_A12B_tI26_139_fij_a.cpp ANRL/aflow_anrl_A_tI2_139_a.cpp ANRL/aflow_anrl_A_tI8_139_h.cpp ANRL/aflow_anrl_A3B_tI8_139_bd_a.cpp ANRL/aflow_anrl_AB2_tI6_139_a_e.cpp ANRL/aflow_anrl_A4B5_tI18_139_i_ah.cpp ANRL/aflow_anrl_A4B_tI10_139_de_a.cpp ANRL/aflow_anrl_A8B_tI18_139_hi_a.cpp ANRL/aflow_anrl_A2B_tI6_139_d_a.cpp ANRL/aflow_anrl_A_tI2_139_a.cpp ANRL/aflow_anrl_A2B_tI12_140_h_a.cpp ANRL/aflow_anrl_AB3_tI16_140_b_ah.cpp ANRL/aflow_anrl_AB_tI16_140_ab_h.cpp ANRL/aflow_anrl_A4BC_tI24_141_h_b_a.cpp ANRL/aflow_anrl_A_tI4_141_a.cpp ANRL/aflow_anrl_A3B4_tI28_141_ad_h.cpp ANRL/aflow_anrl_A2B_tI12_141_e_a.cpp ANRL/aflow_anrl_AB_tI16_141_e_e.cpp ANRL/aflow_anrl_A2B_tI24_141_2e_e.cpp ANRL/aflow_anrl_AB_tI8_141_a_b.cpp ANRL/aflow_anrl_A2B3_tI80_141_ceh_3h.cpp ANRL/aflow_anrl_ABC4_tI96_142_e_ab_2g.cpp ANRL/aflow_anrl_A2B_hP9_147_g_ad.cpp ANRL/aflow_anrl_AB_hR16_148_cf_cf.cpp ANRL/aflow_anrl_AB3_hR8_148_c_f.cpp ANRL/aflow_anrl_AB_hR26_148_b2f_a2f.cpp ANRL/aflow_anrl_AB3C_hR10_148_c_f_c.cpp ANRL/aflow_anrl_A2B_hP9_150_ef_bd.cpp ANRL/aflow_anrl_A3B_hP24_151_3c_2a.cpp ANRL/aflow_anrl_A2B_hP9_152_c_a.cpp ANRL/aflow_anrl_A_hP3_152_a.cpp ANRL/aflow_anrl_AB_hP6_154_a_b.cpp ANRL/aflow_anrl_AB3_hR8_155_c_de.cpp ANRL/aflow_anrl_A3B2_hR5_155_e_c.cpp ANRL/aflow_anrl_AB_hR6_160_b_b.cpp ANRL/aflow_anrl_AB_hR6_160_3a_3a.cpp ANRL/aflow_anrl_ABC3_hR10_161_a_a_b.cpp ANRL/aflow_anrl_AB2_hP9_162_ad_k.cpp ANRL/aflow_anrl_AB2CD2_hP36_163_h_i_bf_i.cpp ANRL/aflow_anrl_A3B2_hP5_164_ad_d.cpp ANRL/aflow_anrl_AB2_hP3_164_a_d.cpp ANRL/aflow_anrl_A3B_hP24_165_adg_f.cpp ANRL/aflow_anrl_AB_hR2_166_a_b.cpp ANRL/aflow_anrl_A_hR2_166_c.cpp ANRL/aflow_anrl_A_hR1_166_a.cpp ANRL/aflow_anrl_A7B6_hR13_166_ah_3c.cpp ANRL/aflow_anrl_A_hR3_166_ac.cpp ANRL/aflow_anrl_A2B3_hR5_166_c_ac.cpp ANRL/aflow_anrl_A_hR1_166_a.cpp ANRL/aflow_anrl_A5B2_hR7_166_a2c_c.cpp ANRL/aflow_anrl_A_hR2_166_c.cpp ANRL/aflow_anrl_A_hR12_166_2h.cpp ANRL/aflow_anrl_ABC2_hR4_166_a_b_c.cpp ANRL/aflow_anrl_A_hR2_166_c.cpp ANRL/aflow_anrl_A_hR105_166_bc9h4i.cpp ANRL/aflow_anrl_A6B_hR7_166_g_a.cpp ANRL/aflow_anrl_ABC3_hR10_167_a_b_e.cpp ANRL/aflow_anrl_ABC3_hR10_167_a_b_e.cpp ANRL/aflow_anrl_A2B3_hR10_167_c_e.cpp ANRL/aflow_anrl_A2B_hP18_180_fi_bd.cpp ANRL/aflow_anrl_AB2_hP9_180_d_j.cpp ANRL/aflow_anrl_A2B_hP9_180_j_c.cpp ANRL/aflow_anrl_AB3_hP8_182_c_g.cpp ANRL/aflow_anrl_A_hP4_186_ab.cpp ANRL/aflow_anrl_AB_hP8_186_ab_ab.cpp ANRL/aflow_anrl_AB_hP4_186_b_b.cpp ANRL/aflow_anrl_AB_hP12_186_a2b_a2b.cpp ANRL/aflow_anrl_A5B3C_hP18_186_2a3b_2ab_b.cpp ANRL/aflow_anrl_AB_hP4_186_b_a.cpp ANRL/aflow_anrl_ABC_hP3_187_a_d_f.cpp ANRL/aflow_anrl_AB_hP2_187_d_a.cpp ANRL/aflow_anrl_A2B_hP9_189_fg_bc.cpp ANRL/aflow_anrl_AB4C_hP6_191_a_h_b.cpp ANRL/aflow_anrl_AB5_hP6_191_a_cg.cpp ANRL/aflow_anrl_A_hP1_191_a.cpp ANRL/aflow_anrl_A3B_hP4_191_bc_a.cpp ANRL/aflow_anrl_AB2_hP3_191_a_d.cpp ANRL/aflow_anrl_A2B_hP6_191_h_e.cpp ANRL/aflow_anrl_AB_hP6_191_f_ad.cpp ANRL/aflow_anrl_AB_hP8_194_ad_f.cpp ANRL/aflow_anrl_A_hP6_194_h.cpp ANRL/aflow_anrl_AB_hP12_194_af_bf.cpp ANRL/aflow_anrl_A_hP4_194_ac.cpp ANRL/aflow_anrl_AB3_hP8_194_c_bf.cpp ANRL/aflow_anrl_AB2_hP6_194_b_f.cpp ANRL/aflow_anrl_AB_hP4_194_c_d.cpp ANRL/aflow_anrl_ABC2_hP8_194_d_a_f.cpp ANRL/aflow_anrl_A3B_hP8_194_h_c.cpp ANRL/aflow_anrl_A_hP4_194_bc.cpp ANRL/aflow_anrl_AB2_hP6_194_c_f.cpp ANRL/aflow_anrl_A5B2_hP14_194_abdf_f.cpp ANRL/aflow_anrl_AB2_hP12_194_f_ah.cpp ANRL/aflow_anrl_ABC_hP6_194_c_d_a.cpp ANRL/aflow_anrl_A_hP4_194_f.cpp ANRL/aflow_anrl_AB2_hP6_194_c_ad.cpp ANRL/aflow_anrl_AB3C4_hP16_194_c_af_ef.cpp ANRL/aflow_anrl_A_hP2_194_c.cpp ANRL/aflow_anrl_AB2_hP24_194_ef_fgh.cpp ANRL/aflow_anrl_AB_hP12_194_df_ce.cpp ANRL/aflow_anrl_AB_hP4_194_c_a.cpp ANRL/aflow_anrl_A2B_hP12_194_cg_f.cpp ANRL/aflow_anrl_A4B_cI40_197_cde_c.cpp ANRL/aflow_anrl_ABC_cP12_198_a_a_a.cpp ANRL/aflow_anrl_A3B_cP16_198_b_a.cpp ANRL/aflow_anrl_A_cP8_198_2a.cpp ANRL/aflow_anrl_AB_cP8_198_a_a.cpp ANRL/aflow_anrl_AB_cP8_198_a_a.cpp ANRL/aflow_anrl_AB_cI16_199_a_a.cpp ANRL/aflow_anrl_AB32C48_cI162_204_a_2efg_2gh.cpp ANRL/aflow_anrl_A3B_cI32_204_g_c.cpp ANRL/aflow_anrl_A12B_cI26_204_g_a.cpp ANRL/aflow_anrl_A_cP8_205_c.cpp ANRL/aflow_anrl_AB_cP16_205_c_c.cpp ANRL/aflow_anrl_AB2_cP12_205_a_c.cpp ANRL/aflow_anrl_AB3C6_cI80_206_a_d_e.cpp ANRL/aflow_anrl_A_cI16_206_c.cpp ANRL/aflow_anrl_A_cP20_213_cd.cpp ANRL/aflow_anrl_A3B4C_cP8_215_d_e_a.cpp ANRL/aflow_anrl_AB4_cP5_215_a_e.cpp ANRL/aflow_anrl_AB3C4_cP8_215_a_c_e.cpp ANRL/aflow_anrl_AB5_cF24_216_a_ce.cpp ANRL/aflow_anrl_ABC_cF12_216_b_c_a.cpp ANRL/aflow_anrl_AB_cF8_216_c_a.cpp ANRL/aflow_anrl_A4B_cI10_217_c_a.cpp ANRL/aflow_anrl_A_cI58_217_ac2g.cpp ANRL/aflow_anrl_A5B8_cI52_217_ce_cg.cpp ANRL/aflow_anrl_A_cI16_220_c.cpp ANRL/aflow_anrl_A3B2_cI40_220_d_c.cpp ANRL/aflow_anrl_AB_cP2_221_b_a.cpp ANRL/aflow_anrl_AB_cP6_221_c_d.cpp ANRL/aflow_anrl_AB3C_cP5_221_a_c_b.cpp ANRL/aflow_anrl_AB27CD3_cP32_221_a_dij_b_c.cpp ANRL/aflow_anrl_AB3_cP4_221_a_c.cpp ANRL/aflow_anrl_A_cP1_221_a.cpp ANRL/aflow_anrl_AB11_cP36_221_c_agij.cpp ANRL/aflow_anrl_AB11CD3_cP16_221_a_dg_b_c.cpp ANRL/aflow_anrl_A3B_cP4_221_d_a.cpp ANRL/aflow_anrl_A6B_cP7_221_f_a.cpp ANRL/aflow_anrl_A3B_cP8_223_c_a.cpp ANRL/aflow_anrl_A_cP46_223_dik.cpp ANRL/aflow_anrl_A2B_cP6_224_b_a.cpp ANRL/aflow_anrl_A7B_cF32_225_bd_a.cpp ANRL/aflow_anrl_AB3_cF16_225_a_bc.cpp ANRL/aflow_anrl_A9B16C7_cF128_225_acd_2f_be.cpp ANRL/aflow_anrl_A12B_cF52_225_i_a.cpp ANRL/aflow_anrl_AB2_cF12_225_a_c.cpp ANRL/aflow_anrl_A6B23_cF116_225_e_acfh.cpp ANRL/aflow_anrl_AB2C_cF16_225_a_c_b.cpp ANRL/aflow_anrl_A_cF4_225_a.cpp ANRL/aflow_anrl_AB18C8_cF108_225_a_eh_f.cpp ANRL/aflow_anrl_AB_cF8_225_a_b.cpp ANRL/aflow_anrl_A2B_cF24_227_c_a.cpp ANRL/aflow_anrl_AB2_cF96_227_e_cf.cpp ANRL/aflow_anrl_AB_cF16_227_a_b.cpp ANRL/aflow_anrl_A_cF136_227_aeg.cpp ANRL/aflow_anrl_A2B_cF24_227_d_a.cpp ANRL/aflow_anrl_A_cF8_227_a.cpp ANRL/aflow_anrl_A2BC4_cF56_227_d_a_e.cpp ANRL/aflow_anrl_AB2_cF48_227_c_e.cpp ANRL/aflow_anrl_AB3C3_cF112_227_c_de_f.cpp ANRL/aflow_anrl_A_cI2_229_a.cpp ANRL/aflow_anrl_A3B_cI8_229_b_a.cpp ANRL/aflow_anrl_A4B3_cI14_229_c_b.cpp ANRL/aflow_anrl_A2B7_cI54_229_e_afh.cpp ANRL/aflow_anrl_AB12C3_cI32_229_a_h_b.cpp ANRL/aflow_anrl_AB4C3_cI16_229_a_c_b.cpp ANRL/aflow_anrl_A4B3_cI112_230_af_g.cpp ANRL/aflow_anrl_A2B_aP6_2_aei_i.cpp ANRL/aflow_anrl_A8B5_mP13_6_a7b_3a2b.cpp ANRL/aflow_anrl_AB_mP4_6_2b_2a.cpp ANRL/aflow_anrl_A2B_mP12_7_4a_2a.cpp ANRL/aflow_anrl_A2B_mP18_7_6a_3a.cpp ANRL/aflow_anrl_A3B_mP16_7_6a_2a.cpp ANRL/aflow_anrl_A9B2_mP22_7_9a_2a.cpp ANRL/aflow_anrl_A5B3_mC32_9_5a_3a.cpp ANRL/aflow_anrl_AB3_mC16_9_a_3a.cpp ANRL/aflow_anrl_A2B_mP6_10_mn_bg.cpp ANRL/aflow_anrl_AB3_mP16_10_mn_3m3n.cpp ANRL/aflow_anrl_ABC2_mP8_10_ac_eh_mn.cpp ANRL/aflow_anrl_AB_mP6_10_en_am.cpp ANRL/aflow_anrl_A_mP8_10_2m2n.cpp ANRL/aflow_anrl_A7B2C2_mC22_12_aij_h_i.cpp ANRL/aflow_anrl_A_mC16_12_4i.cpp ANRL/aflow_anrl_A2B_mP12_13_2g_ef.cpp ANRL/aflow_anrl_A2B_mP6_14_e_a.cpp ANRL/aflow_anrl_A7B8_mP120_14_14e_16e.cpp ANRL/aflow_anrl_AB3_mC16_15_e_cf.cpp ANRL/aflow_anrl_A_mC24_15_2e2f.cpp ANRL/aflow_anrl_A2B_oP12_17_abe_e.cpp ANRL/aflow_anrl_AB3_oP16_19_a_3a.cpp ANRL/aflow_anrl_AB2_oC6_21_a_k.cpp ANRL/aflow_anrl_A2BC2_oF40_22_fi_ad_gh.cpp ANRL/aflow_anrl_AB_oF8_22_a_c.cpp ANRL/aflow_anrl_A3B_oI32_23_ij2k_k.cpp ANRL/aflow_anrl_A8B2C12D2E_oI50_23_bcfk_i_3k_j_a.cpp ANRL/aflow_anrl_ABC2_oI16_23_ab_i_k.cpp ANRL/aflow_anrl_ABC4_oI12_23_a_b_k.cpp ANRL/aflow_anrl_AB7CD2_oI44_24_a_b3d_c_ac.cpp ANRL/aflow_anrl_A2B_oP12_26_abc_ab.cpp ANRL/aflow_anrl_A2B_oP12_26_abc_ab.cpp ANRL/aflow_anrl_A5B_oP24_26_3a3b2c_ab.cpp ANRL/aflow_anrl_A6B4C16D_oP108_27_abcd4e_4e_16e_e.cpp ANRL/aflow_anrl_A2B_oP12_29_2a_a.cpp ANRL/aflow_anrl_AB2_oP12_29_a_2a.cpp ANRL/aflow_anrl_ABC_oP12_29_a_a_a.cpp ANRL/aflow_anrl_A5B3C15_oP46_30_a2c_bc_a7c.cpp ANRL/aflow_anrl_ABC3_oP20_30_2a_c_3c.cpp ANRL/aflow_anrl_A13B2C2_oP34_32_a6c_c_c.cpp ANRL/aflow_anrl_A2B3_oP40_33_4a_6a.cpp ANRL/aflow_anrl_A2B8C_oP22_34_c_4c_a.cpp ANRL/aflow_anrl_AB2_oP6_34_a_c.cpp ANRL/aflow_anrl_AB8C2_oC22_35_a_ab3e_e.cpp ANRL/aflow_anrl_AB_oC8_36_a_a.cpp ANRL/aflow_anrl_A2B5C2_oC36_37_d_c2d_d.cpp ANRL/aflow_anrl_A2B3_oC40_39_2d_2c2d.cpp ANRL/aflow_anrl_A9BC_oC44_39_3c3d_a_c.cpp ANRL/aflow_anrl_AB2C_oC16_40_a_2b_b.cpp ANRL/aflow_anrl_AB3_oC16_40_b_3b.cpp ANRL/aflow_anrl_A10B3_oF52_42_2abce_ab.cpp ANRL/aflow_anrl_AB_oF8_42_a_a.cpp ANRL/aflow_anrl_A2BC2_oI20_45_c_b_c.cpp ANRL/aflow_anrl_ABC_oI36_46_ac_bc_3b.cpp ANRL/aflow_anrl_A2B8CD_oP24_48_k_2m_d_b.cpp ANRL/aflow_anrl_A5B2_oP14_49_dehq_ab.cpp ANRL/aflow_anrl_AB2C8D_oP24_49_g_q_2qr_e.cpp ANRL/aflow_anrl_A2BC4_oP28_50_ij_ac_ijm.cpp ANRL/aflow_anrl_A3BC2_oP48_50_3m_m_2m.cpp ANRL/aflow_anrl_A2B_oP24_52_2e_cd.cpp ANRL/aflow_anrl_A3B2_oP20_52_de_cd.cpp ANRL/aflow_anrl_ABC2_oP16_53_h_e_gh.cpp ANRL/aflow_anrl_ABC3_oP20_53_e_g_hi.cpp ANRL/aflow_anrl_ABC3_oP20_54_e_d_cf.cpp ANRL/aflow_anrl_A2B_oP24_55_2g2h_gh.cpp ANRL/aflow_anrl_A3B5_oP16_55_ch_agh.cpp ANRL/aflow_anrl_A_oP16_55_2g2h.cpp ANRL/aflow_anrl_A2B_oP6_58_g_a.cpp ANRL/aflow_anrl_ABC_oP6_59_a_b_a.cpp ANRL/aflow_anrl_A2B3_oP20_60_d_cd.cpp ANRL/aflow_anrl_A3B_oP32_60_3d_d.cpp ANRL/aflow_anrl_A7B8_oP120_60_7d_8d.cpp ANRL/aflow_anrl_AB_oP48_61_3c_3c.cpp ANRL/aflow_anrl_A2B3_oP20_62_2c_3c.cpp ANRL/aflow_anrl_A2B4C_oP28_62_ac_2cd_c.cpp ANRL/aflow_anrl_A2B_oP12_62_2c_c.cpp ANRL/aflow_anrl_A3B_oP16_62_cd_c.cpp ANRL/aflow_anrl_AB2C3_oP24_62_c_d_cd.cpp ANRL/aflow_anrl_AB3_oP16_62_c_3c.cpp ANRL/aflow_anrl_AB4C_oP24_62_c_2cd_c.cpp ANRL/aflow_anrl_AB_oP8_62_c_c.cpp ANRL/aflow_anrl_A2BC3_oC24_63_e_c_cg.cpp ANRL/aflow_anrl_A43B5C17_oC260_63_c8fg6h_cfg_ce3f2h.cpp ANRL/aflow_anrl_A6B_oC28_63_efg_c.cpp ANRL/aflow_anrl_AB3C_oC20_63_a_cf_c.cpp ANRL/aflow_anrl_AB4C_oC24_63_a_fg_c.cpp ANRL/aflow_anrl_AB4C_oC24_63_c_fg_c.cpp ANRL/aflow_anrl_A2B_oC24_64_2f_f.cpp ANRL/aflow_anrl_A2B4C_oC28_66_l_kl_a.cpp ANRL/aflow_anrl_A3B_oC64_66_gi2lm_2l.cpp ANRL/aflow_anrl_A3B_oC64_66_kl2m_bdl.cpp ANRL/aflow_anrl_A2BC_oC16_67_ag_b_g.cpp ANRL/aflow_anrl_ABC2_oC16_67_b_g_ag.cpp ANRL/aflow_anrl_AB_oC8_67_a_g.cpp ANRL/aflow_anrl_AB_oC8_67_a_g.cpp ANRL/aflow_anrl_AB4_oC20_68_a_i.cpp ANRL/aflow_anrl_AB2_oF48_70_f_fg.cpp ANRL/aflow_anrl_A4B3_oI14_71_gh_cg.cpp ANRL/aflow_anrl_ABC_oI12_71_h_j_g.cpp ANRL/aflow_anrl_ABCD3_oI48_73_d_e_e_ef.cpp ANRL/aflow_anrl_A2B_oI12_74_h_e.cpp ANRL/aflow_anrl_A4B_oI20_74_beh_e.cpp ANRL/aflow_anrl_AB2C12D4_tP76_75_2a2b_2d_12d_4d.cpp ANRL/aflow_anrl_A2BC_tP16_76_2a_a_a.cpp ANRL/aflow_anrl_A3B7_tP40_76_3a_7a.cpp ANRL/aflow_anrl_A2B6CD7_tP64_77_2d_6d_d_ab6d.cpp ANRL/aflow_anrl_A2B_tP48_77_8d_4d.cpp ANRL/aflow_anrl_A2B7C2_tP88_78_4a_14a_4a.cpp ANRL/aflow_anrl_A2BC2_tI20_79_c_2a_c.cpp ANRL/aflow_anrl_AB2_tI48_80_2b_4b.cpp ANRL/aflow_anrl_AB2_tP12_81_adg_2h.cpp ANRL/aflow_anrl_A3B_tI32_82_3g_g.cpp ANRL/aflow_anrl_A3B2_tP10_83_adk_j.cpp ANRL/aflow_anrl_A2B_tP30_85_ab2g_cg.cpp ANRL/aflow_anrl_AB3_tP32_86_g_3g.cpp ANRL/aflow_anrl_A4B_tI20_88_f_a.cpp ANRL/aflow_anrl_AB2_tI96_88_2f_4f.cpp ANRL/aflow_anrl_A17BC4D_tP184_89_17p_p_4p_io.cpp ANRL/aflow_anrl_A4B2C13D_tP40_90_g_d_cef2g_c.cpp ANRL/aflow_anrl_AB4C17D4E_tP54_90_a_g_c4g_g_c.cpp ANRL/aflow_anrl_ABC_tP24_91_d_d_d.cpp ANRL/aflow_anrl_AB32CD4E8_tP184_93_i_16p_af_2p_4p.cpp ANRL/aflow_anrl_A14B3C5_tP44_94_c3g_ad_bg.cpp ANRL/aflow_anrl_A6B2C_tP18_94_eg_c_a.cpp ANRL/aflow_anrl_ABC_tP24_95_d_d_d.cpp ANRL/aflow_anrl_A2B8CD_tI24_97_d_k_a_b.cpp ANRL/aflow_anrl_AB8C2_tI44_97_e_2k_cd.cpp ANRL/aflow_anrl_A2B_tI12_98_f_a.cpp ANRL/aflow_anrl_A2B8C2D_tP26_100_c_abcd_c_a.cpp ANRL/aflow_anrl_A3B11C6_tP40_100_ac_bc2d_cd.cpp ANRL/aflow_anrl_A7B7C2_tP32_101_bde_ade_d.cpp ANRL/aflow_anrl_A2B3_tP20_102_2c_b2c.cpp ANRL/aflow_anrl_AB4_tP10_103_a_d.cpp ANRL/aflow_anrl_A5B5C4_tP28_104_ac_ac_c.cpp ANRL/aflow_anrl_AB6C4_tP22_104_a_2ac_c.cpp ANRL/aflow_anrl_A2BC2_tP20_105_f_ac_2e.cpp ANRL/aflow_anrl_A3BC3D_tP64_106_3c_c_3c_c.cpp ANRL/aflow_anrl_A5B7_tI24_107_ac_abd.cpp ANRL/aflow_anrl_AB_tI4_107_a_a.cpp ANRL/aflow_anrl_A3B5_tI32_108_ac_a2c.cpp ANRL/aflow_anrl_ABC_tI12_109_a_a_a.cpp ANRL/aflow_anrl_AB_tI8_109_a_a.cpp ANRL/aflow_anrl_A2BC8_tI176_110_2b_b_8b.cpp ANRL/aflow_anrl_A2B_tP12_111_2n_adf.cpp ANRL/aflow_anrl_AB_tP8_111_n_n.cpp ANRL/aflow_anrl_AB4C_tP12_112_b_n_e.cpp ANRL/aflow_anrl_A2BC7D2_tP24_113_e_a_cef_e.cpp ANRL/aflow_anrl_A3B_tP32_114_3e_e.cpp ANRL/aflow_anrl_A4B_tP10_114_e_a.cpp ANRL/aflow_anrl_A2B3_tP5_115_g_ag.cpp ANRL/aflow_anrl_AB2_tP12_115_j_egi.cpp ANRL/aflow_anrl_A2B3_tP20_116_bci_fj.cpp ANRL/aflow_anrl_A2B3_tP20_117_i_adgh.cpp ANRL/aflow_anrl_A3B_tP16_118_ei_f.cpp ANRL/aflow_anrl_A5B3_tP32_118_g2i_aceh.cpp ANRL/aflow_anrl_A3B_tI24_119_b2i_af.cpp ANRL/aflow_anrl_AB_tI4_119_c_a.cpp ANRL/aflow_anrl_A4BC2_tI28_120_i_d_e.cpp ANRL/aflow_anrl_A4BC4D_tP10_123_gh_a_i_d.cpp ANRL/aflow_anrl_AB4C_tP12_124_a_m_c.cpp ANRL/aflow_anrl_AB4_tP10_124_a_m.cpp ANRL/aflow_anrl_A4B_tP10_125_m_a.cpp ANRL/aflow_anrl_ABC4_tP12_125_a_b_m.cpp ANRL/aflow_anrl_A2BC4_tP28_126_cd_e_k.cpp ANRL/aflow_anrl_A4B_tP20_127_ehj_g.cpp ANRL/aflow_anrl_A6B2C_tP18_128_eh_d_a.cpp ANRL/aflow_anrl_A7B2C_tP40_128_egi_h_e.cpp ANRL/aflow_anrl_A2BC4_tP28_130_f_c_g.cpp ANRL/aflow_anrl_A5B3_tP32_130_cg_cf.cpp ANRL/aflow_anrl_A2B2C4D_tP18_132_e_i_o_d.cpp ANRL/aflow_anrl_AB6C_tP16_132_d_io_a.cpp ANRL/aflow_anrl_AB3_tP32_133_h_i2j.cpp ANRL/aflow_anrl_A2B_tP24_135_gh_h.cpp ANRL/aflow_anrl_A4B2C_tP28_135_gh_h_d.cpp ANRL/aflow_anrl_A2B3_tP40_137_cdf_3g.cpp ANRL/aflow_anrl_A2B_tP6_137_d_a.cpp ANRL/aflow_anrl_A4BC4_tP18_137_g_b_g.cpp ANRL/aflow_anrl_AB2_tP6_137_a_d.cpp ANRL/aflow_anrl_A_tP12_138_bi.cpp ANRL/aflow_anrl_AB_tI8_139_e_e.cpp ANRL/aflow_anrl_A3B5_tI32_140_ah_bk.cpp ANRL/aflow_anrl_A3B5_tI32_140_ah_cl.cpp ANRL/aflow_anrl_A2B_tI12_141_e_a.cpp ANRL/aflow_anrl_A_tI16_142_f.cpp ANRL/aflow_anrl_A4B14C3_hP21_143_bd_ac4d_d.cpp ANRL/aflow_anrl_A4B6C_hP11_143_bd_2d_a.cpp ANRL/aflow_anrl_AB2_hP12_143_cd_ab2d.cpp ANRL/aflow_anrl_A4B_hP15_144_4a_a.cpp ANRL/aflow_anrl_AB_hP6_144_a_a.cpp ANRL/aflow_anrl_A2B3C3DE7_hP48_145_2a_3a_3a_a_7a.cpp ANRL/aflow_anrl_A3BC_hR5_146_b_a_a.cpp ANRL/aflow_anrl_ABC3_hR10_146_2a_2a_2b.cpp ANRL/aflow_anrl_A2B4C_hR42_148_2f_4f_f.cpp ANRL/aflow_anrl_A2B_hR18_148_2f_f.cpp ANRL/aflow_anrl_AB3_hP24_149_acgi_3l.cpp ANRL/aflow_anrl_A3B_hP24_153_3c_2b.cpp ANRL/aflow_anrl_A_hP9_154_bc.cpp ANRL/aflow_anrl_AB2_hP9_156_b2c_3a2bc.cpp ANRL/aflow_anrl_AB_hP12_156_2ab3c_2ab3c.cpp ANRL/aflow_anrl_AB_hP4_156_ab_ab.cpp ANRL/aflow_anrl_A5B6C2_hP13_157_2ac_2c_b.cpp ANRL/aflow_anrl_A3B_hP8_158_d_a.cpp ANRL/aflow_anrl_A2B3_hP20_159_bc_2c.cpp ANRL/aflow_anrl_A4B3_hP28_159_ab2c_2c.cpp ANRL/aflow_anrl_AB4C7D_hP26_159_b_ac_a2c_b.cpp ANRL/aflow_anrl_A3B_hR4_160_b_a.cpp ANRL/aflow_anrl_A8B5_hR26_160_a3bc_a3b.cpp ANRL/aflow_anrl_ABC_hR3_160_a_a_a.cpp ANRL/aflow_anrl_AB_hR10_160_5a_5a.cpp ANRL/aflow_anrl_A2B3_hP5_164_d_ad.cpp ANRL/aflow_anrl_AB2_hP9_164_bd_c2d.cpp ANRL/aflow_anrl_ABC2_hP4_164_a_b_d.cpp ANRL/aflow_anrl_A3B_hP24_165_bdg_f.cpp ANRL/aflow_anrl_A4B3_hR7_166_2c_ac.cpp ANRL/aflow_anrl_ABC_hR6_166_c_c_c.cpp ANRL/aflow_anrl_AB3C_hR10_167_b_e_a.cpp ANRL/aflow_anrl_ABC2_hR24_167_e_e_2e.cpp ANRL/aflow_anrl_A2B13C4_hP57_168_d_c6d_2d.cpp ANRL/aflow_anrl_AB4C_hP72_168_2d_8d_2d.cpp ANRL/aflow_anrl_A2B3_hP30_169_2a_3a.cpp ANRL/aflow_anrl_A2B3_hP30_170_2a_3a.cpp ANRL/aflow_anrl_A10B2C_hP39_171_5c_c_a.cpp ANRL/aflow_anrl_A10B2C_hP39_172_5c_c_a.cpp ANRL/aflow_anrl_A3B_hP8_173_c_b.cpp ANRL/aflow_anrl_A4B3_hP14_173_bc_c.cpp ANRL/aflow_anrl_A12B7C2_hP21_174_2j2k_ajk_cf.cpp ANRL/aflow_anrl_ABC_hP12_174_cj_fk_aj.cpp ANRL/aflow_anrl_A8B7C6_hP21_175_ck_aj_k.cpp ANRL/aflow_anrl_ABC_hP36_175_jk_jk_jk.cpp ANRL/aflow_anrl_A3B2_hP10_176_h_bc.cpp ANRL/aflow_anrl_A3B3C_hP14_176_h_h_c.cpp ANRL/aflow_anrl_A3B_hP8_176_h_c.cpp ANRL/aflow_anrl_A2B_hP36_177_j2lm_n.cpp ANRL/aflow_anrl_AB3_hP24_178_b_ac.cpp ANRL/aflow_anrl_A_hP6_178_a.cpp ANRL/aflow_anrl_AB3_hP24_179_b_ac.cpp ANRL/aflow_anrl_A2B_hP9_181_j_c.cpp ANRL/aflow_anrl_ABC_hP3_183_a_a_a.cpp ANRL/aflow_anrl_AB_hP6_183_c_ab.cpp ANRL/aflow_anrl_AB4C_hP72_184_d_4d_d.cpp ANRL/aflow_anrl_A3BC_hP30_185_cd_c_ab.cpp ANRL/aflow_anrl_A3B_hP24_185_ab2c_c.cpp ANRL/aflow_anrl_A3B_hP8_185_c_a.cpp ANRL/aflow_anrl_AB3_hP24_185_c_ab2c.cpp ANRL/aflow_anrl_A3B7_hP20_186_c_b2c.cpp ANRL/aflow_anrl_AB3_hP4_187_e_fh.cpp ANRL/aflow_anrl_A3BC_hP10_188_k_c_a.cpp ANRL/aflow_anrl_AB9C4_hP28_188_e_kl_ak.cpp ANRL/aflow_anrl_A8BC3D6_hP18_189_bfh_a_g_i.cpp ANRL/aflow_anrl_A9BC3D5_hP18_189_fi_a_g_bh.cpp ANRL/aflow_anrl_A2B_hP18_190_gh_bf.cpp ANRL/aflow_anrl_A5B3_hP16_190_bdh_g.cpp ANRL/aflow_anrl_AB_hP24_190_i_afh.cpp ANRL/aflow_anrl_A2B3C18D6_hP58_192_c_f_lm_l.cpp ANRL/aflow_anrl_AB2_hP72_192_m_j2kl.cpp ANRL/aflow_anrl_A5B3_hP16_193_dg_g.cpp ANRL/aflow_anrl_A3B_hP16_194_gh_ac.cpp ANRL/aflow_anrl_A5B2_hP28_194_ahk_ch.cpp ANRL/aflow_anrl_A9B3C_hP26_194_hk_h_a.cpp ANRL/aflow_anrl_A12BC4_cP34_195_2j_ab_2e.cpp ANRL/aflow_anrl_A12B2C_cF60_196_h_bc_a.cpp ANRL/aflow_anrl_ABC3_cP20_198_a_a_b.cpp ANRL/aflow_anrl_A2B11_cP39_200_f_aghij.cpp ANRL/aflow_anrl_AB3C_cP60_201_be_fh_g.cpp ANRL/aflow_anrl_A6B6C_cF104_202_h_h_c.cpp ANRL/aflow_anrl_A_cF240_202_h2i.cpp ANRL/aflow_anrl_A2BCD3E6_cF208_203_e_c_d_f_g.cpp ANRL/aflow_anrl_A4B2C6D16E_cF232_203_e_d_f_eg_a.cpp ANRL/aflow_anrl_AB3C16_cF160_203_a_bc_eg.cpp ANRL/aflow_anrl_A2B3C6_cP264_205_2d_ab2c2d_6d.cpp ANRL/aflow_anrl_A_cP240_205_10d.cpp ANRL/aflow_anrl_AB3C2_cI96_206_c_e_ad.cpp ANRL/aflow_anrl_A17B15_cP64_207_acfk_eij.cpp ANRL/aflow_anrl_A3B_cP16_208_j_b.cpp ANRL/aflow_anrl_A6B2CD6E_cP64_208_m_ad_b_m_c.cpp ANRL/aflow_anrl_A24BC_cF104_209_j_a_b.cpp ANRL/aflow_anrl_A12B36CD12_cF488_210_h_3h_a_fg.cpp ANRL/aflow_anrl_A12B6C_cF608_210_4h_2h_e.cpp ANRL/aflow_anrl_A2B_cI72_211_hi_i.cpp ANRL/aflow_anrl_A2B_cP12_212_c_a.cpp ANRL/aflow_anrl_A3B3C_cI56_214_g_h_a.cpp ANRL/aflow_anrl_A3BC2_cI48_214_f_a_e.cpp ANRL/aflow_anrl_A4B9_cP52_215_ei_3efgi.cpp ANRL/aflow_anrl_ABCD_cF16_216_c_d_b_a.cpp ANRL/aflow_anrl_A3B4C_cP16_218_c_e_a.cpp ANRL/aflow_anrl_A7BC3D13_cF192_219_de_b_c_ah.cpp ANRL/aflow_anrl_A15B4_cI76_220_ae_c.cpp ANRL/aflow_anrl_A4B3_cI28_220_c_a.cpp ANRL/aflow_anrl_A2B3C6_cP33_221_cd_ag_fh.cpp ANRL/aflow_anrl_A5B3C16_cP96_222_ce_d_fi.cpp ANRL/aflow_anrl_A23B6_cF116_225_bd2f_e.cpp ANRL/aflow_anrl_A6B2C_cF36_225_e_c_a.cpp ANRL/aflow_anrl_AB13_cF112_226_a_bi.cpp ANRL/aflow_anrl_A2B2C7_cF88_227_c_d_af.cpp ANRL/aflow_anrl_A3B4_cF56_227_ad_e.cpp ANRL/aflow_anrl_A5BCD6_cF416_228_eg_c_b_h.cpp ANRL/aflow_anrl_A6B_cF224_228_h_c.cpp ANRL/aflow_anrl_A3B10_cI52_229_e_fh.cpp ANRL/aflow_anrl_A4B_cI10_229_c_a.cpp ANRL/aflow_anrl_A7B3_cI40_229_df_e.cpp ANRL/aflow_anrl_A2B3C12D3_cI160_230_a_c_h_d.cpp 

COMPARE_FUNCS_DH_SRC = aflow_compare_structure.cpp aflow_compare_structure_function.cpp
CONTRIB_FUNCS_SW_SRC = aflow_contrib_shidong_main.cpp aflow_contrib_shidong.cpp aflow_contrib_shidong_auxiliary.cpp aflow_contrib_shidong_cluster_expansion.cpp     
#[OBSOLETE] [JUNKAI] CONTRIB_FUNCS_JX_SRC = aflow_contrib_junkai_cproto.cpp aflow_contrib_junkai_protocheck.cpp aflow_contrib_junkai_primitivetrans.cpp aflow_contrib_junkai_basic.cpp aflow_contrib_junkai_magnet.cpp aflow_contrib_junkai_phasediag.cpp
#[OBSOLETE] [JUNKAI] CONTRIB_FUNCS_JX_SRC = aflow_contrib_junkai_basic.cpp aflow_contrib_junkai_phasediag.cpp
#[OBSOLETE] [KESONG] CONTRIB_FUNCS_KY_SRC = aflow_contrib_kesong_main.cpp aflow_contrib_kesong_ipocc.cpp aflow_contrib_kesong_hnfcell.cpp aflow_contrib_kesong_pocc_basic.cpp aflow_contrib_kesong_multienum.cpp
#DX and CO - START
CONTRIB_FUNCS_RT_SRC = aflow_symmetry_spacegroup.cpp aflow_symmetry_spacegroup_functions.cpp aflow_symmetry_spacegroup_ITC_library.cpp
#DX and CO - END
CONTRIB_FUNCS_WS_SRC = aflow_contrib_wahyu.cpp
CONTRIB_FUNCS_SC_SRC = 
#aflow_contrib_stefano_htsecurity.cpp

#DX and CO - START
CONTRIB_FUNCS_KR_SRC = aflow_bader.cpp
#DX and CO - END
CONTRIB_FUNCS_CT_SRC = aflow_contrib_cormac.cpp

CHULL_FUNCS_SRC = aflow_chull.cpp
#README_SRC= aflow_readme.cpp

DATA_SRC=aflow_data_readme.cpp aflow_data_calculated.cpp aflow_data_htqc.cpp aflow_data_libraries.cpp aflow_data_stokes.cpp aflow_data_nist.cpp aflow_data_latex.cpp aflow_data_aflow_potcars.cpp aflow_data.cpp
DATA_OBJ                = $(DATA_FUNCS_SRC:.cpp=.o)

#CO 180409
POCC_FUNCS_SRC = aflow_pocc.cpp aflow_pocc.h

#AUROSTD_OBJ    = $(AUROSTD_CPP:.cpp=.o)
AUROSTD_OBJ = aurostd.o
AFLOW_CLASS_OBJ        = $(AFLOW_CLASS_SRC:.cpp=.o)
AFLOW_FUNCS_OBJ        = $(AFLOW_FUNCS_SRC:.cpp=.o)
PFLOW_FUNCS_OBJ        = $(PFLOW_FUNCS_SRC:.cpp=.o)
ESTRUCTURE_FUNCS_OBJ   = $(ESTRUCTURE_FUNCS_SRC:.cpp=.o)
POCCUPATION_FUNCS_OBJ  = $(POCCUPATION_FUNCS_SRC:.cpp=.o)
APENNSY_FUNCS_OBJ      = $(APENNSY_FUNCS_SRC:.cpp=.o)
LIBRARIES_FUNCS_OBJ    = $(LIBRARIES_FUNCS_SRC:.cpp=.o)
CONTRIB_FUNCS_OBJ      = $(CONTRIB_FUNCS_SW_SRC:.cpp=.o) $(CONTRIB_FUNCS_JX_SRC:.cpp=.o) $(CONTRIB_FUNCS_KY_SRC:.cpp=.o) $(CONTRIB_FUNCS_RT_SRC:.cpp=.o) $(CONTRIB_FUNCS_WS_SRC:.cpp=.o) $(CONTRIB_FUNCS_SC_SRC:.cpp=.o) $(CONTRIB_FUNCS_KR_SRC:.cpp=.o) $(CONTRIB_FUNCS_CT_SRC:.cpp=.o) 
#README_OBJ = $(README_SRC:.cpp=.o)
CALCULATED_OBJ         = $(CALCULATED_SRC:.cpp=.o)
AFLOWLIB_OBJ           = $(AFLOWLIB_SRC:.cpp=.o)
AGL_OBJ                = $(AGL_SRC:.cpp=.o)
AEL_OBJ                = $(AEL_SRC:.cpp=.o)
APL_OBJ                = $(APL_SRC:.cpp=.o)
ANRL_OBJ               = $(ANRL_SRC:.cpp=.o)
COMPARE_FUNCS_DH_OBJ   = ${COMPARE_FUNCS_DH_SRC:.cpp=.o}
CHULL_OBJ              = $(CHULL_FUNCS_SRC:.cpp=.o)
POCC_OBJ               = $(POCC_FUNCS_SRC:.cpp=.o) #CO 180409
SQLITE_OBJ             = $(SQLITE_C:.c=.o)	#CO 170614

OBJS_ALL= aflow.o $(AUROSTD_OBJ) $(AFLOW_CLASS_OBJ) $(AFLOW_FUNCS_OBJ) $(PFLOW_FUNCS_OBJ) $(ESTRUCTURE_FUNCS_OBJ) $(POCCUPATION_FUNCS_OBJ) $(APENNSY_FUNCS_OBJ) $(AUROSTD) $(LIBRARIES_FUNCS_OBJ) aflow_mix_experiments.o $(CONTRIB_FUNCS_OBJ) aflowlib_webapp_entry.cpp aflowlib_webapp_bands.cpp $(AFLOWLIB_OBJ) $(AGL_OBJ) $(AEL_OBJ) $(APL_OBJ) ${COMPARE_FUNCS_DH_OBJ} $(CHULL_OBJ) $(POCC_OBJ) $(SQLITE_OBJ)
#$(ANRL_OBJ) ${COMPARE_FUNCS_DH_OBJ} $(CHULL_OBJ) $(POCC_OBJ) $(SQLITE_OBJ)

all: aflow aflow_data 
	$(MAKE) check
aflow: $(OBJS_ALL) Makefile 
#	$(MAKE) -C ../AFLOW2_FREE/EXTRA/APL2AGR/ && cd ../../
	cp -f aflow_apl*.o  aflow_aapl*.o  aflow_qha*.o ./APL/
	#cp -f aflow_anrl*.o ./ANRL/   #KSYANG, 20190501
	#cp sqlite3.o ./SQLITE/
#	$(CPP) $(LIBS_AFLOW) -o aflow *.o $(LIB_MKL) $(LIBS_AFLOW)
	$(CPP) $(LIBS_AFLOW) -o aflow `ls *.o | grep -v aflow_data` $(LIB_MKL) $(LIBS_AFLOW)
	ln -sf aflow aflowd
	ln -sf aflow apennsy
static_no_pdf: aflow aflow_data
	rm -rf $(SQLITE_OBJ)
	$(CC) -DSQLITE_ENABLE_JSON1 -c $(SQLITE_C) $(LIBS) $(LIBS_STATIC)
	cp sqlite3.o ./SQLITE/
	rm -f aflow
#	$(CPP) $(LIBS_AFLOW_STATIC) -o aflow *.o $(LIB_MKL) $(LIBS_AFLOW_STATIC)
	$(CPP) $(LIBS_AFLOW_STATIC) -o aflow `ls *.o | grep -v aflow_data` $(LIB_MKL) $(LIBS_AFLOW_STATIC)
	rm -f aflow_data
#	$(CPP) $(VERS) $(INCLUDE) $(CCFLAGS) $(OPTS) -o aflow_data aflow_data.cpp $(LIBS_STATIC)
	$(CPP) $(VERS) $(INCLUDE) $(CCFLAGS) $(OPTS) -o aflow_data aflow_data.o aflow_data_*.o $(LIBS_STATIC)
#static: static_no_pdf
#	tar Jcfvp /www/AFLOW/CURRENT/aflow_$(VERSNUMBER)_amd64_debian_static.tar.xz aflow aflow_data	
$(SQLITE_OBJ): $(SQLITE_C)	#CO 1706147
	#$(CC) -DSQLITE_ENABLE_JSON1 -DSQLITE_OMIT_LOAD_EXTENSION -o $@ -c $^ $(LIBS)	#-DSQLITE_OMIT_LOAD_EXTENSION doesn't really work amalgamation
	#$(CC) -DSQLITE_ENABLE_JSON1 -o $@ -c $^ $(LIBS)
	$(CC) -DSQLITE_ENABLE_JSON1 -c $^ $(LIBS)
	cp sqlite3.o ./SQLITE/
extra:	findsym frozsl apl2agr platon subst
apl2agr:
	echo "MAKING APL2AGR"
	if [ -f $(EXTRA_FREE)/APL2AGR/Makefile ]; then $(MAKE) -C $(EXTRA_FREE)/APL2AGR realclean; fi
	if [ -f $(EXTRA_FREE)/APL2AGR/Makefile ]; then $(MAKE) -C $(EXTRA_FREE)/APL2AGR; fi
	if [ -f $(EXTRA_FREE)/APL2AGR/apl2agr ] && [ -d ../AFLOW3_NONFREE/BIN ]; then cp $(EXTRA_FREE)/APL2AGR/apl2agr ../AFLOW3_NONFREE/BIN; fi
subst:
	echo "MAKING SUBST"
	if [ -f $(EXTRA_FREE)/SUBST/Makefile ]; then $(MAKE) -C $(EXTRA_FREE)/SUBST clean; fi
	if [ -f $(EXTRA_FREE)/SUBST/Makefile ]; then $(MAKE) -C $(EXTRA_FREE)/SUBST; fi
	if [ -f $(EXTRA_FREE)/SUBST/subst ] && [ -d ../AFLOW3_NONFREE/BIN ]; then cp $(EXTRA_FREE)/SUBST/subst ../AFLOW3_NONFREE/BIN; fi
findsym:
	echo "MAKING FINDSYM"
	if [ -f ../AFLOW3_NONFREE/EXTRA/FINDSYM/Makefile ]; then $(MAKE) -C ../AFLOW3_NONFREE/EXTRA/FINDSYM clean; fi
	if [ -f ../AFLOW3_NONFREE/EXTRA/FINDSYM/Makefile ]; then $(MAKE) -C ../AFLOW3_NONFREE/EXTRA/FINDSYM; fi
	if [ -f ../AFLOW3_NONFREE/EXTRA/FINDSYM/findsym ]; then cp ../AFLOW3_NONFREE/EXTRA/FINDSYM/findsym ../AFLOW3_NONFREE/BIN; fi
frozsl:
	echo "MAKING FROZSL"
	if [ -f ../AFLOW3_NONFREE/EXTRA/FROZSL/Makefile ]; then $(MAKE) -C ../AFLOW3_NONFREE/EXTRA/FROZSL clean; fi
	if [ -f ../AFLOW3_NONFREE/EXTRA/FROZSL/Makefile ]; then $(MAKE) -C ../AFLOW3_NONFREE/EXTRA/FROZSL; fi
	if [ -f ../AFLOW3_NONFREE/EXTRA/FROZSL/frozsl ]; then cp ../AFLOW3_NONFREE/EXTRA/FROZSL/frozsl ../AFLOW3_NONFREE/BIN; fi
platon:
	echo "MAKING PLATON"
	if [ -f ../AFLOW3_NONFREE/EXTRA/PLATON/51108/Makefile ]; then $(MAKE) -C ../AFLOW3_NONFREE/EXTRA/PLATON/51108 clean; fi
	if [ -f ../AFLOW3_NONFREE/EXTRA/PLATON/51108/Makefile ]; then $(MAKE) -C ../AFLOW3_NONFREE/EXTRA/PLATON/51108; fi
	if [ -f ../AFLOW3_NONFREE/EXTRA/PLATON/51108/platon ]; then cp ../AFLOW3_NONFREE/EXTRA/PLATON/51108/platon ../AFLOW3_NONFREE/BIN; fi
enum:
	echo "MAKING ENUM"
	if [ -f $(EXTRA_FREE)/GUS/Makefile ]; then $(MAKE) -C $(EXTRA_FREE)/GUS; fi
extra_clean:
	if [ -f ../AFLOW3_NONFREE/EXTRA/FINDSYM/Makefile ]; then $(MAKE) -C ../AFLOW3_NONFREE/EXTRA/FINDSYM clean; fi
	if [ -f ../AFLOW3_NONFREE/EXTRA/FROZSL/Makefile ]; then $(MAKE) -C ../AFLOW3_NONFREE/EXTRA/FROZSL clean; fi
	if [ -f ../AFLOW3_NONFREE/EXTRA/PLATON/51108/Makefile ]; then $(MAKE) -C ../AFLOW3_NONFREE/EXTRA/PLATON/51108 clean; fi
	if [ -f $(EXTRA_FREE)/GUS/Makefile ]; then $(MAKE) -C $(EXTRA_FREE)/GUS realclean; fi
	if [ -f $(EXTRA_FREE)/APL2AGR/Makefile ]; then $(MAKE) -C $(EXTRA_FREE)/APL2AGR realclean; fi
mac:
	$(MAKE) aflow_data
	$(MAKE) -e ARCH=-D_MACOSX_ -j 4
	ln -sf aflow aflowd
	ln -sf aflow apennsy
	$(MAKE) check
#mac_static: mac
#	tar Jcfvp aflow_$(VERSNUMBER)_maxOSX_10.X.X.tar.xz aflow aflow_data
aflow.o: aflow.cpp aflow_test.cpp $(AFLOW_HPP) $(AUROSTD_HPP) Makefile
	$(CPP) $(VERS) $(INCLUDE) $(CCFLAGS) $(OPTS) aflow.cpp -c
#aurostd.o: aurostd.cpp aurostd.h Makefile
#	$(CPP) $(VERS) $(INCLUDE) $(CCFLAGS) $(OPTS) aurostd.cpp -c -o aurostd.o
#aurostd.cpp:
#	cd ./AUROSTD && ppp4cpp.py
#aurostd.h:
#	cd ./AUROSTD && ppp4cpp.py

aurostd.o: $(AUROSTD_HPP) $(AUROSTD_CPP) Makefile
	$(CPP) $(VERS) $(INCLUDE) $(CCFLAGS) $(OPTS) $(ARCH) ./AUROSTD/aurostd.cpp -c -o aurostd.o 
.cpp.o: .cpp Makefile $(AFLOW_HPP) $(AUROSTD_HPP) $(AUROSTD_CPP)
	$(CPP) $(VERS) $(INCLUDE) $(CCFLAGS) $(OPTS) $(@:.o=.cpp) -c
	touch $@
aflow_init.o: aflow_init.cpp Makefile $(AFLOW_HPP) $(AUROSTD_HPP) $(AUROSTD_CPP)
	$(CPP) $(VERS) $(INCLUDE) $(CCFLAGS) $(OPTS) $(ARCH) aflow_init.cpp -c
aflow_apennsy_phases.o: aflow_apennsy_phases.cpp $(AFLOW_HPP) $(AUROSTD_HPP) Makefile
	$(CPP) $(VERS) $(INCLUDE) $(CCFLAGS)  aflow_apennsy_phases.cpp -c
#aflow_mix_nomix_table.o: aflow_mix_nomix_table.cpp $(AFLOW_HPP) $(AUROSTD_HPP) Makefile
#	$(CPP) $(VERS) $(INCLUDE) $(CCFLAGS) $(OPTS) aflow_mix_nomix_table.cpp -c
#CO 170622 - START
aflowlib_webapp_entry.cpp: aflowlib_webapp_entry.js
	rm -rf $@
	echo "// $@ automatic generated" >> $@
	echo "std::string AFLOW_WEBAPP_ENTRY_JS=\"\\" >> $@
	cat $^ | perl -p -e 's/\r//g' | perl -p -e 's/\\u/\\\\u/g' | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> $@
	echo "\";" >> $@
#CO 170622 - END
#CO 180305 - START
aflowlib_webapp_bands.cpp: aflowlib_webapp_bands.js
	rm -rf $@
	echo "// $@ automatic generated" >> $@
	echo "std::string AFLOW_WEBAPP_BANDS_JS=\"\\" >> $@
	cat $^ | perl -p -e 's/\r//g' | perl -p -e 's/\\u/\\\\u/g' | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> $@
	echo "\";" >> $@
#CO 180305 - END
aflowlib_web_interface.o: aflowlib_web_interface.cpp aflowlib_webapp_entry.cpp aflowlib_webapp_bands.cpp
	$(CPP) $(VERS) $(INCLUDE) $(CCFLAGS) $(OPTS) $(ARCH) aflowlib_web_interface.cpp -c
aflow_mix_nomix_table.cpp: Makefile $(AFLOW_HPP) $(AUROSTD_HPP)
	rm -rf nomix.tmp
	find $(PROJECT_DIRECTORY_GNDSTATE)/NOMIX/ -name aflow.in | sed "s/\/common\/$(PROJECT_DIRECTORY_GNDSTATE)\/NOMIX\/LIB2\/LIB\///g" | sed "s/aflow.in//g" >> nomix.tmp
	rm -rf aflow_mix_nomix_table.cpp
	echo "// aflow_mix_nomix_table.cpp automatic generated " >> aflow_mix_nomix_table.cpp
	echo "#include <string>" >> aflow_mix_nomix_table.cpp
	echo "std::string NomixTable=\"\\" >> aflow_mix_nomix_table.cpp
	cat nomix.tmp | sort | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_mix_nomix_table.cpp
	echo "\";" >> aflow_mix_nomix_table.cpp
	chmod 755 aflow_mix_nomix_table.cpp
	rm -rf nomix.tmp
aflow_mix_experiments.o: aflow_mix_experiments.cpp $(AFLOW_HPP) $(AUROSTD_HPP) Makefile
	$(CPP) $(VERS) $(INCLUDE) $(CCFLAGS) $(OPTS) aflow_mix_experiments.cpp -c
aflow_xproto.o: aflow_xproto.cpp aflow_xproto_library_default.cpp Makefile
	$(CPP) $(VERS) $(INCLUDE) $(CCFLAGS) $(OPTS) aflow_xproto.cpp -c
aflow_library_icsd.dat.xz: $(ICSD_LIB)/README_LIBRARY_ICSD1.TXT $(ICSD_LIB)/README_LIBRARY_ICSD2.TXT $(ICSD_LIB)/README_LIBRARY_ICSD3.TXT $(ICSD_LIB)/README_LIBRARY_ICSD4.TXT $(ICSD_LIB)/README_LIBRARY_ICSD5.TXT $(ICSD_LIB)/README_LIBRARY_ICSD6.TXT $(ICSD_LIB)/README_LIBRARY_ICSD7.TXT $(ICSD_LIB)/README_LIBRARY_ICSD8.TXT $(ICSD_LIB)/README_LIBRARY_ICSD9.TXT
	$(MAKE) aflow -j 16
	rm -rf aflow_library_icsd.dat
	echo "[AFLOW] aflow_library_icsd.dat  " >> aflow_library_icsd.dat
	echo "[README_LIBRARY_ICSD1.TXT]START " >> aflow_library_icsd.dat
	cat $(ICSD_LIB)/README_LIBRARY_ICSD1.TXT >> aflow_library_icsd.dat
	echo "[README_LIBRARY_ICSD1.TXT]STOP " >> aflow_library_icsd.dat
	echo "[README_LIBRARY_ICSD2.TXT]START " >> aflow_library_icsd.dat
	cat $(ICSD_LIB)/README_LIBRARY_ICSD2.TXT >> aflow_library_icsd.dat
	echo "[README_LIBRARY_ICSD2.TXT]STOP " >> aflow_library_icsd.dat
	echo "[README_LIBRARY_ICSD3.TXT]START " >> aflow_library_icsd.dat
	cat $(ICSD_LIB)/README_LIBRARY_ICSD3.TXT >> aflow_library_icsd.dat
	echo "[README_LIBRARY_ICSD3.TXT]STOP " >> aflow_library_icsd.dat
	echo "[README_LIBRARY_ICSD4.TXT]START " >> aflow_library_icsd.dat
	cat $(ICSD_LIB)/README_LIBRARY_ICSD4.TXT >> aflow_library_icsd.dat
	echo "[README_LIBRARY_ICSD4.TXT]STOP " >> aflow_library_icsd.dat
	echo "[README_LIBRARY_ICSD5.TXT]START " >> aflow_library_icsd.dat
	cat $(ICSD_LIB)/README_LIBRARY_ICSD5.TXT >> aflow_library_icsd.dat
	echo "[README_LIBRARY_ICSD5.TXT]STOP " >> aflow_library_icsd.dat
	echo "[README_LIBRARY_ICSD6.TXT]START " >> aflow_library_icsd.dat
	cat $(ICSD_LIB)/README_LIBRARY_ICSD6.TXT >> aflow_library_icsd.dat
	echo "[README_LIBRARY_ICSD6.TXT]STOP " >> aflow_library_icsd.dat
	echo "[README_LIBRARY_ICSD7.TXT]START " >> aflow_library_icsd.dat
	cat $(ICSD_LIB)/README_LIBRARY_ICSD7.TXT >> aflow_library_icsd.dat
	echo "[README_LIBRARY_ICSD7.TXT]STOP " >> aflow_library_icsd.dat
	echo "[README_LIBRARY_ICSD8.TXT]START " >> aflow_library_icsd.dat
	cat $(ICSD_LIB)/README_LIBRARY_ICSD8.TXT >> aflow_library_icsd.dat
	echo "[README_LIBRARY_ICSD8.TXT]STOP " >> aflow_library_icsd.dat
	echo "[README_LIBRARY_ICSD9.TXT]START " >> aflow_library_icsd.dat
	cat $(ICSD_LIB)/README_LIBRARY_ICSD9.TXT >> aflow_library_icsd.dat
	echo "[README_LIBRARY_ICSD9.TXT]STOP " >> aflow_library_icsd.dat
	cat aflow_library_icsd.dat |sed "s/WYCK_ICSD/WYC_ICSD/g" |sed "s/  / /g"|sed "s/  / /g"|sed "s/  / /g"|sed "s/  / /g"|sed "s/  / /g" | grep -v "//" > aflow_library_icsd.dat2
	mv aflow_library_icsd.dat2 aflow_library_icsd.dat
	xz -9vf aflow_library_icsd.dat
	cp -f aflow_library_icsd.dat.xz ../AFLOW3_NONFREE/LIBS
	cp -f aflow_library_icsd.dat.xz /common/AFLOW/LIBS
	rm -f aflow_data_*
	$(MAKE) aflow_data
#	rm -f $(EXTRA_FREE)/NIST/ICSD_List_2014.txt*
	./aflow --protos_icsd > $(EXTRA_FREE)/NIST/ICSD_List_2014.txt
	rm -f $(EXTRA_FREE)/NIST/ICSD_List_2014.txt.xz
	xz -9vf $(EXTRA_FREE)/NIST/ICSD_List_2014.txt
$(ICSD_LIB)/README_LIBRARY_ICSD1.TXT:
	rm -rf $(ICSD_LIB)/README_LIBRARY_ICSD1.TXT
	echo "// ***************************************************************************">>$(ICSD_LIB)/README_LIBRARY_ICSD1.TXT
	xzcat $(NIST_LIB)/README_LIBRARY_ICSD1.TXT.xz >>$(ICSD_LIB)/README_LIBRARY_ICSD1.TXT
	echo "// ***************************************************************************">>$(ICSD_LIB)/README_LIBRARY_ICSD1.TXT
	echo " STRUCTURE ___PROTO_END___">>$(ICSD_LIB)/README_LIBRARY_ICSD1.TXT
	echo " END">>$(ICSD_LIB)/README_LIBRARY_ICSD1.TXT
	echo "// ***************************************************************************">>$(ICSD_LIB)/README_LIBRARY_ICSD1.TXT
$(ICSD_LIB)/README_LIBRARY_ICSD2.TXT:
	rm -rf $(ICSD_LIB)/README_LIBRARY_ICSD2.TXT
	echo "// ***************************************************************************">>$(ICSD_LIB)/README_LIBRARY_ICSD2.TXT
	xzcat $(NIST_LIB)/README_LIBRARY_ICSD2.TXT.xz >>$(ICSD_LIB)/README_LIBRARY_ICSD2.TXT
	echo "// ***************************************************************************">>$(ICSD_LIB)/README_LIBRARY_ICSD2.TXT
	echo " STRUCTURE ___PROTO_END___">>$(ICSD_LIB)/README_LIBRARY_ICSD2.TXT
	echo " END">>$(ICSD_LIB)/README_LIBRARY_ICSD2.TXT
	echo "// ***************************************************************************">>$(ICSD_LIB)/README_LIBRARY_ICSD2.TXT
$(ICSD_LIB)/README_LIBRARY_ICSD3.TXT:
	rm -rf $(ICSD_LIB)/README_LIBRARY_ICSD3.TXT
	echo "// ***************************************************************************">>$(ICSD_LIB)/README_LIBRARY_ICSD3.TXT
	xzcat $(NIST_LIB)/README_LIBRARY_ICSD3.TXT.xz >>$(ICSD_LIB)/README_LIBRARY_ICSD3.TXT
	echo "// ***************************************************************************">>$(ICSD_LIB)/README_LIBRARY_ICSD3.TXT
	echo " STRUCTURE ___PROTO_END___">>$(ICSD_LIB)/README_LIBRARY_ICSD3.TXT
	echo " END">>$(ICSD_LIB)/README_LIBRARY_ICSD3.TXT
	echo "// ***************************************************************************">>$(ICSD_LIB)/README_LIBRARY_ICSD3.TXT
$(ICSD_LIB)/README_LIBRARY_ICSD4.TXT:
	rm -rf $(ICSD_LIB)/README_LIBRARY_ICSD4.TXT
	echo "// ***************************************************************************">>$(ICSD_LIB)/README_LIBRARY_ICSD4.TXT
	xzcat $(NIST_LIB)/README_LIBRARY_ICSD4.TXT.xz >>$(ICSD_LIB)/README_LIBRARY_ICSD4.TXT
	echo "// ***************************************************************************">>$(ICSD_LIB)/README_LIBRARY_ICSD4.TXT
	echo " STRUCTURE ___PROTO_END___">>$(ICSD_LIB)/README_LIBRARY_ICSD4.TXT
	echo " END">>$(ICSD_LIB)/README_LIBRARY_ICSD4.TXT
	echo "// ***************************************************************************">>$(ICSD_LIB)/README_LIBRARY_ICSD4.TXT
$(ICSD_LIB)/README_LIBRARY_ICSD5.TXT:
	rm -rf $(ICSD_LIB)/README_LIBRARY_ICSD5.TXT
	echo "// ***************************************************************************">>$(ICSD_LIB)/README_LIBRARY_ICSD5.TXT
	xzcat $(NIST_LIB)/README_LIBRARY_ICSD5.TXT.xz >>$(ICSD_LIB)/README_LIBRARY_ICSD5.TXT
	echo "// ***************************************************************************">>$(ICSD_LIB)/README_LIBRARY_ICSD5.TXT
	echo " STRUCTURE ___PROTO_END___">>$(ICSD_LIB)/README_LIBRARY_ICSD5.TXT
	echo " END">>$(ICSD_LIB)/README_LIBRARY_ICSD5.TXT
	echo "// ***************************************************************************">>$(ICSD_LIB)/README_LIBRARY_ICSD5.TXT
$(ICSD_LIB)/README_LIBRARY_ICSD6.TXT:
	rm -rf $(ICSD_LIB)/README_LIBRARY_ICSD6.TXT
	echo "// ***************************************************************************">>$(ICSD_LIB)/README_LIBRARY_ICSD6.TXT
	xzcat $(NIST_LIB)/README_LIBRARY_ICSD6.TXT.xz >>$(ICSD_LIB)/README_LIBRARY_ICSD6.TXT
	echo "// ***************************************************************************">>$(ICSD_LIB)/README_LIBRARY_ICSD6.TXT
	echo " STRUCTURE ___PROTO_END___">>$(ICSD_LIB)/README_LIBRARY_ICSD6.TXT
	echo " END">>$(ICSD_LIB)/README_LIBRARY_ICSD6.TXT
	echo "// ***************************************************************************">>$(ICSD_LIB)/README_LIBRARY_ICSD6.TXT
$(ICSD_LIB)/README_LIBRARY_ICSD7.TXT:
	rm -rf $(ICSD_LIB)/README_LIBRARY_ICSD7.TXT
	echo "// ***************************************************************************">>$(ICSD_LIB)/README_LIBRARY_ICSD7.TXT
	xzcat $(NIST_LIB)/README_LIBRARY_ICSD7.TXT.xz >>$(ICSD_LIB)/README_LIBRARY_ICSD7.TXT
	echo "// ***************************************************************************">>$(ICSD_LIB)/README_LIBRARY_ICSD7.TXT
	echo " STRUCTURE ___PROTO_END___">>$(ICSD_LIB)/README_LIBRARY_ICSD7.TXT
	echo " END">>$(ICSD_LIB)/README_LIBRARY_ICSD7.TXT
	echo "// ***************************************************************************">>$(ICSD_LIB)/README_LIBRARY_ICSD7.TXT
$(ICSD_LIB)/README_LIBRARY_ICSD8.TXT:
	rm -rf $(ICSD_LIB)/README_LIBRARY_ICSD8.TXT
	echo "// ***************************************************************************">>$(ICSD_LIB)/README_LIBRARY_ICSD8.TXT
	xzcat $(NIST_LIB)/README_LIBRARY_ICSD8.TXT.xz >>$(ICSD_LIB)/README_LIBRARY_ICSD8.TXT
	echo "// ***************************************************************************">>$(ICSD_LIB)/README_LIBRARY_ICSD8.TXT
	echo " STRUCTURE ___PROTO_END___">>$(ICSD_LIB)/README_LIBRARY_ICSD8.TXT
	echo " END">>$(ICSD_LIB)/README_LIBRARY_ICSD8.TXT
	echo "// ***************************************************************************">>$(ICSD_LIB)/README_LIBRARY_ICSD8.TXT
$(ICSD_LIB)/README_LIBRARY_ICSD9.TXT:
	rm -rf $(ICSD_LIB)/README_LIBRARY_ICSD9.TXT
	echo "// ***************************************************************************">>$(ICSD_LIB)/README_LIBRARY_ICSD9.TXT
	xzcat $(NIST_LIB)/README_LIBRARY_ICSD9.TXT.xz >>$(ICSD_LIB)/README_LIBRARY_ICSD9.TXT
	echo "// ***************************************************************************">>$(ICSD_LIB)/README_LIBRARY_ICSD9.TXT
	echo " STRUCTURE ___PROTO_END___">>$(ICSD_LIB)/README_LIBRARY_ICSD9.TXT
	echo " END">>$(ICSD_LIB)/README_LIBRARY_ICSD9.TXT
	echo "// ***************************************************************************">>$(ICSD_LIB)/README_LIBRARY_ICSD9.TXT
############################################################################################################
aflow_data_calculated.cpp: Makefile
	rm -rf aflow_data_calculated.cpp
	if [ -f $(EXTRA_FREE)/CALCULATED/vAUID.dat.xz ]; then echo "// aflow_data_calculated.cpp automatic generated from CALCULATED/vAUID.dat" >> aflow_data_calculated.cpp ; fi
	if [ -f $(EXTRA_FREE)/CALCULATED/vAUID.dat.xz ]; then echo "#include <string>" >> aflow_data_calculated.cpp ; fi
	if [ -f $(EXTRA_FREE)/CALCULATED/vAUID.dat.xz ]; then echo "std::string vAUID=\"\\" >> aflow_data_calculated.cpp ; fi
	if [ -f $(EXTRA_FREE)/CALCULATED/vAUID.dat.xz ]; then xzcat $(EXTRA_FREE)/CALCULATED/vAUID.dat.xz | grep -v "//" | perl -p -e 's/\n/\\n\\\n/g' >> aflow_data_calculated.cpp ; fi
	if [ -f $(EXTRA_FREE)/CALCULATED/vAUID.dat.xz ]; then echo "\";" >> aflow_data_calculated.cpp ; fi
	if [ -f $(EXTRA_FREE)/CALCULATED/vAURL.dat.xz ]; then echo "// aflow_data_calculated.cpp automatic generated from CALCULATED/vAURL.dat" >> aflow_data_calculated.cpp ; fi
	if [ -f $(EXTRA_FREE)/CALCULATED/vAURL.dat.xz ]; then echo "#include <string>" >> aflow_data_calculated.cpp ; fi
	if [ -f $(EXTRA_FREE)/CALCULATED/vAURL.dat.xz ]; then echo "std::string vAURL=\"\\" >> aflow_data_calculated.cpp ; fi
	if [ -f $(EXTRA_FREE)/CALCULATED/vAURL.dat.xz ]; then xzcat $(EXTRA_FREE)/CALCULATED/vAURL.dat.xz | grep -v "//" | perl -p -e 's/\n/\\n\\\n/g' >> aflow_data_calculated.cpp ; fi
	if [ -f $(EXTRA_FREE)/CALCULATED/vAURL.dat.xz ]; then echo "\";" >> aflow_data_calculated.cpp ; fi
	if [ -f $(EXTRA_FREE)/CALCULATED/vLOOP.dat.xz ]; then echo "// aflow_data_calculated.cpp automatic generated from CALCULATED/vLOOP.dat" >> aflow_data_calculated.cpp ; fi
	if [ -f $(EXTRA_FREE)/CALCULATED/vLOOP.dat.xz ]; then echo "#include <string>" >> aflow_data_calculated.cpp ; fi
	if [ -f $(EXTRA_FREE)/CALCULATED/vLOOP.dat.xz ]; then echo "std::string vLOOP=\"\\" >> aflow_data_calculated.cpp ; fi
	if [ -f $(EXTRA_FREE)/CALCULATED/vLOOP.dat.xz ]; then xzcat $(EXTRA_FREE)/CALCULATED/vLOOP.dat.xz | grep -v "//" | perl -p -e 's/\n/\\n\\\n/g' >> aflow_data_calculated.cpp ; fi
	if [ -f $(EXTRA_FREE)/CALCULATED/vLOOP.dat.xz ]; then echo "\";" >> aflow_data_calculated.cpp ; fi
	if [ ! -f aflow_data_calculated.cpp ]; then rm -rf aflow_data_calculated.cpp.xz*; wget http://materials.duke.edu/AFLOW/AFLOW3_FREE/aflow_data_calculated.cpp.xz && xz -dv aflow_data_calculated.cpp.xz ; fi
	touch $@
aflow_data_htqc.cpp: README_LIBRARY_HTQC.TXT Makefile
	rm -f aflow_data_htqc.cpp
	echo "// aflow_data_htqc.cpp automatic generated from README_LIBRARY_HTQC.TXT" >> aflow_data_htqc.cpp
	echo "#include <string>" >> aflow_data_htqc.cpp
	echo "std::string Library_HTQC=\"\\" >> aflow_data_htqc.cpp
	cat README_LIBRARY_HTQC.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_htqc.cpp
	cat README_LIBRARY_HTQC_BORIDES.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_htqc.cpp
	cat README_LIBRARY_HTQC_CARBIDES.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_htqc.cpp
	echo "STRUCTURE ___PROTO_END___"  | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_htqc.cpp
	echo " END"  | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_htqc.cpp
	echo "\";" >> aflow_data_htqc.cpp
#aflow_data_binary.cpp:  Makefile
#	rm -rf aflow_data_binary.cpp
#	echo "// aflow_data_binary.cpp automatic generated from PROTO/Binary.exp" >> aflow_data_binary.cpp
#	echo "std::string AFLOW_BinaryRead=\"\\" >> aflow_data_binary.cpp
#	cat $(EXTRA_FREE)/PROTO/Binary.exp | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> aflow_data_binary.cpp
#	echo "\";" >> aflow_data_binary.cpp
#	echo "std::string AFLOW_Binary_Angle_Read=\"\\" >> aflow_data_binary.cpp
#	cat $(EXTRA_FREE)/PROTO/Binary_angle.exp | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> aflow_data_binary.cpp
#	echo "\";" >> aflow_data_binary.cpp
aflow_data_libraries.cpp: Makefile
	rm -rf aflow_data_libraries.cpp
	echo "// aflow_data_libraries.cpp automatic generated from ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_libX.dat" >> aflow_data_libraries.cpp
	echo "#include <string>" >> aflow_data_libraries.cpp
#	echo "std::string aflowlib_icsd=\"\\" >> aflow_data_libraries.cpp && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_icsd.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\\Gamma/\\\\Gamma/g' | perl -p -e 's/\\Sigma/\\\\Sigma/g' >> aflow_data_libraries.cpp && echo "\";" >> aflow_data_libraries.cpp
	echo "std::string aflowlib_icsd=\"\";" >> aflow_data_libraries.cpp
	echo "LIB1"
#	echo "std::string aflowlib_lib1=\"\\" >> aflow_data_libraries.cpp && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib1.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\\Gamma/\\\\Gamma/g' | perl -p -e 's/\\Sigma/\\\\Sigma/g' >> aflow_data_libraries.cpp && echo "\";" >> aflow_data_libraries.cpp	
	echo "std::string aflowlib_lib1=\"\";" >> aflow_data_libraries.cpp
	echo "LIB2"
#	echo "std::string aflowlib_lib2=\"\\" >> aflow_data_libraries.cpp && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib2.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\\Gamma/\\\\Gamma/g' | perl -p -e 's/\\Sigma/\\\\Sigma/g' >> aflow_data_libraries.cpp && echo "\";" >> aflow_data_libraries.cpp	
	echo "std::string aflowlib_lib2=\"\";" >> aflow_data_libraries.cpp
	echo "LIB3"
#	echo "std::string aflowlib_lib3=\"\\" >> aflow_data_libraries.cpp && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib3.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\\Gamma/\\\\Gamma/g' | perl -p -e 's/\\Sigma/\\\\Sigma/g' >> aflow_data_libraries.cpp && echo "\";" >> aflow_data_libraries.cpp
	echo "std::string aflowlib_lib3=\"\";" >> aflow_data_libraries.cpp
	echo "LIB4"
#	echo "std::string aflowlib_lib4=\"\\" >> aflow_data_libraries.cpp && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib4.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\\Gamma/\\\\Gamma/g' | perl -p -e 's/\\Sigma/\\\\Sigma/g' >> aflow_data_libraries.cpp && echo "\";" >> aflow_data_libraries.cpp
	echo "std::string aflowlib_lib4=\"\";" >> aflow_data_libraries.cpp
	echo "LIB5"
#	echo "std::string aflowlib_lib5=\"\\" >> aflow_data_libraries.cpp && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib5.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\\Gamma/\\\\Gamma/g' | perl -p -e 's/\\Sigma/\\\\Sigma/g' >> aflow_data_libraries.cpp && echo "\";" >> aflow_data_libraries.cpp
	echo "std::string aflowlib_lib5=\"\";" >> aflow_data_libraries.cpp
	echo "LIB6"
#	echo "std::string aflowlib_lib6=\"\\" >> aflow_data_libraries.cpp && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib6.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\\Gamma/\\\\Gamma/g' | perl -p -e 's/\\Sigma/\\\\Sigma/g' >> aflow_data_libraries.cpp && echo "\";" >> aflow_data_libraries.cpp
	echo "std::string aflowlib_lib6=\"\";" >> aflow_data_libraries.cpp
	echo "LIB7"
#	echo "std::string aflowlib_lib7=\"\\" >> aflow_data_libraries.cpp && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib7.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\\Gamma/\\\\Gamma/g' | perl -p -e 's/\\Sigma/\\\\Sigma/g' >> aflow_data_libraries.cpp && echo "\";" >> aflow_data_libraries.cpp
	echo "std::string aflowlib_lib7=\"\";" >> aflow_data_libraries.cpp
	echo "LIB8"
#	echo "std::string aflowlib_lib8=\"\\" >> aflow_data_libraries.cpp && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib8.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\\Gamma/\\\\Gamma/g' | perl -p -e 's/\\Sigma/\\\\Sigma/g' >> aflow_data_libraries.cpp && echo "\";" >> aflow_data_libraries.cpp
	echo "std::string aflowlib_lib8=\"\";" >> aflow_data_libraries.cpp
	echo "LIB9"
#	echo "std::string aflowlib_lib9=\"\\" >> aflow_data_libraries.cpp && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib9.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\\Gamma/\\\\Gamma/g' | perl -p -e 's/\\Sigma/\\\\Sigma/g' >> aflow_data_libraries.cpp && echo "\";" >> aflow_data_libraries.cpp
	echo "std::string aflowlib_lib9=\"\";" >> aflow_data_libraries.cpp
aflow_data_stokes.cpp: Makefile 
	rm -rf aflow_data_stokes.cpp
	if [ -f $(EXTRA_FREE)/DATA_FINDSYM/data_space.txt.xz ]; then echo "// aflow_data_stokes.cpp automatic generated " >> aflow_data_stokes.cpp ; fi
	if [ -f $(EXTRA_FREE)/DATA_FINDSYM/data_space.txt.xz ]; then echo "#include <string>" >> aflow_data_stokes.cpp ; fi
	if [ -f $(EXTRA_FREE)/DATA_FINDSYM/data_space.txt.xz ]; then echo "std::string FINDSYM_data_space_txt=\"\\" >> aflow_data_stokes.cpp ; fi
	if [ -f $(EXTRA_FREE)/DATA_FINDSYM/data_space.txt.xz ]; then xzcat $(EXTRA_FREE)/DATA_FINDSYM/data_space.txt.xz | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\0//g' >> aflow_data_stokes.cpp ; fi
	if [ -f $(EXTRA_FREE)/DATA_FINDSYM/data_space.txt.xz ]; then echo "\";" >> aflow_data_stokes.cpp ; fi
	if [ -f $(EXTRA_FREE)/DATA_FINDSYM/data_wyckoff.txt.xz ]; then echo "std::string FINDSYM_data_wyckoff_txt=\"\\" >> aflow_data_stokes.cpp ; fi
	if [ -f $(EXTRA_FREE)/DATA_FINDSYM/data_wyckoff.txt.xz ]; then xzcat $(EXTRA_FREE)/DATA_FINDSYM/data_wyckoff.txt.xz | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\0//g' >> aflow_data_stokes.cpp ; fi
	if [ -f $(EXTRA_FREE)/DATA_FINDSYM/data_wyckoff.txt.xz ]; then echo "\";" >> aflow_data_stokes.cpp ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/data_space.txt.xz ]; then echo "std::string FROZSL_data_space_txt=\"\\" >> aflow_data_stokes.cpp ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/data_space.txt.xz ]; then xzcat $(EXTRA_FREE)/DATA_FROZSL/data_space.txt.xz | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\0//g' >> aflow_data_stokes.cpp ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/data_space.txt.xz ]; then echo "\";" >> aflow_data_stokes.cpp ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/data_wyckoff.txt.xz ]; then echo "std::string FROZSL_data_wyckoff_txt=\"\\" >> aflow_data_stokes.cpp ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/data_wyckoff.txt.xz ]; then xzcat $(EXTRA_FREE)/DATA_FROZSL/data_wyckoff.txt.xz | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\0//g' >> aflow_data_stokes.cpp ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/data_wyckoff.txt.xz ]; then echo "\";" >> aflow_data_stokes.cpp ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/data_images.txt.xz ]; then echo "std::string FROZSL_data_images_txt=\"\\" >> aflow_data_stokes.cpp ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/data_images.txt.xz ]; then xzcat $(EXTRA_FREE)/DATA_FROZSL/data_images.txt.xz | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g'  | perl -p -e 's/\0//g' >> aflow_data_stokes.cpp ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/data_images.txt.xz ]; then echo "\";" >> aflow_data_stokes.cpp ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/data_irreps.txt.xz ]; then echo "std::string FROZSL_data_irreps_txt=\"\\" >> aflow_data_stokes.cpp ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/data_irreps.txt.xz ]; then xzcat $(EXTRA_FREE)/DATA_FROZSL/data_irreps.txt.xz | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\0//g' >> aflow_data_stokes.cpp ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/data_irreps.txt.xz ]; then echo "\";" >> aflow_data_stokes.cpp ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/data_isotropy.txt.xz ]; then echo "std::string FROZSL_data_isotropy_txt=\"\\" >> aflow_data_stokes.cpp ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/data_isotropy.txt.xz ]; then xzcat $(EXTRA_FREE)/DATA_FROZSL/data_isotropy.txt.xz | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\0//g' >> aflow_data_stokes.cpp ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/data_isotropy.txt.xz ]; then echo "\";" >> aflow_data_stokes.cpp ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/data_little.txt.xz ]; then echo "std::string FROZSL_data_little_txt=\"\\" >> aflow_data_stokes.cpp ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/data_little.txt.xz ]; then xzcat $(EXTRA_FREE)/DATA_FROZSL/data_little.txt.xz | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\0//g' >> aflow_data_stokes.cpp ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/data_little.txt.xz ]; then echo "\";" >> aflow_data_stokes.cpp ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/symmetry2.dat.xz ]; then echo "std::string FROZSL_symmetry2_dat=\"\\" >> aflow_data_stokes.cpp ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/symmetry2.dat.xz ]; then xzcat $(EXTRA_FREE)/DATA_FROZSL/symmetry2.dat.xz | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\0//g' >> aflow_data_stokes.cpp ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/symmetry2.dat.xz ]; then echo "\";" >> aflow_data_stokes.cpp ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/const.dat.xz ]; then echo "std::string FROZSL_const_dat=\"\\" >> aflow_data_stokes.cpp ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/const.dat.xz ]; then xzcat $(EXTRA_FREE)/DATA_FROZSL/const.dat.xz | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\0//g' >> aflow_data_stokes.cpp ; fi
	if [ -f $(EXTRA_FREE)/DATA_FROZSL/const.dat.xz ]; then echo "\";" >> aflow_data_stokes.cpp ; fi
	if [ -f $(EXTRA_FREE)/DATA_PHVASP/phvaspsetup_AFLOW.xz ]; then echo "std::string FROZSL_phvaspsetup_AFLOW=\"\\" >> aflow_data_stokes.cpp ; fi
	if [ -f $(EXTRA_FREE)/DATA_PHVASP/phvaspsetup_AFLOW.xz ]; then xzcat $(EXTRA_FREE)/DATA_PHVASP/phvaspsetup_AFLOW.xz | perl -p -e 's/\#/_ASCII_23_/g' | perl -p -e 's/\"/_ASCII_22_/g' | perl -p -e 's/\\/_ASCII_5C_/g' | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_stokes.cpp ; fi
	if [ -f $(EXTRA_FREE)/DATA_PHVASP/phvaspsetup_AFLOW.xz ]; then echo "\";" >> aflow_data_stokes.cpp ; fi
	if [ -f $(EXTRA_FREE)/DATA_PHVASP/phvaspsetup_POSCAR.xz ]; then echo "std::string FROZSL_phvaspsetup_POSCAR=\"\\" >> aflow_data_stokes.cpp ; fi
	if [ -f $(EXTRA_FREE)/DATA_PHVASP/phvaspsetup_POSCAR.xz ]; then xzcat $(EXTRA_FREE)/DATA_PHVASP/phvaspsetup_POSCAR.xz | perl -p -e 's/\#/_ASCII_23_/g' | perl -p -e 's/\"/_ASCII_22_/g' | perl -p -e 's/\\/_ASCII_5C_/g' | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_stokes.cpp ; fi
	if [ -f $(EXTRA_FREE)/DATA_PHVASP/phvaspsetup_POSCAR.xz ]; then echo "\";" >> aflow_data_stokes.cpp ; fi
	if [ ! -f aflow_data_stokes.cpp ]; then rm -rf aflow_data_stokes.cpp.xz*; wget http://materials.duke.edu/AFLOW/AFLOW3_FREE/aflow_data_stokes.cpp.xz && xz -dv aflow_data_stokes.cpp.xz ; fi
	touch $@
aflow_data_nist.cpp: Makefile 
	rm -rf aflow_data_nist.cpp
	if [ -f $(EXTRA_FREE)/NIST/ElectronStoppingPower.txt.xz ]; then echo "// aflow_data_nist.cpp automatic generated from NIST" >> aflow_data_nist.cpp ; fi
	if [ -f $(EXTRA_FREE)/NIST/ElectronStoppingPower.txt.xz ]; then echo "#include <string>" >> aflow_data_nist.cpp ; fi
	if [ -f $(EXTRA_FREE)/NIST/ElectronStoppingPower.txt.xz ]; then echo "std::string ElectronStoppingPower_txt=\"\\" >> aflow_data_nist.cpp ; fi
	if [ -f $(EXTRA_FREE)/NIST/ElectronStoppingPower.txt.xz ]; then xzcat $(EXTRA_FREE)/NIST/ElectronStoppingPower.txt.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_nist.cpp ; fi
	if [ -f $(EXTRA_FREE)/NIST/ElectronStoppingPower.txt.xz ]; then echo "\";" >> aflow_data_nist.cpp ; fi
	if [ -f $(EXTRA_FREE)/NIST/PhotonCrossSection.txt.xz ]; then echo "std::string PhotonCrossSection_txt=\"\\" >> aflow_data_nist.cpp ; fi
	if [ -f $(EXTRA_FREE)/NIST/PhotonCrossSection.txt.xz ]; then xzcat $(EXTRA_FREE)/NIST/PhotonCrossSection.txt.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_nist.cpp ; fi
	if [ -f $(EXTRA_FREE)/NIST/PhotonCrossSection.txt.xz ]; then echo "\";" >> aflow_data_nist.cpp ; fi
	if [ -f $(EXTRA_FREE)/NIST/PhotonStoppingPower.txt.xz ]; then echo "std::string PhotonStoppingPower_txt=\"\\" >> aflow_data_nist.cpp ; fi
	if [ -f $(EXTRA_FREE)/NIST/PhotonStoppingPower.txt.xz ]; then xzcat $(EXTRA_FREE)/NIST/PhotonStoppingPower.txt.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_nist.cpp ; fi
	if [ -f $(EXTRA_FREE)/NIST/PhotonStoppingPower.txt.xz ]; then echo "\";" >> aflow_data_nist.cpp ; fi
	if [ -f $(EXTRA_FREE)/NIST/ICSD_List_2014.txt.xz ]; then echo "std::string ICSD_List_txt=\"\\" >> aflow_data_nist.cpp ; fi
	if [ -f $(EXTRA_FREE)/NIST/ICSD_List_2014.txt.xz ]; then xzcat $(EXTRA_FREE)/NIST/ICSD_List_2014.txt.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_nist.cpp ; fi
	if [ -f $(EXTRA_FREE)/NIST/ICSD_List_2014.txt.xz ]; then echo "\";" >> aflow_data_nist.cpp ; fi
	if [ ! -f aflow_data_nist.cpp ]; then rm -rf aflow_data_nist.cpp.xz*; wget http://materials.duke.edu/AFLOW/AFLOW3_FREE/aflow_data_nist.cpp.xz && xz -dv aflow_data_nist.cpp.xz ; fi
	touch $@
aflow_data_aflow_potcars.cpp: Makefile 
	rm -rf aflow_data_aflow_potcars.cpp
	echo "// aflow_data_aflow_potcars.cpp automatic generated from AFLOW_PSEUDOPOTENTIALS.TXT" >> aflow_data_aflow_potcars.cpp
	echo "#include <string>" >> aflow_data_aflow_potcars.cpp
	if [ -f ../AFLOW3_NONFREE/EXTRA/NONE_AFLOW_PSEUDOPOTENTIALS/AFLOW_PSEUDOPOTENTIALS.TXT.xz ]; then echo "bool AFLOW_PSEUDOPOTENTIALS=1;" >> aflow_data_aflow_potcars.cpp ; else echo "bool AFLOW_PSEUDOPOTENTIALS=0;" >> aflow_data_aflow_potcars.cpp ; fi
	echo "std::string AFLOW_PSEUDOPOTENTIALS_TXT=\"\\" >> aflow_data_aflow_potcars.cpp
	#if[ -f ../AFLOW3_NONFREE/EXTRA/NONE_AFLOW_PSEUDOPOTENTIALS/AFLOW_PSEUDOPOTENTIALS.TXT.xz ]; then xzcat ../AFLOW3_NONFREE/EXTRA/NONE_AFLOW_PSEUDOPOTENTIALS/AFLOW_PSEUDOPOTENTIALS.TXT.xz | perl -p -e 's/\n/\\n\\\n/g' >> aflow_data_aflow_potcars.cpp ; fi
	if [ -f ../AFLOW3_NONFREE/EXTRA/NONE_AFLOW_PSEUDOPOTENTIALS/AFLOW_PSEUDOPOTENTIALS.TXT_XZ_BASE64.xz ]; then xzcat ../AFLOW3_NONFREE/EXTRA/NONE_AFLOW_PSEUDOPOTENTIALS/AFLOW_PSEUDOPOTENTIALS.TXT_XZ_BASE64.xz | perl -p -e 's/\n/\\n\\\n/g' >> aflow_data_aflow_potcars.cpp ; fi
	echo "\";" >> aflow_data_aflow_potcars.cpp
	echo "std::string AFLOW_PSEUDOPOTENTIALS_LIST_TXT=\"\\" >> aflow_data_aflow_potcars.cpp
	#if [ -f ../AFLOW3_NONFREE/EXTRA/NONE_AFLOW_PSEUDOPOTENTIALS/AFLOW_PSEUDOPOTENTIALS.TXT.xz ]; then xzcat ../AFLOW3_NONFREE/EXTRA/NONE_AFLOW_PSEUDOPOTENTIALS/AFLOW_PSEUDOPOTENTIALS.TXT.xz | grep "POTCAR=" | perl -p -e 's/\n/\\n\\\n/g' >> aflow_data_aflow_potcars.cpp ; fi
	if [ -f ../AFLOW3_NONFREE/EXTRA/NONE_AFLOW_PSEUDOPOTENTIALS/AFLOW_PSEUDOPOTENTIALS_LIST.TXT.xz ]; then xzcat ../AFLOW3_NONFREE/EXTRA/NONE_AFLOW_PSEUDOPOTENTIALS/AFLOW_PSEUDOPOTENTIALS_LIST.TXT.xz | grep "POTCAR=" | perl -p -e 's/\n/\\n\\\n/g' >> aflow_data_aflow_potcars.cpp ; fi
	echo "\";" >> aflow_data_aflow_potcars.cpp
aflow_data_latex.cpp: Makefile
	rm -rf aflow_data_latex.cpp
	echo "// aflow_data_latex.cpp automatic generated from LATEX" >> aflow_data_latex.cpp
	echo "#include <string>" >> aflow_data_latex.cpp
	echo "std::string f144468a7ccc2d3a72ba44000715efdb=\"\\" >> aflow_data_latex.cpp
	if [ -f ../AFLOW3_NONFREE/EXTRA/LATEX/f144468a7ccc2d3a72ba44000715efdb ]; then cat ../AFLOW3_NONFREE/EXTRA/LATEX/f144468a7ccc2d3a72ba44000715efdb* | perl -p -e 's/\\/\\\\/g' | perl -p -e 's/\"/\\\"/g'  | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_latex.cpp ; fi
	echo "\";" >> aflow_data_latex.cpp
	echo "std::string d0f1b0e47f178ae627a388d3bf65d2d2=\"\\" >> aflow_data_latex.cpp
	if [ -f ../AFLOW3_NONFREE/EXTRA/LATEX/d0f1b0e47f178ae627a388d3bf65d2d2 ]; then cat ../AFLOW3_NONFREE/EXTRA/LATEX/d0f1b0e47f178ae627a388d3bf65d2d2 | perl -p -e 's/\\/\\\\/g' | perl -p -e 's/\"/\\\"/g'  | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_latex.cpp ; fi
	echo "\";" >> aflow_data_latex.cpp
	echo "std::string decf00ca3ad2fe494eea8e543e929068=\"\\" >> aflow_data_latex.cpp
	if [ -f ../AFLOW3_NONFREE/EXTRA/LATEX/decf00ca3ad2fe494eea8e543e929068 ]; then cat ../AFLOW3_NONFREE/EXTRA/LATEX/decf00ca3ad2fe494eea8e543e929068 | perl -p -e 's/\\/\\\\/g' | perl -p -e 's/\"/\\\"/g'  | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_latex.cpp ; fi
	echo "\";" >> aflow_data_latex.cpp
aflow_data_readme.cpp: README_AFLOW_LICENSE_GPL3.TXT README_AFLOW_ACONVASP.TXT README_AFLOW_APL.TXT README_AFLOW_QHA_SCQHA_QHA3P.TXT README_AFLOW_AGL.TXT README_AFLOW_AEL.TXT README_AFLOW_POCC.TXT README_AFLOW_ANRL.TXT README_AFLOW_COMPARE.TXT README_AFLOW.TXT README_AFLOW_HTRESOURCES.TXT README_PROTO.TXT README_AFLOW_APENNSY.TXT README_AFLOW_FROZSL.TXT README_AFLOW_SCRIPTING.TXT README_AFLOW_SYM.TXT README_AFLOW_CHULL.TXT README_AFLOW_XAFLOW.TXT README_AFLOW_AFLOWRC.TXT README_AFLOW_EXCEPTIONS.TXT README_CONTRIBS.TXT README_AFLOW_VERSIONS_HISTORY.TXT Makefile
	rm -rf aflow_data_readme.cpp
	echo "// aflow_data_readme.cpp automatic generated" >> aflow_data_readme.cpp
	echo "#include <string>" >> aflow_data_readme.cpp
	echo "std::string README_AFLOW_LICENSE_GPL3_TXT=\"\\" >> aflow_data_readme.cpp
	cat README_AFLOW_LICENSE_GPL3.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> aflow_data_readme.cpp
	echo "\";" >> aflow_data_readme.cpp
	echo "std::string README_AFLOW_TXT=\"\\" >> aflow_data_readme.cpp
	cat README_AFLOW.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> aflow_data_readme.cpp
	echo "\";" >> aflow_data_readme.cpp
	echo "std::string README_AFLOW_VERSIONS_HISTORY_TXT=\"\\" >> aflow_data_readme.cpp
	cat README_AFLOW_VERSIONS_HISTORY.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> aflow_data_readme.cpp
	echo "\";" >> aflow_data_readme.cpp
	echo "std::string README_AFLOW_PFLOW_TXT=\"\\" >> aflow_data_readme.cpp
	cat README_AFLOW_ACONVASP.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> aflow_data_readme.cpp
	echo "\";" >> aflow_data_readme.cpp
	echo "std::string README_AFLOW_APENNSY_TXT=\"\\" >> aflow_data_readme.cpp
	cat README_AFLOW_APENNSY.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> aflow_data_readme.cpp
	echo "\";" >> aflow_data_readme.cpp
	echo "std::string README_AFLOW_SCRIPTING_TXT=\"\\" >> aflow_data_readme.cpp
	cat README_AFLOW_SCRIPTING.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> aflow_data_readme.cpp
	echo "\";" >> aflow_data_readme.cpp
	echo "std::string README_AFLOW_FROZSL_TXT=\"\\" >> aflow_data_readme.cpp
	cat README_AFLOW_FROZSL.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> aflow_data_readme.cpp
	echo "\";" >> aflow_data_readme.cpp
	echo "std::string README_AFLOW_POCC_TXT=\"\\" >> aflow_data_readme.cpp
	cat README_AFLOW_POCC.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> aflow_data_readme.cpp
	echo "\";" >> aflow_data_readme.cpp
	echo "std::string README_AFLOW_APL_TXT=\"\\" >> aflow_data_readme.cpp
	cat README_AFLOW_APL.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> aflow_data_readme.cpp
	echo "\";" >> aflow_data_readme.cpp
	echo "std::string README_AFLOW_QHA_SCQHA_QHA3P_TXT=\"\\" >> aflow_data_readme.cpp
	cat README_AFLOW_QHA_SCQHA_QHA3P.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> aflow_data_readme.cpp
	echo "\";" >> aflow_data_readme.cpp
	echo "std::string README_AFLOW_AGL_TXT=\"\\" >> aflow_data_readme.cpp
	cat README_AFLOW_AGL.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> aflow_data_readme.cpp
	echo "\";" >> aflow_data_readme.cpp
	echo "std::string README_AFLOW_AEL_TXT=\"\\" >> aflow_data_readme.cpp
	cat README_AFLOW_AEL.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> aflow_data_readme.cpp
	echo "\";" >> aflow_data_readme.cpp
	echo "std::string README_AFLOW_ANRL_TXT=\"\\" >> aflow_data_readme.cpp
	cat README_AFLOW_ANRL.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> aflow_data_readme.cpp
	echo "\";" >> aflow_data_readme.cpp
	echo "std::string README_AFLOW_COMPARE_TXT=\"\\" >> aflow_data_readme.cpp
	cat README_AFLOW_COMPARE.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> aflow_data_readme.cpp
	echo "\";" >> aflow_data_readme.cpp
	echo "std::string README_AFLOW_SYM_TXT=\"\\" >> aflow_data_readme.cpp
	cat README_AFLOW_SYM.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> aflow_data_readme.cpp
	echo "\";" >> aflow_data_readme.cpp
	echo "std::string README_AFLOW_CHULL_TXT=\"\\" >> aflow_data_readme.cpp
	cat README_AFLOW_CHULL.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> aflow_data_readme.cpp
	echo "\";" >> aflow_data_readme.cpp
	echo "std::string README_AFLOW_EXCEPTIONS_TXT=\"\\" >> aflow_data_readme.cpp
	cat README_AFLOW_EXCEPTIONS.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> aflow_data_readme.cpp
	echo "\";" >> aflow_data_readme.cpp
	echo "std::string README_AFLOW_HTRESOURCES_TXT=\"\\" >> aflow_data_readme.cpp
	cat README_AFLOW_HTRESOURCES.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> aflow_data_readme.cpp
	echo "\";" >> aflow_data_readme.cpp
	echo "std::string README_PROTO_TXT=\"\\" >> aflow_data_readme.cpp
	cat README_PROTO.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> aflow_data_readme.cpp
	echo "\";" >> aflow_data_readme.cpp
	echo "std::string README_AFLOW_XAFLOW_TXT=\"\\" >> aflow_data_readme.cpp
	cat README_AFLOW_XAFLOW.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> aflow_data_readme.cpp
	echo "\";" >> aflow_data_readme.cpp
	echo "std::string README_AFLOW_AFLOWRC_TXT=\"\\" >> aflow_data_readme.cpp
	cat README_AFLOW_AFLOWRC.TXT | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' | perl -p -e 's/\"/\\"/g' >> aflow_data_readme.cpp
	echo "\";" >> aflow_data_readme.cpp
aflow_data_extra.cpp: Makefile aflow_library_icsd.dat.xz
	rm -rf aflow_data_extra.cpp aflow_data_libraries.cpp
	echo "#include <string>" >> aflow_data_libraries.cpp
	echo "std::string icsd_1_ary=\"\\" >> aflow_data_extra.cpp && xzcat $(NIST_LIB)/1-ary.icsd.xz | grep -v "//" | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\?/\? /g' | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_extra.cpp && echo "\";" >> aflow_data_extra.cpp	
	echo "std::string icsd_2_ary=\"\\" >> aflow_data_extra.cpp && xzcat $(NIST_LIB)/2-ary.icsd.xz | grep -v "//" | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\?/\? /g' | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_extra.cpp && echo "\";" >> aflow_data_extra.cpp	
	echo "std::string icsd_3_ary=\"\\" >> aflow_data_extra.cpp && xzcat $(NIST_LIB)/3-ary.icsd.xz | grep -v "//" | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\?/\? /g' | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_extra.cpp && echo "\";" >> aflow_data_extra.cpp	
	echo "std::string icsd_4_ary=\"\\" >> aflow_data_extra.cpp && xzcat $(NIST_LIB)/4-ary.icsd.xz | grep -v "//" | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\?/\? /g' | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_extra.cpp && echo "\";" >> aflow_data_extra.cpp	
	echo "std::string icsd_5_ary=\"\\" >> aflow_data_extra.cpp && xzcat $(NIST_LIB)/5-ary.icsd.xz | grep -v "//" | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\?/\? /g' | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_extra.cpp && echo "\";" >> aflow_data_extra.cpp	
	echo "std::string icsd_6_ary=\"\\" >> aflow_data_extra.cpp && xzcat $(NIST_LIB)/6-ary.icsd.xz | grep -v "//" | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\?/\? /g' | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_extra.cpp && echo "\";" >> aflow_data_extra.cpp	
	echo "std::string icsd_7_ary=\"\\" >> aflow_data_extra.cpp && xzcat $(NIST_LIB)/7-ary.icsd.xz | grep -v "//" | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\?/\? /g' | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_extra.cpp && echo "\";" >> aflow_data_extra.cpp	
	echo "std::string icsd_8_ary=\"\\" >> aflow_data_extra.cpp && xzcat $(NIST_LIB)/8-ary.icsd.xz | grep -v "//" | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\?/\? /g' | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_extra.cpp && echo "\";" >> aflow_data_extra.cpp	
	echo "std::string icsd_9_ary=\"\\" >> aflow_data_extra.cpp && xzcat $(NIST_LIB)/9-ary.icsd.xz | grep -v "//" | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\?/\? /g' | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_extra.cpp && echo "\";" >> aflow_data_extra.cpp	
	echo "std::string icsd_10_ary=\"\\" >> aflow_data_extra.cpp && xzcat $(NIST_LIB)/10-ary.icsd.xz | grep -v "//" | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\?/\? /g' | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_extra.cpp && echo "\";" >> aflow_data_extra.cpp
	echo "std::string icsd_11_ary=\"\\" >> aflow_data_extra.cpp && xzcat $(NIST_LIB)/11-ary.icsd.xz | grep -v "//" | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\?/\? /g' | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_extra.cpp && echo "\";" >> aflow_data_extra.cpp
	echo "std::string icsd_12_ary=\"\\" >> aflow_data_extra.cpp && xzcat $(NIST_LIB)/12-ary.icsd.xz | grep -v "//" | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\?/\? /g' | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_extra.cpp && echo "\";" >> aflow_data_extra.cpp
	echo "std::string icsd_13_ary=\"\\" >> aflow_data_extra.cpp && xzcat $(NIST_LIB)/13-ary.icsd.xz | grep -v "//" | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\?/\? /g' | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_extra.cpp && echo "\";" >> aflow_data_extra.cpp
	echo "std::string icsd_14_ary=\"\\" >> aflow_data_extra.cpp && xzcat $(NIST_LIB)/14-ary.icsd.xz | grep -v "//" | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\?/\? /g' | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_extra.cpp && echo "\";" >> aflow_data_extra.cpp
	echo "std::string icsd_15_ary=\"\\" >> aflow_data_extra.cpp && xzcat $(NIST_LIB)/15-ary.icsd.xz | grep -v "//" | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\?/\? /g' | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_extra.cpp && echo "\";" >> aflow_data_extra.cpp
#	echo "std::string Library_ICSD=\"\\" >> aflow_data_extra.cpp && xzcat  ../AFLOW3_NONFREE/LIBS/aflow_library_icsd.dat.xz | grep -v "//" | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\?/\? /g' | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_extra.cpp && echo "\";" >> aflow_data_extra.cpp
	echo "std::string aflowlib_lib1=\"\\" >> aflow_data_libraries.cpp && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib1.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_libraries.cpp && echo "\";" >> aflow_data_libraries.cpp	
	echo "std::string aflowlib_lib2=\"\\" >> aflow_data_libraries.cpp && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib2.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_libraries.cpp && echo "\";" >> aflow_data_libraries.cpp	
	echo "std::string aflowlib_lib3=\"\\" >> aflow_data_libraries.cpp && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib3.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_libraries.cpp && echo "\";" >> aflow_data_libraries.cpp
	echo "std::string aflowlib_lib4=\"\\" >> aflow_data_libraries.cpp && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib4.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_libraries.cpp && echo "\";" >> aflow_data_libraries.cpp
	echo "std::string aflowlib_lib5=\"\\" >> aflow_data_libraries.cpp && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib5.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_libraries.cpp && echo "\";" >> aflow_data_libraries.cpp
	echo "std::string aflowlib_lib6=\"\\" >> aflow_data_libraries.cpp && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib6.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_libraries.cpp && echo "\";" >> aflow_data_libraries.cpp
	echo "std::string aflowlib_lib7=\"\\" >> aflow_data_libraries.cpp && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib7.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_libraries.cpp && echo "\";" >> aflow_data_libraries.cpp
	echo "std::string aflowlib_lib8=\"\\" >> aflow_data_libraries.cpp && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib8.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_libraries.cpp && echo "\";" >> aflow_data_libraries.cpp
	echo "std::string aflowlib_lib9=\"\\" >> aflow_data_libraries.cpp && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_lib9.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_libraries.cpp && echo "\";" >> aflow_data_libraries.cpp
	echo "std::string aflowlib_icsd=\"\\" >> aflow_data_libraries.cpp && xzcat ../AFLOW3_NONFREE/EXTRA/CALCULATED/aflowlib_icsd.dat.xz | grep -v "//" | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_libraries.cpp && echo "\";" >> aflow_data_libraries.cpp
#aflow_data: aflow_data_readme.cpp aflow_data_calculated.cpp aflow_data_htqc.cpp aflow_data_libraries.cpp aflow_data_stokes.cpp aflow_data_nist.cpp aflow_data_latex.cpp aflow_data_aflow_potcars.cpp aflow_data.cpp
aflow_data: aflow_data_readme.o aflow_data_calculated.o aflow_data_htqc.o aflow_data_libraries.o aflow_data_stokes.o aflow_data_nist.o aflow_data_latex.o aflow_data_aflow_potcars.o aflow_data.cpp
	rm -f aflow_data_extra.cpp
	echo "" >> aflow_data_extra.cpp
	if [ -f ../AFLOW3_NONFREE/LIBS/aflow_library_icsd.dat.xz ]; then echo "#define defined_Library_ICSD" >> aflow_data_extra.cpp; fi
	if [ -f ../AFLOW3_NONFREE/LIBS/aflow_library_icsd.dat.xz ]; then echo "std::string Library_ICSD=\"\\" >> aflow_data_extra.cpp; fi
	if [ -f ../AFLOW3_NONFREE/LIBS/aflow_library_icsd.dat.xz ]; then xzcat  ../AFLOW3_NONFREE/LIBS/aflow_library_icsd.dat.xz | grep -v "//" | perl -p -e 's/\"/\\"/g' | perl -p -e 's/\?/\? /g' | perl -p -e 's/\n/ \\n\\\n/g' >> aflow_data_extra.cpp; fi
	if [ -f ../AFLOW3_NONFREE/LIBS/aflow_library_icsd.dat.xz ]; then echo "\";" >> aflow_data_extra.cpp; fi
	$(CPP) $(VERS) $(INCLUDE) $(CCFLAGS) $(OPTS) -c aflow_data.cpp
	$(CPP) $(VERS) $(INCLUDE) $(CCFLAGS) $(OPTS) -o aflow_data aflow_data.o aflow_data_*.o 
amir:	aflow_data_readme.cpp aflow_data_calculated.cpp aflow_data_htqc.cpp aflow_data_libraries.cpp aflow_data_stokes.cpp aflow_data_nist.cpp aflow_data_latex.cpp aflow_data_aflow_potcars.cpp aflow_data.cpp aflow_data_extra.cpp
	$(CPP) $(VERS) $(INCLUDE) $(CCFLAGS) $(OPTS) -DICSD -o aflow_data aflow_data.cpp
############################################################################################################
aflow_phonons.o: aflow_phonons.cpp $(AFLOW_HPP) $(AUROSTD_HPP)  Makefile
	$(CPP) $(VERS) $(INCLUDE) $(CCFLAGS) $(OPTS) aflow_phonons.cpp -c
aflow_matlab_funcs.cpp: Makefile
	rm -rf aflow_matlab_funcs.cpp
	if [ -f $(EXTRA_FREE)/MATLAB/param.m ]; then echo "// aflow_matlab_funcs.cpp automatic generated from MATLAB/*" >>  aflow_matlab_funcs.cpp ; fi
	if [ -f $(EXTRA_FREE)/MATLAB/param.m ]; then echo "#include <sstream>" >> aflow_matlab_funcs.cpp ; fi
	if [ -f $(EXTRA_FREE)/MATLAB/param.m ]; then echo "// aflow_matlab_funcs.cpp automatic generated from MATLAB/param.m" >> aflow_matlab_funcs.cpp ; fi
	if [ -f $(EXTRA_FREE)/MATLAB/param.m ]; then echo "std::string MATLAB_FUNCS_param(void){" >> aflow_matlab_funcs.cpp ; fi
	if [ -f $(EXTRA_FREE)/MATLAB/param.m ]; then echo "   std::stringstream strstream; strstream << \"\\" >> aflow_matlab_funcs.cpp ; fi
	if [ -f $(EXTRA_FREE)/MATLAB/param.m ]; then cat $(EXTRA_FREE)/MATLAB/param.m | perl -p -e 's/\"/\\\"/g' | perl -p -e 's/\\G/\\\\G/g' | perl -p -e 's/\\S/\\\\S/g' | perl -p -e 's/\n/\" << std::endl;\n   strstream << \"/g' >> aflow_matlab_funcs.cpp ; fi
	if [ -f $(EXTRA_FREE)/MATLAB/param.m ]; then echo "\" << std::endl; return strstream.str();};" >> aflow_matlab_funcs.cpp ; fi
	if [ -f $(EXTRA_FREE)/MATLAB/plotband.m ]; then echo "// aflow_matlab_funcs.cpp automatic generated from MATLAB/plotband.m" >> aflow_matlab_funcs.cpp ; fi
	if [ -f $(EXTRA_FREE)/MATLAB/plotband.m ]; then echo "// std::string MATLAB_FUNCS_plotband(std::string DIRECTORY,std::string OPTION1){" >> aflow_matlab_funcs.cpp ; fi
	if [ -f $(EXTRA_FREE)/MATLAB/plotband.m ]; then echo "std::string MATLAB_FUNCS_plotband(void){" >> aflow_matlab_funcs.cpp ; fi
	if [ -f $(EXTRA_FREE)/MATLAB/plotband.m ]; then echo "   std::stringstream strstream; strstream << \"\\" >> aflow_matlab_funcs.cpp ; fi
	if [ -f $(EXTRA_FREE)/MATLAB/plotband.m ]; then cat $(EXTRA_FREE)/MATLAB/plotband.m | perl -p -e 's/\"/\\\"/g' | perl -p -e 's/\\G/\\\\G/g' | perl -p -e 's/\\S/\\\\S/g' | perl -p -e 's/\n/\" << std::endl;\n   strstream << \"/g' | perl -p -e 's/MATLAB_FUNCS_DIRECTORY/\" << DIRECTORY << \"/g' | perl -p -e 's/DOS_SCALE_MODE/\" << OPTION1 << \"/g' >> aflow_matlab_funcs.cpp ; fi
	if [ -f $(EXTRA_FREE)/MATLAB/plotband.m ]; then echo "\" << std::endl; return strstream.str();};" >> aflow_matlab_funcs.cpp ; fi
	if [ ! -f aflow_matlab_funcs.cpp ]; then rm -rf aflow_matlab_funcs.cpp.xz*; wget http://materials.duke.edu/AFLOW/AFLOW3_FREE/aflow_matlab_funcs.cpp.xz && xz -dv aflow_matlab_funcs.cpp.xz ; fi
	touch $@
aflow_gnuplot_funcs.cpp: Makefile
	rm -rf aflow_gnuplot_funcs.cpp
	if [ -f $(EXTRA_FREE)/GNUPLOT/plotbz.sh ]; then echo "// aflow_gnuplot_funcs.cpp automatic generated from GNUPLOT/*" >> aflow_gnuplot_funcs.cpp ; fi
	if [ -f $(EXTRA_FREE)/GNUPLOT/plotbz.sh ]; then echo "#include <sstream>" >> aflow_gnuplot_funcs.cpp ; fi
	if [ -f $(EXTRA_FREE)/GNUPLOT/plotbz.sh ]; then echo "// aflow_gnuplot_funcs.cpp automatic generated from GNUPLOT/plotbz.sh" >> aflow_gnuplot_funcs.cpp ; fi
	if [ -f $(EXTRA_FREE)/GNUPLOT/plotbz.sh ]; then echo "//std::string GNUPLOT_FUNCS_plotbz(std::string DIRECTORY,std::string OPTION1){" >> aflow_gnuplot_funcs.cpp ; fi
	if [ -f $(EXTRA_FREE)/GNUPLOT/plotbz.sh ]; then echo "std::string GNUPLOT_FUNCS_plotbz(void){" >> aflow_gnuplot_funcs.cpp ; fi
	if [ -f $(EXTRA_FREE)/GNUPLOT/plotbz.sh ]; then echo "   std::stringstream strstream; strstream << \"\\" >> aflow_gnuplot_funcs.cpp ; fi
	if [ -f $(EXTRA_FREE)/GNUPLOT/plotbz.sh ]; then cat $(EXTRA_FREE)/GNUPLOT/plotbz.sh | perl -p -e 's/\\/\\\\/g' | perl -p -e 's/\"/\\\"/g' | perl -p -e 's/\\G/\\\\G/g' | perl -p -e 's/\\S/\\\\S/g' | perl -p -e 's/\n/\" << std::endl;\n   strstream << \"/g' | perl -p -e 's/GNUPLOT_FUNCS_DIRECTORY/\" << DIRECTORY << \"/g' | perl -p -e 's/GNUPLOT_FUNCS_OPTION1/\" << OPTION1 << \"/g' >> aflow_gnuplot_funcs.cpp ; fi
	if [ -f $(EXTRA_FREE)/GNUPLOT/plotbz.sh ]; then echo "\" << std::endl; return strstream.str();};" >> aflow_gnuplot_funcs.cpp ; fi
	if [ ! -f aflow_gnuplot_funcs.cpp ]; then rm -rf aflow_gnuplot_funcs.cpp.xz*; wget http://materials.duke.edu/AFLOW/AFLOW3_FREE/aflow_gnuplot_funcs.cpp.xz && xz -dv aflow_gnuplot_funcs.cpp.xz ; fi
	touch $@
aurostd_xtemplates.o: aurostd_xtemplates.cpp $(AFLOW_HPP) $(AUROSTD_HPP) Makefile
	$(CPP) $(VERS) $(INCLUDE) $(CCFLAGS) $(OPTS) aurostd_xtemplates.cpp -c
#
###### HTRESOURCES #################################################
aflow_contrib_stefano_htcurriculum.o: $(AFLOWLIB_HPP) $(AFLOW_HPP) $(AUROSTD_HPP) Makefile
	$(CPP) $(VERS) $(INCLUDE) $(CCFLAGS) $(OPTS)  -c
aflow_contrib_stefano_htsecurity.o: aflow_contrib_stefano_htsecurity.cpp $(AFLOWLIB_HPP) $(AFLOW_HPP) $(AUROSTD_HPP) Makefile
	$(CPP) $(VERS) $(INCLUDE) $(CCFLAGS) $(OPTS) aflow_contrib_stefano_htsecurity.cpp -c
aflowlib_web_outreach.o: aflowlib_web_outreach.cpp $(AFLOWLIB_HPP) $(AFLOW_HPP) $(AUROSTD_HPP) Makefile
	$(CPP) $(VERS) $(INCLUDE) $(CCFLAGS) $(OPTS) aflowlib_web_outreach.cpp -c

###### SHIDONG #################################################
aflow_contrib_shidong_main.o: aflow_contrib_shidong_main.cpp aflow_contrib_shidong_main.h aflow_contrib_shidong.h aflow_contrib_shidong_cluster_expansion.h aflow_contrib_shidong_auxiliary.h 
aflow_contrib_shidong.o: aflow_contrib_shidong.cpp aflow_contrib_shidong.h aflow_contrib_shidong_auxiliary.h aflow_contrib_shidong_cluster_expansion.h 
aflow_contrib_shidong_auxiliary.o: aflow_contrib_shidong_auxiliary.cpp aflow_contrib_shidong_auxiliary.h aflow_contrib_shidong_cluster_expansion.h
aflow_contrib_shidong_cexstructure.o: aflow_contrib_shidong_cluster_expansion.cpp aflow_contrib_shidong_cluster_expansion.h
aflow_contrib_shidong_mc.o: aflow_contrib_shidong_mc.cpp aflow_contrib_shidong_mc.h aflow_contrib_shidong_statistic.h aflow_contrib_shidong_cluster_expansion.h 
aflow_contrib_shidong_statistic.o: aflow_contrib_shidong_statistic.cpp aflow_contrib_shidong_statistic.h

#####COMPARE###################################################
#aflow_compare_structure_function.o: aflow_compare_structure_function.cpp Makefile $(AFLOW_HPP) $(AUROSTD_HPP) $(AUROSTD_CPP)
#	${CPP} $(VERS) ${INCLUDE} ${CFLAGS} ${LIBS} -std=gnu++11 aflow_compare_structure_function.cpp -c
#aflow_compare_structure.o: aflow_compare_structure.cpp Makefile $(AFLOW_HPP) $(AUROSTD_HPP) $(AUROSTD_CPP)
#	${CPP} $(VERS) ${INCLUDE} ${CFLAGS} ${LIBS} -std=gnu++11 aflow_compare_structure.cpp -c

#[OBSOLETE] [JUNKAI]###### JUNKAI #################################################
#[OBSOLETE] [JUNKAI] aflow_contrib_junkai_cproto.o:aflow_contrib_junkai_cproto.cpp aflow_contrib_junkai_basic.h
#[OBSOLETE] [JUNKAI] aflow_contrib_junkai_primitivetrans.o:aflow_contrib_junkai_primitivetrans.cpp aflow_contrib_junkai_basic.h
#[OBSOLETE] [JUNKAI] aflow_contrib_junkai_basic.o:aflow_contrib_junkai_basic.cpp aflow_contrib_junkai_basic.h
#[OBSOLETE] [JUNKAI] aflow_contrib_junkai_protocheck.o:aflow_contrib_junkai_protocheck.cpp aflow_contrib_junkai_basic.h
#[OBSOLETE] [JUNKAI] aflow_contrib_junkai_magnet.o:aflow_contrib_junkai_magnet.cpp aflow_contrib_junkai_basic.h
#[OBSOLETE] [JUNKAI] aflow_contrib_junkai_phasediag.o:aflow_contrib_junkai_phasediag.cpp aflow_contrib_junkai_basic.h
###### ACONVASP #################################################
aflow_aconvasp_main.o: aflow_aconvasp_main.cpp $(AFLOW_HPP) $(AUROSTD_HPP) Makefile
	$(CPP) $(VERS) $(INCLUDE) $(CCFLAGS) $(OPTS) aflow_aconvasp_main.cpp -c
aflow_apennsy.o: aflow_apennsy.cpp $(AFLOW_HPP) $(AUROSTD_HPP) aflow_apennsy.h Makefile
	$(CPP) $(VERS) $(INCLUDE) $(CCFLAGS) $(OPTS) aflow_apennsy.cpp -c
aflow_apennsy_main.o: aflow_apennsy_main.cpp $(AFLOW_HPP) $(AUROSTD_HPP) aflow_apennsy.h Makefile
	$(CPP) $(VERS) $(INCLUDE) $(CCFLAGS) $(OPTS) aflow_apennsy_main.cpp -c
aflow_wyckoff.o: aflow_wyckoff.cpp $(AFLOW_HPP) $(AUROSTD_HPP)
	$(CPP) $(VERS) -D_WICKOFF_PRE_ aflow_wyckoff.cpp -c
####### APL STUFF ####################################################
#DX and CO - START
APL/aflow_apl_doscalc.o: APL/aflow_apl_doscalc.cpp $(AFLOW_HPP) $(AUROSTD_HPP) Makefile
#	${CPP} $(VERS) $(INCLUDE) $(CCFLAGS) ${LIBS} -o $@ $<
	${CPP} $(VERS) $(INCLUDE) $(CCFLAGS_MT) -std=gnu++0x APL/aflow_apl_doscalc.cpp -c
	touch $@
APL/aflow_apl_pdisc.o: APL/aflow_apl_pdisc.cpp $(AFLOW_HPP) $(AUROSTD_HPP) Makefile
#	${CPP} $(VERS) $(INCLUDE) $(CCFLAGS) ${LIBS} -o $@ $<
	${CPP} $(VERS) $(INCLUDE) $(CCFLAGS_MT) -std=gnu++0x APL/aflow_apl_pdisc.cpp -c
	touch $@
APL/aflow_aapl_tcond.o: APL/aflow_aapl_tcond.cpp $(AFLOW_HPP) $(AUROSTD_HPP) Makefile
#	${CPP} $(VERS) $(INCLUDE) $(CCFLAGS) ${LIBS} -o $@ $<
	${CPP} $(VERS) $(INCLUDE) $(CCFLAGS_MT) -std=gnu++0x APL/aflow_aapl_tcond.cpp -c
	touch $@
APL/aflow_qha_dm_pdos_save.o: APL/aflow_qha_dm_pdos_save.cpp $(AFLOW_HPP) $(AUROSTD_HPP) Makefile
#	${CPP} $(VERS) $(INCLUDE) $(CCFLAGS) ${LIBS} -o $@ $<
	${CPP} $(VERS) $(INCLUDE) $(CCFLAGS_MT) -std=gnu++0x APL/aflow_qha_dm_pdos_save.cpp -c
	touch $@
APL/aflow_apl_group_velocity.o: APL/aflow_apl_group_velocity.cpp $(AFLOW_HPP) $(AUROSTD_HPP) Makefile
#	${CPP} $(VERS) $(INCLUDE) $(CCFLAGS) ${LIBS} -o $@ $<
	${CPP} $(VERS) $(INCLUDE) $(CCFLAGS_MT) -std=gnu++0x APL/aflow_apl_group_velocity.cpp -c
	touch $@
APL/aflow_apl_atomic_disp.o: APL/aflow_apl_atomic_disp.cpp $(AFLOW_HPP) $(AUROSTD_HPP) Makefile
#	${CPP} $(VERS) $(INCLUDE) $(CCFLAGS) ${LIBS} -o $@ $<
	${CPP} $(VERS) $(INCLUDE) $(CCFLAGS_MT) -std=gnu++0x APL/aflow_apl_atomic_disp.cpp -c
	touch $@
APL/aflow_qha_gruneisen.o: APL/aflow_qha_gruneisen.cpp $(AFLOW_HPP) $(AUROSTD_HPP) Makefile
#	${CPP} $(VERS) $(INCLUDE) $(CCFLAGS) ${LIBS} -o $@ $<
	${CPP} $(VERS) $(INCLUDE) $(CCFLAGS_MT) -std=gnu++0x APL/aflow_qha_gruneisen.cpp -c
	touch $@
APL/aflow_scqha_gruneisen.o: APL/aflow_scqha_gruneisen.cpp $(AFLOW_HPP) $(AUROSTD_HPP) Makefile
#	${CPP} $(VERS) $(INCLUDE) $(CCFLAGS) ${LIBS} -o $@ $<
	${CPP} $(VERS) $(INCLUDE) $(CCFLAGS_MT) -std=gnu++0x APL/aflow_scqha_gruneisen.cpp -c
	touch $@
APL/aflow_scqha_eos.o: APL/aflow_scqha_eos.cpp $(AFLOW_HPP) $(AUROSTD_HPP) Makefile
#	${CPP} $(VERS) $(INCLUDE) $(CCFLAGS) ${LIBS} -o $@ $<
	${CPP} $(VERS) $(INCLUDE) $(CCFLAGS_MT) -std=gnu++0x APL/aflow_scqha_eos.cpp -c
	touch $@
APL/aflow_qha3phonons_eos.o: APL/aflow_qha3phonons_eos.cpp $(AFLOW_HPP) $(AUROSTD_HPP) Makefile
#	${CPP} $(VERS) $(INCLUDE) $(CCFLAGS) ${LIBS} -o $@ $<
	${CPP} $(VERS) $(INCLUDE) $(CCFLAGS_MT) -std=gnu++0x APL/aflow_qha3phonons_eos.cpp -c
	touch $@
APL/aflow_scqha_T_freqs.o: APL/aflow_scqha_T_freqs.cpp $(AFLOW_HPP) $(AUROSTD_HPP) Makefile
#	${CPP} $(VERS) $(INCLUDE) $(CCFLAGS) ${LIBS} -o $@ $<
	${CPP} $(VERS) $(INCLUDE) $(CCFLAGS_MT) -std=gnu++0x APL/aflow_scqha_T_freqs.cpp -c
	touch $@
#DX - CO - END
#APL/aflow_apl_gsa.o: APL/aflow_apl_gsa.cpp $(AFLOW_HPP) $(AUROSTD_HPP) Makefile
##	${CPP} $(VERS) $(INCLUDE) $(CCFLAGS) ${LIBS} -o $@ $<
#	${CPP} $(VERS) $(INCLUDE) $(CCFLAGS_GSA) APL/aflow_apl_gsa.cpp -c
#	touch $@

BIN_tscc = /projects/kyang-group/bin
################################################################################################
install: all
	cp -f aflow $(BIN)/aflow.$(VERSNUMBER).$(DAY)
	cp -f aflow_data  $(BIN)/
	ls -las $(BIN)/aflow*
	rm -f $(BIN)/aflow 
	ln -s $(BIN)/aflow.$(VERSNUMBER).$(DAY)	 $(BIN)/aflow
	rm -f $(BIN)/aconvasp 
	ln -s $(BIN)/aflow  $(BIN)/aconvasp
	ls -las $(BIN)/aconvasp
	$(MAKE) check
#Add by Kesong, 2019-07-01
install_tscc: all
	cp -f aflow $(BIN_tscc)/aflow.$(VERSNUMBER).$(DAY)
	cp -f aflow_data  $(BIN_tscc)/
	ls -las $(BIN_tscc)/aflow*
	rm -f $(BIN_tscc)/aflow 
	ln -s $(BIN_tscc)/aflow.$(VERSNUMBER).$(DAY)	 $(BIN_tscc)/aflow
	rm -f $(BIN_tscc)/aconvasp 
	ln -s $(BIN_tscc)/aflow  $(BIN_tscc)/aconvasp
	ls -las $(BIN_tscc)/aconvasp
	$(MAKE) check
install_mac: mac        
	cp -f aflow $(BIN)/aflow.$(VERSNUMBER).$(DAY)
	cp -f aflow_data  $(BIN)/
	ls -las $(BIN)/aflow*
	rm -f $(BIN)/aflow 
	ln -s $(BIN)/aflow.$(VERSNUMBER).$(DAY)	 $(BIN)/aflow
	rm -f $(BIN)/aconvasp 
	ln -s $(BIN)/aflow  $(BIN)/aconvasp
	ls -las $(BIN)/aconvasp
	$(MAKE) check
install_auro: install_user
	./aflow --version | grep AFLOW >> $(BIN)/etc_aflow_status.txt
	chmod 644 $(BIN)/etc_aflow_status.txt
	./aflow --calculated | grep Calcs >> $(BIN)/etc_aflowlib_status.txt
	chmod 644 $(BIN)/etc_aflowlib_status.txt
	./aflow --apool > $(BIN)/etc_aflowlib_apool.txt
	chmod 644 $(BIN)/etc_aflowlib_apool.txt
install_root: all
	cp -f aflow /usr/local/bin/aflow.$(VERSNUMBER)
	cp -f aflow_data /usr/local/bin/aflow_data
	cp -f /usr/local/bin/aflow.$(VERSNUMBER) /usr/local/bin/aflow
	ls -las /usr/local/bin/aflow*
	$(MAKE) check
niet:	
	dropbox stop
	$(MAKE) clean
	$(MAKE) icsdclean
	$(MAKE) extra_clean
	if [ -f $(EXTRA_FREE)/GUS/Makefile ]; then $(MAKE) -C $(EXTRA_FREE)/GUS realclean; fi
	echo "" > $(ICSD_LIB)/README_LIBRARY_ICSD1.TXT
	echo "" > $(ICSD_LIB)/README_LIBRARY_ICSD2.TXT
	echo "" > $(ICSD_LIB)/README_LIBRARY_ICSD3.TXT
	echo "" > $(ICSD_LIB)/README_LIBRARY_ICSD4.TXT
	echo "" > $(ICSD_LIB)/README_LIBRARY_ICSD5.TXT
	echo "" > $(ICSD_LIB)/README_LIBRARY_ICSD6.TXT
	echo "" > $(ICSD_LIB)/README_LIBRARY_ICSD7.TXT
	echo "" > $(ICSD_LIB)/README_LIBRARY_ICSD8.TXT
	$(MAKE) web
	$(MAKE) backup
	$(MAKE) icsdclean
	$(MAKE) -j 32
	rm -fv aflow_data_calculated.cpp* && $(MAKE) aflow_data_calculated.cpp && cp -f aflow_data_calculated.cpp /www/AFLOW/AFLOW3_FREE/ && xz -9fv /www/AFLOW/AFLOW3_FREE/aflow_data_calculated.cpp
	rm -fv aflow_data_stokes.cpp* && $(MAKE) aflow_data_stokes.cpp && cp -f aflow_data_stokes.cpp /www/AFLOW/AFLOW3_FREE/ && xz -9fv /www/AFLOW/AFLOW3_FREE/aflow_data_stokes.cpp
	rm -fv aflow_data_nist.cpp* && $(MAKE) aflow_data_nist.cpp && cp -f aflow_data_nist.cpp /www/AFLOW/AFLOW3_FREE/ && xz -9fv /www/AFLOW/AFLOW3_FREE/aflow_data_nist.cpp
	rm -fv aflow_matlab_funcs.cpp* && $(MAKE) aflow_matlab_funcs.cpp && cp -f aflow_matlab_funcs.cpp /www/AFLOW/AFLOW3_FREE/ && xz -9fv /www/AFLOW/AFLOW3_FREE/aflow_matlab_funcs.cpp
	rm -fv aflow_gnuplot_funcs.cpp* && $(MAKE) aflow_gnuplot_funcs.cpp && cp -f aflow_gnuplot_funcs.cpp /www/AFLOW/AFLOW3_FREE/ && xz -9fv /www/AFLOW/AFLOW3_FREE/aflow_gnuplot_funcs.cpp
	dropbox start
backup: realclean clean icsdclean 
	if [ -f $(EXTRA_FREE)/GUS/Makefile ]; then $(MAKE) -C $(EXTRA_FREE)/GUS realclean ; fi
	mkdir -p aflow.$(VERSNUMBER)
#	mkdir -p aflow.$(VERSNUMBER)/SOURCES
	cp -drvf a*cpp a*.h *.js README* Makefile APL ANRL* AUROSTD SQLITE aflow.$(VERSNUMBER)
# no EXTRA in SRC	if [ -d $(EXTRA_FREE) ]; then cp -drvf $(EXTRA_FREE) aflow.$(VERSNUMBER) ; fi
	if [ -d SOURCE ]; then cp -drvf SOURCE aflow.$(VERSNUMBER) ; fi
	tar Jcfvp aflow.$(VERSNUMBER).tar.xz aflow.$(VERSNUMBER)
	rm -rfv $(BACKUP)/aflow.$(VERSNUMBER) $(BACKUP)/aflow.$(VERSNUMBER).tar.xz
	if [ -d /www/AFLOW ]; then cp -f aflow.$(VERSNUMBER).tar.xz /www/AFLOW ; fi
	mv aflow.$(VERSNUMBER).tar.xz $(BACKUP)
	rm -rfv aflow.$(VERSNUMBER)
web:	realclean
	rm -rvf aflow3.tar.xz aflow_test.tar.xz aflow3 ~/aflow3.tar.xz ~/aflow_test.tar.xz ~/aflow3
	tar Jcfvp ~/aflow3.tar.xz Makefile *cpp *.h *.js READ* APL ANRL* AUROSTD SQLITE
# no EXTRA IN SRC	if [ -d $(EXTRA_FREE) ]; then tar Jcfvp ~/aflow3.tar.xz Makefile *cpp *.h *.js READ* APL ANRL* AUROSTD SQLITE $(EXTRA_FREE) ; fi
# no other SOURCES in SRC	if [ -d SOURCES ]; then tar Jcfvp ~/aflow3.tar.xz Makefile *cpp *.h *.js READ* APL ANRL* AUROSTD SQLITE $(EXTRA_FREE) SOURCES ; fi
	if [ -d /www/AFLOW ]; then cp README_AFLOW.TXT /www/AFLOW ; fi
	if [ -d /www/AFLOW ]; then cp README_AFLOW_AEL.TXT /www/AFLOW ; fi
	if [ -d /www/AFLOW ]; then cp README_AFLOW_AGL.TXT /www/AFLOW ; fi
	if [ -d /www/AFLOW ]; then cp README_AFLOW_ANRL.TXT /www/AFLOW ; fi
	if [ -d /www/AFLOW ]; then cp README_AFLOW_APENNSY.TXT /www/AFLOW ; fi
	if [ -d /www/AFLOW ]; then cp README_AFLOW_APL.TXT /www/AFLOW ; fi
	if [ -d /www/AFLOW ]; then cp README_AFLOW_QHA_SCQHA_QHA3P.TXT /www/AFLOW ; fi
	if [ -d /www/AFLOW ]; then cp README_AFLOW_ACONVASP.TXT /www/AFLOW ; fi
	if [ -d /www/AFLOW ]; then cp README_AFLOW_COMPARE.TXT /www/AFLOW ; fi
	if [ -d /www/AFLOW ]; then cp README_AFLOW_FROZSL.TXT /www/AFLOW ; fi
	if [ -d /www/AFLOW ]; then cp README_AFLOW_LICENSE_GPL3.TXT /www/AFLOW ; fi
	if [ -d /www/AFLOW ]; then cp README_AFLOW_POCC.TXT /www/AFLOW ; fi
	if [ -d /www/AFLOW ]; then cp README_AFLOW_SCRIPTING.TXT /www/AFLOW ; fi
	if [ -d /www/AFLOW ]; then cp README_AFLOW_SYM.TXT /www/AFLOW ; fi
	if [ -d /www/AFLOW ]; then cp README_AFLOW_CHULL.TXT /www/AFLOW ; fi
	if [ -d /www/AFLOW ]; then cp README_AFLOW_EXCEPTIONS.TXT /www/AFLOW ; fi
	if [ -d /www/AFLOW ]; then cp README_AFLOW_VERSIONS_HISTORY.TXT /www/AFLOW ; fi
	if [ -d /www/AFLOW ]; then cp README_AFLOW_XAFLOW.TXT /www/AFLOW ; fi
	if [ -d /www/AFLOW ]; then cp README_AFLOW_AFLOWRC.TXT /www/AFLOW ; fi
	if [ -d /www/AFLOW ]; then cp README_CONTRIBS.TXT /www/AFLOW ; fi
	if [ -d /www/AFLOW ]; then cp README_PROTO.TXT /www/AFLOW ; fi
	if [ -f ../AFLOW3_NONFREE/EXTRA/FROZSL/Makefile ]; then $(MAKE) -C ../AFLOW3_NONFREE/EXTRA/FROZSL; fi
	if [ -f ../AFLOW3_NONFREE/EXTRA/FROZSL/frozsl ]; then cp ../AFLOW3_NONFREE/EXTRA/FROZSL/frozsl /www/AFLOW ; fi
	if [ -f ../AFLOW3_NONFREE/EXTRA/FROZSL/frozsl_init ]; then cp ../AFLOW3_NONFREE/EXTRA/FROZSL/frozsl_init /www/AFLOW ; fi
	if [ -f ../AFLOW3_NONFREE/EXTRA/FINDSYM/Makefile ]; then $(MAKE) -C ../AFLOW3_NONFREE/EXTRA/FINDSYM; fi
	if [ -f ../AFLOW3_NONFREE/EXTRA/FINDSYM/findsym ]; then cp ../AFLOW3_NONFREE/EXTRA/FINDSYM/findsym /www/AFLOW ; fi
	if [ -d /www/auro ]; then mv -f  ~/aflow3.tar.xz /www/auro ; fi
wget:	realclean
	rm -rf aflow3.tar.xz* aflow_test.tar.xz*
	wget http://materials.duke.edu/auro/aflow3.tar.xz
	tar xfvp aflow3.tar.xz
	rm -rf aflow3.tar.xz 
distroget: wget
	rm -rf TESTS PLATON* FINDSYM* FROZSL*
	cat Makefile | sed "s/Makefile /Makefile /g" | sed "s/\*\.a//g" > M
	mv M Makefile
patch_makefile1: Makefile
	cat Makefile | sed "s/ \-Wno\-unused\-but\-set\-variable//g" > M	
	mv M Makefile
patch_makefile2: Makefile
	cat Makefile | sed "s/\-O1//g" | sed "s/\-O2//g" | sed "s/\-O3//g" > M	
	mv M Makefile
check:
	-./aflow --proto=4:Bi:Tc | ./aflow --sprim | md5sum
	#174c76b2b2928dcdf2f3b39069a8b59
#	./aflow --checki | grep FAIL
www:	web
clean:
	rm -rf *~ *.o *.o-* *.rpo AUROSTD/*~ AUROSTD/*.o AUROSTD/*.rpo ./APL/*~ ./APL/*o ./APL/*a ./ANRL/*~ ./ANRL/*o ./ANRL/*a ./ANRL*/*~ ./ANRL*/*o ./ANRL*/*a ./SQLITE/*o aflow aflow_data aflowd aflow2 aflowd2 aflow1 aflowd1 aflow.exe aconvasp apennsy aflow_convasp* *.eps *.pdf *ps *aux *dvi *log *m *tex nohup* core* VASP_TEST* GRND_TEST* xa xb xc xd xe xf aflow.out* core* *~ \#* *.o aflow_pennsy *.eps *.pdf *ps *aux *dvi *log *m *tex *exe aurostd_xtensor_extra.cpp aflow_data_*cpp xaflowauto boo* xlist* _xdat `find . -name *~`
	rm -rfv `find . -name "LOCK*" | grep -v TESTS | grep -v SOURCES` *.out aflow.rasmol.xyz* aflow_potcar.cpp aflow_xpotcar.cpp aflow3.tar.xz ./AFLOWDATA ./LIBRARY* *html *toc *m *.jpg *.gif *.png
	rm -rfv ./APL/*~ ./ANRL/*~ ./ANRL*/*~ ./SQLITE/*~ ../AFLOW3_NONFREE/BIN/* ./EXTRA/GUS/*~ ./EXTRA/GUS/*.o 
realclean: clean
	rm -rf aurostd.cpp aurostd.h aflowlib_webapp_entry.cpp aflowlib_webapp_bands.cpp aflow_data_readme.cpp aflow_readme_aflow.cpp aflow_readme_aconvasp.cpp aflow_readme_apennsy.cpp aflow_readme_proto.cpp aflow_xproto_icsd_lib.cpp aflow_matlab_funcs.cpp aflow_gnuplot_funcs.cpp aflow_data_*cpp *.dat *conflic* 
	if [ -f $(EXTRA_FREE)/GUS/Makefile ]; then $(MAKE) -C $(EXTRA_FREE)/GUS realclean; fi
	if [ -f ../AFLOW3_NONFREE/EXTRA/FINDSYM/Makefile ]; then $(MAKE) -C ../AFLOW3_NONFREE/EXTRA/FINDSYM realclean; fi
	if [ -f ../AFLOW3_NONFREE/EXTRA/FROZSL/Makefile ]; then $(MAKE) -C ../AFLOW3_NONFREE/EXTRA/FROZSL clean; fi
	if [ -f ../AFLOW3_NONFREE/EXTRA/FROZSL/Makefile ]; then rm -fv ../AFLOW3_NONFREE/EXTRA/FROZSL/*conflic*; fi
	if [ -f $(EXTRA_FREE)/APL2AGR/Makefile ]; then $(MAKE) -C $(EXTRA_FREE)/APL2AGR realclean; fi
# ****************************************************************************
# icsd
icsdclean: clean realclean
	echo "dummy" > $(ICSD_LIB)//dummy.txt
	rm $(ICSD_LIB)//*
icsd:	$(ICSD_LIB)/README_LIBRARY_ICSD1.TXT $(ICSD_LIB)/README_LIBRARY_ICSD2.TXT $(ICSD_LIB)/README_LIBRARY_ICSD3.TXT $(ICSD_LIB)/README_LIBRARY_ICSD4.TXT $(ICSD_LIB)/README_LIBRARY_ICSD5.TXT $(ICSD_LIB)/README_LIBRARY_ICSD6.TXT $(ICSD_LIB)/README_LIBRARY_ICSD7.TXT $(ICSD_LIB)/README_LIBRARY_ICSD8.TXT
icsd_test: aflow
	cat $(ICSD_LIB)/README_LIBRARY_ICSD2.TXT | grep PROTOTYPE | sed "s/\[/\n/g" | grep PROTOTYPE | sed "s/PROTOTYPE/\.\/aflow --proto_icsd /g" > xgo
	sh ./xgo > /dev/null
# ****************************************************************************
# PENNSY PENNSY PENNSY PENNSY PENNSY PENNSY PENNSY PENNSY PENNSY PENNSY PENNSY
energy:	all
	./aflow --energy > apennsy.tex
	latex apennsy.tex | grep -v Overfull | grep -v cmtt | grep -v type | grep -v "\[\]" | grep "\["
psenergy: all
	./aflow --psenergy > psapennsy.tex
	latex psapennsy.tex
	dvips psapennsy.dvi -o psapennsy.eps
	ps2pdf psapennsy.eps
	cp psapennsy*ps EPS
	cp psapennsy*pdf EPS
MIX:	aflow
	./aflow --libraryX --mix > aflow_nomix.$(TIME).cpp
	chmod 755 aflow_nomix.$(TIME).cpp
	rm -rf aflow_nomix.cpp
	ln -sf aflow_nomix.$(TIME).cpp aflow_nomix.cpp
	grep -c UNKNOWN aflow_nomix.2*
mix:	MIX
NOMIX:	aflow
	cat aflow_nomix.cpp | grep NOMIX | sed "s/\/\//\n/g" | grep calcs | sed "s/\=/\n/g" | grep calcs | sed "s/ calcs/xxxx/g" | sed "s/ /cd \.\//g" | sed "s/xxxx/ \&\& aflow --multi \&\& cd \.\.\//g" > xnomix
	scp xnomix fslcollab8@mary:LIBS/LIB2/LIB
	rm -rf xnomix
nomix: NOMIX
EXPERIMENTS:aflow
	./aflow --libraryX --experiments > aflow_mix_experiments.$(TIME).cpp
experiments:EXPERIMENTS
TABLE:aflow
	./aflow --libraryX --table > aflow_table.$(TIME).tex
table:TABLE
xaflowauto: aflow
	./aflow --libraryX --xaflowauto > xaflowauto
	scp xaflowauto fslcollab8@marylou5.byu.edu:LIBS
boo: aflow
	rm -rf boo
	./aflow --calculated=icsd --remove | sed "s/\.\//\.\/ICSD\//g" >> boo
	./aflow --calculated=lib1 --remove | sed "s/\.\//\.\/LIB1\/LIB\//g" >> boo
	./aflow --calculated=lib2 --remove | sed "s/\.\//\.\/LIB2\/LIB\//g" >> boo
	./aflow --calculated=lib3 --remove | sed "s/\.\//\.\/LIB3\/LIB\//g" >> boo
	./aflow --calculated=lib4 --remove | sed "s/\.\//\.\/LIB4\/LIB\//g" >> boo
	./aflow --calculated=lib5 --remove | sed "s/\.\//\.\/LIB5\/LIB\//g" >> boo
	./aflow --calculated=lib6 --remove | sed "s/\.\//\.\/LIB6\/LIB\//g" >> boo
	./aflow --calculated=lib7 --remove | sed "s/\.\//\.\/LIB7\/LIB\//g" >> boo
	./aflow --calculated=lib8 --remove | sed "s/\.\//\.\/LIB8\/LIB\//g" >> boo
	./aflow --calculated=lib9 --remove | sed "s/\.\//\.\/LIB9\/LIB\//g" >> boo
	#scp boo aflow@lonsdale.tchpc.tcd.ie:LIB3/
	scp boo fslcollab8@mary:LIBS/
WEB:	aflow
	rm -rf *eps *pdf apennsy.m *html
	#	./aflow --libraryX --QUIET --print=pdf --hull > apennsy.m
	./aflow --libraryX --QUIET --print=gif --print=jpg --hull --oss=cout > apennsy.m
	./aflow --libraryX --QUIET --print=html > apennsy.html
	/usr/local/bin/matlab -nodesktop -r apennsy
	mv *jpg *gif /www/AFLOW
	rm -rf *eps
APOOL_OLD: aflow aflow_apennsy_main.cpp .cpp
	rm -rf /www/php/apool.html
	cat /www/php/apool_pre.html >> /www/php/apool.html
#	./aflow --apool_test >> /www/php/apool.html
	./aflow --apool >> /www/php/apool.html
	cat /www/php/apool_post.html >> /www/php/apool.html	
	./aflow --apool | grep available | sed "s/ checked//g" | sed "s/></XXXX\n/g" | grep radio | sed "s/<input type=\"radio\" name=\"alloy\" value=//g" | sed "s/XXXX/,/g" > /home/auro/work/AWRAPPER/awrapper_tables_public.cpp
	#
	#
#	rm -rf /www/php/apool2.html
#	cat /www/php/apool_pre.html >> /www/php/apool2.html
#	./aflow --apool_private >> /www/php/apool2.html
#	cat /www/php/apool_post.html >> /www/php/apool2.html	
#	#
shull:	all
	./aflow --library2 --shull --oss=cout > apennsys.m
vasp:	all
	./aflow --VASPIN | sed "s/\t/ /g" | sed "s/  / /g" | sed "s/  / /g" > structures.tex
