// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *                                                                         *
// ***************************************************************************
// Written by SC 2017

#include "aflow.h"

// ***************************************************************************
#ifdef _AFLOW_AFLOWRC_H_

#define AFLOWRC_FILENAME_LOCAL   XHOST.Home+"/.aflow.rc"
#define AFLOWRC_FILENAME_GLOBAL  "/etc/aflow.conf"

// DEFAULT DEFINITIONS

// GENERIC MPI   // DONE
#define AFLOWRC_MPI_START_DEFAULT                     string("")
#define         MPI_START_DEFAULT                     XHOST.adefault.getattachedscheme("MPI_START_DEFAULT")
#define AFLOWRC_MPI_STOP_DEFAULT                      string("")
#define         MPI_STOP_DEFAULT                      XHOST.adefault.getattachedscheme("MPI_STOP_DEFAULT")
#define AFLOWRC_MPI_COMMAND_DEFAULT                   string("mpirun -np")
#define         MPI_COMMAND_DEFAULT                   XHOST.adefault.getattachedscheme("MPI_COMMAND_DEFAULT")
#define AFLOWRC_MPI_NCPUS_DEFAULT                     4
#define         MPI_NCPUS_DEFAULT                     XHOST.adefault.getattachedutype<int>("MPI_NCPUS_DEFAULT")
#define AFLOWRC_MPI_NCPUS_MAX                         init::GetCPUCores() //16  //CO 180124
#define         MPI_NCPUS_MAX                         XHOST.adefault.getattachedutype<int>("MPI_NCPUS_MAX")

// BINARY    // DONE
#define AFLOWRC_DEFAULT_VASP_GAMMA_BIN                string("vasp46s_gamma")
#define         DEFAULT_VASP_GAMMA_BIN                XHOST.adefault.getattachedscheme("DEFAULT_VASP_GAMMA_BIN")
#define AFLOWRC_DEFAULT_VASP_GAMMA_MPI_BIN            string("mpivasp46s_gamma")
#define         DEFAULT_VASP_GAMMA_MPI_BIN            XHOST.adefault.getattachedscheme("DEFAULT_VASP_GAMMA_MPI_BIN")
#define AFLOWRC_DEFAULT_VASP_BIN                      string("vasp46s")
#define         DEFAULT_VASP_BIN                      XHOST.adefault.getattachedscheme("DEFAULT_VASP_BIN")
#define AFLOWRC_DEFAULT_VASP_MPI_BIN                  string("mpivasp46s")
#define         DEFAULT_VASP_MPI_BIN                  XHOST.adefault.getattachedscheme("DEFAULT_VASP_MPI_BIN")
#define AFLOWRC_DEFAULT_VASP5_BIN                     string("vasp53s")
#define         DEFAULT_VASP5_BIN                     XHOST.adefault.getattachedscheme("DEFAULT_VASP5_BIN")
#define AFLOWRC_DEFAULT_VASP5_MPI_BIN                 string("mpivasp53s")
#define         DEFAULT_VASP5_MPI_BIN                 XHOST.adefault.getattachedscheme("DEFAULT_VASP5_MPI_BIN")
//aims
#define AFLOWRC_DEFAULT_AIMS_BIN                      string("aims")
#define         DEFAULT_AIMS_BIN                      XHOST.adefault.getattachedscheme("DEFAULT_AIMS_BIN")

// POTCARS // DONE
#define AFLOWRC_DEFAULT_VASP_POTCAR_DIRECTORIES       string("/common/VASP,/common/AFLOW/VASP,/home/aflow/common/AFLOW/VASP,/fslhome/fslcollab8/group/VASP,/fslhome/glh43/src/,/share/home/00470/tg457283/common/AFLOW/VASP/,/share/home/00457/tg457357/common/AFLOW/VASP/,/home/mehl/bin/AFLOW/VASP/,~/common/VASP/,~/common/AFLOW/VASP/,/nics/a/proj/aflow/common/AFLOW/VASP/,/home/users/aflow/common/VASP,/share/apps/AFLOW3/VASP,/projects/kyang-group/common/VASP,/somewhere/")  // first is default, tokenized with ","
#define         DEFAULT_VASP_POTCAR_DIRECTORIES       XHOST.adefault.getattachedscheme("DEFAULT_VASP_POTCAR_DIRECTORIES")
#define AFLOWRC_DEFAULT_VASP_POTCAR_DATE              string("current")
#define         DEFAULT_VASP_POTCAR_DATE              XHOST.adefault.getattachedscheme("DEFAULT_VASP_POTCAR_DATE")
#define AFLOWRC_DEFAULT_VASP_POTCAR_SUFFIX            string("/POTCAR")
#define         DEFAULT_VASP_POTCAR_SUFFIX            XHOST.adefault.getattachedscheme("DEFAULT_VASP_POTCAR_SUFFIX")
#define AFLOWRC_DEFAULT_VASP_POTCAR_DATE_POT_LDA      string("01Apr2000")   // when no date is given for pot_LDA
#define         DEFAULT_VASP_POTCAR_DATE_POT_LDA      XHOST.adefault.getattachedscheme("DEFAULT_VASP_POTCAR_DATE_POT_LDA")
#define AFLOWRC_DEFAULT_VASP_POTCAR_DATE_POT_GGA      string("01Apr2000")   // when no date is given for pot_GGA
#define         DEFAULT_VASP_POTCAR_DATE_POT_GGA      XHOST.adefault.getattachedscheme("DEFAULT_VASP_POTCAR_DATE_POT_GGA")


// KPOINTS/DOS // DONE
#define AFLOWRC_DEFAULT_BANDS_GRID                    20
#define         DEFAULT_BANDS_GRID                    XHOST.adefault.getattachedutype<int>("DEFAULT_BANDS_GRID") 
#define AFLOWRC_DEFAULT_BANDS_LATTICE                 string("AUTO")
#define         DEFAULT_BANDS_LATTICE                 XHOST.adefault.getattachedscheme("DEFAULT_BANDS_LATTICE")
#define AFLOWRC_DEFAULT_KSCHEME                       string("M")
#define         DEFAULT_KSCHEME                       XHOST.adefault.getattachedscheme("DEFAULT_KSCHEME")
#define AFLOWRC_DEFAULT_KPPRA                         6000
#define         DEFAULT_KPPRA                         XHOST.adefault.getattachedutype<int>("DEFAULT_KPPRA")
#define AFLOWRC_DEFAULT_KPPRA_STATIC                  10000
#define         DEFAULT_KPPRA_STATIC                  XHOST.adefault.getattachedutype<int>("DEFAULT_KPPRA_STATIC")
#define AFLOWRC_DEFAULT_STATIC_KSCHEME                string("M") // WAHYU DEFAULT
#define         DEFAULT_STATIC_KSCHEME                XHOST.adefault.getattachedscheme("DEFAULT_STATIC_KSCHEME")
#define AFLOWRC_DEFAULT_KPPRA_ICSD                    8000
#define         DEFAULT_KPPRA_ICSD                    XHOST.adefault.getattachedutype<int>("DEFAULT_KPPRA_ICSD")
#define AFLOWRC_DEFAULT_UNARY_BANDS_GRID              128
#define         DEFAULT_UNARY_BANDS_GRID              XHOST.adefault.getattachedutype<int>("DEFAULT_UNARY_BANDS_GRID")
#define AFLOWRC_DEFAULT_UNARY_KPPRA                   8000 // 32768 // 27000
#define         DEFAULT_UNARY_KPPRA                   XHOST.adefault.getattachedutype<int>("DEFAULT_UNARY_KPPRA")
#define AFLOWRC_DEFAULT_UNARY_KPPRA_STATIC            8000 // 32768 // 27000
#define         DEFAULT_UNARY_KPPRA_STATIC            XHOST.adefault.getattachedutype<int>("DEFAULT_UNARY_KPPRA_STATIC")
#define AFLOWRC_DEFAULT_PHONONS_KSCHEME               string("G")
#define         DEFAULT_PHONONS_KSCHEME               XHOST.adefault.getattachedscheme("DEFAULT_PHONONS_KSCHEME")
#define AFLOWRC_DEFAULT_DOS_EMIN                      -10.0
#define         DEFAULT_DOS_EMIN                      XHOST.adefault.getattachedutype<double>("DEFAULT_DOS_EMIN")
#define AFLOWRC_DEFAULT_DOS_EMAX                      10.0
#define         DEFAULT_DOS_EMAX                      XHOST.adefault.getattachedutype<double>("DEFAULT_DOS_EMAX")
#define AFLOWRC_DEFAULT_DOS_SCALE                     1.2
#define         DEFAULT_DOS_SCALE                     XHOST.adefault.getattachedutype<double>("DEFAULT_DOS_SCALE")

// PRECISION // DONE
#define AFLOWRC_DEFAULT_VASP_PREC_ENMAX_LOW                   1.0
#define         DEFAULT_VASP_PREC_ENMAX_LOW                   XHOST.adefault.getattachedutype<double>("DEFAULT_VASP_PREC_ENMAX_LOW")
#define AFLOWRC_DEFAULT_VASP_PREC_ENMAX_MEDIUM                1.3
#define         DEFAULT_VASP_PREC_ENMAX_MEDIUM                XHOST.adefault.getattachedutype<double>("DEFAULT_VASP_PREC_ENMAX_MEDIUM")
#define AFLOWRC_DEFAULT_VASP_PREC_ENMAX_NORMAL                1.3
#define         DEFAULT_VASP_PREC_ENMAX_NORMAL                XHOST.adefault.getattachedutype<double>("DEFAULT_VASP_PREC_ENMAX_NORMAL")
#define AFLOWRC_DEFAULT_VASP_PREC_ENMAX_HIGH                  1.4
#define         DEFAULT_VASP_PREC_ENMAX_HIGH                  XHOST.adefault.getattachedutype<double>("DEFAULT_VASP_PREC_ENMAX_HIGH")
#define AFLOWRC_DEFAULT_VASP_PREC_ENMAX_ACCURATE              1.4
#define         DEFAULT_VASP_PREC_ENMAX_ACCURATE              XHOST.adefault.getattachedutype<double>("DEFAULT_VASP_PREC_ENMAX_ACCURATE")
#define AFLOWRC_DEFAULT_VASP_SPIN_REMOVE_CUTOFF               0.05
#define         DEFAULT_VASP_SPIN_REMOVE_CUTOFF               XHOST.adefault.getattachedutype<double>("DEFAULT_VASP_SPIN_REMOVE_CUTOFF")
#define AFLOWRC_DEFAULT_VASP_PREC_POTIM                       0.5
#define         DEFAULT_VASP_PREC_POTIM                       XHOST.adefault.getattachedutype<double>("DEFAULT_VASP_PREC_POTIM")
#define AFLOWRC_DEFAULT_VASP_PREC_EDIFFG                      -1E-3
#define         DEFAULT_VASP_PREC_EDIFFG                      XHOST.adefault.getattachedutype<double>("DEFAULT_VASP_PREC_EDIFFG")

// OPTIONS // DONE
#define AFLOWRC_DEFAULT_VASP_EXTERNAL_INCAR                   string("./INCAR")
#define         DEFAULT_VASP_EXTERNAL_INCAR                   XHOST.adefault.getattachedscheme("DEFAULT_VASP_EXTERNAL_INCAR")
#define AFLOWRC_DEFAULT_AIMS_EXTERNAL_CONTROL                 string("./control.in")
#define         DEFAULT_AIMS_EXTERNAL_CONTROL                 XHOST.adefault.getattachedscheme("DEFAULT_AIMS_EXTERNAL_CONTROL")
#define AFLOWRC_DEFAULT_VASP_EXTERNAL_POSCAR                  string("./POSCAR")
#define         DEFAULT_VASP_EXTERNAL_POSCAR                  XHOST.adefault.getattachedscheme("DEFAULT_VASP_EXTERNAL_POSCAR")
#define AFLOWRC_DEFAULT_AIMS_EXTERNAL_GEOM                    string("./geometry.in")
#define         DEFAULT_AIMS_EXTERNAL_GEOM                    XHOST.adefault.getattachedscheme("DEFAULT_AIMS_EXTERNAL_GEOM")
#define AFLOWRC_DEFAULT_VASP_EXTERNAL_POTCAR                  string("./POTCAR")
#define         DEFAULT_VASP_EXTERNAL_POTCAR                  XHOST.adefault.getattachedscheme("DEFAULT_VASP_EXTERNAL_POTCAR")
#define AFLOWRC_DEFAULT_VASP_EXTERNAL_KPOINTS                 string("./KPOINTS")
#define         DEFAULT_VASP_EXTERNAL_KPOINTS                 XHOST.adefault.getattachedscheme("DEFAULT_VASP_EXTERNAL_KPOINTS")
#define AFLOWRC_DEFAULT_VASP_PSEUDOPOTENTIAL_TYPE             string("potpaw_PBE")
#define         DEFAULT_VASP_PSEUDOPOTENTIAL_TYPE             XHOST.adefault.getattachedscheme("DEFAULT_VASP_PSEUDOPOTENTIAL_TYPE")
#define AFLOWRC_DEFAULT_VASP_FORCE_OPTION_RELAX_MODE_SCHEME   string("ENERGY")
#define         DEFAULT_VASP_FORCE_OPTION_RELAX_MODE_SCHEME   XHOST.adefault.getattachedscheme("DEFAULT_VASP_FORCE_OPTION_RELAX_MODE_SCHEME")
#define AFLOWRC_DEFAULT_VASP_FORCE_OPTION_PREC_SCHEME         string("ACCURATE")
#define         DEFAULT_VASP_FORCE_OPTION_PREC_SCHEME         XHOST.adefault.getattachedscheme("DEFAULT_VASP_FORCE_OPTION_PREC_SCHEME")
#define AFLOWRC_DEFAULT_VASP_FORCE_OPTION_ALGO_SCHEME         string("NORMAL")
#define         DEFAULT_VASP_FORCE_OPTION_ALGO_SCHEME         XHOST.adefault.getattachedscheme("DEFAULT_VASP_FORCE_OPTION_ALGO_SCHEME")
#define AFLOWRC_DEFAULT_VASP_FORCE_OPTION_METAGGA_SCHEME      string("NONE")
#define         DEFAULT_VASP_FORCE_OPTION_METAGGA_SCHEME      XHOST.adefault.getattachedscheme("DEFAULT_VASP_FORCE_OPTION_METAGGA_SCHEME")
#define AFLOWRC_DEFAULT_VASP_FORCE_OPTION_IVDW_SCHEME         string("0")  
#define         DEFAULT_VASP_FORCE_OPTION_IVDW_SCHEME         XHOST.adefault.getattachedscheme("DEFAULT_VASP_FORCE_OPTION_IVDW_SCHEME")
#define AFLOWRC_DEFAULT_VASP_FORCE_OPTION_TYPE_SCHEME         string("DEFAULT")
#define         DEFAULT_VASP_FORCE_OPTION_TYPE_SCHEME         XHOST.adefault.getattachedscheme("DEFAULT_VASP_FORCE_OPTION_TYPE_SCHEME")
#define AFLOWRC_DEFAULT_VASP_FORCE_OPTION_ABMIX_SCHEME        string("DEFAULT")
#define         DEFAULT_VASP_FORCE_OPTION_ABMIX_SCHEME        XHOST.adefault.getattachedscheme("DEFAULT_VASP_FORCE_OPTION_ABMIX_SCHEME")
#define AFLOWRC_DEFAULT_VASP_FORCE_OPTION_SYM                 TRUE
#define         DEFAULT_VASP_FORCE_OPTION_SYM                 XHOST.adefault.getattachedutype<bool>("DEFAULT_VASP_FORCE_OPTION_SYM") 
#define AFLOWRC_DEFAULT_VASP_FORCE_OPTION_SPIN                TRUE
#define         DEFAULT_VASP_FORCE_OPTION_SPIN                XHOST.adefault.getattachedutype<bool>("DEFAULT_VASP_FORCE_OPTION_SPIN") 
#define AFLOWRC_DEFAULT_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_1 FALSE
#define         DEFAULT_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_1 XHOST.adefault.getattachedutype<bool>("DEFAULT_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_1") 
#define AFLOWRC_DEFAULT_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_2 FALSE
#define         DEFAULT_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_2 XHOST.adefault.getattachedutype<bool>("DEFAULT_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_2") 
#define AFLOWRC_DEFAULT_VASP_FORCE_OPTION_BADER               FALSE
#define         DEFAULT_VASP_FORCE_OPTION_BADER               XHOST.adefault.getattachedutype<bool>("DEFAULT_VASP_FORCE_OPTION_BADER") 
#define AFLOWRC_DEFAULT_VASP_FORCE_OPTION_ELF                 FALSE
#define         DEFAULT_VASP_FORCE_OPTION_ELF                 XHOST.adefault.getattachedutype<bool>("DEFAULT_VASP_FORCE_OPTION_ELF") 
#define AFLOWRC_DEFAULT_VASP_FORCE_OPTION_AUTO_MAGMOM         FALSE   // TRUE
#define         DEFAULT_VASP_FORCE_OPTION_AUTO_MAGMOM         XHOST.adefault.getattachedutype<bool>("DEFAULT_VASP_FORCE_OPTION_AUTO_MAGMOM") 
#define AFLOWRC_DEFAULT_VASP_FORCE_OPTION_WAVECAR             FALSE
#define         DEFAULT_VASP_FORCE_OPTION_WAVECAR             XHOST.adefault.getattachedutype<bool>("DEFAULT_VASP_FORCE_OPTION_WAVECAR") 
#define AFLOWRC_DEFAULT_VASP_FORCE_OPTION_CHGCAR              TRUE
#define         DEFAULT_VASP_FORCE_OPTION_CHGCAR              XHOST.adefault.getattachedutype<bool>("DEFAULT_VASP_FORCE_OPTION_CHGCAR") 
#define AFLOWRC_DEFAULT_VASP_FORCE_OPTION_LSCOUPLING          FALSE
#define         DEFAULT_VASP_FORCE_OPTION_LSCOUPLING          XHOST.adefault.getattachedutype<bool>("DEFAULT_VASP_FORCE_OPTION_LSCOUPLING") 

// AFLOW_LIBRARY AFLOW_PROJECT // DONE
#define AFLOWRC_DEFAULT_AFLOW_LIBRARY_DIRECTORIES             string("/common/AFLOW/LIBS/,/home/aflow/common/AFLOW/LIBS/,/fslhome/glh43/src/,/usr/local/bin/,/fslhome/fslcollab8/group/bin/,/home/auro/work/AFLOW3/,~/common/AFLOW/LIBS/,./,/nics/a/proj/aflow/common/AFLOW/LIBS/,/home/users/aflow/common/AFLOW/LIBS,/home/junkai/PROTO_DATABASE/,/projects/kyang-group/common/LIBS,/somewhere/")  // first is default, tokenized with ","
#define         DEFAULT_AFLOW_LIBRARY_DIRECTORIES             XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_LIBRARY_DIRECTORIES")
#define AFLOWRC_DEFAULT_AFLOW_PROJECTS_DIRECTORIES            string("/common/ICSD,/common/LIB1,/common/LIB2,/common/LIB3,/common/LIB4,/common/LIB5,/common/LIB6,/common/LIB7,/common/LIB8,/common/LIB9")  // first is default, tokenized with ","
#define         DEFAULT_AFLOW_PROJECTS_DIRECTORIES            XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_PROJECTS_DIRECTORIES")

// PLATON/FINDSYM // DONE
#define AFLOWRC_DEFAULT_PLATON_P_EQUAL                        FALSE
#define         DEFAULT_PLATON_P_EQUAL                        XHOST.adefault.getattachedutype<bool>("DEFAULT_PLATON_P_EQUAL") 
#define AFLOWRC_DEFAULT_PLATON_P_EXACT                        FALSE
#define         DEFAULT_PLATON_P_EXACT                        XHOST.adefault.getattachedutype<bool>("DEFAULT_PLATON_P_EXACT") 
#define AFLOWRC_DEFAULT_PLATON_P_ANG                          1.0
#define         DEFAULT_PLATON_P_ANG                          XHOST.adefault.getattachedutype<double>("DEFAULT_PLATON_P_ANG") 
#define AFLOWRC_DEFAULT_PLATON_P_D1                           0.25
#define         DEFAULT_PLATON_P_D1                           XHOST.adefault.getattachedutype<double>("DEFAULT_PLATON_P_D1") 
#define AFLOWRC_DEFAULT_PLATON_P_D2                           0.25
#define         DEFAULT_PLATON_P_D2                           XHOST.adefault.getattachedutype<double>("DEFAULT_PLATON_P_D2") 
#define AFLOWRC_DEFAULT_PLATON_P_D3                           0.25
#define         DEFAULT_PLATON_P_D3                           XHOST.adefault.getattachedutype<double>("DEFAULT_PLATON_P_D3") 
#define AFLOWRC_DEFAULT_FINDSYM_TOL                           1.0e-3
#define         DEFAULT_FINDSYM_TOL                           XHOST.adefault.getattachedutype<double>("DEFAULT_FINDSYM_TOL") 

// GNUPLOT // DONE
#define AFLOWRC_DEFAULT_GNUPLOT_EPS_FONT                      string("Helvetica")              // we do not have time to waste with ttf
#define         DEFAULT_GNUPLOT_EPS_FONT                      XHOST.adefault.getattachedscheme("DEFAULT_GNUPLOT_EPS_FONT")
#define AFLOWRC_DEFAULT_GNUPLOT_EPS_FONT_BOLD                 string("Helvetica-Bold")         // we do not have time to waste with ttf
#define         DEFAULT_GNUPLOT_EPS_FONT_BOLD                 XHOST.adefault.getattachedscheme("DEFAULT_GNUPLOT_EPS_FONT_BOLD")
#define AFLOWRC_DEFAULT_GNUPLOT_EPS_FONT_ITALICS              string("Helvetica-Oblique")      // we do not have time to waste with ttf
#define         DEFAULT_GNUPLOT_EPS_FONT_ITALICS              XHOST.adefault.getattachedscheme("DEFAULT_GNUPLOT_EPS_FONT_ITALICS")
#define AFLOWRC_DEFAULT_GNUPLOT_EPS_FONT_BOLD_ITALICS         string("Helvetica-BoldOblique")  // we do not have time to waste with ttf
#define         DEFAULT_GNUPLOT_EPS_FONT_BOLD_ITALICS         XHOST.adefault.getattachedscheme("DEFAULT_GNUPLOT_EPS_FONT_BOLD_ITALICS")
#define AFLOWRC_DEFAULT_GNUPLOT_PNG_FONT                      string("Arial")                  // we do not have time to waste with ttf
#define         DEFAULT_GNUPLOT_PNG_FONT                      XHOST.adefault.getattachedscheme("DEFAULT_GNUPLOT_PNG_FONT")
#define AFLOWRC_DEFAULT_GNUPLOT_PNG_FONT_BOLD                 string("Arial_Bold")             // we do not have time to waste with ttf
#define         DEFAULT_GNUPLOT_PNG_FONT_BOLD                 XHOST.adefault.getattachedscheme("DEFAULT_GNUPLOT_PNG_FONT_BOLD")
#define AFLOWRC_DEFAULT_GNUPLOT_PNG_FONT_ITALICS              string("Arial_Italic")           // we do not have time to waste with ttf
#define         DEFAULT_GNUPLOT_PNG_FONT_ITALICS              XHOST.adefault.getattachedscheme("DEFAULT_GNUPLOT_PNG_FONT_ITALICS")
#define AFLOWRC_DEFAULT_GNUPLOT_PNG_FONT_BOLD_ITALICS         string("Arial_BoldItalic")       // we do not have time to waste with ttf
#define         DEFAULT_GNUPLOT_PNG_FONT_BOLD_ITALICS         XHOST.adefault.getattachedscheme("DEFAULT_GNUPLOT_PNG_FONT_BOLD_ITALICS")

// DEFAULT CHULL
#define AFLOWRC_DEFAULT_CHULL_OUTLIER_MULTIPLIER                        3.25
#define         DEFAULT_CHULL_OUTLIER_MULTIPLIER                        XHOST.adefault.getattachedutype<double>("DEFAULT_CHULL_OUTLIER_MULTIPLIER")
#define AFLOWRC_DEFAULT_CHULL_LATEX_BANNER                              1
#define         DEFAULT_CHULL_LATEX_BANNER                              XHOST.adefault.getattachedutype<int>("DEFAULT_CHULL_LATEX_BANNER")
#define AFLOWRC_DEFAULT_CHULL_LATEX_COMPOUNDS_COLUMN                    FALSE
#define         DEFAULT_CHULL_LATEX_COMPOUNDS_COLUMN                    XHOST.adefault.getattachedutype<bool>("DEFAULT_CHULL_LATEX_COMPOUNDS_COLUMN")
#define AFLOWRC_DEFAULT_CHULL_LATEX_COMPOSITION_HEADER                  FALSE
#define         DEFAULT_CHULL_LATEX_COMPOSITION_HEADER                  XHOST.adefault.getattachedutype<bool>("DEFAULT_CHULL_LATEX_COMPOSITION_HEADER")
#define AFLOWRC_DEFAULT_CHULL_LATEX_PLOT_ENDPOINTS                      FALSE
#define         DEFAULT_CHULL_LATEX_PLOT_ENDPOINTS                      XHOST.adefault.getattachedutype<bool>("DEFAULT_CHULL_LATEX_PLOT_ENDPOINTS")
#define AFLOWRC_DEFAULT_CHULL_LATEX_PLOT_OFF_HULL                       -1
#define         DEFAULT_CHULL_LATEX_PLOT_OFF_HULL                       XHOST.adefault.getattachedutype<int>("DEFAULT_CHULL_LATEX_PLOT_OFF_HULL")
#define AFLOWRC_DEFAULT_CHULL_LATEX_PLOT_UNSTABLE                       FALSE
#define         DEFAULT_CHULL_LATEX_PLOT_UNSTABLE                       XHOST.adefault.getattachedutype<bool>("DEFAULT_CHULL_LATEX_PLOT_UNSTABLE")
#define AFLOWRC_DEFAULT_CHULL_LATEX_FILTER_SCHEME                       string("")
#define         DEFAULT_CHULL_LATEX_FILTER_SCHEME                       XHOST.adefault.getattachedscheme("DEFAULT_CHULL_LATEX_FILTER_SCHEME")
#define AFLOWRC_DEFAULT_CHULL_LATEX_FILTER_VALUE                        0.0
#define         DEFAULT_CHULL_LATEX_FILTER_VALUE                        XHOST.adefault.getattachedutype<double>("DEFAULT_CHULL_LATEX_FILTER_VALUE")
#define AFLOWRC_DEFAULT_CHULL_LATEX_COLOR_BAR                           TRUE
#define         DEFAULT_CHULL_LATEX_COLOR_BAR                           XHOST.adefault.getattachedutype<bool>("DEFAULT_CHULL_LATEX_COLOR_BAR")
#define AFLOWRC_DEFAULT_CHULL_LATEX_HEAT_MAP                            TRUE
#define         DEFAULT_CHULL_LATEX_HEAT_MAP                            XHOST.adefault.getattachedutype<bool>("DEFAULT_CHULL_LATEX_HEAT_MAP")
#define AFLOWRC_DEFAULT_CHULL_LATEX_COLOR_GRADIENT                      TRUE
#define         DEFAULT_CHULL_LATEX_COLOR_GRADIENT                      XHOST.adefault.getattachedutype<bool>("DEFAULT_CHULL_LATEX_COLOR_GRADIENT")
#define AFLOWRC_DEFAULT_CHULL_LATEX_COLOR_MAP                           string("")
#define         DEFAULT_CHULL_LATEX_COLOR_MAP                           XHOST.adefault.getattachedscheme("DEFAULT_CHULL_LATEX_COLOR_MAP")
#define AFLOWRC_DEFAULT_CHULL_LATEX_TERNARY_LABEL_COLOR                 string("")
#define         DEFAULT_CHULL_LATEX_TERNARY_LABEL_COLOR                 XHOST.adefault.getattachedscheme("DEFAULT_CHULL_LATEX_TERNARY_LABEL_COLOR")
#define AFLOWRC_DEFAULT_CHULL_LATEX_REVERSE_AXIS                        FALSE
#define         DEFAULT_CHULL_LATEX_REVERSE_AXIS                        XHOST.adefault.getattachedutype<bool>("DEFAULT_CHULL_LATEX_REVERSE_AXIS")
#define AFLOWRC_DEFAULT_CHULL_LATEX_FACET_LINE_DROP_SHADOW              FALSE
#define         DEFAULT_CHULL_LATEX_FACET_LINE_DROP_SHADOW              XHOST.adefault.getattachedutype<bool>("DEFAULT_CHULL_LATEX_FACET_LINE_DROP_SHADOW")
#define AFLOWRC_DEFAULT_CHULL_LATEX_LINKS                               1
#define         DEFAULT_CHULL_LATEX_LINKS                               XHOST.adefault.getattachedutype<int>("DEFAULT_CHULL_LATEX_LINKS")
#define AFLOWRC_DEFAULT_CHULL_LATEX_LABEL_NAME                          string("")
#define         DEFAULT_CHULL_LATEX_LABEL_NAME                          XHOST.adefault.getattachedscheme("DEFAULT_CHULL_LATEX_LABEL_NAME")
#define AFLOWRC_DEFAULT_CHULL_LATEX_META_LABELS                         FALSE
#define         DEFAULT_CHULL_LATEX_META_LABELS                         XHOST.adefault.getattachedutype<bool>("DEFAULT_CHULL_LATEX_META_LABELS")
#define AFLOWRC_DEFAULT_CHULL_LATEX_LABELS_OFF_HULL                     FALSE
#define         DEFAULT_CHULL_LATEX_LABELS_OFF_HULL                     XHOST.adefault.getattachedutype<bool>("DEFAULT_CHULL_LATEX_LABELS_OFF_HULL")
#define AFLOWRC_DEFAULT_CHULL_LATEX_PLOT_REDUCE_COMPOSITION             -1
#define         DEFAULT_CHULL_LATEX_PLOT_REDUCE_COMPOSITION             XHOST.adefault.getattachedutype<int>("DEFAULT_CHULL_LATEX_PLOT_REDUCE_COMPOSITION")
#define AFLOWRC_DEFAULT_CHULL_LATEX_HELVETICA_FONT                      FALSE
#define         DEFAULT_CHULL_LATEX_HELVETICA_FONT                      XHOST.adefault.getattachedutype<bool>("DEFAULT_CHULL_LATEX_HELVETICA_FONT")
#define AFLOWRC_DEFAULT_CHULL_LATEX_FONT_SIZE                           string("")
#define         DEFAULT_CHULL_LATEX_FONT_SIZE                           XHOST.adefault.getattachedscheme("DEFAULT_CHULL_LATEX_FONT_SIZE")
#define AFLOWRC_DEFAULT_CHULL_LATEX_ROTATE_LABELS                       TRUE
#define         DEFAULT_CHULL_LATEX_ROTATE_LABELS                       XHOST.adefault.getattachedutype<bool>("DEFAULT_CHULL_LATEX_ROTATE_LABELS")
#define AFLOWRC_DEFAULT_CHULL_LATEX_BOLD_LABELS                         TRUE
#define         DEFAULT_CHULL_LATEX_BOLD_LABELS                         XHOST.adefault.getattachedutype<bool>("DEFAULT_CHULL_LATEX_BOLD_LABELS")

// CORES // DONE
#define AFLOWRC_AFLOW_CORE_TEMPERATURE_BEEP                   56.0    // Celsius
#define         AFLOW_CORE_TEMPERATURE_BEEP                   XHOST.adefault.getattachedutype<double>("AFLOW_CORE_TEMPERATURE_BEEP") 
#define AFLOWRC_AFLOW_CORE_TEMPERATURE_HALT                   65.0    // Celsius, you need to run aflow as root to shutdown
#define         AFLOW_CORE_TEMPERATURE_HALT                   XHOST.adefault.getattachedutype<double>("AFLOW_CORE_TEMPERATURE_HALT") 
#define AFLOWRC_AFLOW_CORE_TEMPERATURE_REFRESH                5.0    // seconds
#define         AFLOW_CORE_TEMPERATURE_REFRESH                XHOST.adefault.getattachedutype<double>("AFLOW_CORE_TEMPERATURE_REFRESH") 

// MACHINE DEPENDENT MPI
#define AFLOWRC_MPI_OPTIONS_DUKE_BETA_MPICH                   string("ulimit -s unlimited ") // DUKE_BETA_MPICH
#define         MPI_OPTIONS_DUKE_BETA_MPICH                   XHOST.adefault.getattachedscheme("MPI_OPTIONS_DUKE_BETA_MPICH")
#define AFLOWRC_MPI_COMMAND_DUKE_BETA_MPICH                   string("/usr/bin/mpiexec -np") // DUKE_BETA_MPICH
#define         MPI_COMMAND_DUKE_BETA_MPICH                   XHOST.adefault.getattachedscheme("MPI_COMMAND_DUKE_BETA_MPICH")
#define AFLOWRC_MPI_BINARY_DIR_DUKE_BETA_MPICH                string("/usr/local/bin/") // DUKE_BETA_MPICH
#define         MPI_BINARY_DIR_DUKE_BETA_MPICH                XHOST.adefault.getattachedscheme("MPI_BINARY_DIR_DUKE_BETA_MPICH")

#define AFLOWRC_MPI_OPTIONS_DUKE_BETA_OPENMPI                 string("ulimit -s unlimited ") // DUKE_BETA_OPENMPI
#define         MPI_OPTIONS_DUKE_BETA_OPENMPI                 XHOST.adefault.getattachedscheme("MPI_OPTIONS_DUKE_BETA_OPENMPI")
#define AFLOWRC_MPI_COMMAND_DUKE_BETA_OPENMPI                 string("/usr/bin/mpirun.openmpi -np") // DUKE_BETA_OPENMPI
#define         MPI_COMMAND_DUKE_BETA_OPENMPI                 XHOST.adefault.getattachedscheme("MPI_COMMAND_DUKE_BETA_OPENMPI")
#define AFLOWRC_MPI_BINARY_DIR_DUKE_BETA_OPENMPI              string("/usr/local/bin/") // DUKE_BETA_OPENMPI
#define         MPI_BINARY_DIR_DUKE_BETA_OPENMPI              XHOST.adefault.getattachedscheme("MPI_BINARY_DIR_DUKE_BETA_OPENMPI")

#define AFLOWRC_MPI_OPTIONS_DUKE_MATERIALS                    string("ulimit -s unlimited ") // DUKE_MATERIALS
#define         MPI_OPTIONS_DUKE_MATERIALS                    XHOST.adefault.getattachedscheme("MPI_OPTIONS_DUKE_MATERIALS")
#define AFLOWRC_MPI_COMMAND_DUKE_MATERIALS                    string("/usr/bin/mpiexec -np") // DUKE_MATERIALS
#define         MPI_COMMAND_DUKE_MATERIALS                    XHOST.adefault.getattachedscheme("MPI_COMMAND_DUKE_MATERIALS")
#define AFLOWRC_MPI_BINARY_DIR_DUKE_MATERIALS                 string("/usr/local/bin/")  // DUKE_MATERIALS
#define         MPI_BINARY_DIR_DUKE_MATERIALS                 XHOST.adefault.getattachedscheme("MPI_BINARY_DIR_DUKE_MATERIALS")

#define AFLOWRC_MPI_OPTIONS_DUKE_AFLOWLIB                     string("ulimit -s unlimited ") // DUKE_AFLOWLIB
#define         MPI_OPTIONS_DUKE_AFLOWLIB                     XHOST.adefault.getattachedscheme("MPI_OPTIONS_DUKE_AFLOWLIB")
#define AFLOWRC_MPI_COMMAND_DUKE_AFLOWLIB                     string("/usr/bin/mpiexec -np") // DUKE_AFLOWLIB
#define         MPI_COMMAND_DUKE_AFLOWLIB                     XHOST.adefault.getattachedscheme("MPI_COMMAND_DUKE_AFLOWLIB")
#define AFLOWRC_MPI_BINARY_DIR_DUKE_AFLOWLIB                  string("/usr/local/bin/") // DUKE_AFLOWLIB
#define         MPI_BINARY_DIR_DUKE_AFLOWLIB                  XHOST.adefault.getattachedscheme("MPI_BINARY_DIR_DUKE_AFLOWLIB")

#define AFLOWRC_MPI_OPTIONS_DUKE_HABANA                       string("ulimit -s unlimited ") // DUKE_HABANA
#define         MPI_OPTIONS_DUKE_HABANA                       XHOST.adefault.getattachedscheme("MPI_OPTIONS_DUKE_HABANA")
#define AFLOWRC_MPI_COMMAND_DUKE_HABANA                       string("/usr/bin/mpiexec -np") // DUKE_HABANA
#define         MPI_COMMAND_DUKE_HABANA                       XHOST.adefault.getattachedscheme("MPI_COMMAND_DUKE_HABANA")
#define AFLOWRC_MPI_BINARY_DIR_DUKE_HABANA                    string("/usr/local/bin/") // DUKE_HABANA
#define         MPI_BINARY_DIR_DUKE_HABANA                    XHOST.adefault.getattachedscheme("MPI_BINARY_DIR_DUKE_HABANA")

#define AFLOWRC_MPI_OPTIONS_DUKE_QRATS_MPICH                  string("ulimit -s unlimited ") // DUKE_QRATS_MPICH
#define         MPI_OPTIONS_DUKE_QRATS_MPICH                  XHOST.adefault.getattachedscheme("MPI_OPTIONS_DUKE_QRATS_MPICH")
#define AFLOWRC_MPI_COMMAND_DUKE_QRATS_MPICH                  string("/MAIN/bin/MPICH/bin/mpirun -np ") // DUKE_QRATS_MPICH
#define         MPI_COMMAND_DUKE_QRATS_MPICH                  XHOST.adefault.getattachedscheme("MPI_COMMAND_DUKE_QRATS_MPICH")
#define AFLOWRC_MPI_BINARY_DIR_DUKE_QRATS_MPICH               string("/usr/local/bin/") // DUKE_QRATS_MPICH
#define         MPI_BINARY_DIR_DUKE_QRATS_MPICH               XHOST.adefault.getattachedscheme("MPI_BINARY_DIR_DUKE_QRATS_MPICH")

#define AFLOWRC_MPI_OPTIONS_DUKE_QUSER_OPENMPI                string("ulimit -s unlimited ") // DUKE_QUSER_MPICH
#define         MPI_OPTIONS_DUKE_QUSER_OPENMPI                XHOST.adefault.getattachedscheme("MPI_OPTIONS_DUKE_QUSER_OPENMPI")
#define AFLOWRC_MPI_COMMAND_DUKE_QUSER_OPENMPI                string("/usr/bin/mpirun -n ") // DUKE_QUSER_MPICH
#define         MPI_COMMAND_DUKE_QUSER_OPENMPI                XHOST.adefault.getattachedscheme("MPI_COMMAND_DUKE_QUSER_OPENMPI")
#define AFLOWRC_MPI_BINARY_DIR_DUKE_QUSER_OPENMPI             string("/home/bin/") // DUKE_QUSER_MPICH
#define         MPI_BINARY_DIR_DUKE_QUSER_OPENMPI             XHOST.adefault.getattachedscheme("MPI_BINARY_DIR_DUKE_QUSER_OPENMPI")

#define AFLOWRC_MPI_OPTIONS_MPCDF_EOS_MPIIFORT                string("ulimit -s unlimited ") // MPCDF_EOS_MPICH
#define         MPI_OPTIONS_MPCDF_EOS_MPIIFORT                XHOST.adefault.getattachedscheme("MPI_OPTIONS_MPCDF_EOS_MPIIFORT")
#define AFLOWRC_MPI_COMMAND_MPCDF_EOS_MPIIFORT                string("/usr/bin/srun -n ") // MPCDF_EOS_MPICH
#define         MPI_COMMAND_MPCDF_EOS_MPIIFORT                XHOST.adefault.getattachedscheme("MPI_COMMAND_MPCDF_EOS_MPIIFORT")
#define AFLOWRC_MPI_BINARY_DIR_MPCDF_EOS_MPIIFORT             string("~/bin/") // MPCDF_EOS_MPICH
#define         MPI_BINARY_DIR_MPCDF_EOS_MPIIFORT             XHOST.adefault.getattachedscheme("MPI_BINARY_DIR_MPCDF_EOS_MPIIFORT")

#define AFLOWRC_MPI_OPTIONS_MPCDF_DRACO_MPIIFORT              string("ulimit -s unlimited ") // MPCDF_DRACO_MPICH  // FIX_DRACO
#define         MPI_OPTIONS_MPCDF_DRACO_MPIIFORT              XHOST.adefault.getattachedscheme("MPI_OPTIONS_MPCDF_DRACO_MPIIFORT")  // FIX_DRACO
#define AFLOWRC_MPI_COMMAND_MPCDF_DRACO_MPIIFORT              string("/usr/bin/srun -n ") // MPCDF_DRACO_MPICH // FIX_DRACO
#define         MPI_COMMAND_MPCDF_DRACO_MPIIFORT              XHOST.adefault.getattachedscheme("MPI_COMMAND_MPCDF_DRACO_MPIIFORT") // FIX_DRACO
#define AFLOWRC_MPI_BINARY_DIR_MPCDF_DRACO_MPIIFORT           string("~/bin/") // MPCDF_DRACO_MPICH // FIX_DRACO
#define         MPI_BINARY_DIR_MPCDF_DRACO_MPIIFORT           XHOST.adefault.getattachedscheme("MPI_BINARY_DIR_MPCDF_DRACO_MPIIFORT") // FIX_DRACO

#define AFLOWRC_MPI_OPTIONS_MPCDF_HYDRA_MPIIFORT              string("ulimit -s unlimited ") // MPCDF_HYDRA_MPICH
#define         MPI_OPTIONS_MPCDF_HYDRA_MPIIFORT              XHOST.adefault.getattachedscheme("MPI_OPTIONS_MPCDF_HYDRA_MPIIFORT")
#define AFLOWRC_MPI_COMMAND_MPCDF_HYDRA_MPIIFORT              string("/usr/bin/srun -n ") // MPCDF_HYDRA_MPICH
#define         MPI_COMMAND_MPCDF_HYDRA_MPIIFORT              XHOST.adefault.getattachedscheme("MPI_COMMAND_MPCDF_HYDRA_MPIIFORT")
#define AFLOWRC_MPI_BINARY_DIR_MPCDF_HYDRA_MPIIFORT           string("~/bin/") // MPCDF_HYDRA_MPICH
#define         MPI_BINARY_DIR_MPCDF_HYDRA_MPIIFORT           XHOST.adefault.getattachedscheme("MPI_BINARY_DIR_MPCDF_HYDRA_MPIIFORT")

#define AFLOWRC_MPI_OPTIONS_TERAGRID_RANGER                   string("") // TERAGRID_RANGER
#define         MPI_OPTIONS_TERAGRID_RANGER                   XHOST.adefault.getattachedscheme("MPI_OPTIONS_TERAGRID_RANGER")
#define AFLOWRC_MPI_COMMAND_TERAGRID_RANGER                   string("/share/sge6.2/default/pe_scripts/ibrun") // TERAGRID_RANGER
#define         MPI_COMMAND_TERAGRID_RANGER                   XHOST.adefault.getattachedscheme("MPI_COMMAND_TERAGRID_RANGER")
#define AFLOWRC_MPI_BINARY_DIR_TERAGRID_RANGER                string("~/bin/") // TERAGRID_RANGER
#define         MPI_BINARY_DIR_TERAGRID_RANGER                XHOST.adefault.getattachedscheme("MPI_BINARY_DIR_TERAGRID_RANGER")

#define AFLOWRC_MPI_OPTIONS_TERAGRID_KRAKEN                   string("") // TERAGRID_KRAKEN
#define         MPI_OPTIONS_TERAGRID_KRAKEN                   XHOST.adefault.getattachedscheme("MPI_OPTIONS_TERAGRID_KRAKEN")
#define AFLOWRC_MPI_COMMAND_TERAGRID_KRAKEN                   string("aprun -n") // TERAGRID_KRAKEN
#define         MPI_COMMAND_TERAGRID_KRAKEN                   XHOST.adefault.getattachedscheme("MPI_COMMAND_TERAGRID_KRAKEN")
#define AFLOWRC_MPI_BINARY_DIR_TERAGRID_KRAKEN                string("/nics/a/proj/aflow/bin/") // TERAGRID_KRAKEN
#define         MPI_BINARY_DIR_TERAGRID_KRAKEN                XHOST.adefault.getattachedscheme("MPI_BINARY_DIR_TERAGRID_KRAKEN")

#define AFLOWRC_MPI_OPTIONS_FULTON_MARYLOU                    string("export OMP_NUM_THREADS=1") // FULTON_MARYLOU
#define         MPI_OPTIONS_FULTON_MARYLOU                    XHOST.adefault.getattachedscheme("MPI_OPTIONS_FULTON_MARYLOU")
#define AFLOWRC_MPI_COMMAND_FULTON_MARYLOU                    string("mpiexec") // FULTON_MARYLOU
#define         MPI_COMMAND_FULTON_MARYLOU                    XHOST.adefault.getattachedscheme("MPI_COMMAND_FULTON_MARYLOU")
#define AFLOWRC_MPI_BINARY_DIR_FULTON_MARYLOU                 string("/fslgroup/fslg_datamining/bin/") // FULTON_MARYLOU
#define         MPI_BINARY_DIR_FULTON_MARYLOU                 XHOST.adefault.getattachedscheme("MPI_BINARY_DIR_FULTON_MARYLOU")

#define AFLOWRC_MPI_OPTIONS_TRINITY_PARSONS                   string("") // TRINITY COLLEGE IRELAND
#define         MPI_OPTIONS_TRINITY_PARSONS                   XHOST.adefault.getattachedscheme("MPI_OPTIONS_TRINITY_PARSONS")
#define AFLOWRC_MPI_COMMAND_TRINITY_PARSONS                   string("mpirun -np ") // TRINITY COLLEGE IRELAND
#define         MPI_COMMAND_TRINITY_PARSONS                   XHOST.adefault.getattachedscheme("MPI_COMMAND_TRINITY_PARSONS")
#define AFLOWRC_MPI_BINARY_DIR_TRINITY_PARSONS                string("/home/users/aflow/bin/") // TRINITY COLLEGE IRELAND
#define         MPI_BINARY_DIR_TRINITY_PARSONS                XHOST.adefault.getattachedscheme("MPI_BINARY_DIR_TRINITY_PARSONS")

#define AFLOWRC_MPI_OPTIONS_MACHINE1                          string("") // future expansions
#define         MPI_OPTIONS_MACHINE1                          XHOST.adefault.getattachedscheme("MPI_OPTIONS_MACHINE1")
#define AFLOWRC_MPI_COMMAND_MACHINE1                          string("...something ...")  // future expansions
#define         MPI_COMMAND_MACHINE1                          XHOST.adefault.getattachedscheme("MPI_COMMAND_MACHINE1")
#define AFLOWRC_MPI_BINARY_DIR_MACHINE1                       string("/somewhere/")  // future expansions
#define         MPI_BINARY_DIR_MACHINE1                       XHOST.adefault.getattachedscheme("MPI_BINARY_DIR_MACHINE1")

#define AFLOWRC_MPI_OPTIONS_MACHINE2                          string("") // future expansions
#define         MPI_OPTIONS_MACHINE2                          XHOST.adefault.getattachedscheme("MPI_OPTIONS_MACHINE2")
#define AFLOWRC_MPI_COMMAND_MACHINE2                          string("stub not used")  // future expansions
#define         MPI_COMMAND_MACHINE2                          XHOST.adefault.getattachedscheme("MPI_COMMAND_MACHINE2")
#define AFLOWRC_MPI_BINARY_DIR_MACHINE2                       string("/home/aflow/bin/")  // future expansions
#define         MPI_BINARY_DIR_MACHINE2                       XHOST.adefault.getattachedscheme("MPI_BINARY_DIR_MACHINE2")
 
#endif // _AFLOW_AFLOWRC_H_

// POCC STUFF
// defaults go in 4 positions; Here with #define AFLOWRC_DEFAULT, in read(), in write_default(), and in print_aflowrc()...
// I coded strings (without spaces), <int>.. you can do <doubles> just like the <int>
// for strings with spaces I need to fix the code. Dont add them now. Then you need to go around the whole code and fix the use of DEFAULTS, possibly also in the READMEs if they are specified.
// STRING     string blablabla=DEFAULT_STRING =>  string blablabla=XHOST.adefault.getattachedscheme("DEFAULT_STRING")
// INT        int blablabla=DEFAULT_INT       =>  int blablabla=XHOST.adefault.getattachedscheme<int>("DEFAULT_INT")
// DOUBLE     double blablabla=DEFAULT_DOUBLE =>  double blablabla=XHOST.adefault.getattachedscheme<double>("DEFAULT_DOUBLE")
// ./aflow --machine to check them out


#ifndef _AFLOW_AFLOWRC_CPP_
#define _AFLOW_AFLOWRC_CPP_

// ***************************************************************************
// aflowrc::load_default
// ***************************************************************************
namespace aflowrc {
  bool load_default(string schema,string schema_default) {
    bool found=FALSE;
    string aus,string_to_add=schema_default;
    vector<string> tokens;
    for(uint i=0;i<XHOST.vaflowrc.size()&&!found;i++) {
      aurostd::string2tokens(XHOST.vaflowrc.at(i),tokens,"=");
      if(tokens.size()>0&&!found) {
	if(aurostd::RemoveWhiteSpaces(tokens.at(0))==schema&&!found) {
	  found=TRUE; aus=tokens.at(1);
	  aurostd::string2tokens(aus,tokens,"\""); //	    if(tokens.size()>0) cerr << tokens.at(0) << endl;
	  if(tokens.size()>0) string_to_add=tokens.at(0);
	}
      }
    }
    // fix ~/ with XHOST.User
    if(aurostd::substring2bool(string_to_add,"~/")) aurostd::StringSubst(string_to_add,"~/",XHOST.Home+"/");
    XHOST.adefault.addattachedscheme(schema,string_to_add,TRUE); // add what is present or the default if not present
    return found;
  }
  template<class utype> 
  bool load_default(string schema,utype schema_default) {
    bool found=XHOST.adefault.args2addattachedscheme(XHOST.vaflowrc,schema,string(schema+"="),""); // add what is present
    if(!found) XHOST.adefault.addattachedscheme(schema,aurostd::utype2string<utype>(schema_default),TRUE);  // add default if not present
    return found;
  }
}
  
// ***************************************************************************
// aflowrc::is_available
// ***************************************************************************
namespace aflowrc {
  bool is_available(std::ostream& oss,bool AFLOWRC_VERBOSE) {
    bool LDEBUG=(FALSE || XHOST.DEBUG || AFLOWRC_VERBOSE);   
    bool aflowrc_local=FALSE;
    bool aflowrc_global=FALSE;
    if(LDEBUG) oss << "aflowrc::is_available: BEGIN" << endl;
    if(LDEBUG) oss << "aflowrc::is_available: XHOST.Home=" << XHOST.Home << endl;
    // TESTING LOCAL OR USER BASED
    if(XHOST.aflowrc_filename.empty()) XHOST.aflowrc_filename=AFLOWRC_FILENAME_LOCAL;
    aflowrc_local=aurostd::FileExist(AFLOWRC_FILENAME_LOCAL);
    aflowrc_global=aurostd::FileExist(AFLOWRC_FILENAME_GLOBAL);

    // LOCAL=TRUE && GLOBAL=TRUE => take LOCAL
    if(aflowrc_local && aflowrc_global) {
      if(LDEBUG) oss << "aflowrc::is_available: LOCAL=TRUE && GLOBAL=TRUE => LOCAL " << endl;
      XHOST.aflowrc_filename=AFLOWRC_FILENAME_LOCAL;
      if(LDEBUG) oss << "aflowrc::is_available: XHOST.aflowrc_filename=" << XHOST.aflowrc_filename << endl;
      if(LDEBUG) oss << "aflowrc::is_available: END" << endl;
      return TRUE;
    }
    // LOCAL=TRUE && GLOBAL=FALSE => take LOCAL
    if(aflowrc_local && !aflowrc_global) {
      if(LDEBUG) oss << "aflowrc::is_available: LOCAL=TRUE && GLOBAL=FALSE => LOCAL " << endl;
      XHOST.aflowrc_filename=AFLOWRC_FILENAME_LOCAL; 
      if(LDEBUG) oss << "aflowrc::is_available: XHOST.aflowrc_filename=" << XHOST.aflowrc_filename << endl;
      if(LDEBUG) oss << "aflowrc::is_available: END" << endl;
      return TRUE;
    }
    // LOCAL=FALSE && GLOBAL=TRUE => take GLOBAL
    if(!aflowrc_local && aflowrc_global) {
      if(LDEBUG) oss << "aflowrc::is_available: LOCAL=FALSE && GLOBAL=TRUE => GLOBAL " << endl;
      XHOST.aflowrc_filename=AFLOWRC_FILENAME_GLOBAL;
      if(LDEBUG) oss << "aflowrc::is_available: XHOST.aflowrc_filename=" << XHOST.aflowrc_filename << endl;
      if(LDEBUG) oss << "aflowrc::is_available: END" << endl;
      return TRUE;
    }
    // LOCAL=FALSE && GLOBAL=FALSE => take NOTHING AND REWRITE
    if(!aflowrc_local && !aflowrc_global) {
      if(LDEBUG) oss << "aflowrc::is_available: LOCAL=FALSE && GLOBAL=FALSE => NOTHING " << endl;
      XHOST.aflowrc_filename=AFLOWRC_FILENAME_LOCAL; // because it is going to write it
      if(LDEBUG) oss << "aflowrc::is_available: XHOST.aflowrc_filename=" << XHOST.aflowrc_filename << endl;
      if(LDEBUG) oss << "aflowrc::is_available: END" << endl;
      return FALSE;
    }

    if(LDEBUG) oss << "aflowrc::is_available: END" << endl;
    return FALSE;
  }
} // namespace aflowrc


// ***************************************************************************
// aflowrc::read
// ***************************************************************************
namespace aflowrc {
  bool read(std::ostream& oss,bool AFLOWRC_VERBOSE) {
    bool LDEBUG=(FALSE || XHOST.DEBUG || AFLOWRC_VERBOSE);   
    if(LDEBUG) oss << "aflowrc::read: BEGIN" << endl;
    if(LDEBUG) oss << "aflowrc::read: XHOST.Home=" << XHOST.Home << endl;
    if(XHOST.aflowrc_filename.empty()) XHOST.aflowrc_filename=AFLOWRC_FILENAME_LOCAL;
    if(LDEBUG) oss << "aflowrc::read: XHOST.aflowrc_filename=" << XHOST.aflowrc_filename << endl;

    if(!aflowrc::is_available(oss,AFLOWRC_VERBOSE)) 
      if(!aurostd::substring2bool(XHOST.aflowrc_filename,"/mnt/MAIN"))
	 cout << "WARNING aflowrc::read: " << XHOST.aflowrc_filename << " not found, loading DEFAULT values" << endl;
	 
    aurostd::file2string(XHOST.aflowrc_filename,XHOST.aflowrc_content);
    // oss << "BEGIN" << endl << XHOST.aflowrc_content << "END" << endl;
    // XHOST.aflowrc_content=aurostd::RemoveComments(XHOST.aflowrc_content); // NOW Clean XHOST.aflowrc_content
    //   XHOST.aflowrc_content=aurostd::RemoveWhiteSpaces(XHOST.aflowrc_content); // NOW Clean XHOST.aflowrc_content
    XHOST.aflowrc_content=aurostd::RemoveComments(XHOST.aflowrc_content); // NOW Clean XHOST.aflowrc_content
    // oss << "BEGIN" << endl << XHOST.aflowrc_content << "END" << endl;
    aurostd::string2vectorstring(XHOST.aflowrc_content,XHOST.vaflowrc); // vectorize

    // DEFAULT GENERIC MPI
    aflowrc::load_default("MPI_START_DEFAULT",AFLOWRC_MPI_START_DEFAULT); 
    aflowrc::load_default("MPI_STOP_DEFAULT",AFLOWRC_MPI_STOP_DEFAULT); 
    aflowrc::load_default("MPI_COMMAND_DEFAULT",AFLOWRC_MPI_COMMAND_DEFAULT); 
    aflowrc::load_default("MPI_NCPUS_DEFAULT",AFLOWRC_MPI_NCPUS_DEFAULT); 
    aflowrc::load_default("MPI_NCPUS_MAX",AFLOWRC_MPI_NCPUS_MAX); 

    // BINARY VASP
    aflowrc::load_default("DEFAULT_VASP_GAMMA_BIN",AFLOWRC_DEFAULT_VASP_GAMMA_BIN); 
    aflowrc::load_default("DEFAULT_VASP_GAMMA_MPI_BIN",AFLOWRC_DEFAULT_VASP_GAMMA_MPI_BIN); 
    aflowrc::load_default("DEFAULT_VASP_BIN",AFLOWRC_DEFAULT_VASP_BIN); 
    aflowrc::load_default("DEFAULT_VASP_MPI_BIN",AFLOWRC_DEFAULT_VASP_MPI_BIN); 
    aflowrc::load_default("DEFAULT_VASP5_BIN",AFLOWRC_DEFAULT_VASP5_BIN); 
    aflowrc::load_default("DEFAULT_VASP5_MPI_BIN",AFLOWRC_DEFAULT_VASP5_MPI_BIN); 
    // BINARY AIMS
    aflowrc::load_default("DEFAULT_AIMS_BIN",AFLOWRC_DEFAULT_AIMS_BIN); 
    
    // POTCARS
    aflowrc::load_default("DEFAULT_VASP_POTCAR_DIRECTORIES",AFLOWRC_DEFAULT_VASP_POTCAR_DIRECTORIES); 
    aflowrc::load_default("DEFAULT_VASP_POTCAR_DATE",AFLOWRC_DEFAULT_VASP_POTCAR_DATE); 
    aflowrc::load_default("DEFAULT_VASP_POTCAR_SUFFIX",AFLOWRC_DEFAULT_VASP_POTCAR_SUFFIX); 
    aflowrc::load_default("DEFAULT_VASP_POTCAR_DATE_POT_LDA",AFLOWRC_DEFAULT_VASP_POTCAR_DATE_POT_LDA); 
    aflowrc::load_default("DEFAULT_VASP_POTCAR_DATE_POT_GGA",AFLOWRC_DEFAULT_VASP_POTCAR_DATE_POT_GGA); 

    // DEFAULT KPOINTS/DOS
    aflowrc::load_default("DEFAULT_BANDS_GRID",AFLOWRC_DEFAULT_BANDS_GRID); 
    aflowrc::load_default("DEFAULT_BANDS_LATTICE",AFLOWRC_DEFAULT_BANDS_LATTICE); 
    aflowrc::load_default("DEFAULT_KSCHEME",AFLOWRC_DEFAULT_KSCHEME); 
    aflowrc::load_default("DEFAULT_KPPRA",AFLOWRC_DEFAULT_KPPRA); 
    aflowrc::load_default("DEFAULT_STATIC_KSCHEME",AFLOWRC_DEFAULT_STATIC_KSCHEME); 
    aflowrc::load_default("DEFAULT_KPPRA_STATIC",AFLOWRC_DEFAULT_KPPRA_STATIC); 
    aflowrc::load_default("DEFAULT_KPPRA_ICSD",AFLOWRC_DEFAULT_KPPRA_ICSD); 
    aflowrc::load_default("DEFAULT_UNARY_BANDS_GRID",AFLOWRC_DEFAULT_UNARY_BANDS_GRID); 
    aflowrc::load_default("DEFAULT_UNARY_KPPRA",AFLOWRC_DEFAULT_UNARY_KPPRA); 
    aflowrc::load_default("DEFAULT_UNARY_KPPRA_STATIC",AFLOWRC_DEFAULT_UNARY_KPPRA_STATIC); 
    aflowrc::load_default("DEFAULT_PHONONS_KSCHEME",AFLOWRC_DEFAULT_PHONONS_KSCHEME); 
    aflowrc::load_default("DEFAULT_DOS_EMIN",AFLOWRC_DEFAULT_DOS_EMIN); 
    aflowrc::load_default("DEFAULT_DOS_EMAX",AFLOWRC_DEFAULT_DOS_EMAX); 
    aflowrc::load_default("DEFAULT_DOS_SCALE",AFLOWRC_DEFAULT_DOS_SCALE); 

    // PRECISION
    aflowrc::load_default("DEFAULT_VASP_PREC_ENMAX_LOW",AFLOWRC_DEFAULT_VASP_PREC_ENMAX_LOW);
    aflowrc::load_default("DEFAULT_VASP_PREC_ENMAX_MEDIUM",AFLOWRC_DEFAULT_VASP_PREC_ENMAX_MEDIUM);
    aflowrc::load_default("DEFAULT_VASP_PREC_ENMAX_NORMAL",AFLOWRC_DEFAULT_VASP_PREC_ENMAX_NORMAL);
    aflowrc::load_default("DEFAULT_VASP_PREC_ENMAX_HIGH",AFLOWRC_DEFAULT_VASP_PREC_ENMAX_HIGH);
    aflowrc::load_default("DEFAULT_VASP_PREC_ENMAX_ACCURATE",AFLOWRC_DEFAULT_VASP_PREC_ENMAX_ACCURATE);
    aflowrc::load_default("DEFAULT_VASP_SPIN_REMOVE_CUTOFF",AFLOWRC_DEFAULT_VASP_SPIN_REMOVE_CUTOFF);
    aflowrc::load_default("DEFAULT_VASP_PREC_POTIM",AFLOWRC_DEFAULT_VASP_PREC_POTIM);
    aflowrc::load_default("DEFAULT_VASP_PREC_EDIFFG",AFLOWRC_DEFAULT_VASP_PREC_EDIFFG);

    // OPTIONS
    aflowrc::load_default("DEFAULT_VASP_EXTERNAL_INCAR",AFLOWRC_DEFAULT_VASP_EXTERNAL_INCAR);
    aflowrc::load_default("DEFAULT_VASP_EXTERNAL_POSCAR",AFLOWRC_DEFAULT_VASP_EXTERNAL_POSCAR);
    aflowrc::load_default("DEFAULT_VASP_EXTERNAL_POTCAR",AFLOWRC_DEFAULT_VASP_EXTERNAL_POTCAR);
    aflowrc::load_default("DEFAULT_VASP_EXTERNAL_KPOINTS",AFLOWRC_DEFAULT_VASP_EXTERNAL_KPOINTS);
    aflowrc::load_default("DEFAULT_VASP_PSEUDOPOTENTIAL_TYPE",AFLOWRC_DEFAULT_VASP_PSEUDOPOTENTIAL_TYPE);
    aflowrc::load_default("DEFAULT_VASP_FORCE_OPTION_RELAX_MODE_SCHEME",AFLOWRC_DEFAULT_VASP_FORCE_OPTION_RELAX_MODE_SCHEME);
    aflowrc::load_default("DEFAULT_VASP_FORCE_OPTION_PREC_SCHEME",AFLOWRC_DEFAULT_VASP_FORCE_OPTION_PREC_SCHEME);
    aflowrc::load_default("DEFAULT_VASP_FORCE_OPTION_ALGO_SCHEME",AFLOWRC_DEFAULT_VASP_FORCE_OPTION_ALGO_SCHEME);
    aflowrc::load_default("DEFAULT_VASP_FORCE_OPTION_METAGGA_SCHEME",AFLOWRC_DEFAULT_VASP_FORCE_OPTION_METAGGA_SCHEME);
    aflowrc::load_default("DEFAULT_VASP_FORCE_OPTION_IVDW_SCHEME",AFLOWRC_DEFAULT_VASP_FORCE_OPTION_IVDW_SCHEME);
    aflowrc::load_default("DEFAULT_VASP_FORCE_OPTION_TYPE_SCHEME",AFLOWRC_DEFAULT_VASP_FORCE_OPTION_TYPE_SCHEME);
    aflowrc::load_default("DEFAULT_VASP_FORCE_OPTION_ABMIX_SCHEME",AFLOWRC_DEFAULT_VASP_FORCE_OPTION_ABMIX_SCHEME);
    aflowrc::load_default("DEFAULT_VASP_FORCE_OPTION_SYM",AFLOWRC_DEFAULT_VASP_FORCE_OPTION_SYM);
    aflowrc::load_default("DEFAULT_VASP_FORCE_OPTION_SPIN",AFLOWRC_DEFAULT_VASP_FORCE_OPTION_SPIN);
    aflowrc::load_default("DEFAULT_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_1",AFLOWRC_DEFAULT_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_1);
    aflowrc::load_default("DEFAULT_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_2",AFLOWRC_DEFAULT_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_2);
    aflowrc::load_default("DEFAULT_VASP_FORCE_OPTION_BADER",AFLOWRC_DEFAULT_VASP_FORCE_OPTION_BADER);
    aflowrc::load_default("DEFAULT_VASP_FORCE_OPTION_ELF",AFLOWRC_DEFAULT_VASP_FORCE_OPTION_ELF);
    aflowrc::load_default("DEFAULT_VASP_FORCE_OPTION_AUTO_MAGMOM",AFLOWRC_DEFAULT_VASP_FORCE_OPTION_AUTO_MAGMOM);
    aflowrc::load_default("DEFAULT_VASP_FORCE_OPTION_WAVECAR",AFLOWRC_DEFAULT_VASP_FORCE_OPTION_WAVECAR);
    aflowrc::load_default("DEFAULT_VASP_FORCE_OPTION_CHGCAR",AFLOWRC_DEFAULT_VASP_FORCE_OPTION_CHGCAR);
    aflowrc::load_default("DEFAULT_VASP_FORCE_OPTION_LSCOUPLING",AFLOWRC_DEFAULT_VASP_FORCE_OPTION_LSCOUPLING);

    // AFLOW_LIBRARY AFLOW_PROJECT
    aflowrc::load_default("DEFAULT_AFLOW_LIBRARY_DIRECTORIES",AFLOWRC_DEFAULT_AFLOW_LIBRARY_DIRECTORIES);
    aflowrc::load_default("DEFAULT_AFLOW_PROJECTS_DIRECTORIES",AFLOWRC_DEFAULT_AFLOW_PROJECTS_DIRECTORIES);
    
    // DEFAULT PLATON/FINDSYM
    aflowrc::load_default("DEFAULT_PLATON_P_EQUAL",AFLOWRC_DEFAULT_PLATON_P_EQUAL);
    aflowrc::load_default("DEFAULT_PLATON_P_EXACT",AFLOWRC_DEFAULT_PLATON_P_EXACT);
    aflowrc::load_default("DEFAULT_PLATON_P_ANG",AFLOWRC_DEFAULT_PLATON_P_ANG);
    aflowrc::load_default("DEFAULT_PLATON_P_D1",AFLOWRC_DEFAULT_PLATON_P_D1);
    aflowrc::load_default("DEFAULT_PLATON_P_D2",AFLOWRC_DEFAULT_PLATON_P_D2);
    aflowrc::load_default("DEFAULT_PLATON_P_D3",AFLOWRC_DEFAULT_PLATON_P_D3);
    aflowrc::load_default("DEFAULT_FINDSYM_TOL",AFLOWRC_DEFAULT_FINDSYM_TOL);

    // DEFAULT GNUPLOT
    aflowrc::load_default("DEFAULT_GNUPLOT_EPS_FONT",AFLOWRC_DEFAULT_GNUPLOT_EPS_FONT);
    aflowrc::load_default("DEFAULT_GNUPLOT_EPS_FONT_BOLD",AFLOWRC_DEFAULT_GNUPLOT_EPS_FONT_BOLD); 
    aflowrc::load_default("DEFAULT_GNUPLOT_EPS_FONT_ITALICS",AFLOWRC_DEFAULT_GNUPLOT_EPS_FONT_ITALICS); 
    aflowrc::load_default("DEFAULT_GNUPLOT_EPS_FONT_BOLD_ITALICS",AFLOWRC_DEFAULT_GNUPLOT_EPS_FONT_BOLD_ITALICS); 
    aflowrc::load_default("DEFAULT_GNUPLOT_PNG_FONT",AFLOWRC_DEFAULT_GNUPLOT_PNG_FONT); 
    aflowrc::load_default("DEFAULT_GNUPLOT_PNG_FONT_BOLD",AFLOWRC_DEFAULT_GNUPLOT_PNG_FONT_BOLD); 
    aflowrc::load_default("DEFAULT_GNUPLOT_PNG_FONT_ITALICS",AFLOWRC_DEFAULT_GNUPLOT_PNG_FONT_ITALICS); 
    aflowrc::load_default("DEFAULT_GNUPLOT_PNG_FONT_BOLD_ITALICS",AFLOWRC_DEFAULT_GNUPLOT_PNG_FONT_BOLD_ITALICS); 
 
    // DEFAULT CHULL
    aflowrc::load_default("DEFAULT_CHULL_OUTLIER_MULTIPLIER",AFLOWRC_DEFAULT_CHULL_OUTLIER_MULTIPLIER); 
    aflowrc::load_default("DEFAULT_CHULL_LATEX_BANNER",AFLOWRC_DEFAULT_CHULL_LATEX_BANNER); 
    aflowrc::load_default("DEFAULT_CHULL_LATEX_COMPOUNDS_COLUMN",AFLOWRC_DEFAULT_CHULL_LATEX_COMPOUNDS_COLUMN); 
    aflowrc::load_default("DEFAULT_CHULL_LATEX_COMPOSITION_HEADER",AFLOWRC_DEFAULT_CHULL_LATEX_COMPOSITION_HEADER); 
    aflowrc::load_default("DEFAULT_CHULL_LATEX_PLOT_ENDPOINTS",AFLOWRC_DEFAULT_CHULL_LATEX_PLOT_ENDPOINTS); 
    aflowrc::load_default("DEFAULT_CHULL_LATEX_PLOT_OFF_HULL",AFLOWRC_DEFAULT_CHULL_LATEX_PLOT_OFF_HULL); 
    aflowrc::load_default("DEFAULT_CHULL_LATEX_PLOT_UNSTABLE",AFLOWRC_DEFAULT_CHULL_LATEX_PLOT_UNSTABLE); 
    aflowrc::load_default("DEFAULT_CHULL_LATEX_FILTER_SCHEME",AFLOWRC_DEFAULT_CHULL_LATEX_FILTER_SCHEME); 
    aflowrc::load_default("DEFAULT_CHULL_LATEX_FILTER_VALUE",AFLOWRC_DEFAULT_CHULL_LATEX_FILTER_VALUE); 
    aflowrc::load_default("DEFAULT_CHULL_LATEX_COLOR_BAR",AFLOWRC_DEFAULT_CHULL_LATEX_COLOR_BAR); 
    aflowrc::load_default("DEFAULT_CHULL_LATEX_HEAT_MAP",AFLOWRC_DEFAULT_CHULL_LATEX_HEAT_MAP); 
    aflowrc::load_default("DEFAULT_CHULL_LATEX_COLOR_GRADIENT",AFLOWRC_DEFAULT_CHULL_LATEX_COLOR_GRADIENT); 
    aflowrc::load_default("DEFAULT_CHULL_LATEX_COLOR_MAP",AFLOWRC_DEFAULT_CHULL_LATEX_COLOR_MAP); 
    aflowrc::load_default("DEFAULT_CHULL_LATEX_TERNARY_LABEL_COLOR",AFLOWRC_DEFAULT_CHULL_LATEX_TERNARY_LABEL_COLOR); 
    aflowrc::load_default("DEFAULT_CHULL_LATEX_REVERSE_AXIS",AFLOWRC_DEFAULT_CHULL_LATEX_REVERSE_AXIS); 
    aflowrc::load_default("DEFAULT_CHULL_LATEX_FACET_LINE_DROP_SHADOW",AFLOWRC_DEFAULT_CHULL_LATEX_FACET_LINE_DROP_SHADOW); 
    aflowrc::load_default("DEFAULT_CHULL_LATEX_LINKS",AFLOWRC_DEFAULT_CHULL_LATEX_LINKS); 
    aflowrc::load_default("DEFAULT_CHULL_LATEX_LABEL_NAME",AFLOWRC_DEFAULT_CHULL_LATEX_LABEL_NAME); 
    aflowrc::load_default("DEFAULT_CHULL_LATEX_META_LABELS",AFLOWRC_DEFAULT_CHULL_LATEX_META_LABELS); 
    aflowrc::load_default("DEFAULT_CHULL_LATEX_LABELS_OFF_HULL",AFLOWRC_DEFAULT_CHULL_LATEX_LABELS_OFF_HULL); 
    aflowrc::load_default("DEFAULT_CHULL_LATEX_PLOT_REDUCE_COMPOSITION",AFLOWRC_DEFAULT_CHULL_LATEX_PLOT_REDUCE_COMPOSITION); 
    aflowrc::load_default("DEFAULT_CHULL_LATEX_HELVETICA_FONT",AFLOWRC_DEFAULT_CHULL_LATEX_HELVETICA_FONT); 
    aflowrc::load_default("DEFAULT_CHULL_LATEX_FONT_SIZE",AFLOWRC_DEFAULT_CHULL_LATEX_FONT_SIZE); 
    aflowrc::load_default("DEFAULT_CHULL_LATEX_ROTATE_LABELS",AFLOWRC_DEFAULT_CHULL_LATEX_ROTATE_LABELS); 
    aflowrc::load_default("DEFAULT_CHULL_LATEX_BOLD_LABELS",AFLOWRC_DEFAULT_CHULL_LATEX_BOLD_LABELS); 

 
    // DEFAULT CORE
    aflowrc::load_default("AFLOW_CORE_TEMPERATURE_BEEP",AFLOWRC_AFLOW_CORE_TEMPERATURE_BEEP);
    aflowrc::load_default("AFLOW_CORE_TEMPERATURE_HALT",AFLOWRC_AFLOW_CORE_TEMPERATURE_HALT);
    aflowrc::load_default("AFLOW_CORE_TEMPERATURE_REFRESH",AFLOWRC_AFLOW_CORE_TEMPERATURE_REFRESH);

    // DEFAULT MACHINE DEPENDENT MPI
    aflowrc::load_default("MPI_OPTIONS_DUKE_BETA_MPICH",AFLOWRC_MPI_OPTIONS_DUKE_BETA_MPICH); 
    aflowrc::load_default("MPI_COMMAND_DUKE_BETA_MPICH",AFLOWRC_MPI_COMMAND_DUKE_BETA_MPICH); 
    aflowrc::load_default("MPI_BINARY_DIR_DUKE_BETA_MPICH",AFLOWRC_MPI_BINARY_DIR_DUKE_BETA_MPICH); 

    aflowrc::load_default("MPI_OPTIONS_DUKE_BETA_OPENMPI",AFLOWRC_MPI_OPTIONS_DUKE_BETA_OPENMPI); 
    aflowrc::load_default("MPI_COMMAND_DUKE_BETA_OPENMPI",AFLOWRC_MPI_COMMAND_DUKE_BETA_OPENMPI); 
    aflowrc::load_default("MPI_BINARY_DIR_DUKE_BETA_OPENMPI",AFLOWRC_MPI_BINARY_DIR_DUKE_BETA_OPENMPI); 

    aflowrc::load_default("MPI_OPTIONS_DUKE_MATERIALS",AFLOWRC_MPI_OPTIONS_DUKE_MATERIALS); 
    aflowrc::load_default("MPI_COMMAND_DUKE_MATERIALS",AFLOWRC_MPI_COMMAND_DUKE_MATERIALS); 
    aflowrc::load_default("MPI_BINARY_DIR_DUKE_MATERIALS",AFLOWRC_MPI_BINARY_DIR_DUKE_MATERIALS); 

    aflowrc::load_default("MPI_OPTIONS_DUKE_AFLOWLIB",AFLOWRC_MPI_OPTIONS_DUKE_AFLOWLIB); 
    aflowrc::load_default("MPI_COMMAND_DUKE_AFLOWLIB",AFLOWRC_MPI_COMMAND_DUKE_AFLOWLIB); 
    aflowrc::load_default("MPI_BINARY_DIR_DUKE_AFLOWLIB",AFLOWRC_MPI_BINARY_DIR_DUKE_AFLOWLIB); 

    aflowrc::load_default("MPI_OPTIONS_DUKE_HABANA",AFLOWRC_MPI_OPTIONS_DUKE_HABANA); 
    aflowrc::load_default("MPI_COMMAND_DUKE_HABANA",AFLOWRC_MPI_COMMAND_DUKE_HABANA); 
    aflowrc::load_default("MPI_BINARY_DIR_DUKE_HABANA",AFLOWRC_MPI_BINARY_DIR_DUKE_HABANA); 

    aflowrc::load_default("MPI_OPTIONS_DUKE_QRATS_MPICH",AFLOWRC_MPI_OPTIONS_DUKE_QRATS_MPICH); 
    aflowrc::load_default("MPI_COMMAND_DUKE_QRATS_MPICH",AFLOWRC_MPI_COMMAND_DUKE_QRATS_MPICH); 
    aflowrc::load_default("MPI_BINARY_DIR_DUKE_QRATS_MPICH",AFLOWRC_MPI_BINARY_DIR_DUKE_QRATS_MPICH); 

    aflowrc::load_default("MPI_OPTIONS_DUKE_QUSER_OPENMPI",AFLOWRC_MPI_OPTIONS_DUKE_QUSER_OPENMPI); 
    aflowrc::load_default("MPI_COMMAND_DUKE_QUSER_OPENMPI",AFLOWRC_MPI_COMMAND_DUKE_QUSER_OPENMPI); 
    aflowrc::load_default("MPI_BINARY_DIR_DUKE_QUSER_OPENMPI",AFLOWRC_MPI_BINARY_DIR_DUKE_QUSER_OPENMPI); 

    aflowrc::load_default("MPI_OPTIONS_MPCDF_EOS_MPIIFORT",AFLOWRC_MPI_OPTIONS_MPCDF_EOS_MPIIFORT); 
    aflowrc::load_default("MPI_COMMAND_MPCDF_EOS_MPIIFORT",AFLOWRC_MPI_COMMAND_MPCDF_EOS_MPIIFORT); 
    aflowrc::load_default("MPI_BINARY_DIR_MPCDF_EOS_MPIIFORT",AFLOWRC_MPI_BINARY_DIR_MPCDF_EOS_MPIIFORT); 

    aflowrc::load_default("MPI_OPTIONS_MPCDF_DRACO_MPIIFORT",AFLOWRC_MPI_OPTIONS_MPCDF_DRACO_MPIIFORT); 
    aflowrc::load_default("MPI_COMMAND_MPCDF_DRACO_MPIIFORT",AFLOWRC_MPI_COMMAND_MPCDF_DRACO_MPIIFORT); 
    aflowrc::load_default("MPI_BINARY_DIR_MPCDF_DRACO_MPIIFORT",AFLOWRC_MPI_BINARY_DIR_MPCDF_DRACO_MPIIFORT); 

    aflowrc::load_default("MPI_OPTIONS_MPCDF_HYDRA_MPIIFORT",AFLOWRC_MPI_OPTIONS_MPCDF_HYDRA_MPIIFORT); 
    aflowrc::load_default("MPI_COMMAND_MPCDF_HYDRA_MPIIFORT",AFLOWRC_MPI_COMMAND_MPCDF_HYDRA_MPIIFORT); 
    aflowrc::load_default("MPI_BINARY_DIR_MPCDF_HYDRA_MPIIFORT",AFLOWRC_MPI_BINARY_DIR_MPCDF_HYDRA_MPIIFORT); 

    aflowrc::load_default("MPI_OPTIONS_TERAGRID_RANGER",AFLOWRC_MPI_OPTIONS_TERAGRID_RANGER); 
    aflowrc::load_default("MPI_COMMAND_TERAGRID_RANGER",AFLOWRC_MPI_COMMAND_TERAGRID_RANGER); 
    aflowrc::load_default("MPI_BINARY_DIR_TERAGRID_RANGER",AFLOWRC_MPI_BINARY_DIR_TERAGRID_RANGER); 

    aflowrc::load_default("MPI_OPTIONS_TERAGRID_KRAKEN",AFLOWRC_MPI_OPTIONS_TERAGRID_KRAKEN); 
    aflowrc::load_default("MPI_COMMAND_TERAGRID_KRAKEN",AFLOWRC_MPI_COMMAND_TERAGRID_KRAKEN); 
    aflowrc::load_default("MPI_BINARY_DIR_TERAGRID_KRAKEN",AFLOWRC_MPI_BINARY_DIR_TERAGRID_KRAKEN); 

    aflowrc::load_default("MPI_OPTIONS_FULTON_MARYLOU",AFLOWRC_MPI_OPTIONS_FULTON_MARYLOU); 
    aflowrc::load_default("MPI_COMMAND_FULTON_MARYLOU",AFLOWRC_MPI_COMMAND_FULTON_MARYLOU); 
    aflowrc::load_default("MPI_BINARY_DIR_FULTON_MARYLOU",AFLOWRC_MPI_BINARY_DIR_FULTON_MARYLOU); 

    aflowrc::load_default("MPI_OPTIONS_TRINITY_PARSONS",AFLOWRC_MPI_OPTIONS_TRINITY_PARSONS); 
    aflowrc::load_default("MPI_COMMAND_TRINITY_PARSONS",AFLOWRC_MPI_COMMAND_TRINITY_PARSONS); 
    aflowrc::load_default("MPI_BINARY_DIR_TRINITY_PARSONS",AFLOWRC_MPI_BINARY_DIR_TRINITY_PARSONS); 

    aflowrc::load_default("MPI_OPTIONS_MACHINE1",AFLOWRC_MPI_OPTIONS_MACHINE1); 
    aflowrc::load_default("MPI_COMMAND_MACHINE1",AFLOWRC_MPI_COMMAND_MACHINE1); 
    aflowrc::load_default("MPI_BINARY_DIR_MACHINE1",AFLOWRC_MPI_BINARY_DIR_MACHINE1); 

    aflowrc::load_default("MPI_OPTIONS_MACHINE2",AFLOWRC_MPI_OPTIONS_MACHINE2); 
    aflowrc::load_default("MPI_COMMAND_MACHINE2",AFLOWRC_MPI_COMMAND_MACHINE2); 
    aflowrc::load_default("MPI_BINARY_DIR_MACHINE2",AFLOWRC_MPI_BINARY_DIR_MACHINE2); 
    
    if(LDEBUG) oss << "aflowrc::read: END" << endl;

    return TRUE;
  }
} // namespace aflowrc


// ***************************************************************************
// aflowrc::write_default
// ***************************************************************************
namespace aflowrc {
  bool write_default(std::ostream& oss,bool AFLOWRC_VERBOSE) {
    bool LDEBUG=(FALSE || XHOST.DEBUG || AFLOWRC_VERBOSE);   
    if(LDEBUG) oss << "aflowrc::write_default: BEGIN" << endl;
    if(LDEBUG) oss << "aflowrc::write_default: XHOST.Home=" << XHOST.Home << endl;
    if(XHOST.aflowrc_filename.empty()) XHOST.aflowrc_filename=AFLOWRC_FILENAME_LOCAL;
    if(LDEBUG) oss << "aflowrc::write_default: XHOST.aflowrc_filename=" << XHOST.aflowrc_filename << endl;

    stringstream aflowrc("");
    aflowrc << "// ****************************************************************************************************" << endl;
    aflowrc << "// *                                                                                                  *" << endl;
    aflowrc << "// *                          aflow - Automatic-FLOW for materials discovery                          *" << endl;
    aflowrc << "// *                aflow.org consortium - High-Throughput ab-initio Computing Project                *" << endl;
    aflowrc << "// *                                                                                                  *" << endl;
    aflowrc << "// ****************************************************************************************************" << endl;
    aflowrc << "// DEFAULT .aflow.rc generated by AFLOW V" << string(AFLOW_VERSION) << endl;
    aflowrc << "// comments with // ignored... " << endl;
    aflowrc << "// strings are with=\"...\" " << endl;

    aflowrc << " " << endl;
    aflowrc << "// DEFAULT GENERIC MPI " << endl;
    aflowrc << "MPI_START_DEFAULT=\"" << AFLOWRC_MPI_START_DEFAULT << "\"" << endl;
    aflowrc << "MPI_STOP_DEFAULT=\"" << AFLOWRC_MPI_STOP_DEFAULT << "\"" << endl; 
    aflowrc << "MPI_COMMAND_DEFAULT=\"" << AFLOWRC_MPI_COMMAND_DEFAULT << "\"" << endl;
    aflowrc << "MPI_NCPUS_DEFAULT=" << AFLOWRC_MPI_NCPUS_DEFAULT << endl;
    aflowrc << "MPI_NCPUS_MAX=" << AFLOWRC_MPI_NCPUS_MAX << endl;
    
    aflowrc << " " << endl;
    aflowrc << "// DEFAULTS BINARY" << endl;
    aflowrc << "DEFAULT_VASP_GAMMA_BIN=\"" << AFLOWRC_DEFAULT_VASP_GAMMA_BIN << "\"" << endl;
    aflowrc << "DEFAULT_VASP_GAMMA_MPI_BIN=\"" << AFLOWRC_DEFAULT_VASP_GAMMA_MPI_BIN << "\"" << endl;
    aflowrc << "DEFAULT_VASP_BIN=\"" << AFLOWRC_DEFAULT_VASP_BIN << "\"" << endl;
    aflowrc << "DEFAULT_VASP_MPI_BIN=\"" << AFLOWRC_DEFAULT_VASP_MPI_BIN << "\"" << endl;
    aflowrc << "DEFAULT_VASP5_BIN=\"" << AFLOWRC_DEFAULT_VASP5_BIN << "\"" << endl;
    aflowrc << "DEFAULT_VASP5_MPI_BIN=\"" << AFLOWRC_DEFAULT_VASP5_MPI_BIN << "\"" << endl;
    aflowrc << "DEFAULT_AIMS_BIN=\"" << AFLOWRC_DEFAULT_AIMS_BIN << "\"" << endl;
    
    aflowrc << " " << endl;
    aflowrc << "// DEFAULTS POTCARS" << endl;
    aflowrc << "DEFAULT_VASP_POTCAR_DIRECTORIES=\"" << AFLOWRC_DEFAULT_VASP_POTCAR_DIRECTORIES << "\"" << endl;
    aflowrc << "DEFAULT_VASP_POTCAR_DATE=\"" << AFLOWRC_DEFAULT_VASP_POTCAR_DATE << "\"" << endl;
    aflowrc << "DEFAULT_VASP_POTCAR_SUFFIX=\"" << AFLOWRC_DEFAULT_VASP_POTCAR_SUFFIX << "\"" << endl;
    aflowrc << "DEFAULT_VASP_POTCAR_DATE_POT_LDA=\"" << AFLOWRC_DEFAULT_VASP_POTCAR_DATE_POT_LDA << "\"" << endl;
    aflowrc << "DEFAULT_VASP_POTCAR_DATE_POT_GGA=\"" << AFLOWRC_DEFAULT_VASP_POTCAR_DATE_POT_GGA << "\"" << endl;

    aflowrc << " " << endl;
    aflowrc << "// DEFAULTS KPOINTS/DOS" << endl;
    aflowrc << "DEFAULT_BANDS_GRID=" << AFLOWRC_DEFAULT_BANDS_GRID << endl;
    aflowrc << "DEFAULT_BANDS_LATTICE=\"" << AFLOWRC_DEFAULT_BANDS_LATTICE << "\"" << endl;
    aflowrc << "DEFAULT_KSCHEME=\"" << AFLOWRC_DEFAULT_KSCHEME << "\"" << endl;
    aflowrc << "DEFAULT_KPPRA=" << AFLOWRC_DEFAULT_KPPRA << endl;
    aflowrc << "DEFAULT_STATIC_KSCHEME=\"" << AFLOWRC_DEFAULT_STATIC_KSCHEME << "\"" << endl;
    aflowrc << "DEFAULT_KPPRA_STATIC=" << AFLOWRC_DEFAULT_KPPRA_STATIC << endl;
    aflowrc << "DEFAULT_KPPRA_ICSD=" << AFLOWRC_DEFAULT_KPPRA_ICSD << endl;
    aflowrc << "DEFAULT_UNARY_BANDS_GRID=" << AFLOWRC_DEFAULT_UNARY_BANDS_GRID << endl;
    aflowrc << "DEFAULT_UNARY_KPPRA=" << AFLOWRC_DEFAULT_UNARY_KPPRA << " // 32768 // 27000" << endl;
    aflowrc << "DEFAULT_UNARY_KPPRA_STATIC=" << AFLOWRC_DEFAULT_UNARY_KPPRA_STATIC << "// 32768 // 27000" << endl;
    aflowrc << "DEFAULT_PHONONS_KSCHEME=\"" << AFLOWRC_DEFAULT_PHONONS_KSCHEME << "\"" << endl;
    aflowrc << "DEFAULT_DOS_EMIN=" << AFLOWRC_DEFAULT_DOS_EMIN << endl;
    aflowrc << "DEFAULT_DOS_EMAX=" << AFLOWRC_DEFAULT_DOS_EMAX << endl;
    aflowrc << "DEFAULT_DOS_SCALE=" << AFLOWRC_DEFAULT_DOS_SCALE << endl;

    aflowrc << " " << endl;
    aflowrc << "// DEFAULTS PRECISION" << endl;
    aflowrc << "DEFAULT_VASP_PREC_ENMAX_LOW=" << AFLOWRC_DEFAULT_VASP_PREC_ENMAX_LOW << endl;
    aflowrc << "DEFAULT_VASP_PREC_ENMAX_MEDIUM=" << AFLOWRC_DEFAULT_VASP_PREC_ENMAX_MEDIUM << endl;
    aflowrc << "DEFAULT_VASP_PREC_ENMAX_NORMAL=" << AFLOWRC_DEFAULT_VASP_PREC_ENMAX_NORMAL << endl;
    aflowrc << "DEFAULT_VASP_PREC_ENMAX_HIGH=" << AFLOWRC_DEFAULT_VASP_PREC_ENMAX_HIGH << endl;
    aflowrc << "DEFAULT_VASP_PREC_ENMAX_ACCURATE=" << AFLOWRC_DEFAULT_VASP_PREC_ENMAX_ACCURATE << endl;
    aflowrc << "DEFAULT_VASP_SPIN_REMOVE_CUTOFF=" << AFLOWRC_DEFAULT_VASP_SPIN_REMOVE_CUTOFF << endl;
    aflowrc << "DEFAULT_VASP_PREC_POTIM=" << AFLOWRC_DEFAULT_VASP_PREC_POTIM << endl;
    aflowrc << "DEFAULT_VASP_PREC_EDIFFG=" << AFLOWRC_DEFAULT_VASP_PREC_EDIFFG << endl;

    aflowrc << " " << endl;
    aflowrc << "// DEFAULTS OPTIONS " << endl;
    aflowrc << "DEFAULT_VASP_EXTERNAL_INCAR=" << AFLOWRC_DEFAULT_VASP_EXTERNAL_INCAR << endl;
    aflowrc << "DEFAULT_VASP_EXTERNAL_POSCAR=" << AFLOWRC_DEFAULT_VASP_EXTERNAL_POSCAR << endl;
    aflowrc << "DEFAULT_VASP_EXTERNAL_POTCAR=" << AFLOWRC_DEFAULT_VASP_EXTERNAL_POTCAR << endl;
    aflowrc << "DEFAULT_VASP_EXTERNAL_KPOINTS=" << AFLOWRC_DEFAULT_VASP_EXTERNAL_KPOINTS << endl;
    aflowrc << "DEFAULT_VASP_PSEUDOPOTENTIAL_TYPE=" << AFLOWRC_DEFAULT_VASP_PSEUDOPOTENTIAL_TYPE << endl;
    aflowrc << "DEFAULT_VASP_FORCE_OPTION_RELAX_MODE_SCHEME=" << AFLOWRC_DEFAULT_VASP_FORCE_OPTION_RELAX_MODE_SCHEME << endl;
    aflowrc << "DEFAULT_VASP_FORCE_OPTION_PREC_SCHEME=" << AFLOWRC_DEFAULT_VASP_FORCE_OPTION_PREC_SCHEME << endl;
    aflowrc << "DEFAULT_VASP_FORCE_OPTION_ALGO_SCHEME=" << AFLOWRC_DEFAULT_VASP_FORCE_OPTION_ALGO_SCHEME << endl;
    aflowrc << "DEFAULT_VASP_FORCE_OPTION_METAGGA_SCHEME=" << AFLOWRC_DEFAULT_VASP_FORCE_OPTION_METAGGA_SCHEME << endl;
    aflowrc << "DEFAULT_VASP_FORCE_OPTION_IVDW_SCHEME=" << AFLOWRC_DEFAULT_VASP_FORCE_OPTION_IVDW_SCHEME << endl;
    aflowrc << "DEFAULT_VASP_FORCE_OPTION_TYPE_SCHEME=" << AFLOWRC_DEFAULT_VASP_FORCE_OPTION_TYPE_SCHEME << endl;
    aflowrc << "DEFAULT_VASP_FORCE_OPTION_ABMIX_SCHEME=" << AFLOWRC_DEFAULT_VASP_FORCE_OPTION_ABMIX_SCHEME << endl;
    aflowrc << "DEFAULT_VASP_FORCE_OPTION_SYM=" << AFLOWRC_DEFAULT_VASP_FORCE_OPTION_SYM << endl;
    aflowrc << "DEFAULT_VASP_FORCE_OPTION_SPIN=" << AFLOWRC_DEFAULT_VASP_FORCE_OPTION_SPIN << endl;
    aflowrc << "DEFAULT_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_1=" << AFLOWRC_DEFAULT_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_1 << endl;
    aflowrc << "DEFAULT_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_2=" << AFLOWRC_DEFAULT_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_2 << endl;
    aflowrc << "DEFAULT_VASP_FORCE_OPTION_BADER=" << AFLOWRC_DEFAULT_VASP_FORCE_OPTION_BADER << endl;
    aflowrc << "DEFAULT_VASP_FORCE_OPTION_ELF=" << AFLOWRC_DEFAULT_VASP_FORCE_OPTION_ELF << endl;
    aflowrc << "DEFAULT_VASP_FORCE_OPTION_AUTO_MAGMOM=" << AFLOWRC_DEFAULT_VASP_FORCE_OPTION_AUTO_MAGMOM << endl;
    aflowrc << "DEFAULT_VASP_FORCE_OPTION_WAVECAR=" << AFLOWRC_DEFAULT_VASP_FORCE_OPTION_WAVECAR << endl;
    aflowrc << "DEFAULT_VASP_FORCE_OPTION_CHGCAR=" << AFLOWRC_DEFAULT_VASP_FORCE_OPTION_CHGCAR << endl;
    aflowrc << "DEFAULT_VASP_FORCE_OPTION_LSCOUPLING=" << AFLOWRC_DEFAULT_VASP_FORCE_OPTION_LSCOUPLING << endl;

    aflowrc << " " << endl;
    aflowrc << "// AFLOW_LIBRARY AFLOW_PROJECT" << endl;
    aflowrc << "DEFAULT_AFLOW_LIBRARY_DIRECTORIES=\"" << AFLOWRC_DEFAULT_AFLOW_LIBRARY_DIRECTORIES << "\"" << endl;
    aflowrc << "DEFAULT_AFLOW_PROJECTS_DIRECTORIES=\"" << AFLOWRC_DEFAULT_AFLOW_PROJECTS_DIRECTORIES << "\"" << endl;
    
    aflowrc << " " << endl;
    aflowrc << "// DEFAULT PLATON/FINDSYM" << endl;
    aflowrc << "DEFAULT_PLATON_P_EQUAL=" << AFLOWRC_DEFAULT_PLATON_P_EQUAL << endl;
    aflowrc << "DEFAULT_PLATON_P_EXACT=" << AFLOWRC_DEFAULT_PLATON_P_EXACT << endl;
    aflowrc << "DEFAULT_PLATON_P_ANG=" << AFLOWRC_DEFAULT_PLATON_P_ANG << endl;
    aflowrc << "DEFAULT_PLATON_P_D1=" << AFLOWRC_DEFAULT_PLATON_P_D1 << endl;
    aflowrc << "DEFAULT_PLATON_P_D2=" << AFLOWRC_DEFAULT_PLATON_P_D2 << endl;
    aflowrc << "DEFAULT_PLATON_P_D3=" << AFLOWRC_DEFAULT_PLATON_P_D3 << endl;
    aflowrc << "DEFAULT_FINDSYM_TOL=" << AFLOWRC_DEFAULT_FINDSYM_TOL << endl;

    aflowrc << " " << endl;
    aflowrc << "// DEFAULTS GNUPLOT" << endl;
    aflowrc << "DEFAULT_GNUPLOT_EPS_FONT=\"" << AFLOWRC_DEFAULT_GNUPLOT_EPS_FONT << "\"" << endl;
    aflowrc << "DEFAULT_GNUPLOT_EPS_FONT_BOLD=\"" << AFLOWRC_DEFAULT_GNUPLOT_EPS_FONT_BOLD << "\"" << endl;
    aflowrc << "DEFAULT_GNUPLOT_EPS_FONT_ITALICS=\"" << AFLOWRC_DEFAULT_GNUPLOT_EPS_FONT_ITALICS << "\"" << endl;
    aflowrc << "DEFAULT_GNUPLOT_EPS_FONT_BOLD_ITALICS=\"" << AFLOWRC_DEFAULT_GNUPLOT_EPS_FONT_BOLD_ITALICS << "\"" << endl;
    aflowrc << "DEFAULT_GNUPLOT_PNG_FONT=\"" << AFLOWRC_DEFAULT_GNUPLOT_PNG_FONT << "\"" << endl;
    aflowrc << "DEFAULT_GNUPLOT_PNG_FONT_BOLD=\"" << AFLOWRC_DEFAULT_GNUPLOT_PNG_FONT_BOLD << "\"" << endl;
    aflowrc << "DEFAULT_GNUPLOT_PNG_FONT_ITALICS=\"" << AFLOWRC_DEFAULT_GNUPLOT_PNG_FONT_ITALICS << "\"" << endl;
    aflowrc << "DEFAULT_GNUPLOT_PNG_FONT_BOLD_ITALICS=\"" << AFLOWRC_DEFAULT_GNUPLOT_PNG_FONT_BOLD_ITALICS << "\"" << endl;
    
    aflowrc << " " << endl;
    aflowrc << "// DEFAULTS CHULL" << endl;
    aflowrc << "DEFAULT_CHULL_OUTLIER_MULTIPLIER=" << AFLOWRC_DEFAULT_CHULL_OUTLIER_MULTIPLIER << endl;
    aflowrc << "DEFAULT_CHULL_LATEX_BANNER=" << AFLOWRC_DEFAULT_CHULL_LATEX_BANNER << endl;
    aflowrc << "DEFAULT_CHULL_LATEX_COMPOUNDS_COLUMN=" << AFLOWRC_DEFAULT_CHULL_LATEX_COMPOUNDS_COLUMN << endl;
    aflowrc << "DEFAULT_CHULL_LATEX_COMPOSITION_HEADER=" << AFLOWRC_DEFAULT_CHULL_LATEX_COMPOSITION_HEADER << endl;
    aflowrc << "DEFAULT_CHULL_LATEX_PLOT_ENDPOINTS=" << AFLOWRC_DEFAULT_CHULL_LATEX_PLOT_ENDPOINTS << endl;
    aflowrc << "DEFAULT_CHULL_LATEX_PLOT_OFF_HULL=" << AFLOWRC_DEFAULT_CHULL_LATEX_PLOT_OFF_HULL << endl;
    aflowrc << "DEFAULT_CHULL_LATEX_PLOT_UNSTABLE=" << AFLOWRC_DEFAULT_CHULL_LATEX_PLOT_UNSTABLE << endl;
    aflowrc << "DEFAULT_CHULL_LATEX_FILTER_SCHEME=\"" << AFLOWRC_DEFAULT_CHULL_LATEX_FILTER_SCHEME << "\"" << endl;
    aflowrc << "DEFAULT_CHULL_LATEX_FILTER_VALUE=" << AFLOWRC_DEFAULT_CHULL_LATEX_FILTER_VALUE << endl;
    aflowrc << "DEFAULT_CHULL_LATEX_COLOR_BAR=" << AFLOWRC_DEFAULT_CHULL_LATEX_COLOR_BAR << endl;
    aflowrc << "DEFAULT_CHULL_LATEX_HEAT_MAP=" << AFLOWRC_DEFAULT_CHULL_LATEX_HEAT_MAP << endl;
    aflowrc << "DEFAULT_CHULL_LATEX_COLOR_GRADIENT=" << AFLOWRC_DEFAULT_CHULL_LATEX_COLOR_GRADIENT << endl;
    aflowrc << "DEFAULT_CHULL_LATEX_COLOR_MAP=\"" << AFLOWRC_DEFAULT_CHULL_LATEX_COLOR_MAP << "\"" << endl;
    aflowrc << "DEFAULT_CHULL_LATEX_TERNARY_LABEL_COLOR=\"" << AFLOWRC_DEFAULT_CHULL_LATEX_TERNARY_LABEL_COLOR << "\"" << endl;
    aflowrc << "DEFAULT_CHULL_LATEX_REVERSE_AXIS=" << AFLOWRC_DEFAULT_CHULL_LATEX_REVERSE_AXIS << endl;
    aflowrc << "DEFAULT_CHULL_LATEX_FACET_LINE_DROP_SHADOW=" << AFLOWRC_DEFAULT_CHULL_LATEX_FACET_LINE_DROP_SHADOW << endl;
    aflowrc << "DEFAULT_CHULL_LATEX_LINKS=" << AFLOWRC_DEFAULT_CHULL_LATEX_LINKS << endl;
    aflowrc << "DEFAULT_CHULL_LATEX_LABEL_NAME=\"" << AFLOWRC_DEFAULT_CHULL_LATEX_LABEL_NAME << "\"" << endl;
    aflowrc << "DEFAULT_CHULL_LATEX_META_LABELS=" << AFLOWRC_DEFAULT_CHULL_LATEX_META_LABELS << endl;
    aflowrc << "DEFAULT_CHULL_LATEX_LABELS_OFF_HULL=" << AFLOWRC_DEFAULT_CHULL_LATEX_LABELS_OFF_HULL << endl;
    aflowrc << "DEFAULT_CHULL_LATEX_PLOT_REDUCE_COMPOSITION=" << AFLOWRC_DEFAULT_CHULL_LATEX_PLOT_REDUCE_COMPOSITION << endl;
    aflowrc << "DEFAULT_CHULL_LATEX_HELVETICA_FONT=" << AFLOWRC_DEFAULT_CHULL_LATEX_HELVETICA_FONT << endl;
    aflowrc << "DEFAULT_CHULL_LATEX_FONT_SIZE=\"" << AFLOWRC_DEFAULT_CHULL_LATEX_FONT_SIZE << "\"" << endl;
    aflowrc << "DEFAULT_CHULL_LATEX_ROTATE_LABELS=" << AFLOWRC_DEFAULT_CHULL_LATEX_ROTATE_LABELS << endl;
    aflowrc << "DEFAULT_CHULL_LATEX_BOLD_LABELS=" << AFLOWRC_DEFAULT_CHULL_LATEX_BOLD_LABELS << endl;

    aflowrc << " " << endl;
    aflowrc << "// DEFAULTS CORE" << endl;
    aflowrc << "AFLOW_CORE_TEMPERATURE_BEEP=" << AFLOWRC_AFLOW_CORE_TEMPERATURE_BEEP << " // Celsius" << endl;
    aflowrc << "AFLOW_CORE_TEMPERATURE_HALT=" << AFLOWRC_AFLOW_CORE_TEMPERATURE_HALT << " // Celsius" << endl;
    aflowrc << "AFLOW_CORE_TEMPERATURE_REFRESH=" << AFLOWRC_AFLOW_CORE_TEMPERATURE_REFRESH << " // seconds"  << endl;

    aflowrc << " " << endl;
    aflowrc << "// DEFAULTS MACHINE DEPENDENT MPI" << endl;
    aflowrc << "MPI_OPTIONS_DUKE_BETA_MPICH=\"" << AFLOWRC_MPI_OPTIONS_DUKE_BETA_MPICH << "\"" << "// DUKE_BETA_MPICH" << endl;
    aflowrc << "MPI_COMMAND_DUKE_BETA_MPICH=\"" << AFLOWRC_MPI_COMMAND_DUKE_BETA_MPICH << "\"" << "// DUKE_BETA_MPICH" << endl;
    aflowrc << "MPI_BINARY_DIR_DUKE_BETA_MPICH=\"" << AFLOWRC_MPI_BINARY_DIR_DUKE_BETA_MPICH << "\"" << "// DUKE_BETA_MPICH" << endl; 

    aflowrc << "MPI_OPTIONS_DUKE_BETA_OPENMPI=\"" << AFLOWRC_MPI_OPTIONS_DUKE_BETA_OPENMPI << "\"" << "// DUKE_BETA_OPENMPI" << endl;
    aflowrc << "MPI_COMMAND_DUKE_BETA_OPENMPI=\"" << AFLOWRC_MPI_COMMAND_DUKE_BETA_OPENMPI << "\"" << "// DUKE_BETA_OPENMPI" << endl;
    aflowrc << "MPI_BINARY_DIR_DUKE_BETA_OPENMPI=\"" << AFLOWRC_MPI_BINARY_DIR_DUKE_BETA_OPENMPI << "\"" << "// DUKE_BETA_OPENMPI" << endl; 

    aflowrc << "MPI_OPTIONS_DUKE_MATERIALS=\"" << AFLOWRC_MPI_OPTIONS_DUKE_MATERIALS << "\"" << "// DUKE_MATERIALS" << endl;
    aflowrc << "MPI_COMMAND_DUKE_MATERIALS=\"" << AFLOWRC_MPI_COMMAND_DUKE_MATERIALS << "\"" << "// DUKE_MATERIALS" << endl;
    aflowrc << "MPI_BINARY_DIR_DUKE_MATERIALS=\"" << AFLOWRC_MPI_BINARY_DIR_DUKE_MATERIALS << "\"" << "// DUKE_MATERIALS" << endl; 

    aflowrc << "MPI_OPTIONS_DUKE_AFLOWLIB=\"" << AFLOWRC_MPI_OPTIONS_DUKE_AFLOWLIB << "\"" << "// DUKE_AFLOWLIB" << endl;
    aflowrc << "MPI_COMMAND_DUKE_AFLOWLIB=\"" << AFLOWRC_MPI_COMMAND_DUKE_AFLOWLIB << "\"" << "// DUKE_AFLOWLIB" << endl;
    aflowrc << "MPI_BINARY_DIR_DUKE_AFLOWLIB=\"" << AFLOWRC_MPI_BINARY_DIR_DUKE_AFLOWLIB << "\"" << "// DUKE_AFLOWLIB" << endl; 

    aflowrc << "MPI_OPTIONS_DUKE_HABANA=\"" << AFLOWRC_MPI_OPTIONS_DUKE_HABANA << "\"" << "// DUKE_HABANA" << endl;
    aflowrc << "MPI_COMMAND_DUKE_HABANA=\"" << AFLOWRC_MPI_COMMAND_DUKE_HABANA << "\"" << "// DUKE_HABANA" << endl;
    aflowrc << "MPI_BINARY_DIR_DUKE_HABANA=\"" << AFLOWRC_MPI_BINARY_DIR_DUKE_HABANA << "\"" << "// DUKE_HABANA" << endl; 

    aflowrc << "MPI_OPTIONS_DUKE_QRATS_MPICH=\"" << AFLOWRC_MPI_OPTIONS_DUKE_QRATS_MPICH << "\"" << "// DUKE_QRATS_MPICH" << endl;
    aflowrc << "MPI_COMMAND_DUKE_QRATS_MPICH=\"" << AFLOWRC_MPI_COMMAND_DUKE_QRATS_MPICH << "\"" << "// DUKE_QRATS_MPICH" << endl;
    aflowrc << "MPI_BINARY_DIR_DUKE_QRATS_MPICH=\"" << AFLOWRC_MPI_BINARY_DIR_DUKE_QRATS_MPICH << "\"" << "// DUKE_QRATS_MPICH" << endl; 

    aflowrc << "MPI_OPTIONS_DUKE_QUSER_OPENMPI=\"" << AFLOWRC_MPI_OPTIONS_DUKE_QUSER_OPENMPI << "\"" << "// DUKE_QUSER_OPENMPI" << endl;
    aflowrc << "MPI_COMMAND_DUKE_QUSER_OPENMPI=\"" << AFLOWRC_MPI_COMMAND_DUKE_QUSER_OPENMPI << "\"" << "// DUKE_QUSER_OPENMPI" << endl;
    aflowrc << "MPI_BINARY_DIR_DUKE_QUSER_OPENMPI=\"" << AFLOWRC_MPI_BINARY_DIR_DUKE_QUSER_OPENMPI << "\"" << "// DUKE_QUSER_OPENMPI" << endl; 

    aflowrc << "MPI_OPTIONS_MPCDF_EOS_MPIIFORT=\"" << AFLOWRC_MPI_OPTIONS_MPCDF_EOS_MPIIFORT << "\"" << "// MPCDF_EOS_MPIIFORT" << endl;
    aflowrc << "MPI_COMMAND_MPCDF_EOS_MPIIFORT=\"" << AFLOWRC_MPI_COMMAND_MPCDF_EOS_MPIIFORT << "\"" << "// MPCDF_EOS_MPIIFORT" << endl;
    aflowrc << "MPI_BINARY_DIR_MPCDF_EOS_MPIIFORT=\"" << AFLOWRC_MPI_BINARY_DIR_MPCDF_EOS_MPIIFORT << "\"" << "// MPCDF_EOS_MPIIFORT" << endl; 

    aflowrc << "MPI_OPTIONS_MPCDF_DRACO_MPIIFORT=\"" << AFLOWRC_MPI_OPTIONS_MPCDF_DRACO_MPIIFORT << "\"" << "// MPCDF_DRACO_MPIIFORT" << endl;
    aflowrc << "MPI_COMMAND_MPCDF_DRACO_MPIIFORT=\"" << AFLOWRC_MPI_COMMAND_MPCDF_DRACO_MPIIFORT << "\"" << "// MPCDF_DRACO_MPIIFORT" << endl;
    aflowrc << "MPI_BINARY_DIR_MPCDF_DRACO_MPIIFORT=\"" << AFLOWRC_MPI_BINARY_DIR_MPCDF_DRACO_MPIIFORT << "\"" << "// MPCDF_DRACO_MPIIFORT" << endl; 

    aflowrc << "MPI_OPTIONS_MPCDF_HYDRA_MPIIFORT=\"" << AFLOWRC_MPI_OPTIONS_MPCDF_HYDRA_MPIIFORT << "\"" << "// MPCDF_HYDRA_MPIIFORT" << endl;
    aflowrc << "MPI_COMMAND_MPCDF_HYDRA_MPIIFORT=\"" << AFLOWRC_MPI_COMMAND_MPCDF_HYDRA_MPIIFORT << "\"" << "// MPCDF_HYDRA_MPIIFORT" << endl;
    aflowrc << "MPI_BINARY_DIR_MPCDF_HYDRA_MPIIFORT=\"" << AFLOWRC_MPI_BINARY_DIR_MPCDF_HYDRA_MPIIFORT << "\"" << "// MPCDF_HYDRA_MPIIFORT" << endl; 

    aflowrc << "MPI_OPTIONS_TERAGRID_RANGER=\"" << AFLOWRC_MPI_OPTIONS_TERAGRID_RANGER << "\"" << "// TERAGRID_RANGER" << endl;
    aflowrc << "MPI_COMMAND_TERAGRID_RANGER=\"" << AFLOWRC_MPI_COMMAND_TERAGRID_RANGER << "\"" << "// TERAGRID_RANGER" << endl;
    aflowrc << "MPI_BINARY_DIR_TERAGRID_RANGER=\"" << AFLOWRC_MPI_BINARY_DIR_TERAGRID_RANGER << "\"" << "// TERAGRID_RANGER" << endl; 

    aflowrc << "MPI_OPTIONS_TERAGRID_KRAKEN=\"" << AFLOWRC_MPI_OPTIONS_TERAGRID_KRAKEN << "\"" << "// TERAGRID_KRAKEN" << endl;
    aflowrc << "MPI_COMMAND_TERAGRID_KRAKEN=\"" << AFLOWRC_MPI_COMMAND_TERAGRID_KRAKEN << "\"" << "// TERAGRID_KRAKEN" << endl;
    aflowrc << "MPI_BINARY_DIR_TERAGRID_KRAKEN=\"" << AFLOWRC_MPI_BINARY_DIR_TERAGRID_KRAKEN << "\"" << "// TERAGRID_KRAKEN" << endl; 

    aflowrc << "MPI_OPTIONS_FULTON_MARYLOU=\"" << AFLOWRC_MPI_OPTIONS_FULTON_MARYLOU << "\"" << "// FULTON_MARYLOU" << endl;
    aflowrc << "MPI_COMMAND_FULTON_MARYLOU=\"" << AFLOWRC_MPI_COMMAND_FULTON_MARYLOU << "\"" << "// FULTON_MARYLOU" << endl;
    aflowrc << "MPI_BINARY_DIR_FULTON_MARYLOU=\"" << AFLOWRC_MPI_BINARY_DIR_FULTON_MARYLOU << "\"" << "// FULTON_MARYLOU" << endl; 

    aflowrc << "MPI_OPTIONS_TRINITY_PARSONS=\"" << AFLOWRC_MPI_OPTIONS_TRINITY_PARSONS << "\"" << "// TRINITY_PARSONS" << endl;
    aflowrc << "MPI_COMMAND_TRINITY_PARSONS=\"" << AFLOWRC_MPI_COMMAND_TRINITY_PARSONS << "\"" << "// TRINITY_PARSONS" << endl;
    aflowrc << "MPI_BINARY_DIR_TRINITY_PARSONS=\"" << AFLOWRC_MPI_BINARY_DIR_TRINITY_PARSONS << "\"" << "// TRINITY_PARSONS" << endl; 

    aflowrc << "MPI_OPTIONS_MACHINE1=\"" << AFLOWRC_MPI_OPTIONS_MACHINE1 << "\"" << "// MACHINE1" << endl;
    aflowrc << "MPI_COMMAND_MACHINE1=\"" << AFLOWRC_MPI_COMMAND_MACHINE1 << "\"" << "// MACHINE1" << endl;
    aflowrc << "MPI_BINARY_DIR_MACHINE1=\"" << AFLOWRC_MPI_BINARY_DIR_MACHINE1 << "\"" << "// MACHINE1" << endl; 

    aflowrc << "MPI_OPTIONS_MACHINE2=\"" << AFLOWRC_MPI_OPTIONS_MACHINE2 << "\"" << "// MACHINE2" << endl;
    aflowrc << "MPI_COMMAND_MACHINE2=\"" << AFLOWRC_MPI_COMMAND_MACHINE2 << "\"" << "// MACHINE2" << endl;
    aflowrc << "MPI_BINARY_DIR_MACHINE2=\"" << AFLOWRC_MPI_BINARY_DIR_MACHINE2 << "\"" << "// MACHINE2" << endl; 

   
    aflowrc << " " << endl;
    aflowrc << "// ****************************************************************************************************" << endl;
  
    //   XHOST.DEBUG=TRUE;
    cerr << "WARNING: aflowrc::write_default: WRITING default " << XHOST.aflowrc_filename << endl;
    aurostd::stringstream2file(aflowrc,XHOST.aflowrc_filename);
    if(LDEBUG) oss << "aflowrc::write_default: END" << endl;
    //    exit(0);
    return TRUE;
  }
} // namespace aflowrc

// ***************************************************************************
// aflowrc::print_aflowrc
// ***************************************************************************
namespace aflowrc {
  bool print_aflowrc(std::ostream& oss,bool AFLOWRC_VERBOSE) {
    bool LDEBUG=(FALSE || XHOST.DEBUG || AFLOWRC_VERBOSE);   
    if(LDEBUG) oss << "aflowrc::print_aflowrc: BEGIN" << endl;
    if(LDEBUG) oss << "aflowrc::print_aflowrc: XHOST.Home=" << XHOST.Home << endl;
    if(LDEBUG) oss << "aflowrc::print_aflowrc: XHOST.aflowrc_filename=" << XHOST.aflowrc_filename << endl;

    if(LDEBUG) oss << "// DEFAULT GENERIC MPI" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_START_DEFAULT\")=\"" << MPI_START_DEFAULT <<"\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_STOP_DEFAULT\")=\"" << MPI_STOP_DEFAULT <<"\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_COMMAND_DEFAULT\")=\"" << MPI_COMMAND_DEFAULT <<"\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<int>(\"MPI_NCPUS_DEFAULT\")=" << MPI_NCPUS_DEFAULT << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<int>(\"MPI_NCPUS_MAX\")=" << MPI_NCPUS_MAX << endl;

    if(LDEBUG) oss << "// DEFAULTS BINARY" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_VASP_GAMMA_BIN\")=\"" << DEFAULT_VASP_GAMMA_BIN <<"\""  << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_VASP_GAMMA_MPI_BIN\")=\"" << DEFAULT_VASP_GAMMA_MPI_BIN <<"\""  << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_VASP_BIN\")=\"" << DEFAULT_VASP_BIN <<"\""  << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_VASP_MPI_BIN\")=\"" << DEFAULT_VASP_MPI_BIN <<"\""  << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_VASP5_BIN\")=\"" << DEFAULT_VASP5_BIN <<"\""  << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_VASP5_MPI_BIN\")=\"" << DEFAULT_VASP5_MPI_BIN <<"\""  << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AIMS_BIN\")=\"" << DEFAULT_AIMS_BIN <<"\""  << endl;
 
    if(LDEBUG) oss << "// DEFAULTS POTCARS" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_VASP_POTCAR_DIRECTORIES\")=\"" << DEFAULT_VASP_POTCAR_DIRECTORIES <<"\""  << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_VASP_POTCAR_DATE\")=\"" << DEFAULT_VASP_POTCAR_DATE <<"\""  << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_VASP_POTCAR_SUFFIX\")=\"" << DEFAULT_VASP_POTCAR_SUFFIX <<"\""  << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_VASP_POTCAR_DATE_POT_LDA\")=\"" << DEFAULT_VASP_POTCAR_DATE_POT_LDA <<"\""  << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_VASP_POTCAR_DATE_POT_GGA\")=\"" << DEFAULT_VASP_POTCAR_DATE_POT_GGA <<"\""  << endl;

    if(LDEBUG) oss << "// DEFAULT KPOINTS/DOS" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<int>(\"DEFAULT_BANDS_GRID\")=" << DEFAULT_BANDS_GRID << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_BANDS_LATTICE\")=\"" << DEFAULT_BANDS_LATTICE <<"\""  << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_KSCHEME\")=\"" << DEFAULT_KSCHEME <<"\""  << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<int>(\"DEFAULT_KPPRA\")=" << DEFAULT_KPPRA << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_STATIC_KSCHEME\")=\"" << DEFAULT_STATIC_KSCHEME <<"\""  << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<int>(\"DEFAULT_KPPRA_STATIC\")=" << DEFAULT_KPPRA_STATIC << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<int>(\"DEFAULT_KPPRA_ICSD\")=" << DEFAULT_KPPRA_ICSD << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<int>(\"DEFAULT_UNARY_BANDS_GRID\")=" << DEFAULT_UNARY_BANDS_GRID << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<int>(\"DEFAULT_UNARY_KPPRA\")=" << DEFAULT_UNARY_KPPRA << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<int>(\"DEFAULT_UNARY_KPPRA_STATIC\")=" << DEFAULT_UNARY_KPPRA_STATIC << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_PHONONS_KSCHEME\")=\"" << DEFAULT_PHONONS_KSCHEME <<"\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<double>(\"DEFAULT_DOS_EMIN\")=" << DEFAULT_DOS_EMIN << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<double>(\"DEFAULT_DOS_EMAX\")=" << DEFAULT_DOS_EMAX << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<double>(\"DEFAULT_DOS_SCALE\")=" << DEFAULT_DOS_SCALE << endl;

    if(LDEBUG) oss << "// DEFAULT PRECISION" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<double>(\"DEFAULT_VASP_PREC_ENMAX_LOW\")=" << DEFAULT_VASP_PREC_ENMAX_LOW << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<double>(\"DEFAULT_VASP_PREC_ENMAX_MEDIUM\")=" << DEFAULT_VASP_PREC_ENMAX_MEDIUM << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<double>(\"DEFAULT_VASP_PREC_ENMAX_NORMAL\")=" << DEFAULT_VASP_PREC_ENMAX_NORMAL << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<double>(\"DEFAULT_VASP_PREC_ENMAX_HIGH\")=" << DEFAULT_VASP_PREC_ENMAX_HIGH << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<double>(\"DEFAULT_VASP_PREC_ENMAX_ACCURATE\")=" << DEFAULT_VASP_PREC_ENMAX_ACCURATE << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<double>(\"DEFAULT_VASP_SPIN_REMOVE_CUTOFF\")=" << DEFAULT_VASP_SPIN_REMOVE_CUTOFF << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<double>(\"DEFAULT_VASP_PREC_POTIM\")=" << DEFAULT_VASP_PREC_POTIM << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<double>(\"DEFAULT_VASP_PREC_EDIFFG\")=" << DEFAULT_VASP_PREC_EDIFFG << endl;

    if(LDEBUG) oss << "// DEFAULTS OPTIONS " << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_VASP_EXTERNAL_INCAR\")=\"" << DEFAULT_VASP_EXTERNAL_INCAR <<"\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_VASP_EXTERNAL_POSCAR\")=\"" << DEFAULT_VASP_EXTERNAL_POSCAR <<"\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_VASP_EXTERNAL_POTCAR\")=\"" << DEFAULT_VASP_EXTERNAL_POTCAR <<"\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_VASP_EXTERNAL_KPOINT\")=\"" << DEFAULT_VASP_EXTERNAL_KPOINTS <<"\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_VASP_PSEUDOPOTENTIAL_TYPE\")=\"" << DEFAULT_VASP_PSEUDOPOTENTIAL_TYPE <<"\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_VASP_FORCE_OPTION_RELAX_MODE_SCHEME\")=\"" << DEFAULT_VASP_FORCE_OPTION_RELAX_MODE_SCHEME <<"\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_VASP_FORCE_OPTION_PREC_SCHEME\")=\"" << DEFAULT_VASP_FORCE_OPTION_PREC_SCHEME <<"\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_VASP_FORCE_OPTION_ALGO_SCHEME\")=\"" << DEFAULT_VASP_FORCE_OPTION_ALGO_SCHEME <<"\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_VASP_FORCE_OPTION_METAGGA_SCHEME\")=\"" << DEFAULT_VASP_FORCE_OPTION_METAGGA_SCHEME <<"\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_VASP_FORCE_OPTION_IVDW_SCHEME\")=\"" << DEFAULT_VASP_FORCE_OPTION_IVDW_SCHEME <<"\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_VASP_FORCE_OPTION_TYPE_SCHEME\")=\"" << DEFAULT_VASP_FORCE_OPTION_TYPE_SCHEME <<"\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_VASP_FORCE_OPTION_ABMIX_SCHEME\")=\"" << DEFAULT_VASP_FORCE_OPTION_ABMIX_SCHEME <<"\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<bool>(\"DEFAULT_VASP_FORCE_OPTION_SYM\")=" << DEFAULT_VASP_FORCE_OPTION_SYM << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<bool>(\"DEFAULT_VASP_FORCE_OPTION_SPIN\")=" << DEFAULT_VASP_FORCE_OPTION_SPIN << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<bool>(\"DEFAULT_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_1\")=" << DEFAULT_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_1 << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<bool>(\"DEFAULT_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_2\")=" << DEFAULT_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_2 << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<bool>(\"DEFAULT_VASP_FORCE_OPTION_BADER\")=" << DEFAULT_VASP_FORCE_OPTION_BADER << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<bool>(\"DEFAULT_VASP_FORCE_OPTION_ELF\")=" << DEFAULT_VASP_FORCE_OPTION_ELF << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<bool>(\"DEFAULT_VASP_FORCE_OPTION_AUTO_MAGMOM\")=" << DEFAULT_VASP_FORCE_OPTION_AUTO_MAGMOM << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<bool>(\"DEFAULT_VASP_FORCE_OPTION_WAVECAR\")=" << DEFAULT_VASP_FORCE_OPTION_WAVECAR << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<bool>(\"DEFAULT_VASP_FORCE_OPTION_CHGCAR\")=" << DEFAULT_VASP_FORCE_OPTION_CHGCAR << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<bool>(\"DEFAULT_VASP_FORCE_OPTION_LSCOUPLING\")=" << DEFAULT_VASP_FORCE_OPTION_LSCOUPLING << endl;

    if(LDEBUG) oss << "// AFLOW_LIBRARY AFLOW_PROJECT" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AFLOW_LIBRARY_DIRECTORIES\")=\"" << DEFAULT_AFLOW_LIBRARY_DIRECTORIES <<"\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AFLOW_PROJECTS_DIRECTORIES\")=\"" << DEFAULT_AFLOW_PROJECTS_DIRECTORIES <<"\"" << endl;
    
    if(LDEBUG) oss << "// DEFAULT PLATON/FINDSYM" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<bool>(\"DEFAULT_PLATON_P_EQUAL\")=" << DEFAULT_PLATON_P_EQUAL << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<bool>(\"DEFAULT_PLATON_P_EXACT\")=" << DEFAULT_PLATON_P_EXACT << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<double>(\"DEFAULT_PLATON_P_ANG\")=" << DEFAULT_PLATON_P_ANG << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<double>(\"DEFAULT_PLATON_P_D1\")=" << DEFAULT_PLATON_P_D1 << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<double>(\"DEFAULT_PLATON_P_D2\")=" << DEFAULT_PLATON_P_D2 << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<double>(\"DEFAULT_PLATON_P_D3\")=" << DEFAULT_PLATON_P_D3 << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<double>(\"DEFAULT_FINDSYM_TOL\")=" << DEFAULT_FINDSYM_TOL << endl;

    if(LDEBUG) oss << "// DEFAULT GNUPLOT" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_GNUPLOT_EPS_FONT\")=\"" << DEFAULT_GNUPLOT_EPS_FONT <<"\""  << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_GNUPLOT_EPS_FONT_BOLD\")=\"" << DEFAULT_GNUPLOT_EPS_FONT_BOLD <<"\""  << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_GNUPLOT_EPS_FONT_ITALICS\")=\"" << DEFAULT_GNUPLOT_EPS_FONT_ITALICS <<"\""  << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_GNUPLOT_EPS_FONT_BOLD_ITALICS\")=\"" << DEFAULT_GNUPLOT_EPS_FONT_BOLD_ITALICS <<"\""  << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_GNUPLOT_PNG_FONT\")=\"" << DEFAULT_GNUPLOT_PNG_FONT <<"\""  << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_GNUPLOT_PNG_FONT_BOLD\")=\"" << DEFAULT_GNUPLOT_PNG_FONT_BOLD <<"\""  << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_GNUPLOT_PNG_FONT_ITALICS\")=\"" << DEFAULT_GNUPLOT_PNG_FONT_ITALICS <<"\""  << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_GNUPLOT_PNG_FONT_BOLD_ITALICS\")=\"" << DEFAULT_GNUPLOT_PNG_FONT_BOLD_ITALICS <<"\""  << endl;

    if(LDEBUG) oss << "// DEFAULT CHULL" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<double>(\"DEFAULT_CHULL_OUTLIER_MULTIPLIER\")=" << DEFAULT_CHULL_OUTLIER_MULTIPLIER << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<double>(\"DEFAULT_CHULL_LATEX_BANNER\")=" << DEFAULT_CHULL_LATEX_BANNER << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<double>(\"DEFAULT_CHULL_LATEX_COMPOUNDS_COLUMN\")=" << DEFAULT_CHULL_LATEX_COMPOUNDS_COLUMN << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<double>(\"DEFAULT_CHULL_LATEX_COMPOSITION_HEADER\")=" << DEFAULT_CHULL_LATEX_COMPOSITION_HEADER << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<double>(\"DEFAULT_CHULL_LATEX_PLOT_ENDPOINTS\")=" << DEFAULT_CHULL_LATEX_PLOT_ENDPOINTS << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<double>(\"DEFAULT_CHULL_LATEX_PLOT_OFF_HULL\")=" << DEFAULT_CHULL_LATEX_PLOT_OFF_HULL << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<double>(\"DEFAULT_CHULL_LATEX_PLOT_UNSTABLE\")=" << DEFAULT_CHULL_LATEX_PLOT_UNSTABLE << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<double>(\"DEFAULT_CHULL_LATEX_FILTER_SCHEME\")=\"" << DEFAULT_CHULL_LATEX_FILTER_SCHEME << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<double>(\"DEFAULT_CHULL_LATEX_COLOR_BAR\")=" << DEFAULT_CHULL_LATEX_COLOR_BAR << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<double>(\"DEFAULT_CHULL_LATEX_HEAT_MAP\")=" << DEFAULT_CHULL_LATEX_HEAT_MAP << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<double>(\"DEFAULT_CHULL_LATEX_COLOR_GRADIENT\")=" << DEFAULT_CHULL_LATEX_COLOR_GRADIENT << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<double>(\"DEFAULT_CHULL_LATEX_COLOR_MAP\")=\"" << DEFAULT_CHULL_LATEX_COLOR_MAP << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<double>(\"DEFAULT_CHULL_LATEX_TERNARY_LABEL_COLOR\")=\"" << DEFAULT_CHULL_LATEX_TERNARY_LABEL_COLOR << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<double>(\"DEFAULT_CHULL_LATEX_REVERSE_AXIS\")=" << DEFAULT_CHULL_LATEX_REVERSE_AXIS << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<double>(\"DEFAULT_CHULL_LATEX_FACET_LINE_DROP_SHADOW\")=" << DEFAULT_CHULL_LATEX_FACET_LINE_DROP_SHADOW << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<double>(\"DEFAULT_CHULL_LATEX_LINKS\")=" << DEFAULT_CHULL_LATEX_LINKS << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<double>(\"DEFAULT_CHULL_LATEX_LABEL_NAME\")=\"" << DEFAULT_CHULL_LATEX_LABEL_NAME << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<double>(\"DEFAULT_CHULL_LATEX_META_LABELS\")=" << DEFAULT_CHULL_LATEX_META_LABELS << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<double>(\"DEFAULT_CHULL_LATEX_LABELS_OFF_HULL\")=" << DEFAULT_CHULL_LATEX_LABELS_OFF_HULL << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<double>(\"DEFAULT_CHULL_LATEX_PLOT_REDUCE_COMPOSITION\")=" << DEFAULT_CHULL_LATEX_PLOT_REDUCE_COMPOSITION << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<double>(\"DEFAULT_CHULL_LATEX_HELVETICA_FONT\")=" << DEFAULT_CHULL_LATEX_HELVETICA_FONT << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<double>(\"DEFAULT_CHULL_LATEX_FONT_SIZE\")=\"" << DEFAULT_CHULL_LATEX_FONT_SIZE << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<double>(\"DEFAULT_CHULL_LATEX_ROTATE_LABELS\")=" << DEFAULT_CHULL_LATEX_ROTATE_LABELS << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<double>(\"DEFAULT_CHULL_LATEX_BOLD_LABELS\")=" << DEFAULT_CHULL_LATEX_BOLD_LABELS << endl;

    if(LDEBUG) oss << "// DEFAULT CORE" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<double>(\"AFLOW_CORE_TEMPERATURE_BEEP\")=" << AFLOW_CORE_TEMPERATURE_BEEP << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<double>(\"AFLOW_CORE_TEMPERATURE_HALT\")=" << AFLOW_CORE_TEMPERATURE_HALT << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<double>(\"AFLOW_CORE_TEMPERATURE_REFRESH\")=" << AFLOW_CORE_TEMPERATURE_REFRESH << endl;

    if(LDEBUG) oss << "// DEFAULT MACHINE DEPENDENT MPI" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_OPTIONS_DUKE_BETA_MPICH\")=\"" << MPI_OPTIONS_DUKE_BETA_MPICH << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_COMMAND_DUKE_BETA_MPICH\")=\"" << MPI_COMMAND_DUKE_BETA_MPICH <<"\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_BINARY_DIR_DUKE_BETA_MPICH\")=\"" << MPI_BINARY_DIR_DUKE_BETA_MPICH <<"\"" << endl;

    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_OPTIONS_DUKE_BETA_OPENMPI\")=\"" << MPI_OPTIONS_DUKE_BETA_OPENMPI <<"\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_COMMAND_DUKE_BETA_OPENMPI\")=\"" << MPI_COMMAND_DUKE_BETA_OPENMPI <<"\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_BINARY_DIR_DUKE_BETA_OPENMPI\")=\"" << MPI_BINARY_DIR_DUKE_BETA_OPENMPI <<"\"" << endl;

    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_OPTIONS_DUKE_MATERIALS\")=\"" << MPI_OPTIONS_DUKE_MATERIALS <<"\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_COMMAND_DUKE_MATERIALS\")=\"" << MPI_COMMAND_DUKE_MATERIALS <<"\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_BINARY_DIR_DUKE_MATERIALS\")=\"" << MPI_BINARY_DIR_DUKE_MATERIALS <<"\"" << endl;

    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_OPTIONS_DUKE_AFLOWLIB\")=\"" << MPI_OPTIONS_DUKE_AFLOWLIB <<"\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_COMMAND_DUKE_AFLOWLIB\")=\"" << MPI_COMMAND_DUKE_AFLOWLIB <<"\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_BINARY_DIR_DUKE_AFLOWLIB\")=\"" << MPI_BINARY_DIR_DUKE_AFLOWLIB <<"\"" << endl;

    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_OPTIONS_DUKE_HABANA\")=\"" << MPI_OPTIONS_DUKE_HABANA <<"\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_COMMAND_DUKE_HABANA\")=\"" << MPI_COMMAND_DUKE_HABANA <<"\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_BINARY_DIR_DUKE_HABANA\")=\"" << MPI_BINARY_DIR_DUKE_HABANA <<"\"" << endl;

    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_OPTIONS_DUKE_QRATS_MPICH\")=\"" << MPI_OPTIONS_DUKE_QRATS_MPICH <<"\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_COMMAND_DUKE_QRATS_MPICH\")=\"" << MPI_COMMAND_DUKE_QRATS_MPICH <<"\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_BINARY_DIR_DUKE_QRATS_MPICH\")=\"" << MPI_BINARY_DIR_DUKE_QRATS_MPICH <<"\"" << endl;

    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_OPTIONS_DUKE_QUSER_OPENMPI\")=\"" << MPI_OPTIONS_DUKE_QUSER_OPENMPI <<"\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_COMMAND_DUKE_QUSER_OPENMPI\")=\"" << MPI_COMMAND_DUKE_QUSER_OPENMPI <<"\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_BINARY_DIR_DUKE_QUSER_OPENMPI\")=\"" << MPI_BINARY_DIR_DUKE_QUSER_OPENMPI <<"\"" << endl;

    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_OPTIONS_MPCDF_EOS_MPIIFORT\")=\"" << MPI_OPTIONS_MPCDF_EOS_MPIIFORT <<"\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_COMMAND_MPCDF_EOS_MPIIFORT\")=\"" << MPI_COMMAND_MPCDF_EOS_MPIIFORT <<"\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_BINARY_DIR_MPCDF_EOS_MPIIFORT\")=\"" << MPI_BINARY_DIR_MPCDF_EOS_MPIIFORT <<"\"" << endl;

    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_OPTIONS_MPCDF_DRACO_MPIIFORT\")=\"" << MPI_OPTIONS_MPCDF_DRACO_MPIIFORT <<"\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_COMMAND_MPCDF_DRACO_MPIIFORT\")=\"" << MPI_COMMAND_MPCDF_DRACO_MPIIFORT <<"\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_BINARY_DIR_MPCDF_DRACO_MPIIFORT\")=\"" << MPI_BINARY_DIR_MPCDF_DRACO_MPIIFORT <<"\"" << endl;

    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_OPTIONS_MPCDF_HYDRA_MPIIFORT\")=\"" << MPI_OPTIONS_MPCDF_HYDRA_MPIIFORT <<"\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_COMMAND_MPCDF_HYDRA_MPIIFORT\")=\"" << MPI_COMMAND_MPCDF_HYDRA_MPIIFORT <<"\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_BINARY_DIR_MPCDF_HYDRA_MPIIFORT\")=\"" << MPI_BINARY_DIR_MPCDF_HYDRA_MPIIFORT <<"\"" << endl;

    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_OPTIONS_TERAGRID_RANGER\")=\"" << MPI_OPTIONS_TERAGRID_RANGER <<"\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_COMMAND_TERAGRID_RANGER\")=\"" << MPI_COMMAND_TERAGRID_RANGER <<"\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_BINARY_DIR_TERAGRID_RANGER\")=\"" << MPI_BINARY_DIR_TERAGRID_RANGER <<"\"" << endl;

    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_OPTIONS_TERAGRID_KRAKEN\")=\"" << MPI_OPTIONS_TERAGRID_KRAKEN <<"\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_COMMAND_TERAGRID_KRAKEN\")=\"" << MPI_COMMAND_TERAGRID_KRAKEN <<"\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_BINARY_DIR_TERAGRID_KRAKEN\")=\"" << MPI_BINARY_DIR_TERAGRID_KRAKEN <<"\"" << endl;

    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_OPTIONS_FULTON_MARYLOU\")=\"" << MPI_OPTIONS_FULTON_MARYLOU <<"\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_COMMAND_FULTON_MARYLOU\")=\"" << MPI_COMMAND_FULTON_MARYLOU <<"\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_BINARY_DIR_FULTON_MARYLOU\")=\"" << MPI_BINARY_DIR_FULTON_MARYLOU <<"\"" << endl;

    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_OPTIONS_TRINITY_PARSONS\")=\"" << MPI_OPTIONS_TRINITY_PARSONS <<"\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_COMMAND_TRINITY_PARSONS\")=\"" << MPI_COMMAND_TRINITY_PARSONS <<"\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_BINARY_DIR_TRINITY_PARSONS\")=\"" << MPI_BINARY_DIR_TRINITY_PARSONS <<"\"" << endl;

    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_OPTIONS_MACHINE1\")=\"" << MPI_OPTIONS_MACHINE1 <<"\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_COMMAND_MACHINE1\")=\"" << MPI_COMMAND_MACHINE1 <<"\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_BINARY_DIR_MACHINE1\")=\"" << MPI_BINARY_DIR_MACHINE1 <<"\"" << endl;

    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_OPTIONS_MACHINE2\")=\"" << MPI_OPTIONS_MACHINE2 <<"\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_COMMAND_MACHINE2\")=\"" << MPI_COMMAND_MACHINE2 <<"\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_BINARY_DIR_MACHINE2\")=\"" << MPI_BINARY_DIR_MACHINE2 <<"\"" << endl;

    //   if(LDEBUG) oss << "XHOST.adefault.content=" << XHOST.adefault.content_string << endl;
    
    if(LDEBUG) oss << "aflowrc::print_aflowrc: END" << endl;

    oss.flush();
    //    exit(0);
    return FALSE;
  }
} // namespace aflowrc

// **************************************************************************
// **************************************************************************

#endif // _AFLOW_AFLOWRC_CPP_

// **************************************************************************
// *                                                                        *
// *             STEFANO CURTAROLO - Duke University 2003-2018              *
// *                                                                        *
// **************************************************************************
