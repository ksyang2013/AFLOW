// ***************************************************************************
// *                                                                         *
// *              AFlow KEVIN RASCH - Duke University 2013                   *
// *                                                                         *
// ***************************************************************************
// aflow_bader.h
// functions written by KEVIN RASCH
// 2013: kevin.rasch@duke.edu

#ifndef _AFLOW_BADER_H_
#define _AFLOW_BADER_H_

// ***************************************************************************
#include <string>
#include <vector>

using std::vector;
using std::string;

// perform Bader analysis using the program from Henkelman Group at UT, Austin
namespace bader_functions {
string BaderCalc(aurostd::xoption vpflow);
bool BaderCalc(aurostd::xoption& vpflow, const string& bader_options, string& directory, ostream& oss);
bool BaderCalc(aurostd::xoption& vpflow, const string& bader_options, const string& prototype,
               const deque<string>& vspecies, const deque<int>& num_each_species, const vector<double>& vZVAL,
               const vector<double>& cutoffs, const vector<int>& downsample_ratios, string& directory, ostream& oss);
bool BaderCalc(const string& bader_options, const string& prototype, const deque<string>& vspecies,
               const deque<int>& num_each_species, const vector<double>& vZVAL, const vector<double>& cutoffs,
               const vector<int>& downsample_ratios, string& directory, ostream& oss);
void FixDirectory(string& directory);
bool Flags2BaderCommands(aurostd::xoption& vpflow, string& bader_options, ostream& oss);
bool getPushCommand(const string& misc_option, string& push_command, ostream& oss);
bool listORrange2vec(const string& misc_option, vector<int>& vout, ostream& oss);
bool BaderExtensionFound(const string& FileNameIn, string& FileNameOut, const string& directory);
bool BaderExtensionFound(const string& FileNameIn, const string& directory);
void adjust_header(string& new_header, stringstream& FileIN_ss);
string prepare_CHGCAR_4_Jmol(aurostd::xoption vpflow);
bool get_species_string(string& outcar_file, string& species_string, const string& dir_to_look, const string& file, ostream& oss);
bool prepare_CHGCAR_4_Jmol(const string& _chgcar_file, string& species_header, bool zip_file, ostream& oss);
bool prepare_CHGCAR_4_Jmol(string& _chgcar_file, string& species_header, ostream& oss);
bool prepare_CHGCAR_4_Jmol(string& _chgcar_file, string& species_header);
//[CO 180220 - moved to aurostd]bool efile2tempfile(string _FileNameIn, string& FileNameOut);
}

//CHGCAR2JVXL function
namespace pflow {
string CHGCAR2JVXL(aurostd::xoption& vpflow);
bool CHGCAR2JVXL(vector<string>& chgcar_files, const vector<double>& cutoffs, const vector<int>& downsample_ratios,
                 const bool& cyclic, ostream& oss);
bool CHGCAR2JVXL(vector<string>& chgcar_files, const vector<double>& cutoffs, const vector<int>& downsample_ratios,
                 vector<string>& output_files, const bool& cyclic, ostream& oss);
bool CHGCAR2JVXL(vector<string>& chgcar_files, const vector<double>& cutoffs, const vector<int>& downsample_ratios,
                 ostream& oss);
bool CHGCAR2JVXL(vector<string>& chgcar_files, const vector<double>& cutoffs, const bool& cyclic, ostream& oss);
bool CHGCAR2JVXL(vector<string>& chgcar_files, const vector<double>& cutoffs, ostream& oss);
bool CHGCAR2JVXL(vector<string>& chgcar_files, const vector<double>& cutoffs, const vector<int>& downsample_ratios,
                 vector<string>& output_files, ostringstream& oss);
bool CHGCAR2JVXL(vector<string>& chgcar_files, const vector<double>& cutoffs, vector<string>& output_files,
                 ostringstream& oss);
bool CHGCAR2JVXL(string chgcar_file, const double cutoff, const int downsample_ratio, const string output_file,
                 ostream& oss);
bool CHGCAR2JVXL(string chgcar_file, const double& cutoff, const int& downsample_ratio, ostream& oss);
bool CHGCAR2JVXL(string chgcar_file, const double& cutoff, string& output_file, ostream& oss);
bool CHGCAR2JVXL(string chgcar_file, const double& cutoff, ostream& oss);
string CHGCAR2JVXL_get_output_filename(string chgcar_file, const double& cutoff, const int& downsample_ratio);
string CHGCAR2JVXL_get_output_filename(string chgcar_file, const double& cutoff);
}

// Phonon graphs
namespace pflow {
bool PLOT_PHDISP(vector<string>& argv);
}
// read phonon calculation output and produce graphs
string gnuplotPhononHeader(string directory, string name, string unit, string plot_type, string title);
// return header info for gnuplot .eps graph
string gnuplotPhononDisp(string dis_path, string ylabel, double conversion, vector<string> specialLabels,
                         vector<double> specialPositions, int nbranches, string rightmargin, string botmargin);
// return info for making a plot of the phonon dispersion
string gnuplotPhononDos(string pdos_path, string leftmargin, string botmargin);
// return infor for making a plot of the phonon density of states
string gnuplotPhononThermo(string thermo_path, string topmargin);
// return info for making a plot of the thermodynamics properties F_vib, S_vib, C_v
void callOutputConversion(bool PNG, bool PDF, bool JPG, bool GIF, string directory, string name, string unit, string plot);
// executes `convert' to produce requested formats from .eps
void cleanGnuplotScript(bool file1, bool file2, bool file3, string script);
// removes the gnuplot scripts because they refer to temporary files

#endif

// ***************************************************************************
// *                                                                         *
// *              AFlow KEVIN RASCH - Duke University 2013                   *
// *                                                                         *
// ***************************************************************************
