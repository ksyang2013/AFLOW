#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <stdexcept>

#include "xmlib.h"
#include <unistd.h>

using namespace xmlib;
using namespace std;

#define numberX 1.267381
//#define numberX 0.90

enum IPCFreqFlags {
    NONE           = 0L,
    ALLOW_NEGATIVE = 1L << 1,
    OMEGA          = 1L << 2,
    RAW            = 1L << 3, // eV/A/A/atomic_mass_unit
    HERTZ          = 1L << 4,
    THZ            = 1L << 5,
    RECIPROCAL_CM  = 1L << 6,
    MEV            = 1L << 7
};

inline IPCFreqFlags operator&(const IPCFreqFlags& __a, const IPCFreqFlags& __b) {
    return IPCFreqFlags( static_cast<int>(__a) & static_cast<int>(__b) );
}

void tokenize(const string&,vector<string>&, string);

int main(int argc, char **argv) {    
  ifstream infile;
  
  double factorTHz2MeV = 4.1356673310;
  double factorRCM2MeV = 4.1356673310 * 2.99792458 * 1E-2;
  
  // Arguments from command line ////////////////////////////////////////

    bool showExactQPoints = false;
    bool userSetGlobalYMin = false;
    double globalYMin = 0.0;
    bool userSetGlobalYMax = false;
    double globalYMax = 0.0;
    bool userPrintPS = false;
    bool userPrintEPS = false;
    bool userPrintPDF = false;
    bool userPrintPNG = false;
 
    if( (argc == 2) && ( ( string(argv[1]) == string("-h") ) ||
                         ( string(argv[1]) == string("--help") ) ) ) {
      cout << endl << "AFlow - Automatic-Flow " << endl << "High-Throughput ab-initio Computin" << endl << " Usage: apl2agr [options] " << endl << endl;

        cout << " Possible option switches:" << endl;
        cout << "   -se      or --exact        Will show exact qpoint positions in final agr file." << endl;
        cout << "   -min val or --min val      Set the minimum y value." << endl;
        cout << "   -max val or --max val      Set the maximum y value." << endl;
        cout << "   -ps      or --ps           Will call xmgrace to generate PostScript file from the final agr file." << endl;
        cout << "   -eps     or --eps          Will call xmgrace to generate EPS file from the final agr file." << endl;
        cout << "   -pdf     or --pdf          Will call xmgrace to generate PDF file from the final agr file." << endl;
        cout << "   -png     or --png          Will call xmgrace to generate PNG file from the final agr file." << endl;
        cout << endl;
        return EXIT_SUCCESS;
    }

    for(int i = 1; i < argc; i++) {
        if( ( string(argv[i]) == string("-se") ) ||
            ( string(argv[i]) == string("--exact") ) ) {
            showExactQPoints = true;
            continue;
        }

        if( ( string(argv[i]) == string("-min") ) ||
            ( string(argv[i]) == string("--min") ) ) {
            globalYMin = atof(argv[++i]);
            userSetGlobalYMin = true;
            continue;
        }

        if( ( string(argv[i]) == string("-max") ) ||
            ( string(argv[i]) == string("--max") ) ) {
            globalYMax = atof(argv[++i]);
            userSetGlobalYMax = true;
            continue;
        }

        if( ( string(argv[i]) == string("-ps") ) ||
            ( string(argv[i]) == string("--ps") ) ) {
            userPrintPS = true;
            continue;
        }

        if( ( string(argv[i]) == string("-eps") ) ||
            ( string(argv[i]) == string("--eps") ) ) {
            userPrintEPS = true;
            continue;
        }

        if( ( string(argv[i]) == string("-pdf") ) ||
            ( string(argv[i]) == string("--pdf") ) ) {
            userPrintPDF = true;
            continue;
        }

        if( ( string(argv[i]) == string("-png") ) ||
            ( string(argv[i]) == string("--png") ) ) {
            userPrintPNG = true;
            continue;
        }
    }

    // Prepare phonon dispertion graph ////////////////////////////////////

    XmGraph pdGraph;
    //  pdGraph.setWorldXMax(0.5);

    bool isPDISavalaible = true;
    IPCFreqFlags freqUnits = NONE;

    bool isBZIP=false;
    
    isBZIP=false;
    infile.open("PDIS.bz2");if(infile.is_open()) isBZIP=true;
    if(isBZIP) system("bzip2 -d PDIS.bz2"); 
    infile.close();

    infile.open("PDIS");

    if( infile.is_open() ) {
        string line;
        vector<string> tokens;
        vector<XmSet> xsets;

        vector<double> specialTickPositions;
        vector<string> specialTickLabels;
        vector<double> exactQPointPositions;

        // Read the header
        while( true ) {
            getline(infile,line);
            if( infile.eof() ) break;
            if( line.empty() ) continue;
            if( line[0] != '#' ) break;

            // Read the name of the studied system
            if( line.find("<system>") != string::npos ) {
	      tokenize(line, tokens, string("\"")); //          pdGraph.setSubTitle(tokens[1].c_str());
	      string title=tokens[1];tokens.clear();
	      if(title.find("[Standard_Primitive Unit Cell Form]")!=string::npos) title=title.substr(0,title.find("[Standard_Primitive Unit Cell Form]"));
	      if(title.find("[Standard_Conventional Unit Cell Form]")!=string::npos) title=title.substr(0,title.find("[Standard_Conventional Unit Cell Form]"));
	      if(title.find("(STD_PRIM doi:10.1016/j.commatsci.2010.05.010)")!=string::npos) title=title.substr(0,title.find("(STD_PRIM doi:10.1016/j.commatsci.2010.05.010)"));
	      if(title.find("(STD_CONV doi:10.1016/j.commatsci.2010.05.010)")!=string::npos) title=title.substr(0,title.find("(STD_CONV doi:10.1016/j.commatsci.2010.05.010)"));
	      pdGraph.setSubTitle(title.c_str());
	      continue;
            }

            // Read the units
            if( line.find("<units>") != string::npos ) {
                tokenize(line, tokens, string(" "));
                freqUnits = (IPCFreqFlags)atoi(tokens[2].c_str());
                tokens.clear();
                continue;
            }

            // Read the name of the studied system
            if( line.find("<nbranches>") != string::npos ){
                // Get value
                tokenize(line, tokens, string(" "));
                int nbranches = atoi( tokens[2].c_str() );
                tokens.clear();

                // Prepare the array
                XmSet set;
                for(int i = 0; i < nbranches; i++)
                    xsets.push_back(set);

                //
                continue;
            }

            // Read the label
            if( line.find("<label>") != string::npos )  {
                tokenize(line, tokens, string(" "));
                specialTickPositions.push_back( atof( tokens[2].c_str() ) );
                string label = tokens[3];
                if( label == string("G") || label == string("Gamma") || 
                    label == string("\\Gamma"))
                    label = string("\\xG");
                specialTickLabels.push_back(label);
                tokens.clear();
                continue;
            }        

            // Read positions of exact qpoints
            if( line.find("<exact>") != string::npos ) {
                tokenize(line, tokens, string(" "));
                exactQPointPositions.push_back( atof( tokens[2].c_str() ) );
                tokens.clear();
                continue;
            }
        }

        // Make the fist line
        tokenize(line, tokens, string(" \t"));
        double x = atof( tokens[1].c_str() );
        for(unsigned int i = 0; i < xsets.size(); i++)
            xsets[i].addXY(x, atof( tokens[2+i].c_str() ) );
        tokens.clear();

        // Read sets
        while( true ) {
            getline(infile,line);

            if( infile.eof() ) break;
            if( line.empty() ) continue;
            if( line[0] == '#' ) continue;

            tokenize(line, tokens, string(" \t"));
            double x = atof( tokens[1].c_str() );
            double conversionFactor = 1.0;
            //if( freqUnits & THZ )
            for(unsigned int i = 0; i < xsets.size(); i++)
                xsets[i].addXY(x, conversionFactor * atof( tokens[2+i].c_str() ) );
            tokens.clear();
        }

        // Add these sets to graph...
        pdGraph.addSets(xsets);
        xsets.clear();

        // Setup graph axis...

        // Alternative axis
        XmAxis xaxis;
        xaxis.setActiveOFF();
        if( !showExactQPoints )
            pdGraph.setAltXAxis(xaxis);
        pdGraph.setAltYAxis(xaxis);

        // Normal axis - X
        xaxis.setActiveON();
        xaxis.setTicksPointingOut();
        xaxis.setSpecTickTypeBoth();
        for(unsigned int j = 0; j < specialTickLabels.size(); j++)
            xaxis.addSpecTick(1,j,specialTickPositions[j],specialTickLabels[j].c_str());
        xaxis.setTicksLabelFont(4);
        xaxis.setTicksLabelSize(70.0);
        xaxis.setTickMajorLineStyle(4);
        xaxis.setTickMajorGridON();
        pdGraph.setXAxis(xaxis);

        // Alternative x axis - positoins of exact qpoints
        if( showExactQPoints ) {
            XmAxis axaxis;
            axaxis.setActiveON();
            axaxis.setTicksPointingIn();
            axaxis.setTickMajorColor(2);
            axaxis.setTickMajorSize(0.25);
            axaxis.setTickMajorLineWidth(2);
            axaxis.setSpecTickTypeBoth();
            for(unsigned int j = 0; j < exactQPointPositions.size(); j++)
                axaxis.addSpecTick(1,j,exactQPointPositions[j],"");
            pdGraph.setAltXAxis(axaxis);
        }

        // Normal axis - Y
        XmAxis yaxis;
        yaxis.setTicksPointingOut();
        yaxis.setTicksDrawOnNormalSide();
        yaxis.setTicksLabelFont(4);
        yaxis.setTicksLabelSize(70.0);
        if( freqUnits & THZ )
            yaxis.setLabel("Frequency (THz)");
        else if( freqUnits & RECIPROCAL_CM )
            yaxis.setLabel("Wavenumber (cm\\S-1\\N)");
        else if( freqUnits & MEV)
            yaxis.setLabel("Energy (meV)");
        else
            yaxis.setLabel("Frequency (unknown units)");
        yaxis.setLabelSize(80.0);
        yaxis.setLabelFont(4);
        pdGraph.setYAxis(yaxis);

        // Autoscale Y min
        double min = pdGraph.getWorldYMin();
        if( userSetGlobalYMin ) {
            min = globalYMin;
        } else {
            if( fabs(min) > 0.1 ) {
              min = floor( pdGraph.getWorldYMin() );
              if( fabs( min - (int)min ) > 1E-3 )
                  min = (int)min - 1.0;
	    } else {
	      min = 0.0;
            }
        }
        pdGraph.setWorldYMin( min );

        // Autoscale Y max
        double max = ceil( pdGraph.getWorldYMax() );
        if( userSetGlobalYMax ) {
            max = globalYMax;
        } else {
            if( max - (int)max > 1E-3 )
                max = (int)max + 1.0;
        }
        pdGraph.setWorldYMax( max );
        pdGraph.autoTicksY();
        //pdGraph.autoscale();
        pdGraph.setWorldXMin(0.0);
        pdGraph.setWorldXMax(1.0);
        
        // Setup title and subtitles
        //pdGraph.setTitle("Phonon dispersion");
        //pdGraph.setTitleFont(4);
        //pdGraph.setTitleSize(95);
        pdGraph.setSubTitleFont(4);
        pdGraph.setSubTitleSize(80.0);

        // Add zero line if there are also negative frequencies
        if( pdGraph.getWorldYMin() < 0.1 ) {
            XmLine line;
            line.setLockTypeWorld();
            line.setStartXY(0.0,0.0);
            line.setEndXY(1.0,0.0);
            line.setStyle(2);
            pdGraph.addLine(line);
        }
    }
    else
    {
        isPDISavalaible = false;
    }

    infile.close();
    infile.clear();
    if(isBZIP) system("bzip2 PDIS");

    // Prepare pDos graph //////////////////////////////////////////////////

    bool isPDOSavalaible = true;
    XmGraph dosGraph;

    // Units of pDos will be still in meV, this factor  we need for
    // transformation fo the y axis min and max
    double factor = 1.0;
    if( freqUnits & THZ )
      factor = factorTHz2MeV;
    else if( freqUnits & RECIPROCAL_CM )
      factor = factorRCM2MeV;

    //
    isBZIP=false;
    infile.open("PDOS.bz2");if(infile.is_open()) isBZIP=true;
    if(isBZIP) system("bzip2 -d PDOS.bz2"); 
    infile.close();

    infile.open("PDOS");

    if( infile.is_open() ) {
        string line;
        vector<string> tokens;
        XmSet xset;

        // Read set
        while( true ) {
            getline(infile,line);

            if( infile.eof() ) break;
            if( line.empty() ) continue;
            if( line[0] == '#' ) continue;

            tokenize(line, tokens, string(" \t"));
            xset.addXY( atof( tokens[3].c_str() ), atof( tokens[2].c_str() ) );
            tokens.clear();
        }

        //
        xset.setFillType(1);
        xset.setFillPattern(1);
        xset.setFillColor(1);
        dosGraph.addSet(xset);

        // Alternative axis
        XmAxis xaxis;
        xaxis.setActiveOFF();
        dosGraph.setAltXAxis(xaxis);
        dosGraph.setAltYAxis(xaxis);

        // Normal axis - X
        xaxis.setActiveON();
        xaxis.setTicksOFF();
        xaxis.setTicksLabelOFF();
        dosGraph.setXAxis(xaxis);
  
        // Normal axis - Y
        XmAxis yaxis;
        yaxis.setTicksPointingOut();
        yaxis.setTicksDrawOnOppositeSide();
        yaxis.setTicksLabelSideOpposite();
        yaxis.setTicksLabelFont(4);
        yaxis.setTicksLabelSize(70.0);
        yaxis.setLabel("\\r{180}Energy (meV)");
        yaxis.setLabelSize(80.0);
        yaxis.setLabelFont(4);
        yaxis.setLabelPlaceOpposite();
        dosGraph.setYAxis(yaxis);

        // Autoscale
        dosGraph.autoscale();
        dosGraph.setWorldXMin(0.0);
        dosGraph.setWorldXMax( xset.getXmax() + 0.1 );
        dosGraph.setWorldYMin( pdGraph.getWorldYMin() * factor );
        dosGraph.setWorldYMax( pdGraph.getWorldYMax() * factor );

        // Set subtitle
        dosGraph.setSubTitle("pDOS");
        dosGraph.setSubTitleFont(4);
        dosGraph.setSubTitleSize(80.0);

        // Add zero line if there are also negative frequencies
        if( dosGraph.getWorldYMin() < 0.1 ) {
            XmLine line;
            line.setLockTypeWorld();
            line.setStartXY(0.0,0.0);
            line.setEndXY(dosGraph.getWorldXMax(),0.0);
            line.setStyle(2);
            dosGraph.addLine(line);
        }
    }
    else if( isPDISavalaible ) {        
        isPDOSavalaible = false;

        // Alternative axis
        XmAxis xaxis;
        xaxis.setActiveOFF();
        dosGraph.setAltXAxis(xaxis);
        dosGraph.setAltYAxis(xaxis);

        // Normal axis - X
        xaxis.setActiveON();
        xaxis.setTicksOFF();
        xaxis.setTicksLabelOFF();
        dosGraph.setXAxis(xaxis);
  
        // Normal axis - Y
        XmAxis yaxis;
        yaxis.setTicksPointingOut();
        yaxis.setTicksDrawOnOppositeSide();
        yaxis.setTicksLabelSideOpposite();
        yaxis.setTicksLabelFont(4);
        yaxis.setTicksLabelSize(70.0);
        yaxis.setLabel("\\r{180}Energy (meV)");
        yaxis.setLabelSize(80.0);
        yaxis.setLabelFont(4);
        yaxis.setLabelPlaceOpposite();
        dosGraph.setYAxis(yaxis);

        // Autoscale
        dosGraph.setWorldXMin(0.0);
        dosGraph.setWorldXMax(1.0);
        dosGraph.autoTicksX();
        dosGraph.setWorldYMin( pdGraph.getWorldYMin() * factor );
        dosGraph.setWorldYMax( pdGraph.getWorldYMax() * factor );
        dosGraph.autoTicksY();
    }
    else
    {
        isPDOSavalaible = false;
    }

    infile.close();
    infile.clear();
    if(isBZIP) system("bzip2 PDOS");

    // Prepare thermo graphs ///////////////////////////////////////////////

    bool isTHERMOavalaible = true;

    XmGraph freeEnergyGraph;
    XmGraph entropyGraph;
    XmGraph heatCapacityGraph;

    isBZIP=false;
    infile.open("THERMO.bz2");if(infile.is_open()) isBZIP=true;
    if(isBZIP) system("bzip2 -d THERMO.bz2"); 
    infile.close();

    infile.open("THERMO");

    if( infile.is_open() ) {
        string line;
        vector<string> tokens;
        XmSet xsetFreeEnergy;
        XmSet xsetEntropy;
        XmSet xsetHeatCapacity;

        // Read set
        while( true ) {
            getline(infile,line);

            if( infile.eof() ) break;
            if( line.empty() ) continue;
            if( line[0] == '#' ) continue;

            tokenize(line, tokens, string(" \t"));
            double t = atof( tokens[0].c_str() );
            xsetFreeEnergy.addXY( t, atof( tokens[3].c_str() ) );
            xsetEntropy.addXY( t, atof( tokens[4].c_str() ) );
            xsetHeatCapacity.addXY( t, atof( tokens[5].c_str() ) );
            tokens.clear();
        }

        //
        freeEnergyGraph.addSet(xsetFreeEnergy);
        entropyGraph.addSet(xsetEntropy);
        heatCapacityGraph.addSet(xsetHeatCapacity);

        // Alternative axis
        XmAxis xaxis;
        xaxis.setActiveOFF();
        freeEnergyGraph.setAltXAxis(xaxis);
        freeEnergyGraph.setAltYAxis(xaxis);
        entropyGraph.setAltXAxis(xaxis);
        entropyGraph.setAltYAxis(xaxis);
        heatCapacityGraph.setAltXAxis(xaxis);
        heatCapacityGraph.setAltYAxis(xaxis);

        // Normal axis - X
        xaxis.setActiveON();
        xaxis.setTicksPointingOut();
        xaxis.setTicksLabelFont(4);
        xaxis.setTicksLabelSize(50.0);
        xaxis.setLabel("Temperature (K)");
        xaxis.setLabelSize(60.0);
        xaxis.setLabelFont(4);
        xaxis.setTickMajorLineStyle(2);
        xaxis.setTickMajorGridON();
        freeEnergyGraph.setXAxis(xaxis);
        entropyGraph.setXAxis(xaxis);
        heatCapacityGraph.setXAxis(xaxis);

        // Normal axis - Y
        XmAxis yaxis;
        yaxis.setTicksPointingOut();
        yaxis.setTicksLabelFont(4);
        yaxis.setTicksLabelSize(50.0);
        yaxis.setLabelSize(60.0);
        yaxis.setLabelFont(4);

        yaxis.setLabel("F\\svib\\N (meV/cell)");
        freeEnergyGraph.setYAxis(yaxis);

        yaxis.setLabel("S\\svib\\N (k\\sB\\N/cell)");
        entropyGraph.setYAxis(yaxis);

        yaxis.setTicksLabelSideOpposite();
        yaxis.setLabelPlaceOpposite();
        yaxis.setLabel("\\r{180}c\\sV\\N (k\\sB\\N/cell)");
        heatCapacityGraph.setYAxis(yaxis);

        // Autoscale
        freeEnergyGraph.autoscaleX();
        entropyGraph.autoscaleX();
        heatCapacityGraph.autoscaleX();

        // Refine Y based on special purposes of ploted quantities
        double min = floor( freeEnergyGraph.getWorldYMin() );
        if( fabs( min - (int)min ) > 1E-3 )
            min = (int)min - 1.0;
        freeEnergyGraph.autoscaleY();
        freeEnergyGraph.setWorldYMin(min);
        freeEnergyGraph.autoTicksY();

        double max = ceil( entropyGraph.getWorldYMax() );
        if( max - (int)max > 1E-3 )
            max = (int)max + 1.0;
        entropyGraph.autoscaleY();
        entropyGraph.setWorldYMax(max);
        entropyGraph.autoTicksY();

        max = ceil( heatCapacityGraph.getWorldYMax() );
        if( max - (int)max < 1E-3 )
            max = (int)max + 1.0;
        heatCapacityGraph.autoscaleY();
        heatCapacityGraph.setWorldYMax(max);
        heatCapacityGraph.autoTicksY();

    } else {
        isTHERMOavalaible = false;
    }

    infile.close();
    infile.clear();
    if(isBZIP) system("bzip2 THERMO");

    // Pack it all to the xmgrace file and save it ////////////////////////

    XmGraceFile xgf;

    // Set page
    xgf.setPageSize(LETTER);

    // Set scaling factor depending of the selected page format
    double fx = 1.0;
    double fy = 1.0;
    double fscale=1.0,origx=0.0,origy=0.0;
    if(userPrintPDF || userPrintEPS) {fscale=0.8;origx=-0.03;origy=0.5;}

    if( xgf.getPageSize() == LETTER ) {
        fx = numberX/1.14334917;
	fy = 1.0/fx;
	//   fy = 1.0/fx/fscale;
    }

    // Add dispertion graph and set its location and dimension
    if( isPDISavalaible ) {
        if( !isPDOSavalaible )
            pdGraph.setView(0.15,numberX/fx,0.42/fy,0.85);
        else
            pdGraph.setView(0.15,1.15/fx,0.42/fy,0.85);

	pdGraph._view_xmin*=fscale;pdGraph._view_xmax*=fscale;pdGraph._view_ymin*=fscale;pdGraph._view_ymax*=fscale;
	pdGraph._view_xmin+=origx;pdGraph._view_xmax+=origx;pdGraph._view_ymin+=origy;pdGraph._view_ymax+=origy;

        xgf.addGraph(pdGraph);
    }

    // Add phonon density of states graph and set its lacation and dimension
    if( isPDOSavalaible ) {
        dosGraph.setView(1.16/fx,numberX/fx,0.42/fy,0.85);
	dosGraph._view_xmin*=fscale;dosGraph._view_xmax*=fscale;dosGraph._view_ymin*=fscale;dosGraph._view_ymax*=fscale;
	dosGraph._view_xmin+=origx;dosGraph._view_xmax+=origx;dosGraph._view_ymin+=origy;dosGraph._view_ymax+=origy;
        xgf.addGraph(dosGraph);
    }
    else if( isPDISavalaible ) {
      dosGraph.setView(0.15,numberX/fx,0.42/fy,0.85);
      dosGraph._view_xmin*=fscale;dosGraph._view_xmax*=fscale;dosGraph._view_ymin*=fscale;dosGraph._view_ymax*=fscale;
      dosGraph._view_xmin+=origx;dosGraph._view_xmax+=origx;dosGraph._view_ymin+=origy;dosGraph._view_ymax+=origy;
      xgf.addGraph(dosGraph);
    }    
    
    // Add thermal properties
    if( isTHERMOavalaible ) {
      double width = 0.3 / fx;
      double x0 = 0.15;
      double spacex = 0.12 / fx;
      freeEnergyGraph.setView(x0,x0+width,0.15,0.355882/fy);
      entropyGraph.setView(x0+width+spacex,x0+2.0*width+spacex,0.15,0.355882/fy);
      heatCapacityGraph.setView(numberX/fx-width,numberX/fx,0.15,0.355882/fy);
      
      freeEnergyGraph._view_xmin*=fscale;freeEnergyGraph._view_xmax*=fscale;freeEnergyGraph._view_ymin*=fscale;freeEnergyGraph._view_ymax*=fscale;
      freeEnergyGraph._view_xmin+=origx;freeEnergyGraph._view_xmax+=origx;freeEnergyGraph._view_ymin+=origy;freeEnergyGraph._view_ymax+=origy;
      entropyGraph._view_xmin*=fscale;entropyGraph._view_xmax*=fscale;entropyGraph._view_ymin*=fscale;entropyGraph._view_ymax*=fscale;
      entropyGraph._view_xmin+=origx;entropyGraph._view_xmax+=origx;entropyGraph._view_ymin+=origy;entropyGraph._view_ymax+=origy;
      heatCapacityGraph._view_xmin*=fscale;heatCapacityGraph._view_xmax*=fscale;heatCapacityGraph._view_ymin*=fscale;heatCapacityGraph._view_ymax*=fscale;
      heatCapacityGraph._view_xmin+=origx;heatCapacityGraph._view_xmax+=origx;heatCapacityGraph._view_ymin+=origy;heatCapacityGraph._view_ymax+=origy;
      
      xgf.addGraph(freeEnergyGraph);
      xgf.addGraph(entropyGraph);
      xgf.addGraph(heatCapacityGraph);
    }

    // Save it
    xgf.saveAs("apl.agr");

    // ////////////////////////////////////////////////////////////////////

    if( userPrintPS ) {
        // Print file throught xmgrace
        // xmgrace -hdevice PostScript -hardcopy apl.agr
        char* args[6];
        args[0] = (char*)("xmgrace");
        args[1] = (char*)("-hdevice");
        args[2] = (char*)("PostScript");
        args[3] = (char*)("-hardcopy");
        args[4] = (char*)("apl.agr");
        args[5] = NULL;
        pid_t pid = fork();
        if( pid == 0 ) {
            (void)execvp("xmgrace", args);
            return EXIT_SUCCESS;
        }
        if( pid < 0 )
          cout << "Problem to fork." << std::endl;  
    }

    if( userPrintEPS ) {
        // Print file throught xmgrace
        // xmgrace -hdevice EPS -hardcopy apl.agr
        char* args[6];
        args[0] = (char*)("xmgrace");
        args[1] = (char*)("-hdevice");
        args[2] = (char*)("EPS");
        args[3] = (char*)("-hardcopy");
        args[4] = (char*)("apl.agr");
        args[5] = NULL;
        pid_t pid = fork();
        if( pid == 0 ) {
            (void)execvp("xmgrace", args);
           return EXIT_SUCCESS;
        }
        if( pid < 0 )
          cout << "Problem to fork." << std::endl;  
    }

    if( userPrintPDF ) {
        // Print file throught xmgrace
        // xmgrace -hdevice PDF -hardcopy apl.agr
        char* args[6];
        args[0] = (char*)("xmgrace");
        args[1] = (char*)("-hdevice");
        args[2] = (char*)("EPS"); // PDF
        args[3] = (char*)("-hardcopy");
        args[4] = (char*)("apl.agr");
        args[5] = NULL;
        pid_t pid = fork();
        if( pid == 0 ) {
	  //	  (void)execvp("xmgrace", args);
	  (void)system(string("xmgrace "+string(args[1])+" "+string(args[2])+" "+string(args[3])+" "+string(args[4])+" ").c_str());
	  //	  (void)system("convert apl.eps -rotate \"+90>\" aplP.eps");
	  (void)system("ps2pdf apl.eps apl.pdf");
	  return EXIT_SUCCESS;
        }
        if( pid < 0 )
          cout << "Problem to fork." << std::endl;  
    }

    if( userPrintPNG ) {
        // Print file throught xmgrace
        // xmgrace -hdevice PNG -hardcopy apl.agr
        char* args[6];
        args[0] = (char*)("xmgrace");
        args[1] = (char*)("-hdevice");
        args[2] = (char*)("PNG");
        args[3] = (char*)("-hardcopy");
        args[4] = (char*)("apl.agr");
        args[5] = NULL;
        pid_t pid = fork();
        if( pid == 0 ) {
            (void)execvp("xmgrace", args);
            return EXIT_SUCCESS;
        }
        if( pid < 0 )
          cout << "Problem to fork." << std::endl;  
    }

    // ////////////////////////////////////////////////////////////////////
    
    return EXIT_SUCCESS;
}

// //////////////////////////////////////////////////////////////////////////

void tokenize(const string& str,vector<string>& tokens, string del)
{
    string delimiters = del;

    // Skip delimiters at beginning.
    string::size_type lastPos = str.find_first_not_of(delimiters, 0);

    // Find first "non-delimiter".
    string::size_type pos = str.find_first_of(delimiters, lastPos);

    while (string::npos != pos || string::npos != lastPos) {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));

        // Skip delimiters. Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);

        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }
}
