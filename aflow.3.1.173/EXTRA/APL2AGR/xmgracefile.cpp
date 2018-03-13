#include <fstream>
#include "xmlib.h"

using namespace std;

namespace xmlib {

// ///////////////////////////////////////////////////////////////////////////

XmGraceFile::XmGraceFile()
{
    setPageSize(A4);
}

// ///////////////////////////////////////////////////////////////////////////

XmGraceFile::~XmGraceFile()
{
    clear();
}

// ///////////////////////////////////////////////////////////////////////////

void XmGraceFile::clear()
{
    for(unsigned int i = 0; i < _graphs.size(); i++)
        _graphs[i].clear();
    _graphs.clear();
}

// ///////////////////////////////////////////////////////////////////////////

void XmGraceFile::addGraph(const XmGraph& g)
{
    _graphs.push_back(g);
}

// ///////////////////////////////////////////////////////////////////////////

void XmGraceFile::setPageSize(const PageSize& ps)
{
    _pageSize = ps;
}

const PageSize& XmGraceFile::getPageSize()
{
    return _pageSize;
}

// ///////////////////////////////////////////////////////////////////////////

void XmGraceFile::saveAs(const char* filename)
{
    // Create file
    ofstream outfile(filename,ios_base::out);

    if( !outfile.is_open() )
        throw XMLIBRuntimeError("xmlib::XmGraceFile::write(); Can not open output file.");

    // Write header of file
    outfile << "# Grace project file" << endl;
    outfile << "#" << endl;
    outfile << "@version 50121" << endl;
    outfile << "@page size ";
    if( _pageSize == A4 )
        outfile << "842, 594";
    if( _pageSize == LETTER )
        outfile << "792, 612";
    outfile << endl;
    outfile << "@page scroll 5%" << endl;
    outfile << "@page inout 5%" << endl;
    outfile << "@link page off" << endl;
    outfile << "@map font 0 to \"Times-Roman\", \"Times-Roman\"" << endl;
    outfile << "@map font 1 to \"Times-Italic\", \"Times-Italic\"" << endl;
    outfile << "@map font 2 to \"Times-Bold\", \"Times-Bold\"" << endl;
    outfile << "@map font 3 to \"Times-BoldItalic\", \"Times-BoldItalic\"" << endl;
    outfile << "@map font 4 to \"Helvetica\", \"Helvetica\"" << endl;
    outfile << "@map font 5 to \"Helvetica-Oblique\", \"Helvetica-Oblique\"" << endl;
    outfile << "@map font 6 to \"Helvetica-Bold\", \"Helvetica-Bold\"" << endl;
    outfile << "@map font 7 to \"Helvetica-BoldOblique\", \"Helvetica-BoldOblique\"" << endl;
    outfile << "@map font 8 to \"Courier\", \"Courier\"" << endl;
    outfile << "@map font 9 to \"Courier-Oblique\", \"Courier-Oblique\"" << endl;
    outfile << "@map font 10 to \"Courier-Bold\", \"Courier-Bold\"" << endl;
    outfile << "@map font 11 to \"Courier-BoldOblique\", \"Courier-BoldOblique\"" << endl;
    outfile << "@map font 12 to \"Symbol\", \"Symbol\"" << endl;
    outfile << "@map font 13 to \"ZapfDingbats\", \"ZapfDingbats\"" << endl;
    outfile << "@map color 0 to (255, 255, 255), \"white\"" << endl;
    outfile << "@map color 1 to (0, 0, 0), \"black\"" << endl;
    outfile << "@map color 2 to (255, 0, 0), \"red\"" << endl;
    outfile << "@map color 3 to (0, 255, 0), \"green\"" << endl;
    outfile << "@map color 4 to (0, 0, 255), \"blue\"" << endl;
    outfile << "@map color 5 to (255, 255, 0), \"yellow\"" << endl;
    outfile << "@map color 6 to (188, 143, 143), \"brown\"" << endl;
    outfile << "@map color 7 to (220, 220, 220), \"grey\"" << endl;
    outfile << "@map color 8 to (148, 0, 211), \"violet\"" << endl;
    outfile << "@map color 9 to (0, 255, 255), \"cyan\"" << endl;
    outfile << "@map color 10 to (255, 0, 255), \"magenta\"" << endl;
    outfile << "@map color 11 to (255, 165, 0), \"orange\"" << endl;
    outfile << "@map color 12 to (114, 33, 188), \"indigo\"" << endl;
    outfile << "@map color 13 to (103, 7, 72), \"maroon\"" << endl;
    outfile << "@map color 14 to (64, 224, 208), \"turquoise\"" << endl;
    outfile << "@map color 15 to (0, 139, 0), \"green4\"" << endl;
    outfile << "@reference date 0" << endl;
    outfile << "@date wrap off" << endl;
    outfile << "@date wrap year 1950" << endl;
    outfile << "@default linewidth 1.0" << endl;
    outfile << "@default linestyle 1" << endl;
    outfile << "@default color 1" << endl;
    outfile << "@default pattern 1" << endl;
    outfile << "@default font 0" << endl;
    outfile << "@default char size 1.000000" << endl;
    outfile << "@default symbol size 1.000000" << endl;
    outfile << "@default sformat \"%.8g\"" << endl;
    outfile << "@background color 0" << endl;
    outfile << "@page background fill on" << endl;

    // Write lines
    for(unsigned int i = 0; i < _graphs.size(); i++)
    {
        for(unsigned int j = 0; j < _graphs[i].getLines().size(); j++)
             _graphs[i].getLine(j).write(outfile,i);
    }

    // Write graphs
    for(unsigned int i = 0; i < _graphs.size(); i++)
	     _graphs[i].write(outfile,i);

    // Close file
    outfile.close();
    outfile.clear();
}

// ///////////////////////////////////////////////////////////////////////////

} // end namespace xmlib
