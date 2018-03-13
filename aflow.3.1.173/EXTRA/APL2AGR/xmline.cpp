#include <iostream>
#include <iomanip>
#include "xmlib.h"

using namespace std;

namespace xmlib {

// ///////////////////////////////////////////////////////////////////////////

XmLine::XmLine()
{
    _active = true;
    setLockTypeWorld();
    _startx = _starty = _endx = _endy = 0.0;
    _linewidth = 1.0;
    _linestyle = 1;
    _linecolor = 1;
}

// ///////////////////////////////////////////////////////////////////////////

XmLine::~XmLine()
{
}

// ///////////////////////////////////////////////////////////////////////////

void XmLine::setActiveON()
{
    _active = true;
}

// ///////////////////////////////////////////////////////////////////////////

void XmLine::setActiveOFF()
{
    _active = false;
}

// ///////////////////////////////////////////////////////////////////////////

void XmLine::setLockTypeWorld()
{
    _locktype = string("world");
}

// ///////////////////////////////////////////////////////////////////////////

void XmLine::setLockTypeView()
{
    _locktype = string("view");
}

// ///////////////////////////////////////////////////////////////////////////
    
void XmLine::setStartXY(double x, double y )
{
    _startx = x;
    _starty = y;
}

// ///////////////////////////////////////////////////////////////////////////
    
void XmLine::setEndXY(double x, double y )
{
    _endx = x;
    _endy = y;
}

// ///////////////////////////////////////////////////////////////////////////

void XmLine::setWidth(double w)
{
    _linewidth = w;
}

// ///////////////////////////////////////////////////////////////////////////

void XmLine::setStyle(int s)
{
    _linestyle = s;
}

// ///////////////////////////////////////////////////////////////////////////

void XmLine::setColor(int c)
{
    _linecolor = c;
}

// ///////////////////////////////////////////////////////////////////////////

XmLine& XmLine::operator= (const XmLine& that)
{
    if( this != &that )
    {
        _active = that._active;
        _locktype = that._locktype;
        _startx = that._startx;
        _starty = that._starty;
        _endx = that._endx;
        _endy = that._endy;
        _linewidth = that._linewidth;
        _linestyle = that._linestyle;
        _linecolor = that._linecolor;
    }

    return *this;
}

// ///////////////////////////////////////////////////////////////////////////

void XmLine::write(std::ostream& os, int belong_graph)
{
    os << "@with line" << endl;    
    os << "@    line "  << (_active ? "on":"off") << endl;	
    os << "@    line loctype " << _locktype << endl;
    os << "@    line g" << belong_graph << endl;
    os << "@    line " << _startx << ", " << _starty << ", "
	                     << _endx << ", " << _endy << endl;
    os << "@    line linewidth " << _linewidth << endl;
    os << "@    line linestyle " << _linestyle << endl;
    os << "@    line color " << _linecolor << endl;
    os << "@    line arrow 0" << endl;
    os << "@    line arrow type 0" << endl;
    os << "@    line arrow length 1.000000" << endl;
    os << "@    line arrow layout 1.000000, 1.000000" << endl;
    os << "@line def" << endl;
}

// ///////////////////////////////////////////////////////////////////////////

} // end namespace xmlib
