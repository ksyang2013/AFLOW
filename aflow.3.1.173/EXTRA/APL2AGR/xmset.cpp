#include <iostream>
#include <iomanip>
#include "xmlib.h"

using namespace std;

namespace xmlib {

// ///////////////////////////////////////////////////////////////////////////

XmSet::XmSet()
{
    setDefault();
}

// ///////////////////////////////////////////////////////////////////////////

XmSet::~XmSet()
{
    clear();
}

// ///////////////////////////////////////////////////////////////////////////

void XmSet::clear()
{
    _x.clear();
    _y.clear();
}

// ///////////////////////////////////////////////////////////////////////////

void XmSet::setDefault()
{
    clear();
    _xmin = _xmax = _ymin = _ymax = 0;
    _lineStyle = _lineType = _lineColor = _linePattern = 1;
    _lineWidth = 1.0;
    _fillType = 0;
    _fillColor = 0;
    _fillPattern = 1;
    _comment = string("");
}

// ///////////////////////////////////////////////////////////////////////////

void XmSet::addXY(double x, double y)
{
    _x.push_back(x);
    _y.push_back(y);
    updateExtremas();
}

// ///////////////////////////////////////////////////////////////////////////

void XmSet::setXY(unsigned int n, double x, double y)
{
    if( n > _x.size() )
        throw XMLIBRuntimeError("xmlib::XmSet::setXY(); Out of range.");
    _x[n] = x;
    _y[n] = y;
    updateExtremas(n);
}

// ///////////////////////////////////////////////////////////////////////////

XmSet& XmSet::operator= (const XmSet& that)
{
    if( this != &that )
    {
        _x.clear();
        _y.clear();
        _x = that._x;
        _y = that._y;

        _xmin = that._xmin;
        _xmax = that._xmax;
        _ymin = that._ymin;
        _ymax = that._ymax;

        _lineStyle = that._lineStyle;
        _lineType = that._lineType;
        _lineWidth = that._lineWidth;
        _lineColor = that._lineColor;
        _linePattern = that._linePattern;

        _comment = that._comment;

        _fillType = that._fillType;
        _fillColor = that._fillColor;
        _fillPattern = that._fillPattern;
    }

    return *this;
}

// ///////////////////////////////////////////////////////////////////////////

void XmSet::setLineType(int n)
{
    _lineType = n;
}

void XmSet::setLineStyle(int n)
{
    _lineStyle = n;
}

void XmSet::setLineWidth(double w)
{
    _lineWidth = w;
}

void XmSet::setLineColor(int n)
{
    _lineColor = n;
}

void XmSet::setLinePattern(int n)
{
    _linePattern = n;
}

int XmSet::getLineType()
{
    return _lineType;
}

int XmSet::getLineStyle()
{
    return _lineStyle;
}

double XmSet::getLineWidth()
{
    return _lineWidth;
}

int XmSet::getLineColor()
{
    return _lineColor;
}

int XmSet::getLinePattern()
{
    return _linePattern;
}

// ///////////////////////////////////////////////////////////////////////////

void XmSet::setComment(const char* s)
{
    _comment = string(s);
}

void XmSet::setComment(const string& s)
{
    _comment = s;
}

const string& XmSet::getComment()
{
    return _comment;
}

// ///////////////////////////////////////////////////////////////////////////

void XmSet::setFillType(int n)
{
    _fillType = n;
}

void XmSet::setFillColor(int n)
{
    _fillColor = n;
}

void XmSet::setFillPattern(int n)
{
    _fillPattern = n;
}

int XmSet::getFillType()
{
    return _fillType;
}

int XmSet::getFillColor()
{
    return _fillColor;
}

int XmSet::getFillPattern()
{
    return _fillPattern;
}

// ///////////////////////////////////////////////////////////////////////////

double XmSet::getXmin()
{
    return _x[_xmin];
}

double XmSet::getXmax()
{
    return _x[_xmax];
}

double XmSet::getYmin()
{
    return _y[_ymin];
}

double XmSet::getYmax()
{
    return _y[_ymax];
}

// ///////////////////////////////////////////////////////////////////////////

void XmSet::write(std::ostream& os, int id)
{
    os << "@type xy" << std::endl;
    for(unsigned int i = 0; i < _x.size(); i++)
      	os << setw(15) << setprecision(8) << _x[i]
	         << setw(15) << setprecision(8) << _y[i] << std::endl;
    os << "&" << std::endl;
}

// ///////////////////////////////////////////////////////////////////////////

void XmSet::updateExtremas()
{
    updateExtremas(_x.size()-1);
}

void XmSet::updateExtremas(int i)
{
    if( _x[i] > _x[_xmax] ) _xmax = i;
    if( _y[i] > _y[_ymax] ) _ymax = i;
    if( _x[i] < _x[_xmin] ) _xmin = i;
    if( _y[i] < _y[_ymin] ) _ymin = i;
}

// ///////////////////////////////////////////////////////////////////////////

} // end namespace xmlib
