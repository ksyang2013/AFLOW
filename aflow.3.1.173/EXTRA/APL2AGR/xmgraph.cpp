#include <iostream>
#include <cmath>
#include "xmlib.h"

using namespace std;

namespace xmlib {

// ///////////////////////////////////////////////////////////////////////////

XmGraph::XmGraph()
{
    initAtributes();
}

// ///////////////////////////////////////////////////////////////////////////

XmGraph::~XmGraph()
{
    clear();
}

// ///////////////////////////////////////////////////////////////////////////

void XmGraph::clear()
{
    _sets.clear();
    _lines.clear();
    _xaxis.clear();
    _yaxis.clear();
    _altxaxis.clear();
    _altyaxis.clear();    
}

// ///////////////////////////////////////////////////////////////////////////

void XmGraph::initAtributes()
{
    _world_xmin = _world_ymin =  1000000.0;
    _world_xmax = _world_ymax = -1000000.0;

    _world_xmin = 0.150000;
    _world_xmax = 1.266801;
    _world_ymin = 0.150000;
    _world_ymax = 0.850000;

    _title = string("");
    _subtitle = string("");
    _titleFont = _subtitleFont = 0;
    _titleSize = _subtitleSize = 1.5;
    _titleColor = _subtitleColor = 1;

    _xaxis.setID("xaxis");
    _yaxis.setID("yaxis");
    _altxaxis.setID("altxaxis");
    _altyaxis.setID("altyaxis");

    _hidden = false;
}
// ///////////////////////////////////////////////////////////////////////////

void XmGraph::addSet(const XmSet& inset)
{
    _sets.push_back(inset);
    updateExtremas();
}

void XmGraph::addSets(const std::vector<XmSet>& sets)
{
    for(unsigned int i = 0; i < sets.size(); i++)
        addSet(sets[i]);
}

// ///////////////////////////////////////////////////////////////////////////

XmGraph& XmGraph::operator= (const XmGraph& that)
{
    if( this != &that )
    {
        _sets.clear();
        _sets = that._sets;
        _lines.clear();
        _lines = that._lines;

        _world_xmin = that._world_xmin;
        _world_ymin = that._world_ymin;
        _world_xmax = that._world_xmax;
        _world_ymax = that._world_ymax;

        _view_xmin = that._view_xmin;
        _view_ymin = that._view_ymin;
        _view_xmax = that._view_xmax;
        _view_ymax = that._view_ymax;

        _title = that._title;
        _titleFont = that._titleFont;
        _titleSize = that._titleSize;
        _titleColor = that._titleColor;

        _subtitle = that._subtitle;
        _subtitleFont = that._subtitleFont;
        _subtitleSize = that._subtitleSize;
        _subtitleColor = that._subtitleColor;

        _xaxis = that._xaxis;
        _yaxis = that._yaxis;
        _altxaxis = that._altxaxis;
        _altyaxis = that._altyaxis;

        _hidden = that._hidden;
    }

    return *this;
}

// ///////////////////////////////////////////////////////////////////////////

void XmGraph::addLine(const XmLine& line)
{
    _lines.push_back(line);
}

// ///////////////////////////////////////////////////////////////////////////

std::vector<XmLine>& XmGraph::getLines()
{
    return _lines;
}

XmLine& XmGraph::getLine(int i)
{
    return _lines[i];
}

// ///////////////////////////////////////////////////////////////////////////

void XmGraph::setXAxis(const XmAxis& xa)
{
    _xaxis = xa;
    _xaxis.setID("xaxis");
}

void XmGraph::setYAxis(const XmAxis& ya)
{
    _yaxis = ya;
    _yaxis.setID("yaxis");
}

void XmGraph::setAltXAxis(const XmAxis& xa)
{
    _altxaxis = xa;
    _altxaxis.setID("altxaxis");
}

void XmGraph::setAltYAxis(const XmAxis& ya)
{
    _altyaxis = ya;
    _altyaxis.setID("altyaxis");
}

// ///////////////////////////////////////////////////////////////////////////

double XmGraph::getWorldXMin()
{
    return _world_xmin;
}

double XmGraph::getWorldYMin()
{
    return _world_ymin;
}

double XmGraph::getWorldXMax()
{
    return _world_xmax;
}

double XmGraph::getWorldYMax()
{
    return _world_ymax;
}

void XmGraph::setWorldXMin(double min)
{
    _world_xmin = min;
}

void XmGraph::setWorldXMax(double max)
{
    _world_xmax = max;
}

void XmGraph::setWorldYMin(double min)
{
    _world_ymin = min;
}

void XmGraph::setWorldYMax(double max)
{
    _world_ymax = max;
}

// ///////////////////////////////////////////////////////////////////////////

void XmGraph::setView(double xmin, double xmax, double ymin, double ymax)
{
    _view_xmin = xmin;
    _view_xmax = xmax;
    _view_ymin = ymin;
    _view_ymax = ymax;
}

// ///////////////////////////////////////////////////////////////////////////

bool XmGraph::isHidden()
{
    return ( _hidden ? true : false );
}

void XmGraph::hide()
{
    _hidden = true;
}

void XmGraph::show()
{
    _hidden = false;
}

// ///////////////////////////////////////////////////////////////////////////

void XmGraph::setTitle(const char *str)
{
    _title = string(str);
}

void XmGraph::setTitleFont(int fn)
{
    _titleFont = fn;
}

void XmGraph::setTitleSize(double fs)
{
    _titleSize = fs;
}

// ///////////////////////////////////////////////////////////////////////////

void XmGraph::setSubTitle(const char *str)
{
    _subtitle = string(str);
}

void XmGraph::setSubTitleFont(int fn)
{
    _subtitleFont = fn;
}

void XmGraph::setSubTitleSize(double fs)
{
    _subtitleSize = fs;
}

// ///////////////////////////////////////////////////////////////////////////

const XmAxis& XmGraph::getXAxis()
{
    return _xaxis;
}

const XmAxis& XmGraph::getYAxis()
{
    return _yaxis;
}

// ///////////////////////////////////////////////////////////////////////////

void XmGraph::autoscale()
{
    round_axis_limits(&_world_xmin,&_world_xmax,1);
    round_axis_limits(&_world_ymin,&_world_ymax,1);
    auto_ticks(_world_xmin,_world_xmax,_xaxis);
    auto_ticks(_world_xmin,_world_xmax,_altxaxis);
    auto_ticks(_world_ymin,_world_ymax,_yaxis);
    auto_ticks(_world_ymin,_world_ymax,_altyaxis);
}

void XmGraph::autoscaleX()
{
    round_axis_limits(&_world_xmin,&_world_xmax,1);
    auto_ticks(_world_xmin,_world_xmax,_xaxis);
    auto_ticks(_world_xmin,_world_xmax,_altxaxis);
}

void XmGraph::autoscaleY()
{
    round_axis_limits(&_world_ymin,&_world_ymax,1);
    auto_ticks(_world_ymin,_world_ymax,_yaxis);
    auto_ticks(_world_ymin,_world_ymax,_altyaxis);
}

void XmGraph::autoTicksX()
{
    auto_ticks(_world_xmin,_world_xmax,_xaxis);
    auto_ticks(_world_xmin,_world_xmax,_altxaxis);
}

void XmGraph::autoTicksY()
{
    auto_ticks(_world_ymin,_world_ymax,_yaxis);
    auto_ticks(_world_ymin,_world_ymax,_altyaxis);
}

// ///////////////////////////////////////////////////////////////////////////

void XmGraph::write(std::ostream& os, int id)
{
    os << "@g" << id << " on" << endl;
    os << "@g" << id << " hidden " << ((_hidden) ? "true":"false") << endl;
    os << "@g" << id << " type xy" << endl;
    os << "@g" << id << " stacked false" << endl;
    os << "@g" << id << " bar hgap 0.000000" << endl;
    os << "@g" << id << " fixedpoint off" << endl;
    os << "@g" << id << " fixedpoint type 0" << endl;
    os << "@g" << id << " fixedpoint xy 0.000000, 0.000000" << endl;
    os << "@g" << id << " fixedpoint format general general" << endl;
    os << "@g" << id << " fixedpoint prec 6, 6" << endl;

    os << "@with g" << id << endl;
    os << "@    world xmin " << _world_xmin << endl;
    os << "@    world xmax " << _world_xmax << endl;
    os << "@    world ymin " << _world_ymin << endl;
    os << "@    world ymax " << _world_ymax << endl;
    os << "@    stack world 0, 0, 0, 0" << endl;
    os << "@    znorm 1" << endl;
    os << "@    view xmin " << _view_xmin << endl;
    os << "@    view xmax " << _view_xmax << endl;
    os << "@    view ymin " << _view_ymin << endl;
    os << "@    view ymax " << _view_ymax << endl;
    os << "@    title \"" << _title << "\"" << endl;
    os << "@    title font " << _titleFont << endl;
    os << "@    title size " << _titleSize/100.0 << endl;
    os << "@    title color " << _titleColor << endl;
    os << "@    subtitle \"" << _subtitle << "\"" << endl;
    os << "@    subtitle font " << _subtitleFont << endl;
    os << "@    subtitle size " << _subtitleSize/100.0 << endl;
    os << "@    subtitle color " << _subtitleColor << endl;
    os << "@    xaxes scale Normal" << endl;
    os << "@    yaxes scale Normal" << endl;
    os << "@    xaxes invert off" << endl;
    os << "@    yaxes invert off" << endl;

    _xaxis.write(os);
    _yaxis.write(os);
    _altxaxis.write(os);
    _altyaxis.write(os);

    for(unsigned int i = 0; i < _sets.size(); i++)
    {
	      os << "@    s" << i << " line type " << _sets[i].getLineType() << endl;
	      os << "@    s" << i << " line linestyle " << _sets[i].getLineStyle() << endl;
	      os << "@    s" << i << " line linewidth " << _sets[i].getLineWidth() << endl;
	      os << "@    s" << i << " line color " << _sets[i].getLineColor() << endl;
	      os << "@    s" << i << " line pattern " << _sets[i].getLinePattern() << endl;
        os << "@    s" << i << " comment \"" << _sets[i].getComment() << "\"" << endl;
        os << "@    s" << i << " fill type " << _sets[i].getFillType() << endl;
        os << "@    s" << i << " fill color " << _sets[i].getFillColor() << endl;
        os << "@    s" << i << " fill pattern " << _sets[i].getFillPattern() << endl;
    }

    for(unsigned int i = 0; i < _sets.size(); i++)
    {
        os << "@target G" << id << ".S" << i << endl;
	      _sets[i].write(os,i);
    }
}

// ///////////////////////////////////////////////////////////////////////////

void XmGraph::updateExtremas()
{
    if( _sets.back().getXmin() < _world_xmin )
        _world_xmin = _sets.back().getXmin();
    if( _sets.back().getXmax() > _world_xmax )
        _world_xmax = _sets.back().getXmax();
    if( _sets.back().getYmin() < _world_ymin )
        _world_ymin = _sets.back().getYmin();
    if( _sets.back().getYmax() > _world_ymax )
        _world_ymax = _sets.back().getYmax();
}

// ///////////////////////////////////////////////////////////////////////////

int XmGraph::sign(double a)
{
    if( a > 0.0 )
    {
        return +1;
    }
    else if ( a < 0.0 )
    {
        return -1;
    }
    else
    {
        return 0;
    }
}

// ///////////////////////////////////////////////////////////////////////////

void XmGraph::round_axis_limits(double *amin, double *amax, int scale)
{
    double smin, smax;
    int nrange;

    if (*amin == *amax) {
        switch (sign(*amin)) {
        case 0:
            *amin = -1.0;
            *amax = +1.0;
            break;
        case 1:
            *amin /= 2.0;
            *amax *= 2.0;
            break;
        case -1:
            *amin *= 2.0;
            *amax /= 2.0;
            break;
        }
    }

    smin = *amin;
    smax = *amax;

    if (sign(smin) == sign(smax)) {
        nrange = -rint(log10(fabs(2*(smax - smin)/(smax + smin))));
        nrange = MAX2(0, nrange);
    } else {
        nrange = 0;
    }
    smin = nicenum(smin, nrange, NICE_FLOOR);
    smax = nicenum(smax, nrange, NICE_CEIL);
    if (sign(smin) == sign(smax)) {
        if (smax/smin > 5.0) {
            smin = 0.0;
        } else if (smin/smax > 5.0) {
            smax = 0.0;
        }
    }

    *amin = smin;
    *amax = smax;
}

// ///////////////////////////////////////////////////////////////////////////

double XmGraph::nicenum(double x, int nrange, int round)
{
    int xsign;
    double f, y, fexp, rx, sx;

    if (x == 0.0) {
        return(0.0);
    }

    xsign = sign(x);
    x = fabs(x);

    fexp = floor(log10(x)) - nrange;
    sx = x/pow(10.0, fexp)/10.0;            /* scaled x */
    rx = floor(sx);                         /* rounded x */
    f = 10*(sx - rx);                       /* fraction between 0 and 10 */

    if ((round == NICE_FLOOR && xsign == +1) ||
        (round == NICE_CEIL  && xsign == -1)) {
        y = (int) floor(f);
    } else if ((round == NICE_FLOOR && xsign == -1) ||
               (round == NICE_CEIL  && xsign == +1)) {
	y = (int) ceil(f);
    } else {    /* round == NICE_ROUND */
	if (f < 1.5)
	    y = 1;
	else if (f < 3.)
	    y = 2;
	else if (f < 7.)
	    y = 5;
	else
	    y = 10;
    }

    sx = rx + (double) y/10.0;

    return (xsign*sx*10.0*pow(10.0, fexp));
}

// ///////////////////////////////////////////////////////////////////////////

void XmGraph::auto_ticks(double tmpmin, double tmpmax, XmAxis& axis)
{
    axis.setTickMajor( nicenum((tmpmax-tmpmin)/(axis.getTickDefault()-1), 0, NICE_ROUND) );
    if (axis.getTickMinorTicks() < 0 || axis.getTickMinorTicks() > 10)
	     axis.setTickMinorTicks(1);
}

// ///////////////////////////////////////////////////////////////////////////

} // end namespace xmlib
