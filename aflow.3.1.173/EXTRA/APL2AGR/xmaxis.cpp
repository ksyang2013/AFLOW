#include <iostream>
#include "xmlib.h"

using namespace std;

namespace xmlib {

// ///////////////////////////////////////////////////////////////////////////

XmAxis::XmAxis(const char *s)
{
    setAtributes();
    setID(s);
}

// ///////////////////////////////////////////////////////////////////////////

XmAxis::XmAxis()
{
    setAtributes();
}

// ///////////////////////////////////////////////////////////////////////////

XmAxis::~XmAxis()
{
    clear();
}

// ///////////////////////////////////////////////////////////////////////////

void XmAxis::clear()
{
    _tickspec0.clear();
    _tickspec1.clear();
    _tickspec2.clear();
    _tickspec3.clear();
}

// ///////////////////////////////////////////////////////////////////////////

void XmAxis::setAtributes()
{
    clear();

    _id = string("");
    _active = true;
    _label = string("");
    _labelSize = 1.0;
    _labelFont = 0;
    _labelColor = 1;
    setLabelPlaceNormal();

    setTicksON();
    _tick_major = 0.2;
    _tick_minor_ticks = 1;
    _tick_autonum = 6;
    setTicksPointingIn();
    setTicksDrawOnBothSides();

    _tick_major_size = 0.50;
    _tick_major_color = 1;
    _tick_major_linewidth = 1.0;
    _tick_major_linestyle = 1;
    setTickMajorGridOFF();

    _tick_minor_size = 0.25;
    _tick_minor_color = 1;
    _tick_minor_linewidth = 1.0;
    _tick_minor_linestyle = 1;
    setTickMinorGridOFF();
   
    setTicksLabelON();
    setTicksLabelSideNormal();
    _ticklabelSize = 1.0;
    _ticklabelFont = 0;
    _ticklabelColor = 1;
    
    setSpecTickTypeNone();
}

// ///////////////////////////////////////////////////////////////////////////

XmAxis& XmAxis::operator= (const XmAxis& that)
{
    if( this != &that )
    {
        _id = that._id;
        _active = that._active;
        _label = that._label;
        _labelSize = that._labelSize;
        _labelFont = that._labelFont;
        _labelColor = that._labelColor;
        _label_place = that._label_place;
        
        _activeTicks = that._activeTicks;
        _tick_major = that._tick_major;
        _tick_minor_ticks = that._tick_minor_ticks;
        _tick_pointing = that._tick_pointing;
        _tick_drawon = that._tick_drawon;
        _tick_autonum = that._tick_autonum;
        
        _activeTickLabel = that._activeTickLabel;
        _ticklabel_side = that._ticklabel_side;
        _ticklabelSize = that._ticklabelSize;
        _ticklabelFont = that._ticklabelFont;
        _ticklabelColor = that._ticklabelColor;

        _tick_major_size = that._tick_major_size;
        _tick_major_color = that._tick_major_color;
        _tick_major_linewidth = that._tick_major_linewidth;
        _tick_major_linestyle = that._tick_major_linestyle;
        _activeTickMajorGrid = that._activeTickMajorGrid;

        _tick_minor_size = that._tick_minor_size;
        _tick_minor_color = that._tick_minor_color;
        _tick_minor_linewidth = that._tick_minor_linewidth;
        _tick_minor_linestyle = that._tick_minor_linestyle;
        _activeTickMinorGrid = that._activeTickMinorGrid;

        _tickspec_type = that._tickspec_type;
        _tickspec0.clear();
        _tickspec1.clear();
        _tickspec2.clear();
        _tickspec3.clear();
        _tickspec0 = that._tickspec0;
        _tickspec1 = that._tickspec1;
        _tickspec2 = that._tickspec2;
        _tickspec3 = that._tickspec3;
    }
    
    return *this;
}

// ///////////////////////////////////////////////////////////////////////////

void XmAxis::setActiveON()
{
    _active = true;
}

void XmAxis::setActiveOFF()
{
    _active = false;
}

// ///////////////////////////////////////////////////////////////////////////

void XmAxis::setID(const char* s)
{
    _id = string(s);
}

// ///////////////////////////////////////////////////////////////////////////

void XmAxis::setLabel(const char* s)
{
    _label = string(s);
}

void XmAxis::setLabelSize(double h)
{
    _labelSize = h;
}

void XmAxis::setLabelFont(int f)
{
    _labelFont = f;
}

void XmAxis::setLabelColor(int c)
{
    _labelColor = c;
}

void XmAxis::setLabelPlaceNormal()
{
    _label_place = string("normal");
}

void XmAxis::setLabelPlaceOpposite()
{
    _label_place = string("opposite");
}

void XmAxis::setLabelPlaceBoth()
{
    _label_place = string("both");
}

// ///////////////////////////////////////////////////////////////////////////

void XmAxis::setTicksON()
{
    _activeTicks = true;
}

void XmAxis::setTicksOFF()
{
    _activeTicks = false;
}

// ///////////////////////////////////////////////////////////////////////////

void XmAxis::setTicksPointingOut()
{
    _tick_pointing = string("out");
}

void XmAxis::setTicksPointingIn()
{
    _tick_pointing = string("in");
}

void XmAxis::setTicksPointingBoth()
{
    _tick_pointing = string("both");
}

// ///////////////////////////////////////////////////////////////////////////

void XmAxis::setTicksDrawOnNormalSide()
{
    _tick_drawon = string("normal");
}

void XmAxis::setTicksDrawOnOppositeSide()
{
    _tick_drawon = string("opposite");
}

void XmAxis::setTicksDrawOnBothSides()
{
    _tick_drawon = string("both");
}

// ///////////////////////////////////////////////////////////////////////////

void XmAxis::setTicksLabelON()
{
    _activeTickLabel = true;
}

void XmAxis::setTicksLabelOFF()
{
    _activeTickLabel = false;
}

// ///////////////////////////////////////////////////////////////////////////

void XmAxis::setTicksLabelSideNormal()
{
    _ticklabel_side = string("normal");
}

void XmAxis::setTicksLabelSideOpposite()
{
    _ticklabel_side = string("opposite");
}

void XmAxis::setTicksLabelSideBoth() 
{
    _ticklabel_side = string("both");
}

void XmAxis::setTicksLabelSize(double h)
{
    _ticklabelSize = h;
}

void XmAxis::setTicksLabelFont(int f)
{
    _ticklabelFont = f;
}

void XmAxis::setTicksLabelColor(int c)
{
    _ticklabelColor = c;
}

// ///////////////////////////////////////////////////////////////////////////

void XmAxis::setSpecTickTypeNone()
{
    _tickspec_type = string("none");
}

void XmAxis::setSpecTickTypeTicks()
{
    _tickspec_type = string("ticks");
}
    
void XmAxis::setSpecTickTypeBoth()
{
    _tickspec_type = string("both");
}

// ///////////////////////////////////////////////////////////////////////////

void XmAxis::addSpecTick(int mm, int i, double x ,const char* l)
{
    _tickspec0.push_back(mm);
    _tickspec1.push_back(i);
    _tickspec2.push_back(x);
    _tickspec3.push_back(string(l));
}

// ///////////////////////////////////////////////////////////////////////////

void XmAxis::setTickMajor(double d)
{
    _tick_major = d;
}

double XmAxis::getTickMajor()
{
    return _tick_major;
}

// ///////////////////////////////////////////////////////////////////////////

void XmAxis::setTickMinorTicks(int n)
{
    _tick_minor_ticks = n;
}

int XmAxis::getTickMinorTicks()
{
    return _tick_minor_ticks;
}

// ///////////////////////////////////////////////////////////////////////////

void XmAxis::setTickMajorSize(double size)
{
    _tick_major_size = size;
}

void XmAxis::setTickMajorColor(int colorID)
{
    _tick_major_color = colorID;
}

void XmAxis::setTickMajorLineWidth(double width)
{
    _tick_major_linewidth = width;
}

void XmAxis::setTickMajorLineStyle(int styleID)
{
    _tick_major_linestyle = styleID;
}

void XmAxis::setTickMajorGridON()
{
    _activeTickMajorGrid = true;
}

void XmAxis::setTickMajorGridOFF()
{
    _activeTickMajorGrid = false;
}

// ///////////////////////////////////////////////////////////////////////////

void XmAxis::setTickMinorSize(double size)
{
    _tick_minor_size = size;
}

void XmAxis::setTickMinorColor(int colorID)
{
    _tick_minor_color = colorID;
}

void XmAxis::setTickMinorLineWidth(double width)
{
    _tick_minor_linewidth = width;
}

void XmAxis::setTickMinorLineStyle(int styleID)
{
    _tick_minor_linestyle = styleID;
}

void XmAxis::setTickMinorGridON()
{
    _activeTickMinorGrid = true;
}

void XmAxis::setTickMinorGridOFF()
{
    _activeTickMinorGrid = false;
}

// ///////////////////////////////////////////////////////////////////////////

int XmAxis::getTickDefault()
{
    return _tick_autonum;
}

// ///////////////////////////////////////////////////////////////////////////

void XmAxis::write(std::ostream& os)
{
    os << "@    " << _id;
    if( !_active ) 
	     os << "  off" << endl;
    else
    {
       os << "  on" << endl;
       os << "@    " << _id << "  type zero false" << endl;
       os << "@    " << _id << "  offset 0.000000 , 0.000000" << endl;
       os << "@    " << _id << "  bar on" << endl;
       os << "@    " << _id << "  bar color 1" << endl;
       os << "@    " << _id << "  bar linestyle 1" << endl;
       os << "@    " << _id << "  bar linewidth 1.0" << endl;
        
       os << "@    " << _id << "  label \"" << _label << "\"" << endl;
       os << "@    " << _id << "  label layout para" << endl;
       os << "@    " << _id << "  label place auto" << endl;
       os << "@    " << _id << "  label char size " << _labelSize/100.0 << endl;
       os << "@    " << _id << "  label font " << _labelFont << endl;
       os << "@    " << _id << "  label color " << _labelColor << endl;
       os << "@    " << _id << "  label place " << _label_place << endl;
        
       os << "@    " << _id << "  tick  " << (_activeTicks ? "on":"off") << endl;
       os << "@    " << _id << "  tick major " << _tick_major << endl;
       os << "@    " << _id << "  tick minor ticks " << _tick_minor_ticks << endl;
       os << "@    " << _id << "  tick default " << _tick_autonum << endl;
       os << "@    " << _id << "  tick place rounded true" << endl;
       os << "@    " << _id << "  tick " << _tick_pointing << endl;
       os << "@    " << _id << "  tick place " << _tick_drawon << endl;
        
       os << "@    " << _id << "  tick major size " << _tick_major_size << endl;
       os << "@    " << _id << "  tick major color " << _tick_major_color << endl;
       os << "@    " << _id << "  tick major linewidth " << _tick_major_linewidth << endl;
       os << "@    " << _id << "  tick major linestyle " << _tick_major_linestyle << endl;
       os << "@    " << _id << "  tick major grid " << (_activeTickMajorGrid ? "on":"off") << endl;
        
       os << "@    " << _id << "  tick minor size " << _tick_minor_size << endl;
       os << "@    " << _id << "  tick minor color " << _tick_minor_color << endl;
       os << "@    " << _id << "  tick minor linewidth " << _tick_minor_linewidth << endl;
       os << "@    " << _id << "  tick minor linestyle " << _tick_minor_linestyle << endl;
       os << "@    " << _id << "  tick minor grid " << (_activeTickMinorGrid ? "on":"off") << endl;
        
       os << "@    " << _id << "  ticklabel " << (_activeTickLabel ? "on":"off") << endl;
       os << "@    " << _id << "  ticklabel format general" << endl;
       os << "@    " << _id << "  ticklabel prec 5" << endl;
       os << "@    " << _id << "  ticklabel formula \"\"" << endl;
       os << "@    " << _id << "  ticklabel append \"\"" << endl;
       os << "@    " << _id << "  ticklabel prepend \"\"" << endl;
       os << "@    " << _id << "  ticklabel angle 0" << endl;
       os << "@    " << _id << "  ticklabel skip 0" << endl;
       os << "@    " << _id << "  ticklabel stagger 0" << endl;
       os << "@    " << _id << "  ticklabel place " << _ticklabel_side << endl;
       os << "@    " << _id << "  ticklabel offset auto" << endl;
       os << "@    " << _id << "  ticklabel offset 0.000000 , 0.010000" << endl;
       os << "@    " << _id << "  ticklabel start type auto" << endl;
       os << "@    " << _id << "  ticklabel start 0.000000" << endl;
       os << "@    " << _id << "  ticklabel stop type auto" << endl;
       os << "@    " << _id << "  ticklabel stop 0.000000" << endl;
       os << "@    " << _id << "  ticklabel char size " << _ticklabelSize/100.0 << endl;
       os << "@    " << _id << "  ticklabel font " << _ticklabelFont << endl;
       os << "@    " << _id << "  ticklabel color " << _ticklabelColor << endl;

       os << "@    " << _id << "  tick spec type " << _tickspec_type << endl;
       if( _tickspec_type != string("none") )
       {
          if( _tickspec0.size() != 0 )
          {
              os << "@    " << _id << "  tick spec " << _tickspec0.size() << endl;
              for(unsigned int i = 0; i < _tickspec0.size(); i++)
              {
                  os << "@    " << _id << "  tick ";
                  if( _tickspec0[i] == 0 )
                      os << "minor ";
                  else
                      os << "major ";
              os << _tickspec1[i] << ", " << _tickspec2[i] << endl;
              if( _tickspec0[i] == 1 )
                  os << "@    " << _id << "  ticklabel " << _tickspec1[i] << ", \"" << _tickspec3[i] << "\"" << endl;
              }
          }
       }
    }
}

} // end namespace xmlib
