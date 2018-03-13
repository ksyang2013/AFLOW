#ifndef XMLIB_XMAXIS_H
#define XMLIB_XMAXIS_H

#include <vector>
#include <string>
#include "xmlib.h"

namespace xmlib {

class XmAxis {
    
  private:
    std::string _id;
    bool _active; 
    
    std::string _label;
    double _labelSize;
    int _labelFont;
    int _labelColor;
    std::string _label_place;

    bool _activeTicks;
    double _tick_major;
    int _tick_minor_ticks;
    int _tick_autonum;
    std::string _tick_pointing;
    std::string _tick_drawon;
    
    bool _activeTickLabel;
    std::string _ticklabel_side;
    double _ticklabelSize;
    int _ticklabelFont;
    int _ticklabelColor;

    double _tick_minor_size;
    int _tick_minor_color;
    double _tick_minor_linewidth;
    int _tick_minor_linestyle;
    bool _activeTickMinorGrid;

    double _tick_major_size;
    int _tick_major_color;
    double _tick_major_linewidth;
    int _tick_major_linestyle;
    bool _activeTickMajorGrid;
    
    std::string _tickspec_type;
    std::vector<int> _tickspec0;
    std::vector<int> _tickspec1;
    std::vector<double> _tickspec2;
    std::vector<std::string> _tickspec3;

  private:
    void setAtributes();
        
  public:
    XmAxis();
    XmAxis(const char*);
    ~XmAxis();
    void clear();
    
    void setActiveON();
    void setActiveOFF();
    void setID(const char*);
    
    void setLabel(const char*);
    void setLabelSize(double);
    void setLabelFont(int);
    void setLabelColor(int);
    void setLabelPlaceNormal();
    void setLabelPlaceOpposite();
    void setLabelPlaceBoth();
    
    void setTicksON();
    void setTicksOFF();
    void setTicksPointingOut();
    void setTicksPointingIn();
    void setTicksPointingBoth();
    void setTicksDrawOnNormalSide();
    void setTicksDrawOnOppositeSide();
    void setTicksDrawOnBothSides();
    
    void setTicksLabelON();
    void setTicksLabelOFF();
    void setTicksLabelSideNormal();
    void setTicksLabelSideOpposite();
    void setTicksLabelSideBoth();    
    void setTicksLabelSize(double);
    void setTicksLabelFont(int);
    void setTicksLabelColor(int);
    
    void setSpecTickTypeNone();
    void setSpecTickTypeTicks();
    void setSpecTickTypeBoth();
    void addSpecTick(int,int,double,const char*);
    
    void setTickMajor(double);
    double getTickMajor();
    void setTickMinorTicks(int);
    int getTickMinorTicks();
    int getTickDefault();

    void setTickMajorSize(double);
    void setTickMajorColor(int);
    void setTickMajorLineWidth(double);
    void setTickMajorLineStyle(int);
    void setTickMajorGridON();
    void setTickMajorGridOFF();

    void setTickMinorSize(double);
    void setTickMinorColor(int);
    void setTickMinorLineWidth(double);
    void setTickMinorLineStyle(int);
    void setTickMinorGridON();
    void setTickMinorGridOFF();
    
    XmAxis& operator= (const XmAxis&);
    void write(std::ostream&);        
};

} // end namespace xmlib

#endif
