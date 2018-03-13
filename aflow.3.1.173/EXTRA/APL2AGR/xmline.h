#ifndef XMLIB_XMLINE_H
#define XMLIB_XMLINE_H

#include <string>
#include <vector>
#include "xmlib.h"

namespace xmlib {

class XmLine {
    
 private:
    bool _active;
    std::string _locktype;
    double _startx;
    double _starty;
    double _endx;
    double _endy;
    double _linewidth;
    int _linestyle;
    int _linecolor;
    
 public:
    XmLine();
    ~XmLine();

    void setActiveON();
    void setActiveOFF();
    void setLockTypeWorld();
    void setLockTypeView();
    void setStartXY(double,double);
    void setEndXY(double,double);
    void setWidth(double);
    void setStyle(int);
    void setColor(int);
  
    XmLine& operator= (const XmLine&);
    void write(std::ostream&, int);
};

} // end namespace xmlib

#endif
