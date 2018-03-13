#ifndef XMLIB_XMSET_H
#define XMLIB_XMSET_H

#include <string>
#include <vector>
#include "xmlib.h"

namespace xmlib {

class XmSet {

  friend class XmGraph;

  private:
    std::vector<double> _x;
    std::vector<double> _y;

    int _xmin;
    int _xmax;
    int _ymin;
    int _ymax;

    int _lineType;
    int _lineStyle;
    int _lineColor;
    int _linePattern;

    double _lineWidth;

    int _fillType;
    int _fillColor;
    int _fillPattern;

    std::string _comment;

  private:
    void updateExtremas();
    void updateExtremas(int);
    void setDefault();

  public:
    XmSet();
   ~XmSet();

    void clear();

    void addXY(double,double);
    void setXY(unsigned int,double,double);

    XmSet& operator= (const XmSet&);
    void write(std::ostream&,int);

    void setLineType(int);
    void setLineStyle(int);
    void setLineWidth(double);
    void setLineColor(int);
    void setLinePattern(int);
    void setComment(const char*);
    void setComment(const std::string&);
    void setFillType(int);
    void setFillColor(int);
    void setFillPattern(int);

    int getLineType();
    int getLineStyle();
    double getLineWidth();
    int getLineColor();
    int getLinePattern();
    const std::string& getComment();
    int getFillType();
    int getFillColor();
    int getFillPattern();

    double getXmin();
    double getXmax();
    double getYmin();
    double getYmax();

};

} // end namespace xmlib

#endif
