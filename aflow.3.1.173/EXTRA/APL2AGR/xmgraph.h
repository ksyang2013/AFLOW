#ifndef XMLIB_XMGRAPH_H
#define XMLIB_XMGRAPH_H

#include <string>
#include <vector>
#include "xmlib.h"

#define NICE_FLOOR   0
#define NICE_CEIL    1
#define NICE_ROUND   2

#define MAX2(a, b) (((a) > (b)) ? a : b)

namespace xmlib {

class XmGraph {

  private:
    std::vector<XmSet> _sets;
    std::vector<XmLine> _lines;
    
    XmAxis _xaxis;
    XmAxis _yaxis;
    XmAxis _altxaxis;
    XmAxis _altyaxis;
    
    double _world_xmin;
    double _world_xmax;
    double _world_ymin;
    double _world_ymax;

  public:
    double _view_xmin;
    double _view_xmax;
    double _view_ymin;
    double _view_ymax;
    
  private:
    std::string _title;
    std::string _subtitle;
    int _titleFont;
    int _subtitleFont;
    double _titleSize;
    double _subtitleSize;
    int _titleColor;
    int _subtitleColor;
    
    bool _hidden;
            
  public:
    XmGraph();
    ~XmGraph();
    void clear();
    
    void addSets(const std::vector<XmSet>&);
    void addSet(const XmSet&);

    void addLine(const XmLine&);
    std::vector<XmLine>& getLines();
    XmLine& getLine(int);

    XmGraph& operator= (const XmGraph&);

    void setXAxis(const XmAxis&);
    void setYAxis(const XmAxis&);
    void setAltXAxis(const XmAxis&);
    void setAltYAxis(const XmAxis&);
    void write(std::ostream&, int);

    void setView(double,double,double,double);

    double getWorldXMin();
    double getWorldYMin();
    double getWorldXMax();
    double getWorldYMax();
    void setWorldXMax(double);
    void setWorldXMin(double);
    void setWorldYMax(double);
    void setWorldYMin(double);
    
    bool isHidden();
    void hide();
    void show();
    
    void setTitle(const char*);
    void setTitleFont(int);
    void setTitleSize(double);
    void setSubTitle(const char*);
    void setSubTitleFont(int);
    void setSubTitleSize(double);

    const XmAxis& getXAxis();
    const XmAxis& getYAxis();

    void autoscale();
    void autoscaleX();
    void autoscaleY();
    void autoTicksX();
    void autoTicksY();

  private:
    void updateExtremas();
    void initAtributes();
    int sign(double);
    void round_axis_limits(double*, double*, int);
    double nicenum(double, int, int);
    void auto_ticks(double, double, XmAxis& );
    
};

} // end namespace xmlib

#endif
