#ifndef XMGRACEFILE_H
#define XMGRACEFILE_H

#include <vector>
#include <string>

#include "xmlib.h"

namespace xmlib {

class XmGraceFile {
    
private:
    PageSize _pageSize;
    std::vector<XmGraph> _graphs;
        
public:
    XmGraceFile();
    ~XmGraceFile();

    void clear();
    void addGraph(const XmGraph&);
    void setPageSize(const PageSize&);
    const PageSize& getPageSize();
    void saveAs(const char*);
};

} // end of namespace xmlib

#endif
