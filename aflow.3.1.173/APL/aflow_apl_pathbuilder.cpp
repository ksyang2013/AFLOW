#include "aflow_apl.h"

using namespace std;

namespace apl {

  // ///////////////////////////////////////////////////////////////////////////

  PathBuilder::PathBuilder() {
    this->clear();
  }

  PathBuilder::PathBuilder(ModeEnumType mode) {
    this->clear();
    setMode(mode);
  }

  PathBuilder::~PathBuilder() {
    this->clear();
  }

  // ///////////////////////////////////////////////////////////////////////////

  void PathBuilder::clear() {
    _mode = SINGLE_POINT_MODE;
    _store = CARTESIAN_LATTICE;
    _path.clear();
    _points.clear();
    _labels.clear();
    _pointsVectorDimension = 0;
    _pointsVectorStartingIndex = 0;
    _nPointsPerSubPath = 0;
  }

  // ///////////////////////////////////////////////////////////////////////////

  void PathBuilder::setMode(ModeEnumType mode) {
    _mode = mode;
  }

  // ///////////////////////////////////////////////////////////////////////////

  void PathBuilder::setStore(StoreEnumType store) {
    _store = store;
  }

  // ///////////////////////////////////////////////////////////////////////////

  void PathBuilder::addPoint(const string& l, int dim, ...) {
    va_list arguments;
    xvector<double> point(dim,1);

    va_start(arguments, dim);
    for(int i = 0; i < dim; i++) {
      point(i+1) = va_arg(arguments, double);
    }
    va_end(arguments);
    //
    addPoint(l, point);
  }

  // ///////////////////////////////////////////////////////////////////////////

  void PathBuilder::addPoint(const string& l, const xvector<double>& p) {
    if( _points.empty() ) {
      _pointsVectorDimension = p.rows;
      _pointsVectorStartingIndex = p.lrows;
    }

    if( p.rows != _pointsVectorDimension ) {
      throw APLRuntimeError("apl::PathBuilder::addPoint(); Wrong dimension of the point.");
    }

    if( p.lrows != _pointsVectorStartingIndex ) {
      xvector<double> pp(_pointsVectorDimension+_pointsVectorStartingIndex,_pointsVectorStartingIndex);
      for(int i = 0; i < _pointsVectorDimension; i++)
	pp(_pointsVectorStartingIndex+i) = p(p.lrows+i);
      _points.push_back(pp);
    } else {
      _points.push_back(p);
    }

    _labels.push_back(l);
  }

  // ///////////////////////////////////////////////////////////////////////////

  int PathBuilder::getDensity() {
    return _nPointsPerSubPath;
  }

  // ///////////////////////////////////////////////////////////////////////////

  void PathBuilder::setDensity(int n) {
    if( n < 0 ) {
      throw APLRuntimeError("apl::PathBuilder::setDensity(); The density should be >= 0.");
    }
    _nPointsPerSubPath = n;
  }

  // ///////////////////////////////////////////////////////////////////////////

  void PathBuilder::buildPath() {
    // Remove the old path
    _path.clear();

    // Quick solution...
    if( _points.empty() ) {
      throw APLRuntimeError("apl::PathBuilder::buildPath; There are no points.");
    };

    // Create path in the SINGLE_POINT_MODE
    if( _mode == SINGLE_POINT_MODE ) {
      for(uint i = 1; i < _points.size(); i++) {
	xvector<double> startPoint = _points[i-1];
	_path.push_back(startPoint);
	if( _nPointsPerSubPath != 0 ) {
	  xvector<double> dPoint = ( _points[i] - _points[i-1] ) * ( 1.0 / (_nPointsPerSubPath) );
	  for(int j = 1; j <= (int)_nPointsPerSubPath; j++)
	    {
	      xvector<double> p = startPoint + j*dPoint;
	      _path.push_back(p);
	    }
	}
      }
      _path.push_back(_points[_points.size()-1]);
    }

    // Create path in the COUPLE_POINT_MODE
    if( _mode == COUPLE_POINT_MODE ) {
      // If the number of points is odd -> escape
      if( _points.size() % 2 != 0 ) {
	throw APLRuntimeError("apl::PathBuilder::buildPath(); The number of points is odd.");
      }

      xvector<double> startPoint(3), endPoint(3);
      endPoint = _points[0];
      for(uint i = 0; i < _points.size(); i+=2) {
	startPoint = _points[i];
	//if( aurostd::modulus( startPoint - endPoint ) > _AFLOW_APL_EPS_ )
	//    _path.push_back(endPoint);
	_path.push_back(startPoint);
	endPoint = _points[i+1];

	if( _nPointsPerSubPath != 0 ) {
	  xvector<double> dPoint = ( endPoint - startPoint ) * ( 1.0 / (_nPointsPerSubPath) );
	  for(int j = 1; j <= (int)_nPointsPerSubPath; j++) {
	      xvector<double> p = startPoint + j*dPoint;
	      _path.push_back(p);
	    }
	}
      }
      //_path.push_back(endPoint);
    }
  }

  // ///////////////////////////////////////////////////////////////////////////

  uint PathBuilder::getPathSize() {
    return _path.size();
  }

  // ///////////////////////////////////////////////////////////////////////////

  uint PathBuilder::getPointSize() {
    if( _mode == SINGLE_POINT_MODE ) {
      return( _points.size() );
    }

    if( _mode == COUPLE_POINT_MODE ) {
      return( ( _points.size() / 2 ) + 1 );
    }

    throw APLRuntimeError("apl::PathBuilder::getPointSize(); Unknown mode.");
  }

  // ///////////////////////////////////////////////////////////////////////////

  double PathBuilder::getPathLength() {
    // Quick solution...
    if( _points.empty() ) return 0;

    //
    uint npaths;
    if( _mode == SINGLE_POINT_MODE )
      npaths = _points.size()-1;
    else if( _mode == COUPLE_POINT_MODE )
      npaths = _points.size()/2;
    else
      throw APLRuntimeError("apl::PathBuilder::getPointLength(); Unknown mode.");

    //
    double length = 0.0;

    for(uint i = 1; i <= npaths; i++) {
      length += getPathLength(i);
    }

    return length;
  }

  // ///////////////////////////////////////////////////////////////////////////

  double PathBuilder::getPathLength(uint i) {
    if( i <= 0 ) {
      throw APLRuntimeError("apl::PathBuilder::getPathLength(); Wrong index. The index has to start from 1.");
    }

    // Quick solution 1...
    if( _points.empty() ) return 0;

    double length;
    if( _mode == SINGLE_POINT_MODE ) {
      // Quick solution 2...
      if( i > _points.size() ) {
	throw APLRuntimeError("apl::PathBuilder::getPathLength(); Wrong index.");
      }
      if( _store == RECIPROCAL_LATTICE ) {
	length = aurostd::modulus( F2C(trasp(reciprocalLattice),_points[i]) - F2C(trasp(reciprocalLattice),_points[i-1]) );
      }
      else
        {
	  length = aurostd::modulus( _points[i] - _points[i-1] );
        }
      return length;
    }
    else if( _mode == COUPLE_POINT_MODE ) {
      // Quick solution 2...
      if( i > _points.size()/2 ) {
	throw APLRuntimeError("apl::PathBuilder::getPathLength(); Wrong index.");
      }
      if( _store == RECIPROCAL_LATTICE ) {
	length = aurostd::modulus( F2C(trasp(reciprocalLattice),_points[(i*2)-1]) - F2C(trasp(reciprocalLattice),_points[(i-1)*2]) );
      } else {
	length = aurostd::modulus( _points[(i*2)-1] - _points[(i-1)*2] );
      }
      return length;
    }
    
    throw APLRuntimeError("apl::PathBuilder::getPathLength(); Unknown mode.");
  }
  
  // ///////////////////////////////////////////////////////////////////////////
  
  xvector<double> PathBuilder::getPoint(uint i) {
    if( i <= 0 ) {
      throw APLRuntimeError("apl::PathBuilder::getPoint(); Wrong index. The index has to start from 1.");
    }
    if( _mode == SINGLE_POINT_MODE ) {
      if( i > _points.size() ) {
	throw APLRuntimeError("apl::PathBuilder::getPoint(); Wrong index.");
      }
      return _points[i-1];
    } else if( _mode == COUPLE_POINT_MODE ) {
      if( i > (_points.size()/2)+1 ) {
	throw APLRuntimeError("apl::PathBuilder::getPoint(); Wrong index.");
      }
      if( i == 1 ) return _points[0];
      if( i == (_points.size()/2)+1 ) return _points[_points.size()-1];
      return _points[2*i-3];
    }
    
    throw APLRuntimeError("apl::PathBuilder::getPoint(); Unknown mode.");
  }
  
  // ///////////////////////////////////////////////////////////////////////////
  
  string PathBuilder::getPointLabel(uint i) {
    if( i <= 0 ) {
      throw APLRuntimeError("apl::PathBuilder::getPointLabel(); Wrong index. The index has to start from 1.");
    }
    
    if( _mode == SINGLE_POINT_MODE ) {
      if( i > _labels.size() ) {
	throw APLRuntimeError("apl::PathBuilder::getPointLabel(); Wrong index.");
      }
      return _labels[i-1];
    } else if( _mode == COUPLE_POINT_MODE ) {
      if( i > (_labels.size()/2)+1 ) {
	throw APLRuntimeError("apl::PathBuilder::getPointLabel(); Wrong index.");
      }
      if( i == 1 ) return _labels[0];
      if( i == (_labels.size()/2)+1 ) return _labels[_labels.size()-1];
      
      int h = 2*i-2;
      int l = 2*i-3;
      
      if( _labels[l] == _labels[h] )
	return _labels[l];
      else
	return(_labels[l]+"|"+_labels[h]);
    }
    
    throw APLRuntimeError("apl::PathBuilder::getPointLabel(); Unknown mode.");
  }
  
  // ///////////////////////////////////////////////////////////////////////////
  
  vector< aurostd::xvector<double> > PathBuilder::getPath(ModeEnumType mode, const string& userPath) {
    vector< xvector<double> > new_points;
    vector<string> new_labels;
    vector<string> tokens;

    // Generate users path by points
    tokenize(userPath,tokens,"|-");
    for(uint i = 0; i < tokens.size(); i++) {
      uint j = 0;
      for( ; j < _labels.size(); j++) {
	// This works fine, but notorious changes of the Gammma point label
	// in aflow_kpoints.cpp lead me to do it like: G, Gamma, of \Gamma or
	// or anything like GAMMA \G or GaMmA or I do not know what else is
	// still a gamma point labeled as G!
	//if( _labels[j] == tokens.at(i) ) break;
	if( _labels[j].find(tokens.at(i)) != string::npos ) break;
      }
      if( j == _labels.size() ) {
	throw APLRuntimeError("apl::PathBuilder::getPath(); Undefined label of the point.");
      }

      new_points.push_back(_points[j]);      
      //new_labels.push_back(_labels[j]);
      // I will use my labels...
      new_labels.push_back(tokens.at(i));
    }

    // update
    _points.clear();
    _labels.clear();
    for(uint i = 0; i < new_points.size(); i++) {
      _points.push_back(new_points[i]);
      _labels.push_back(new_labels[i]);
    }
    new_points.clear();
    new_labels.clear();

    // Build full path
    _mode = mode;
    buildPath();

    //
    return _path;
  }

  // ///////////////////////////////////////////////////////////////////////////

  std::vector< aurostd::xvector<double> > PathBuilder::getPath() {
    buildPath();
    return _path;
  }

  // ///////////////////////////////////////////////////////////////////////////

  void PathBuilder::transform(const aurostd::xmatrix<double>& m) {
    for(uint i = 0; i < _points.size(); i++) {
      //cout << "IN: " << _labels[i] << ": " << _points[i] << std::endl;
      _points[i] = m * _points[i];
      //cout << "OUT: " << _labels[i] << ": " << _points[i] << std::endl;
    }
  }

  // ///////////////////////////////////////////////////////////////////////////

  void PathBuilder::tokenize(const string& str,vector<string>& tokens, string del) {
    string delimiters = del;

    // Skip delimiters at beginning.
    string::size_type lastPos = str.find_first_not_of(delimiters, 0);

    // Find first "non-delimiter".
    string::size_type pos = str.find_first_of(delimiters, lastPos);

    while (string::npos != pos || string::npos != lastPos) {
      // Found a token, add it to the vector.
      tokens.push_back(str.substr(lastPos, pos - lastPos));

      // Skip delimiters. Note the "not_of"
      lastPos = str.find_first_not_of(delimiters, pos);

      // Find next "non-delimiter"
      pos = str.find_first_of(delimiters, lastPos);
    }
  }

  // ///////////////////////////////////////////////////////////////////////////

  void PathBuilder::pointsAreDefinedFor(const xstructure& primitiveStructure, StoreEnumType store) {
    // Transform from the reciprocal lattice of the primitive cell to cartesian coordinates
    if( store == RECIPROCAL_LATTICE ) {
      transform( trasp(ReciprocalLattice(primitiveStructure.lattice)) );
    }

    //
    _store = store;
  }

  // ///////////////////////////////////////////////////////////////////////////
  
  /*
    void PathBuilder::transformPointsFor(const xstructure& supercellStructure, StoreEnumType store) {
        if( store == CARTESIAN_LATTICE ) {
      if( _store == RECIPROCAL_LATTICE )
      transform( trasp(ReciprocalLattice(supercellStructure.lattice)) );
      else
      }
      else if( store == RECIPROCAL_LATTICE ) {
      transform( inverse(trasp(ReciprocalLattice(supercellStructure.lattice))) );
      }

      _store = store;
      cartesianLattice = supercellStructure.lattice;
      reciprocalLattice = ReciprocalLattice(supercellStructure.lattice);
  
   }
  */

  // ///////////////////////////////////////////////////////////////////////////

  void PathBuilder::takeAflowElectronicPath(const string& latticeName,
					    const xstructure& pcStructure,
					    const xstructure& scStructure) {
    // Get path from electronic structure...
    stringstream fileKPOINTS;
    // input is not the reciprocal lattice!
    bool foundBZ;
    fileKPOINTS << LATTICE::KPOINTS_Directions(latticeName,pcStructure.lattice,10,pcStructure.iomode,foundBZ);
    fileKPOINTS.flush();

    if( !foundBZ )
      throw APLRuntimeError("apl::PathBuilder::takeAflowElectronicPath(); The BZ not found for this lattice.");

    //   cerr << fileKPOINTS.str() << std::endl;
    string line;
    vector<string> tokens;
    int nLine = 0;
    while( getline( fileKPOINTS, line ) ) {
      if( nLine++ < 4 ) continue;
      if( line.empty() ) continue;
      if( line.size() == 1 ) continue;
      //cout << line << std::endl;
      tokenize(line,tokens,string(" !"));
      addPoint(tokens[3],3,atof(tokens[0].c_str()),atof(tokens[1].c_str()),atof(tokens[2].c_str()));
      tokens.clear();
      line.clear();
    }
    fileKPOINTS.clear();

    // Transform from the reciprocal lattice of the primitive cell to
    // cartesian coordinates (note; in klattice vectors are rows! )
    transform( trasp(ReciprocalLattice(pcStructure.lattice)) );

    //
    cartesianLattice = scStructure.lattice;
    reciprocalLattice = ReciprocalLattice(scStructure.lattice);

    // set mode
    setMode(COUPLE_POINT_MODE);
    _store = CARTESIAN_LATTICE;
  }

  // ///////////////////////////////////////////////////////////////////////////

} // namespace apl
