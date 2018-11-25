// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2015           *
// *                Aflow PINKU NATH - Duke University 2014-2016             *
// *                                                                         *
// ***************************************************************************
// Written by Pinku Nath
// pn49@duke.edu

#include "aflow_apl.h"
#include <iterator>

#define _isnegative(a) (a<MIN_EIGEN_TRESHOLD) ? true : false

#undef AFLOW_APL_MULTITHREADS_ENABLE

#if GCC_VERSION >= 40400   // added two zeros
#define AFLOW_APL_MULTITHREADS_ENABLE 1
#include <thread>
#else
#warning "The multithread parts of APL will be not included, since they need gcc 4.4 and higher (C++0x support)."
#endif
/*
  #  if __GNUC_PREREQ(4,9) // enable thread if gcc version is >= 4.9
  #define AFLOW_APL_MULTITHREADS_ENABLE
  #include <thread>
  #endif 
*/
namespace apl
{
  // ***************************************************************************************
  AtomicDisplacements::AtomicDisplacements(IPhononCalculator& pc, UniformMesh& mp, Logger& l):_pc(pc), _mp(mp), _logger(l)
  {
    _logger<<"Preparing setup for Quasi-harmonic Gruneisen calculation "<<apl::endl;
    clear();
    _nBranches=_pc.getNumberOfBranches();
    _CUTOFF_FREQ=0.00;
    _is_freq_negative=false;
  }
  // ***************************************************************************************
  AtomicDisplacements::~AtomicDisplacements(){this->clear();}
  // ***************************************************************************************
  void AtomicDisplacements::clear()
  {
    _atomic_masses_amu.clear();
    _atomic_species.clear();
    _rlattice.clear();
    _kpoints.clear();
    _weights.clear();
    _freq_Thz.clear();
    _DM.clear(); 
    _eigenvectors.clear();
    _atomic_c_positions.clear();
    _eigenvectors_path.clear();
    _kpoints_path.clear();
  }
  // ***************************************************************************************
  //initializing variables
  void AtomicDisplacements::populate_variables(const xstructure& xs)
  {
    _logger<<"Populationg variables to calculate AtomicDisplacements properties" <<apl::endl;
    xstructure xstr(xs);
    
    _atomic_masses_amu.clear();
    _atomic_species.clear();
    _atomic_masses_amu=_pc.get_ATOMIC_MASSES_AMU();
    for(uint i=0; i!=xstr.atoms.size(); i++){
      _atomic_species.push_back(xstr.atoms[i].name);
    }
    _rlattice.clear();
    _rlattice=_mp.get_rlattice();
    _klattice=_mp.get_klattice();
    _kpoints.clear();
    _kpoints=_mp.get_kpoints();
    _weights.clear();
    _weights=_mp.get_weights();

    _freq_Thz.clear();
    xvector<double> dv(_nBranches,1);
    _freq_Thz.resize(_kpoints.size(), dv);
    _DM.clear(); _eigenvectors.clear();
    xmatrix< xcomplex<double> > dd(_nBranches,_nBranches,1,1);
    _DM.resize(_kpoints.size(), dd);//at +ve volume
    _eigenvectors.resize(_kpoints.size(), dd);

    _atomic_c_positions.clear();
    xvector<double> tmp(3,1);
    _atomic_c_positions.resize(_atomic_masses_amu.size(), tmp);

    for (uint i=0;i<xstr.atoms.size(); i++)
      {
        _atomic_c_positions[i]=xstr.atoms[i].cpos;
      }
  }
  // ***************************************************************************************
  bool AtomicDisplacements::eigen_solver()
  {
    _freq_test.clear();
    bool pass=false;
#ifdef AFLOW_APL_MULTITHREADS_ENABLE
    // Get the number of CPUS
    int ncpus = sysconf(_SC_NPROCESSORS_ONLN);// AFLOW_MachineNCPUs;
    if(ncpus<1) ncpus=1;
//    int qpointsPerCPU = _kpoints.size() / ncpus;  OBSOLETE 180801
    _freq_test.resize(ncpus, false);
    // Show info 
    if( ncpus == 1 )
      {
        _logger.initProgressBar("Calculating eigenvalues for AtomicDisplacements");
      }
    else
      {
	_logger.initProgressBar("Calculating eigenvalues for AtomicDisplacements (" + stringify(ncpus) + " threads)");
      }

    // Distribute the calculation
    int startIndex, endIndex;
    std::vector< std::thread* > threads;
    vector<vector<int> > thread_dist = getThreadDistribution((int) _kpoints.size(), ncpus);
    for (int icpu = 0; icpu < ncpus; icpu++) {
      startIndex = thread_dist[icpu][0];
      endIndex = thread_dist[icpu][1];
      threads.push_back( new std::thread(&AtomicDisplacements::solve_eigenvalues_in_threads,this,startIndex,endIndex, icpu) );
    }

/* OBSOLETE ME 180801
    for(int icpu = 0; icpu < ncpus; icpu++) {
      startIndex = icpu * qpointsPerCPU;
      endIndex = startIndex + qpointsPerCPU;
      if( ( (uint)endIndex > _kpoints.size() ) ||
          ( ( icpu == ncpus-1 ) && ( (uint)endIndex < _kpoints.size() ) ) )
        endIndex = _kpoints.size();
      threads.push_back( new std::thread(&AtomicDisplacements::solve_eigenvalues_in_threads,this,startIndex,endIndex, icpu) );
    }
*/
    // Wait to finish all threads here!
    for(uint i = 0; i < threads.size(); i++) {
      threads[i]->join();
      delete threads[i];
    }

    // Done
    _logger.finishProgressBar();


    uint bool_size=0;
    for(uint id=0; id!=_freq_test.size(); id++)
      {
        bool_size+=(uint)_freq_test[id];
      }
    if(bool_size==(uint)ncpus)pass=true;
    _freq_test.clear();
#else
    _freq_test.resize(1, false);
    solve_eigenvalues_in_threads(0, (int)_kpoints.size(), 0);
    uint bool_size=0;
    for(uint id=0; id!=_freq_test.size(); id++)
      {
        bool_size+=(uint)_freq_test[id];
      }
    if(bool_size==1)pass=true;
#endif

    return pass;
  }
  // ***************************************************************************************
  void AtomicDisplacements::solve_eigenvalues_in_threads(int startIndex, int endIndex, int cpuid)
  {
    for(int In=startIndex;In<endIndex;In++)
      {
        xvector<double> qpoint;
        xmatrix<xcomplex<double> >  DM(_nBranches,_nBranches,1,1);

	qpoint=_kpoints[In];
	DM=_pc.getDynamicalMatrix(qpoint);
	_DM[In]=DM;

	//calculate eigenvalue and eigenvectors of a dynamical matirx
	xvector<double> eigenvalues(_nBranches, 1);
	xmatrix<xcomplex<double> > eigenvectors(_nBranches, _nBranches, 1,1);
	apl::aplEigensystems e;
	e.eigen_calculation(DM, eigenvalues, eigenvectors, APL_MV_EIGEN_SORT_VAL_ASC);

	for(uint j=1; j<=_nBranches; j++)
	  {
	    if(_iszero(eigenvalues[j]))eigenvalues[j]=0.0;
	    if(eigenvalues[j]<0){
	      if(_isnegative(eigenvalues[j]))
		{
		  if(!_is_freq_negative)
		    {
                      _logger<<  apl::warning <<"AtomicDisplacements:: negative eigenvalue: = " <<eigenvalues[j]<<apl::endl;
		    }
		  _is_freq_negative=true;
		  return;
		}else eigenvalues[j]=0.00;
	    }
	    _freq_Thz[In][j]=sqrt(eigenvalues[j])*RAW2Hz;
	  }
	_eigenvectors[In]=eigenvectors;
      }
    _freq_test[cpuid]=true;
  }//fn end
  // ***************************************************************************************
  //thermal displacements calculation for a given range of temperature 
  //formulas can be found in the link http://atztogo.github.io/phonopy/thermal-displacement.html
  void AtomicDisplacements::thermal_displacements(double Ts, double Te, int Tinc)
  {
    _logger<<"calculating mean square displacements "<<apl::endl;
    if(_atomic_masses_amu.size()!=_atomic_species.size())
      {
        _logger<<apl::error <<"error in atomic masses [remove apl.xml LOCK and run again] "<<apl::endl;
        return;
      }
    
    uint qmesh_size=_eigenvectors.size();

    stringstream os_disp;
    os_disp<<"# Cartesian atomic meansquare displacements in Angstrom uint\n";

    std::string STAR = std::string(5*_nBranches, '*');

    os_disp<<"[AFLOW] "<<STAR<<"\n";
    os_disp<< "[DISPLACEMENTS]START" <<"\n";
    os_disp<<"#"<<setw(9)<<"T(K)"<<setw(15)<<"AtomType"<<setw(15)<<"xdir"<<setw(15)<<"ydir"<<setw(15)<<"zdir"<<"\n";
    os_disp << std::fixed << std::showpoint;

    for(double t=Ts; t<=Te; t=t+Tinc)
      {
        os_disp<<setw(10)<<setprecision(2)<<t<<"\n";
        _VEC_ disps(_nBranches,0.0);

        for(uint kvec=0; kvec!=qmesh_size; kvec++)
          {
            _CMAT_ eigenvectors=xmat2mat(_eigenvectors[kvec]);
            _CMAT_ vecs2(_nBranches,
                         _CVEC_(_nBranches, _CD_(0.,0.)) );

            for(uint i=0; i!=_nBranches; i++){
              for(uint j=0; j!=_nBranches; j++){
                vecs2[i][j]= std::abs(eigenvectors[i][j])*
                  std::abs(eigenvectors[i][j]);
              }}

            _MAT_ vecs2_T(_nBranches,
                          _VEC_ (_nBranches, 0.0) );

            for(uint i=0; i!=_nBranches; i++){
              for(uint j=0; j!=_nBranches; j++){
                vecs2_T[i][j]=vecs2[j][i].real();
              }}

            _VEC_ f=xvec2vec(_freq_Thz[kvec]);

            for(uint i=0; i!=_nBranches; i++){
              uint cnt=0;
              //run through all atoms and their directions
              for(uint atom_alpha=0; atom_alpha!=_nBranches; atom_alpha++){
                if(f[i]<_CUTOFF_FREQ)continue;
                if((atom_alpha%3==0) && (atom_alpha!=0))cnt++;
                double c=vecs2_T[i][atom_alpha]/(_atomic_masses_amu[cnt]*AMU2Kg);
                disps[atom_alpha] += get_Q2(f[i], t) * c ;
              }
            }
          }//kvec loop
        os_disp<<setprecision(8);
        size_t atom_cnt=0;
        for(uint atom_alpha=0; atom_alpha!=_nBranches; atom_alpha++)
          {
            if((atom_alpha%3==0)&& atom_alpha!=0)
              {
                os_disp<<"\n";
              }
            if(atom_alpha%3==0)
	      {
		os_disp<<setw(25)<<_atomic_species[atom_cnt];
		atom_cnt++;
	      }
            os_disp<<setw(15)<<disps[atom_alpha]/((double)qmesh_size);
          }
        os_disp<<"\n";
      }//t loop

    os_disp<<"[AFLOW] "<<STAR<<"\n";
    os_disp << "[DISPLACEMENTS]STOP" <<"\n";

    string outfile =  "aflow.apl.displacements.out";
    if(!aurostd::stringstream2file(os_disp, outfile, "WRITE")) {
      throw APLRuntimeError("Cannot write aflow.apl.displacements.out");
    }
    aurostd::StringstreamClean(os_disp);
  }
  // ***************************************************************************************
  // ***************************************************************************************
  //projected displacements along [hkl] direction for a given range of temperature
  //formulas can be found in the link http://atztogo.github.io/phonopy/thermal-displacement.html
  void  AtomicDisplacements::projected_displacement(const _VEC_ &direction, double Ts, double Te, int Tinc)
  {
    _logger<<"calculating projected mean square displacements along ["<<NumToStr<int>((int)direction[0])
           <<" "<<NumToStr<int>((int)direction[1])<<" "<<NumToStr<int>((int)direction[2])<<"]" << apl::endl;

    uint qmesh_size=_eigenvectors.size();

    if(qmesh_size==0){_logger << apl::error <<" eigenvectors are zero " <<apl::endl; return;}


    if(_atomic_species.size()!=_atomic_masses_amu.size())
      {
        _logger << apl::error <<"error in atomic masses [remove apl.xml LOCK and run again]" <<apl::endl;
        return;
      }

    vector<double> projector(3,0.);
    _MAT_ lattice=xmatd2matd(_rlattice);

    //get projected vector
    for(uint i=0; i!=lattice.size(); i++){
      vector<double> a=lattice[i];
      projector[i]=apl_inner_product(a, direction);
    }

    double norm=0.0;
    for(uint i=0; i!=projector.size(); i++)
      norm+=projector[i]*projector[i];
    norm=sqrt(norm);
    for(uint i=0; i!=projector.size(); i++)
      projector[i]/=norm;

    //make transpose of eigenvector
    vector<_CMAT_> p_eigenvectors;
    for(uint i=0; i!=qmesh_size; i++)//sum over qpoints
      {
        //transpose of eigenvector
        _CMAT_ EVE=xmat2mat(_eigenvectors[i]);
        _CMAT_ EVE_T(_nBranches, _CVEC_(_nBranches, _CD_(0.0,0.0)));

        for(uint j=0; j!=_nBranches; j++)
          for(uint k=0; k!=_nBranches; k++)
            EVE_T[j][k]=EVE[k][j];

        _CMAT_ tmp2d;
        for(uint j=0; j!=_nBranches; j++){
          //seperate real and imeginary parts
          _VEC_ vecsR(_nBranches,0.00);
          _VEC_ vecsI(_nBranches,0.00);
          for(uint k=0; k<_nBranches; k++){
            vecsR[k]=EVE_T[j][k].real();
            vecsI[k]=EVE_T[j][k].imag();
          }
          _CVEC_ tmp1d;
          for(uint k=0; k<_nBranches; k++){
            if(k%3==0){
              vector<double> tmp1(vecsR.begin()+k, vecsR.begin()+k+3);
              vector<double> tmp2(vecsI.begin()+k, vecsI.begin()+k+3);
	      double real=apl_inner_product(tmp1, projector);
	      double imag=apl_inner_product(tmp2, projector);
              tmp1d.push_back(_CD_(real,imag));
            }
          }
          tmp2d.push_back(tmp1d);
        }
        p_eigenvectors.push_back(tmp2d);
      }

    uint nATOMs=_atomic_species.size();
    _CVEC_ tmp1d(_nBranches,_CD_(0.,0.));
    _CMAT_ tmp2d(nATOMs,tmp1d);
    vector<_CMAT_> p_eigenvectorsT(qmesh_size, tmp2d);

    for(uint i=0; i!=qmesh_size; i++){
      for(uint j=0; j!=_nBranches; j++){
        for(uint k=0; k!=nATOMs; k++){
          p_eigenvectorsT[i][k][j]=p_eigenvectors[i][j][k];
        }}}
    stringstream os_disp;
    os_disp<<"# Cartesian projected atomic meansquare displacements in Angstrom uint \n";
    std::string STAR = std::string(10*_nBranches, '*');
    os_disp<<"[AFLOW] "<<STAR<<"\n";
    os_disp<< "[PROJECTED_DISPLACEMENTS]START" <<"\n";
    os_disp<<setw(10)<<"# Projection direction set to ["<<(int)direction[0]<<" "<<(int)direction[1]<<" "<<(int)direction[2]<<"]\n";

    os_disp<<"#"<<setw(14)<<"T(K)";
    for(uint i=0; i!=nATOMs; i++)
      os_disp<<setw(15)<<_atomic_species[i];
    os_disp<<"\n";

    //temperature dependent displacements
    for(double t=Ts; t<=Te; t=t+Tinc){
      _VEC_ disps(nATOMs,0.00);
      for(uint kvec=0; kvec!=qmesh_size; kvec++)//qpoint sum
        {
          _CMAT_ eigenvectors=p_eigenvectorsT[kvec];
          _CMAT_ vecs2(nATOMs, _CVEC_(_nBranches, _CD_(0.,0.)));

          for(uint i=0; i!=nATOMs; i++){
            for(uint j=0; j!=_nBranches; j++){
              vecs2[i][j]= std::abs(eigenvectors[i][j])*
                std::abs(eigenvectors[i][j]);
            }}

          //transpose of vecs2 and its complex parts are zero
          _MAT_ vecs2_T(_nBranches, _VEC_(nATOMs, 0.0) );

          for(uint i=0; i!=nATOMs; i++){
            for(uint j=0; j!=_nBranches; j++){
              vecs2_T[j][i]=vecs2[i][j].real();
            }}

          _VEC_ f= xvec2vec(_freq_Thz[kvec]);

          for(uint i=0; i!=_nBranches; i++){
            for(uint j=0; j<nATOMs; j++){
              if(f[i]<_CUTOFF_FREQ)continue;
              double c = vecs2_T[i][j]/(_atomic_masses_amu[j]*AMU2Kg);
              disps[j]+=get_Q2(f[i], t)*c;
            }
          }
        }//qpoint loop
      //PRINT
      for(uint i=0; i!=nATOMs; i++)
        disps[i]/=(double)qmesh_size;

      os_disp << setw(15)<<std::setprecision(2)<<t<<setw(15);

      for(uint i=0; i!=nATOMs; i++){
        os_disp << std::fixed << std::showpoint;
        os_disp <<std::setprecision(8)<<disps[i]<<setw(15);
      }
      os_disp<<"\n";
    }
    os_disp<<"[AFLOW] "<<STAR<<"\n";
    os_disp << "[PROJECTED_DISPLACEMENTS]STOP" <<"\n";

    string outfile =  "aflow.apl.projected_displacements.out";
    if(!aurostd::stringstream2file(os_disp, outfile, "WRITE")) {
      throw APLRuntimeError("Cannot write aflow.apl.projected_displacements.out");
    }
    aurostd::StringstreamClean(os_disp);
  }
  // ***************************************************************************************
  //Zone-center phonon modes with directions indicated by arrows. This file can be visualized by XcrySDen.  
  //bool  AtomicDisplacements::write_normal_mode_direction(string user_kpoints)
  bool  AtomicDisplacements::write_normal_mode_direction(const vector<xvector<double> >& hs_kpoints)
  {
    if(hs_kpoints.size()==0)
      {
	_logger<<  apl::warning <<"AtomicDisplacements:: hs_kpoints.size()=0, "<<apl::endl;
	return false;
      }
    //uint nk = user_kpoints_double.size();
    uint nk = hs_kpoints.size();
    vector<xmatrix<xcomplex<double> > > eigenvectors; eigenvectors.clear();
    xmatrix<xcomplex<double> > tmp_c_m(_nBranches, _nBranches, 1,1);
    eigenvectors.resize(nk, tmp_c_m);
    xvector<double> tmp_d_v(_nBranches, 1);
    vector<xvector<double> > eigenvalues; eigenvalues.clear();
    eigenvalues.resize(nk, tmp_d_v);
 
    for(uint i=0; i!=nk; i++)
      {
	xvector<double> kpoint(3,1);
	//kpoint[1]=user_kpoints_double[i][0];
	//kpoint[2]=user_kpoints_double[i][1];
	//kpoint[3]=user_kpoints_double[i][2];
	kpoint[1]=hs_kpoints[i][1];
	kpoint[2]=hs_kpoints[i][2];
	kpoint[3]=hs_kpoints[i][3];
	//kpoint = trasp(_klattice) * kpoint;

        xmatrix<xcomplex<double> >  DM(_nBranches,_nBranches,1,1);
        DM=_pc.getDynamicalMatrix(kpoint);
        //calculate eigenvalue and eigenvectors of a dynamical matirx
	apl::aplEigensystems e;
	e.eigen_calculation(DM, eigenvalues[i], eigenvectors[i], APL_MV_EIGEN_SORT_VAL_ASC);
      }
    

    std::ofstream ofs_anime("aflow.apl.normal_mode_direction.axsf");
    if (!ofs_anime.is_open())
      {
        _logger<<apl::error<<"aflow.apl.normal_mode_direction.axsf unable to open"<<apl::endl;
        exit(0);
      }
    ofs_anime.setf(std::ios::scientific);

    uint natmin = _atomic_masses_amu.size();
    uint nbands=3 * natmin;

    double force_factor = 100.0;

    vector<xvector<double> > xmod(natmin, xvector<double>(3,1));
    vector<string> kd_tmp(natmin, "");

    ofs_anime << "ANIMSTEPS " << nbands * nk << std::endl;
    ofs_anime << "CRYSTAL" << std::endl;
    ofs_anime << "PRIMVEC" << std::endl;


    for (uint i = 0; i < 3; ++i) {
      for (uint j = 0; j < 3; ++j) {
	ofs_anime << std::setw(15) << _rlattice[i+1][j+1];
      }
      ofs_anime << std::endl;
    }

    for (uint i = 0; i < natmin; ++i) {
      for (uint j = 1; j <= 3; ++j) {
	xmod[i][j] = _atomic_c_positions[i][j];
      }

      kd_tmp[i] = _atomic_species[i];
    }
    double norm=0.0;
    uint i = 0;

    double amu_ry=911.444242;
   
    for (uint ik = 0; ik < nk; ++ik) 
      {
        for (uint imode = 0; imode < nbands; ++imode) 
	  {
            ofs_anime << "PRIMCOORD " << std::setw(10) << i + 1 << std::endl;
            ofs_anime << std::setw(10) << natmin << std::setw(10) << 1 << std::endl;
            norm = 0.0;

            for (uint j = 0; j < 3 * natmin; ++j) {
	      //std::complex<double>  evec_tmp = std::complex<double> (eigenvectors[ik][imode+1][j+1].re, eigenvectors[ik][imode+1][j+1].im);
	      std::complex<double>  evec_tmp = std::complex<double> (eigenvectors[ik][j+1][imode+1].re, eigenvectors[ik][j+1][imode+1].im);
	      norm += std::pow(evec_tmp.real(), 2) + std::pow(evec_tmp.imag(), 2);
            }

            norm *= force_factor / (double)(natmin);

            for (uint j = 0; j < natmin; ++j) {

	      ofs_anime << std::setw(10) << kd_tmp[j];

	      for (uint k = 0; k < 3; ++k) {
		ofs_anime << std::setw(15) << xmod[j][k+1];
	      }
	      for (uint k = 0; k < 3; ++k) {
		ofs_anime << std::setw(15)
			  << eigenvectors[ik][3 * j + (k+1)][imode+1].re
		  / (std::sqrt(_atomic_masses_amu[j]*amu_ry) * norm);
		//<< eigenvectors[ik][imode+1][3 * (j+1) + (k+1)].re
	      }
	      ofs_anime << std::endl;
            }

            ++i;
	  }
      }

    ofs_anime.close();
    return true;
  }
  // ***************************************************************************************
  void AtomicDisplacements::calc_participation_ratio_all()
  {
    ofstream out ("aflow.apl.apr.out");
    if (!out.is_open())
      {
	_logger<<apl::error <<"aflow.apl.apr.out not able to open "<<apl::endl;
	return;
      }
    std::string STAR = std::string(10*_nBranches, '*');
    out<<"[AFLOW] "<<STAR<<"\n";
    out<< "[ATOMIC PARTICIPATION RATIO]START" <<"\n";

    out.setf(std::ios::scientific);

    out << "# Atomic participation ratio of each phonon modes at k points\n";
    out << "# kpoint, mode, atom, frequency[kpoint][mode] (THz), APR[kpoint][mode][atom]\n";
    for(uint ik=0; ik!=_eigenvectors.size(); ik++)
      {
	out << "#" << std::setw(8) << ik + 1;
	out << " k = ";
	out << std::setw(15) << _kpoints[ik][1] << setw(15) << _kpoints[ik][2]<< setw(15) << _kpoints[ik][3]<<"\n";
	_CMAT_ data2D = xmat2mat(_eigenvectors[ik]);
	for(uint nBr=0; nBr!=data2D.size(); nBr++)
	  {
	    vector<double> ret(_atomic_masses_amu.size(), 0.0);
	    for(uint atom=0; atom!=_atomic_masses_amu.size(); atom++)
	      {
		ret[atom] = (std::norm(data2D[3 * atom][nBr])
			     + std::norm(data2D[3 * atom + 1][nBr])
			     + std::norm(data2D[3 * atom + 2][nBr])) / _atomic_masses_amu[atom];
	      }
	    double sum = 0.0;
	    for (uint iat = 0; iat < _atomic_masses_amu.size(); ++iat) sum += ret[iat] * ret[iat];

	    for (uint iat = 0; iat < _atomic_masses_amu.size(); ++iat)
	      ret[iat] /= std::sqrt( (double)(_atomic_masses_amu.size()) * sum);

	    sum = 0.0;
            for (uint iat = 0; iat < _atomic_masses_amu.size(); ++iat) {
	      sum += ret[iat];
            }
	    for(uint i=0; i!=ret.size(); i++)
	      {
		out << std::setw(8) << ik  + 1;
		out << std::setw(5) << nBr + 1;
		out << std::setw(5) << i + 1;
		out<<setw(15)<<_freq_Thz[ik][nBr+1]<<setw(15)<<ret[i];
		out<<"\n";
	      }
	    out<<"#Participation ratio \n";
	    out<<setw(15)<<sum<<"\n";
	  }
	out<<"\n";
      }
    out<< "[ATOMIC PARTICIPATION RATIO]END" <<"\n";
    out<<"[AFLOW] "<<STAR<<"\n";
    out.close();
  }
  // ***************************************************************************************
  bool AtomicDisplacements::calc_participation_ratio_along_path(const vector< xvector<double> > &qpoints)
  {
    if(qpoints.size()==0)
      {
	_logger<<  apl::warning <<"Can't calculate participation ratio along path, qpoints.size()==0 "<<apl::endl;
	return false;
      }
    _kpoints_path.clear();
    _kpoints_path=qpoints;
    xmatrix< xcomplex<double> > dd(_nBranches,_nBranches,1,1);
    _eigenvectors_path.resize(_kpoints_path.size(), dd);

    return eigen_solver_path();
  }
  // ***************************************************************************************
  void AtomicDisplacements::write_participation_ratio_along_path(const vector <double> &path, const vector<int> &path_segment)
  {
    if(path.size()==0)
      {
	_logger<<  apl::warning <<"Can't calculate participation ratio along path, path.size()==0 "<<apl::endl;
	return;
      }
    if(path_segment.size()==0)
      {
	_logger<<  apl::warning <<"Can't calculate participation ratio along path, path_segment.size()==0 "<<apl::endl;
	return;
      }
    ofstream out ("aflow.apl.apr.path.out");
    if (!out.is_open())
      {
        _logger<<apl::error <<"aflow.apl.apr.path.out not able to open "<<apl::endl;
        return;
      }
    //out.setf(std::ios::scientific);

    for(uint ik=0; ik!=_eigenvectors_path.size(); ik++)
      {
	out<< setw(5) << path_segment[ik];
	out<< setprecision(6)<<std::fixed << std::showpoint
           << setw(15) << path[ik];

	out<<setprecision(6)<<std::fixed << std::showpoint;

        _CMAT_ data2D = xmat2mat(_eigenvectors[ik]);
        for(uint nBr=0; nBr!=data2D.size(); nBr++)
          {
            vector<double> ret(_atomic_masses_amu.size(), 0.0);
            for(uint atom=0; atom!=_atomic_masses_amu.size(); atom++)
              {
                ret[atom] = (std::norm(data2D[3 * atom][nBr])
                             + std::norm(data2D[3 * atom + 1][nBr])
                             + std::norm(data2D[3 * atom + 2][nBr])) / _atomic_masses_amu[atom];
              }
            double sum = 0.0;
            for (uint iat = 0; iat < _atomic_masses_amu.size(); ++iat) sum += ret[iat] * ret[iat];

            for (uint iat = 0; iat < _atomic_masses_amu.size(); ++iat)
              ret[iat] /= std::sqrt( (double)(_atomic_masses_amu.size()) * sum);

            sum = 0.0;
            for (uint iat = 0; iat < _atomic_masses_amu.size(); ++iat) {
              sum += ret[iat];
            }
            for(uint i=0; i!=ret.size(); i++)
              {
                out<<setw(15)<<ret[i];
              }
          }
	out<<"\n";
      }
    out.close();
  }
  // ***************************************************************************************
  void AtomicDisplacements::solve_eigenvalues_in_threads_path(int startIndex, int endIndex, int cpuid)
  {
    for(int In=startIndex;In<endIndex;In++)
      {
        xvector<double> qpoint;
        xmatrix<xcomplex<double> >  DM(_nBranches,_nBranches,1,1);

        qpoint=_kpoints_path[In];
        DM=_pc.getDynamicalMatrix(qpoint);

        //calculate eigenvalue and eigenvectors of a dynamical matirx
        xvector<double> eigenvalues(_nBranches, 1);
        xmatrix<xcomplex<double> > eigenvectors(_nBranches, _nBranches, 1,1);
        apl::aplEigensystems e;
        e.eigen_calculation(DM, eigenvalues, eigenvectors, APL_MV_EIGEN_SORT_VAL_ASC);

        for(uint j=1; j<=_nBranches; j++)
          {
            if(_iszero(eigenvalues[j]))eigenvalues[j]=0.0;
            if(eigenvalues[j]<0){
              if(_isnegative(eigenvalues[j]))
                {
                  if(!_is_freq_negative)
                    {
                      _logger<<  apl::warning <<"AtomicDisplacements:: negative eigenvalue: = " <<eigenvalues[j]<<apl::endl;
                    }
		  _is_freq_negative=true;
		  return;
                }else eigenvalues[j]=0.00;
            }
          }
        _eigenvectors_path[In]=eigenvectors;
      }
    _freq_test[cpuid]=true;
  }//fn end
  // ***************************************************************************************
  bool AtomicDisplacements::eigen_solver_path()
  {
    _freq_test.clear();
    bool pass=false;
#ifdef AFLOW_APL_MULTITHREADS_ENABLE
    // Get the number of CPUS
    int ncpus = sysconf(_SC_NPROCESSORS_ONLN);// AFLOW_MachineNCPUs;
    if(ncpus<1) ncpus=1;
//    int qpointsPerCPU = _kpoints_path.size() / ncpus;
    _freq_test.resize(ncpus, false);
    // Show info 
    if( ncpus == 1 )
      {
        _logger.initProgressBar("Calculating eigenvalues for AtomicDisplacements");
      }
    else
      {
        _logger.initProgressBar("Calculating eigenvalues for AtomicDisplacements (" + stringify(ncpus) + " threads)");
      }

    // Distribute the calculation
    int startIndex, endIndex;
    std::vector< std::thread* > threads;
    vector<vector<int> > thread_dist = getThreadDistribution((int) _kpoints_path.size(), ncpus);
    for (int icpu = 0; icpu < ncpus; icpu++) {
      startIndex = thread_dist[icpu][0];
      endIndex = thread_dist[icpu][1];
      threads.push_back( new std::thread(&AtomicDisplacements::solve_eigenvalues_in_threads_path,this,startIndex,endIndex, icpu) );
    }

/* OBSOLETE ME180801
    for(int icpu = 0; icpu < ncpus; icpu++) {
      startIndex = icpu * qpointsPerCPU;
      endIndex = startIndex + qpointsPerCPU;
      if( ( (uint)endIndex > _kpoints_path.size() ) ||
          ( ( icpu == ncpus-1 ) && ( (uint)endIndex < _kpoints_path.size() ) ) )
        endIndex = _kpoints_path.size();
      threads.push_back( new std::thread(&AtomicDisplacements::solve_eigenvalues_in_threads_path,this,startIndex,endIndex, icpu) );
    }
*/

    // Wait to finish all threads here!
    for(uint i = 0; i < threads.size(); i++) {
      threads[i]->join();
      delete threads[i];
    }

    // Done
    _logger.finishProgressBar();


    uint bool_size=0;
    for(uint id=0; id!=_freq_test.size(); id++)
      {
        bool_size+=(uint)_freq_test[id];
      }
    if(bool_size==(uint)ncpus)pass=true;
    _freq_test.clear();
#else
    _freq_test.resize(1, false);
    solve_eigenvalues_in_threads_path(0, (int)_kpoints_path.size(), 0);
    uint bool_size=0;
    for(uint id=0; id!=_freq_test.size(); id++)
      {
        bool_size+=(uint)_freq_test[id];
      }
    if(bool_size==1)pass=true;
#endif

    return pass;
  }
  // ***************************************************************************************
  double AtomicDisplacements::get_Q2(double freq, double t)
  {
    return (10.545721821978764) * (  //hbar*ev/A^2//
                                   (get_population(freq, t) + 0.5) / (freq * 2. * M_PI));
  }
  // ***************************************************************************************
  double AtomicDisplacements::get_population(double freq, double t)
  {
    if (t < 1.)
      return 0.0;
    else
      {
        return 1.0 / (exp((freq * 0.00413566733) / (8.61733825681e-05 * t)) - 1.);
      }
  }
  // ***************************************************************************************
  _CMAT_ AtomicDisplacements::xmat2mat(const xmatrix<xcomplex<double> > &M)
  {
    _CMAT_ m(M.cols, _CVEC_(M.rows, _CD_(0.0,0.0)));
    for(int i=0; i<M.rows; i++)
      for(int j=0; j<M.cols; j++)
        m[i][j]=_CD_(M[i+1][j+1].re, M[i+1][j+1].im);

    return m;
  }
  // ***************************************************************************************
  _VEC_ AtomicDisplacements::xvec2vec(const xvector<double> &V)
  {
    _VEC_ v(V.rows, 0.0);
    for(int i=0; i<V.rows; i++)
      v[i]=V[i+1];

    return v;
  }
  // ***************************************************************************************
  template <typename T>
  string AtomicDisplacements:: NumToStr ( T Number ){
    ostringstream ss;
    ss << Number;
    return ss.str();
  }
  // ***************************************************************************************
  _MAT_ AtomicDisplacements::xmatd2matd(const xmatrix<double> &M)
  {
    _MAT_ m(M.cols, _VEC_(M.rows, 0.0));
    for(int i=0; i<M.rows; i++)
      for(int j=0; j<M.cols; j++)
        m[i][j]=M[i+1][j+1];

    return m;
  }
  // ***************************************************************************************
  double AtomicDisplacements::apl_inner_product(const vector<double> &a, const vector<double> &b)
  {
    if(a.size()!=b.size())
      _logger<<apl::error<<"apl_inner_product() vector size must be equal" << apl::endl;

    double sum=0.0;
    for(uint i=0; i!=a.size(); i++)sum+=a[i]*b[i];
    return sum;
  }
  // ***************************************************************************************
}//apl namespace end
