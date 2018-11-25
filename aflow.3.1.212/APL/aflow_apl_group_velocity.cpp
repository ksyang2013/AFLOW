// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2015           *
// *                Aflow PINKU NATH - Duke University 2014-2016             *
// *                                                                         *
// ***************************************************************************
// Written by Pinku Nath
// pn49@duke.edu

#include "aflow_apl.h"

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
  GroupVelocity::GroupVelocity(IPhononCalculator& pc,  UniformMesh& umesh, Logger& l):_pc(pc),_umesh(umesh), _logger(l)
  {
    _logger<<"Setting variables to calculate group velocities "<<apl::endl;
    clear();
    _nBranches=_pc.getNumberOfBranches();
    _kshift=1.0e-4;
    _sound_speed=0.0;
    populate_variables();
  }
  // ***************************************************************************************
  GroupVelocity::~GroupVelocity(){this->clear();}
  // ***************************************************************************************
  void GroupVelocity::clear()
  {
    _freq_kp.clear();
    _freq_km.clear();
    _freq.clear();
    _kpoints_kp.clear();
    _kpoints_km.clear();
    _kpoints.clear();
    _weights.clear();
    _gv.clear();
    _phvel.clear();
  }
  // ***************************************************************************************
  void GroupVelocity::populate_variables()
  {
    _logger<<"Populating group velocity auxiliary variables "<<apl::endl;
    _kpoints_kp.clear(); _kpoints_km.clear();_kpoints.clear(); 
    _gv.clear();
    _freq_kp.clear(); _freq_km.clear(); _freq.clear();


    _kpoints=_umesh.get_kpoints();
    _weights=_umesh.get_weights();
    xvector<double> dv(3,1);
    _kpoints_kp.resize(_kpoints.size(), dv);
    _kpoints_km.resize(_kpoints.size(), dv);

    xvector<double> dvp(_nBranches,1);
    _freq_kp.resize(_kpoints.size(), dvp);
    _freq_km.resize(_kpoints.size(), dvp);
    _freq.resize(_kpoints.size(), dvp);
    xmatrix<double> dm(_nBranches,3,1,1);
    _gv.resize(_kpoints.size(), dm);
    _phvel.resize(_kpoints.size(), dvp);

    //delete
    xmatrix<xcomplex<double> > tmp_xmatrix(_nBranches, _nBranches, 1, 1);
    _eigenvectors.clear();
    _eigenvectors.resize(_kpoints.size(), tmp_xmatrix);
    //delete

    if (_kpoints.size()==0)
      {
	_logger<<apl::error<<"group velocity calculation _kpoints.size()==0 "<<apl::endl; 
	exit(0);
      }

  }
  // ***************************************************************************************
  bool GroupVelocity::check_negative_frequencies()
  {
    if(!eigen_solver(1))return false;//at k-0

    if(!eigen_solver())return false;//three different, kx, ky and kz directions
    
    return true;
  }
  // ***************************************************************************************
  bool GroupVelocity::eigen_solver()
  {
  
    for(int dir=1; dir<=3; dir++)
      { 
	for(uint i=0; i!=_kpoints.size(); i++)
	  {
	    _kpoints_kp[i][dir]=0.0;
	    _kpoints_km[i][dir]=0.0;
	    _kpoints_kp[i][dir] = _kpoints[i][dir] + _kshift;
	    _kpoints_km[i][dir] = _kpoints[i][dir] - _kshift;
	  }
	if(!eigen_solver(dir+1))return false;
      }
    for(uint i=0; i!=_gv.size(); i++)
      {
	for(int j=1; j<=_gv[i].rows; j++)
	  {
	    _phvel[i][j] = std::sqrt(std::pow(_gv[i][j][1], 2)+ std::pow(_gv[i][j][2], 2) + std::pow(_gv[i][j][3], 2));
          }
      }
    //calculate sound speed
    sound_speed();
    return true;
  }
  // ***************************************************************************************
  bool GroupVelocity::compute_group_velocities()
  {
    if(!eigen_solver())return false;
    return true;
  }
  // ***************************************************************************************
  void GroupVelocity::sound_speed()
  {
    //finding wave vector near gamma point 
    vector<double> kmod; kmod.clear();
    kmod.resize(_kpoints.size(), 100.0);
    for(uint i=0; i!=_kpoints.size(); i++)
      {
	double k_norm_tmp = std::sqrt(std::pow(_kpoints[i][1], 2)+ std::pow(_kpoints[i][2], 2) + std::pow(_kpoints[i][3], 2));
	if(!_iszero(k_norm_tmp))
	  {
	    kmod[i]=k_norm_tmp;
	  }
      }
    //vector<double>::iterator result = std::min_element(std::begin(kmod), std::end(kmod));
    int min_index= indexofSmallestElement(kmod);//std::distance(std::begin(kmod), result);
    double k_norm = kmod[min_index]; 
    kmod.clear();


    stringstream dm;
    dm << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
    dm << setprecision(8);

    //ofstream dm("aflow.apl.gamma_point.out");
    //dm << std::setprecision(6);
    //dm << std::fixed;

    std::string STAR60 = std::string(60, '*');
    dm <<"[AFLOW] "<<STAR60<<"\n";
    dm <<"[AFLOW_GAMMA_POINT_PROPERTIES]START"<<"\n";
    dm <<"[AFLOW] "<<STAR60<<"\n";
    _sound_speed=0.0;
    dm<<"#"<<setw(15)<<"|kpoint|"<<setw(15)<<"kx"<<setw(15)<<"ky"<<setw(15)<<"kz"<<"\n";
    dm<<setw(15)<<k_norm<<setw(15)<<_kpoints[min_index][2]<<setw(15)<<_kpoints[min_index][2]<<setw(15)<<_kpoints[min_index][3]<<"\n";
    dm<<"#"<<setw(15)<<"Acoustic"<<setw(15)<<"freq(THz)"<<setw(15)<<"freq(THz)"<<setw(15)<<"freq(THz)"<<"\n";
    dm<<setw(30)<<_freq[min_index][1]<<setw(15)<<_freq[min_index][2]<<setw(15)<<_freq[min_index][3]<<"\n";


    for(int k=1; k<=3; k++)
      {
	double tmp=(2.0*M_PI*_freq[min_index][k])/k_norm;
	tmp=1.0/tmp;
	_sound_speed+=std::pow(tmp, 3);
      }
    _sound_speed=(_sound_speed)/3.0;
    _sound_speed=std::pow(_sound_speed, -0.33333);
    dm<<"sound_speed "<<setw(15)<<_sound_speed*100.0<<setw(15)<<"m/sec"<<"\n";
    dm <<"[AFLOW] "<<STAR60<<"\n";
    dm <<"[AFLOW_GAMMA_POINT_PROPERTIES]END"<<"\n";

    string filename = "aflow.apl.gamma_point.out";
    aurostd::stringstream2file(dm, filename);
    if (!aurostd::FileExist(filename)) {
      throw apl::APLRuntimeError("Cannot open output aflow.apl.gamma_point.out.");
    }
    //dm.close();
  }
  // ***************************************************************************************
  void GroupVelocity::write()
  {
    stringstream dm;
    dm << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
    dm << setprecision(8);

    std::string STAR100 = std::string(100, '*');
    dm <<"[AFLOW] "<<STAR100<<"\n";
    dm <<"#freq = frequency THz uint "<<"\n";
    dm <<"#|v| = average group velocity sqrt(v_x^2+v_y^2+v_z^2) (THz Angstrom) uint "<<"\n";
    dm <<"#v_i = group velocities in the i-th direction (THz Angstrom) uint "<<"\n";
    dm <<"#s = speed of sound (THz Angstrom) uint "<<"\n";
    dm <<"[AFLOW] "<<STAR100<<"\n";
    dm <<"[AFLOW_GROUP_VELOCITY]START"<<"\n";
    dm <<"[AFLOW] "<<STAR100<<"\n";
    dm<<"#"<<setw(12)<<"freq"<<setw(15)<<"|v|"<<setw(15)<<"v_x"<<setw(15)<<"v_y"<<setw(15)<<"v_z"<<setw(15)<<"s"<<"\n";
    dm << std::setprecision(6);
    dm << std::fixed;
    for(uint i=0; i!=_gv.size(); i++)
      {
	dm<<"#k-vector"<<setw(15)<<_kpoints[i][1]<<setw(15)<<_kpoints[i][2]<<setw(15)<<_kpoints[i][3]<<"\n";
	for(int j=1; j<=_gv[i].rows; j++)
	  {
	    dm<<setw(15)<<_freq[i][j]<<setw(15)<<_phvel[i][j]<<setw(15)<<_gv[i][j][1]
	      <<setw(15)<<_gv[i][j][2]<<setw(15)<<_gv[i][j][3]<<setw(15)<<_sound_speed<<"\n";
	  }
      }
    dm <<"[AFLOW] "<<STAR100<<"\n";
    dm <<"[AFLOW_GROUP_VELOCITY]END"<<"\n";

    string filename = "aflow.apl.group_velocities.out";
    aurostd::stringstream2file(dm, filename);
    if (!aurostd::FileExist(filename)) {
      throw apl::APLRuntimeError("Cannot open output aflow.apl.group_velocities.out file.");
    }
    //clear unnecessary variables
    clear_auxiliary_variables();
  }
  // ***************************************************************************************
  void GroupVelocity::clear_auxiliary_variables()
  {
    _freq_kp.clear();
    _freq_km.clear();
    _freq_test.clear();
    _kpoints_kp.clear();
    _kpoints_km.clear();
  }
  // ***************************************************************************************
  bool GroupVelocity::eigen_solver(int ktype)
  {
    (void)ktype;
    _freq_test.clear();
    bool pass=false;
#ifdef AFLOW_APL_MULTITHREADS_ENABLE
    // Get the number of CPUS
    int ncpus = sysconf(_SC_NPROCESSORS_ONLN);// AFLOW_MachineNCPUs;
    if(ncpus<1) ncpus=1;
//    int qpointsPerCPU = _kpoints.size() / ncpus;  OBSOLETE ME180801
    _freq_test.resize(ncpus, false);
    // Show info 
    if( ncpus == 1 )
      {
        if(ktype==1)_logger.initProgressBar("Calculating eigenvalues for k-zero");
        else if(ktype==2)_logger.initProgressBar("Calculating eigenvalues for kx-dir");
        else if(ktype==3)_logger.initProgressBar("Calculating eigenvalues for ky-dir");
        else if(ktype==4)_logger.initProgressBar("Calculating eigenvalues for kz-dir");
      }
    else
      {
        if(ktype==1)_logger.initProgressBar("Calculating eigenvalues for k-zero (" + stringify(ncpus) + " threads)");
        else if(ktype==2)_logger.initProgressBar("Calculating eigenvalues for kx-dir (" + stringify(ncpus) + " threads)");
        else if(ktype==3)_logger.initProgressBar("Calculating eigenvalues for ky-dir (" + stringify(ncpus) + " threads)");
        else if(ktype==4)_logger.initProgressBar("Calculating eigenvalues for kz-dir (" + stringify(ncpus) + " threads)");
      }

    // Distribute the calculation
    int startIndex, endIndex;
    std::vector< std::thread* > threads;
    vector<vector<int> > thread_dist = getThreadDistribution((int) _kpoints.size(), ncpus);
    for (int icpu = 0; icpu < ncpus; icpu++) {
      startIndex = thread_dist[icpu][0];
      endIndex = thread_dist[icpu][1];
      threads.push_back( new std::thread(&GroupVelocity::solve_eigenvalues_at_k,this,startIndex,endIndex, icpu, ktype) );
    }

/* OBSOLETE ME180801
    for(int icpu = 0; icpu < ncpus; icpu++) {
      startIndex = icpu * qpointsPerCPU;
      endIndex = startIndex + qpointsPerCPU;
      if( ( (uint)endIndex > _kpoints.size() ) ||
          ( ( icpu == ncpus-1 ) && ( (uint)endIndex < _kpoints.size() ) ) )
        endIndex = _kpoints.size();
      threads.push_back( new std::thread(&GroupVelocity::solve_eigenvalues_at_k,this,startIndex,endIndex, icpu, ktype) );
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
    //calculate_gp_in_path(0, (int)_kpoints.size(), 0, ktype);
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
  void GroupVelocity::solve_eigenvalues_at_k(int startIndex, int endIndex, int cpuid, int ktype)
  {
    for(int In=startIndex;In<endIndex;In++)
      {
        xvector<double> qpoint;
        xmatrix<xcomplex<double> >  DM(_nBranches,_nBranches,1,1);
        xvector<double> qp;
        xmatrix<xcomplex<double> >  DMp(_nBranches,_nBranches,1,1);
        xvector<double> qm;
        xmatrix<xcomplex<double> >  DMm(_nBranches,_nBranches,1,1);

        if(ktype==1)
	  {
	    qpoint=_kpoints[In];
	    DM=_pc.getDynamicalMatrix(qpoint);
	  }else
	  {
	    qp=_kpoints_kp[In];
	    DMp=_pc.getDynamicalMatrix(qp);
	    qm=_kpoints_km[In];
	    DMm=_pc.getDynamicalMatrix(qm);
	  }


	if(ktype==1)
	  {
	    //calculate eigenvalue and eigenvectors of a dynamical matirx
            xvector<double> eigenvalues(_nBranches, 1);
            xmatrix<xcomplex<double> > eigenvectors(_nBranches, _nBranches, 1,1);
	    apl::aplEigensystems e;
	    e.eigen_calculation(DM, eigenvalues, eigenvectors, APL_MV_EIGEN_SORT_VAL_ASC);
	    _eigenvectors[In]=eigenvectors;

	    for(uint j=1; j<=_nBranches; j++)
	      {
		if(_iszero(eigenvalues[j]))eigenvalues[j]=0.0;
		if(eigenvalues[j]<0){
		  if(_isnegative(eigenvalues[j]))
		    {
		      _logger<<  apl::warning <<"GroupVelocity:: negative eigenvalue: = " <<eigenvalues[j]<<apl::endl;
		      return;
		    }else eigenvalues[j]=0.00;
		}
		_freq[In][j]=sqrt(eigenvalues[j])*RAW2Hz;
	      }
	  }else
	  {
            xvector<double> eigenvalues_p(_nBranches, 1);
            xmatrix<xcomplex<double> > eigenvectors_p(_nBranches, _nBranches, 1,1);
            xvector<double> eigenvalues_m(_nBranches, 1);
            xmatrix<xcomplex<double> > eigenvectors_m(_nBranches, _nBranches, 1,1);
	    apl::aplEigensystems e_p;
	    e_p.eigen_calculation(DMp, eigenvalues_p, eigenvectors_p, APL_MV_EIGEN_SORT_VAL_ASC);
	    apl::aplEigensystems e_m;
	    e_p.eigen_calculation(DMm, eigenvalues_m, eigenvectors_m, APL_MV_EIGEN_SORT_VAL_ASC);

	    for(uint j=1; j<=_nBranches; j++)
	      {
		if(_iszero(eigenvalues_p[j]))eigenvalues_p[j]=0.0;
		if(_iszero(eigenvalues_m[j]))eigenvalues_m[j]=0.0;

		if(eigenvalues_p[j]<0)
		  {
		    if(_isnegative(eigenvalues_p[j]))
		      {
			_logger<<  apl::warning <<"GroupVelocity:: negative eigenvalue: = " <<eigenvalues_p[j]<<apl::endl;
			return;
		      }else eigenvalues_p[j]=0.00;
		  }
		if(eigenvalues_m[j]<0)
		  {
		    if(_isnegative(eigenvalues_m[j]))
		      {
			_logger<<  apl::warning <<"GroupVelocity:: negative eigenvalue: = " <<eigenvalues_m[j]<<apl::endl;
			return;
		      }else eigenvalues_m[j]=0.00;
		  }


		_freq_kp[In][j]=sqrt(eigenvalues_p[j])*RAW2Hz;
		_freq_km[In][j]=sqrt(eigenvalues_m[j])*RAW2Hz;

		if(ktype==2)_gv[In][j][1]=(_freq_kp[In][j]-_freq_km[In][j])/(2.0*_kshift);
		else if(ktype==3)_gv[In][j][2]=(_freq_kp[In][j]-_freq_km[In][j])/(2.0*_kshift);
		else if(ktype==4)_gv[In][j][3]=(_freq_kp[In][j]-_freq_km[In][j])/(2.0*_kshift);
	      }
	  }
      }
    _freq_test[cpuid]=true;
  }//fn end
  // ***************************************************************************************
  int  GroupVelocity::indexofSmallestElement(const vector<double> &array)
  {
    int index = 0;
    int size=(int)array.size();

    for(int i = 1; i < size; i++)
    {
        if(array[i] < array[index])
            index = i;              
    }
    return index;
   }
  // ***************************************************************************************
}//apl namespace end
