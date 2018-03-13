****************************************************************************************************
*                                                                                                  *
*                        aflow - STEFANO CURTAROLO Duke University 2003-2018                       *
*                            High-Throughput ab-initio Computing Project                           *
*                                                                                                  *
****************************************************************************************************
LATEST VERSION OF THE FILE:                            materials.duke.edu/AFLOW/README_AFLOW_APL.TXT
****************************************************************************************************

Phonon calculations with deformations written by Michal Jahnatek, Pinku Nath, Jose J. Plata, Corey Oses and Stefano Curtarolo.

The current version of APL can calculate the phonon dispersion curves (output file PDIS),
phonon density of state (file PDOS), and thermal properties (file THERMO). All are
calculated in harmonic approximation. Two methods of calculation can be selected:
direct (supercell) method (DM) or linear response method (LR). The resulted file can
be transformed by next tool apl2agr into agr, ps, eps, or png form (located in
SOURCES/APL2AGR).

APL expects the well relaxed structure on input. There is no additional relaxation
inside the APL. So, you have to prepare a such structure by separate aflow run with
relaxation mode for forces. For phonon calculation, generally, you have to setup
3 things: 1) Tell to aflow you want a phonon calculation. 2) Which method of calculation
to use. 3) What to calculate.

1) To perform a phonon calculation this line has to be present in aflow.in:

   [AFLOW_APL]CALC

   The APL is checking its existence in aflow.in. If it does not exist, the job is not supposed
   to be for phonon calculation and no other phonon settings are read.

2) The user has possibility to switch between two methods of calculation by keyword ENGINE.
   For direct method (DM), also known as supercell method, one has to write:

   [AFLOW_APL]ENGINE=DM

   This is also a default setting. For linear response method it is:

   [AFLOW_APL]ENGINE=LR

   Since this parameter has its own default settings, you do not need to specify it if you are
   fine with DM. Each method has its own, but also common parameters, which can be used to tune
   the whole calculation.

   2.1) Common parameters used by both methods

        The most important setting is the supercell dimension (the supercell is required by both
        methods). Supercell can be set in three ways:
        [AFLOW_APL]SUPERCELL=numberXnumberXnumber
        [AFLOW_APL]MINSHELL=number
        [AFLOW_APL]MINATOMS=number

        EXAMPLES:
        Direct specification, like
        [AFLOW_APL]SUPERCELL=4x4x4

        Minimal nearest neighbor shell, specification like:
        [AFLOW_APL]MINSHELL=6
        The default settings is MINSHELL=6 and APL will construct the appropriate supercell. HOWEVER,
        THIS IS NOT THE RECOMMENDED WAY, BECAUSE INCLUSION OF 6 SHELL IS ENOUGH ONLY FOR VERY SIMPLE
        MATERIALS. There is no way or rule how to find the number of enough shells included into
        calculation. FOR BEST RESULTS YOU HAVE TO SPECIFY THE SUPERCELL DIMENSION DIRECTLY. We found
        the 4x4x4 is usually enough for all fcc systems; 3x3x2 for hcp; bcc are complicated, but
        generally it should be more that 3x3x3. Typically, the supercell with about 100 atoms should
        be fine for 90% systems. You might find the MINATOMS option useful for large unit cells.

        Minimal number of atoms, like:
        [AFLOW_APL]MINATOMS=number
        The code will homogeneously increase the volume of the supercells, until it finds a SUPERCELL(i,j,k) containing
        at least "number" atoms.

        [AFLOW_APL]KPPRA=number
        By default, AFLOW takes the number of KPOINTS (KPPRA) from the "[VASP_KPOINTS_FILE]KPPRA=number" entry.
        You can override the number with "[AFLOW_APL]KPPRA=number". Look at README_AFLOW.TXT for KPPRA.

        [AFLOW_APL]KSCHEME=scheme
        By default, AFLOW takes the KPOINTS scheme from the "[VASP_KPOINTS_FILE]KSCHEME=scheme" entry.
        You can override the scheme with "[AFLOW_APL]KSCHEME=scheme". Look at README_AFLOW.TXT for KSCHEME.

        [AFLOW_APL]KPOINTS=EVEN | ODD
        By default, AFLOW takes the KPOINTS parity from the "[VASP_FORCE_OPTION]KPOINTS=EVEN | ODD" entries, if specified.
        You can override the entries with "[AFLOW_APL]KPOINTS=EVEN | ODD". Look at README_AFLOW.TXT for KPOINTS.
        This is useful for a small number of KPOINTS.

  2.2) Direct method parameters

        The supercell (direct) method can be modified by this parameters:

        [AFLOW_APL]DXYZONLY=n
        [AFLOW_APL]DMAG=0.015
        [AFLOW_APL]DPM=y
        [AFLOW_APL]ZEROSTATE=n

        The DPM, DXYZONLY, and DMAG parameters control how the distortion vectors will be generated.
        Generally, APL will try an atom distortions along the [100],[010],[001],[110], [101],[011],
        [111] direction. It will apply all site symmetry operations, and that direction which produces
        the highest number of new independent vectors will by used in next calculation. For each unique
        atom of the primitive cell we need to know the forces acting along three independent directions.
        Hence, if one distortion direction can produce by the symmetry three linear independent vectors,
        this atom needs to be distorted in only one direction and the rest two are calculated by the
        symmetry. This is the case of fcc lattice.

        if you specify DXYZONLY=y, only the three trial directions ([100],[010],[001]) will be tried.
        This is important for cases when one primitive cell vector is much longer then others. All
        founded directions of distortion are of unit length. The real distortion is applied by some
        specified magnitude. This magnitude should be small enough. It can be set by keyword DMAG.
        The default value is 0.015 Angs., what is generally accepted value also by other software
        packages and works well for all systems.

	Usually, in effort to filter out the anharmonic effects from the forces,  it is important
        to calculate forces in both directions of distortion, i.e., DMAG is applied by positive and
        also negative sign. This is controlled by keyword DPM. If you specify DPM=y (default), twice
        more distortion vectors will be generated, and the final force acting on each atom is calculated
        as 0.5 * ( f(+) - f(-) ).

        Sometimes, it can happen the calculating structure is not well relaxed or can not be in a given
        precision. In this case, even the unstrained state has some non-zero forces. In this case is
        good to subtract these forces from all states. This can be done by ZEROSTATE=y, what means the
        APL will generate the calculation of this unstrained state and the obtained forces for this state
        will be subtracted from all other strained calculations. (!not well tested part)

   2.3) Linear response parameters

        The use of the linear response method as is implemented in VASP is more complicated. Every time
        you specify this method (ENGINE=LR), APL will check you binary file of VASP for version.
        If it is not 5.2 or higher it will refuse to continue. Moreover, vasp52 binary has to
        be run with special settings in you .bashrc, otherwise the calculation will crash after the
        first completely converged electronic loop. Hence, check you have these lines in you .bashrc:

        ulimit -s unlimited

        In case of LR method, APL is not responsible for finding the appropriate set of distortion vectors
        and all stuff for generation and creation of supercells. All is done by VASP. Hence, expect
        quite long run. Please note, the VASP is using the mode XYZONLY, hence finally more distortions
        are generated as are needed.  As result, the VASP will generate file DYNMAT, where all forces acting
        on atoms after each distortion are located. This is the main input file for APL for the next calculation.

   2.4) Linear response parameters

        If you know the studied system is polar, the phonon properties have to be calculated specially,
        where the dipole-dipole interaction is taking into account. This can be specified by keyword POLAR.

        [AFLOW_APL]POLAR=y

        In this case, APL will prepare also next run by VASP52, where Born effective charges and dielectric
        constant matrix will be calculated with Linear Response (DFTP).
        These both quantities will be used later in correction of dynamic matrices for correct TO-LO splittings.
        The default setting is POLAR=n.
	Note that TO-LO splitting is alcualted only with >VASP52.
        Epsilon and the Born charges are calculated in two ways: with LEPSILON or LCALCEPS. Default is the latter
        You can switch the code by changing the lines of aflow_apl.h
          #undef   AFLOW_APL_VASP_USE_LEPSILON
          #define  AFLOW_APL_VASP_USE_LCALCEPS


3) Generally, you can specify what to calculate by these keyword:

   [AFLOW_APL]DC=y
   [AFLOW_APL]DPM=y
   [AFLOW_APL]ZEROSTATE=y
   [AFLOW_APL]DOS=y
   [AFLOW_APL]TP=y

   The default setup is no for all. Each calculation mode has its own set of parameters.

   3.1) Phonon dispersion curves calculation (DC=y)
    
        The most important setting for this module is the path in the reciprocal space for which the
        frequencies should be calculated. The APL can do it automatically, via AFLOW-SYM,
        The procedure is: AFLOW-SYM is called to obtain the spacegroup number of primitive
        cell. Then, this number is used to get the path from the settings used for the calculation of
        the electronic structure in aflow. If you wish to skip this calculation, the space group has to be
        specified in aflow.in, like:

        [AFLOW_APL]DCINITPATH = RHL
        where you can specify the lattice variation as introduced in [doi=10.1016/j.commatsci.2010.05.010]

        or

        [AFLOW_APL]DCINITSG = 166

        If everything is OK, the APL will calculate the dispersion curves along the same path as is
        used for electronic structure. If you need to calculate it by an another path, you can
        specified it as:

        [AFLOW_APL]DCUSERPATH = G-X|X-U|K-G|G-L

        where G,X,U,L points are the same points as are defined inside the aflow (for electronic structure).
        The number of points in which each sub-path (G-X,X-U,...) will be calculated can be set as:

        [AFLOW_APL]DCPOINTS=100

        what is also a default setting. The calculated dispersion curves are saved in file PDIS, where
        the 1st column is ID of points in sub-path, the 2nd is the position in the whole path, ... the next
        column correspond to one phonon branch. There is 3 * natoms_in_primitive_cell branches. The
        frequencies are in units specified by (default)

        [AFLOW_APL]FREQFORMAT = "THZ | ALLOW_NEGATIVE"

        which means the frequency is returned in THz and also negative values are allowed to return. The OR
        operator is important for use. If ALLOW_NEGATIVE is not added, only the positive values are
        returned (negative values are returned as zero).

        Other options are

        [AFLOW_APL]FREQFORMAT = "RECIPROCAL_CM" or "CM-1"
        [AFLOW_APL]FREQFORMAT = "MEV"
        [AFLOW_APL]FREQFORMAT = "OMEGA"
        [AFLOW_APL]FREQFORMAT = "HERTZ"

   3.2) Phonon density of states calculation (DOS=y)

        There are implemented two method for integration in Brillouin zone - Linear Tetrahedron method (LT)
        and root sampling method (RS). As default is selected LT. User can switch between them by as:

        [AFLOW_APL]DOSMETHOD=LT

        or

        [AFLOW_APL]DOSMETHOD=RS

        The mesh of qpoints used for numerical integration should be set as:

        [AFLOW_APL]DOSMESH=21x21x21

        The setting 21x21x21 is also the default setting. The APL will automatically find the irreducible part
        of BZ (the list of irreducible qpoints with appropriate summation weights). Only for these qpoints a
        dynamic matrix is constructed and diagonalized. In case of linear tetrahedron method, the list of
        irreducible tetrahedrons with their weights is constructed and used in summation. The function of
        density of states will be evaluated in finite number of points, which can be set as:

        [AFLOW_APL]DOSPOINTS=2000

        The density of 2000 points is the default setting. The resulted DOS function can be smeared by gaussians
        of specific width. This is quite important for RS method to have a nice looking pictures. The default
        setting of smearing for RS method is 0.05, but can be changed as:

        [AFLOW_APL]DOSSMEAR=0.05

        For LT method the default setting is 0.0, so no smearing.

   3.3) Thermal properties calculation (TP=y)

        Generally, APL is calculating the zero point vibrational energy (U0), the internal energy (U), the
        vibrational free energy (Fv), the vibrational entropy (Sv), and isochoric specific heat (cv) for
        a given range of temperature. This range can be set as

        [AFLOW_APL]TPT=0:2000:10

        what means it will calculate all upper mentioned properties from 0 to 2000 K by step of 10 K. This
        is also a default setting.


4) The APL has integrated hibernate system, which is switched on by default. To turn it down:

   [AFLOW_APL]HIBERNATE=n

   This system means the calculated force constant matrices are stored in external file called
   apl.xml.bz2. Once this file exists and you will run aflow for phonon calculation again, it
   will automatically used this data (awake from it) and speed up calculation.

5) Inputs and inheritance.
   The minimal input is:

   [AFLOW_APL]CALC
   [AFLOW_APL]MINATOMS=100
   #[AFLOW_APL]SUPERCELL=4x4x4
   [AFLOW_APL]DC=y
   [AFLOW_APL]DPM=y
   [AFLOW_APL]ZEROSTATE=y
   [AFLOW_APL]DOS=y
   [AFLOW_APL]TP=y

   5.1) Inheritance
        The generated ARUN.APL***/aflow.in input files inherit some options from the master aflow.in
        If present, the following options are passed to the local aflow.in
        [VASP_FORCE_OPTION]AUTO_PSEUDOPOTENTIALS=...
        [VASP_FORCE_OPTION]SPIN=...
        [VASP_FORCE_OPTION]TYPE=...
        [VASP_FORCE_OPTION]AUTO_MAGMOM=...
        [VASP_FORCE_OPTION]BADER=...
        [VASP_FORCE_OPTION]ELF=...
        [VASP_FORCE_OPTION]ABMIX=...
        [VASP_FORCE_OPTION]KPOINTS=...
        [VASP_FORCE_OPTION]IGNORE_AFIX=...

        Other options are added by default, they are:
        [VASP_FORCE_OPTION]WAVECAR=OFF
        [VASP_FORCE_OPTION]CHGCAR=OFF
        [VASP_FORCE_OPTION]PREC=ACCURATE
        [VASP_FORCE_OPTION]ALGO=NORMAL

        To avoid intercepting NPAR warnings intruduced in vasp52+, the extra option is added
        [VASP_FORCE_OPTION]IGNORE_AFIX=NPARC

6) EXAMPLE:  PbS with Linear Response use the following:
   For Deformation engine, put [AFLOW_APL]ENGINE=DM and uncomment #[AFLOW_APL]DMAG=0.015

    [AFLOW] **********************************************************************************
    [AFLOW]SYSTEM=PbS
    [AFLOW] **********************************************************************************
    [AFLOW_MODE=VASP]
    [AFLOW_MODE_ZIP=bzip2]
    [AFLOW_MODE_BINARY=vasp53s]
    [AFLOW] **********************************************************************************
    #[AFLOW_MODE_MPI]
    [AFLOW_MODE_MPI_MODE]NCPUS=MAX
    [AFLOW_MODE_MPI_MODE]COMMAND ="mpirun -np"
    [AFLOW_MODE_MPI_MODE]AUTOTUNE
    [AFLOW_MODE_MPI_MODE]BINARY="mpivasp53s"
    [AFLOW] **********************************************************************************
    [AFLOW_APL]CALC
    [AFLOW_APL]ENGINE=LR
    #[AFLOW_APL]DMAG=0.015
    [AFLOW_APL]POLAR=Y
    [AFLOW_APL]SUPERCELL=4x4x4
    [AFLOW_APL]DC=Y       
    #[AFLOW_APL]DCINITSG=227
    #[AFLOW_APL]DCUSERPATH=Gamma-X|X-U|K-Gamma|Gamma-L
    [AFLOW_APL]DCINITPATH=FCC
    [AFLOW_APL]DOS=Y                               
    [AFLOW_APL]TP=Y                                
    [AFLOW_APL]TPT=0:2000:10                       
    [AFLOW_APL]KPPRA=1024                          
    [AFLOW_APL]KSCHEME=G                           
    [AFLOW_APL]KPOINTS=EVEN                        
    [AFLOW_APL]DOSMESH=30x30x30
    [AFLOW] **********************************************************************************
    [VASP_RUN]RELAX=2
    #[VASP_FORCE_OPTION]KPOINTS=KEEPK,EVEN,KSHIFT_GAMMA_EVEN,KSCHEME_GAMMA,GAMMA,IBZKPT
    [VASP_FORCE_OPTION]SYM=ON
    [VASP_FORCE_OPTION]AUTO_PSEUDOPOTENTIALS=potpaw_PBE
    [VASP_FORCE_OPTION]NBANDS                       
    [VASP_FORCE_OPTION]SPIN=OFF 
    [VASP_FORCE_OPTION]RELAX_MODE=FORCES
    [VASP_FORCE_OPTION]PREC=PHONONS 
    [VASP_FORCE_OPTION]ALGO=FAST     
    [VASP_FORCE_OPTION]RELAX
    [VASP_FORCE_OPTION]TYPE=DEFAULT      
    [VASP_FORCE_OPTION]CONVERT_UNIT_CELL=SPRIM,MINK
    [AFLOW] **********************************************************************************
    [VASP_INCAR_MODE_EXPLICIT]START
    SYSTEM=PbS
    [VASP_INCAR_MODE_EXPLICIT]STOP
    [AFLOW] **********************************************************************************
    [VASP_KPOINTS_MODE_IMPLICIT]
    [VASP_KPOINTS_FILE]KSCHEME=G
    [VASP_KPOINTS_FILE]KPPRA=400
    [AFLOW] **********************************************************************************
    [VASP_POTCAR_MODE_IMPLICIT]
    [VASP_POTCAR_FILE]Pb
    [VASP_POTCAR_FILE]S
    [AFLOW] **********************************************************************************
    [VASP_POSCAR_MODE_EXPLICIT]START
    PbS   [FCC,FCC,cF8]
    1.224745
       0.00000000000000   2.45291416610433   2.45291416610433
       2.45291416610433   0.00000000000000   2.45291416610433
       2.45291416610433   2.45291416610433   0.00000000000000
    1 1
    Direct(2) [A1B1]
       0.00000000000000   0.00000000000000   0.00000000000000  Pb
       0.50000000000000   0.50000000000000   0.50000000000000  S 
    [VASP_POSCAR_MODE_EXPLICIT]STOP
    [AFLOW] **********************************************************************************

#CO - START 170601
*********************************************************************************************************************
NB:
Any shared tags between APL/QHA/AAPPL can be equivalently specified with [AFLOW_APL], [AFLOW_QHA], [AFLOW_AAPL]
in the aflow.in. 
This allows users the flexibility to first run an APL calculation, then proceed with a subsequent analysis (QHA/AAPL)
without the need to change the tags.
*********************************************************************************************************************
#CO - END 170601

7) Quasi-harmonic approximation (Gruneisen parameter and Equation of state calculations)
   Thermal properties calculated using a quasi-harmonic(QH) model.
   This implementation in AFLOW written by Pinku Nath and Stefano Curtarolo,
   and is described in the paper "Computational Materials Science Volume 125, December 2016, Pages 82–91".

   Fillowing properties can be calculated..
   Gruneisen parameper along high-symmetry path
   Mode Gruneisen in q-mesh
   Average Gruneisen parameter
   Mean-square displacement
   Projected mean-square displacement
   Equation of states (Cp, thermal-expansion, Bulk Modulus etc.)

   The current version of  AFLOW-QH can calculate the Gruneisen parameter, heat capacity, vibrational free energy, 
   and other thermal properties, as a function of temperature. The results are written in a variety of formats in different output files.
   The following keywords can be used to calculate AFLOW-QH properties.

   [AFLOW_QHA]CALC                  // Gruneisen parameter option
   [AFLOW_QHA]DISPLACEMENTS=y       // Phonon mean square displacements option
   [AFLOW_QHA]PROJECTION_DIR=h:k:l  // Phonon mean square displacements option along a particular direction (hkl)
   [AFLOW_QHA]EOS=y                 // Equation of state calculation


   7.1) Gruneisen Parameter calculation ([AFLOW_QHA]CALC)
        Calculate Gruneisen Parameters using the quasi-harmonic approximation.

        [AFLOW_QHA]CALC
        [AFLOW_QHA]GP_VOL_DISTORTION_PERCENTAGE=0.03 [percentage of volume distortion magnitude w.r.t equilibrium volume, default value is 0.03]
        [AFLOW_QHA]CUTOFF_FREQ=0.001 [in AMU uint] [default value is 1e-5] frequency cutoff 

        During the first [AFLOW] execution two directories are going to be created with positive and negative distortions. 
	For example with "GP_VOL_DISTORTION_PERCENTAGE=0.03" will create ARUN.APL_PHONON_P0.03 and ARUN.APL_PHONON_M0.03.
 	Data can be extracted from the second run by performing the following two steps
	
  	(1) REMOVE/RENAME LOCK files from current directory and ARUN.APL_PHONON_* directories e,g. [rm -rf LOCK* APL_PHONON_*/LOCK*]
      	[**WARNING DO NOT REMOVE LOCK from any other directories except current and ARUN.APL_PHONON_ directories**]

  	(2) aflow -run=multi -rsort -D ./ > out &
      	    During run a directory "[ARUN.APL.QH.TMP]" is going to be created which stores some necessary files from other sub-directories.
  	
	(3) If you are using a differently named aflow.in and LOCK files (e.g. apl_aflow.in and apl_LOCK) then use the following commands
       	    (a) rm *.bz2 apl_LOCK  ARUN.APL_PHONON_*/apl_LOCK*
       	    (b) aflow --use_aflow.in=apl_aflow.in --use_LOCK=apl_LOCK -run=multi -rsort -D ./ > out & 

 	OUTPUT file names:

  	[aflow.apl.qh.out]           file storing the average and acoustic Gruneisen parameter for user defined temperature
  	[aflow.apl.gp.pdis.out]      file storing the calculated Gruneisen parameter along high-symmetry k-points, format is like PDIS file.
  	[aflow.apl.gp.mesh.out]      file storing the calculated Gruneisen parameter in Monkhorst-Pack k-points and corresponding k-points file name will be [aflow.apl.qpoints.mesh.out]
  	[aflow.apl.gp.freq.mesh.out] file storing Gruneisen parameter and frequencies for each phonon mode in uniform qpoint-mesh
        [aflow.apl.branch.index.out] file storing phonon mode indices that has generated iafter unfolding phonon overlapping modes
	
	If the phonon frequencies are negative there will be no GP calculations.

   7.2) [AFLOW_QHA]DISPLACEMENTS=y
        Mean square displacements [Angstrom unit] in the user defined temperature range. While calculating displacements, GRUNEISEN key should be ON
        [OUTPUT FILE NAME] aflow.apl.displacements.out

   7.3) [AFLOW_QHA]PROJECTION_DIR=h:k:l
         Projected mean square displacements [Angstrom units] along Cartesian axes [h k l] direction in the user defined temperatures. While calculating displacements, GRUNEISEN key should be ON
         [OUTPUT FILE NAME] aflow.apl.projected_displacements.out 

   7.4) Equation of state calculations (EOS=y)

        Quasi-harmonic equation of states can be calculated using EOS keywords. While calculating EOS, GRUNEISEN keyword should be ON. 
	The following additional keywords are necessary to calculate EOS:

        [AFLOW_QHA]EOS_VOLRANGE_DIST:=-2:4:0.5  // Percentage of distortion magnitude w.r.t equilibrium configuration. 
						    // Format is start_volume:end_volume:increment. For example, -2:4:0.5 represents -2% to 4% with increament 0.5%.

        [AFLOW_QHA]EOS_KPOINTS_MODE:=8000:10000:20:50000 // Default values are KPPRA=8000, STATIC_KPPRA=10000, BANDS_GRID=20, electronic DOS=50000

 	During first [AFLOW] execution, two sets of directories with prefixes "ARUN.APL_PHONON_" and "ARUN.APL_STATIC_"  will be created. 
  	For example, "ARUN.APL_PHONON_P1.5" represents phonon directory with positive 1.5% distortion, and "ARUN.APL_STATIC_P3.5" 
  	represents static directory with positive 3.5% distortion respectively.
 
	Data can be extracted in the second run by using following steps:

  	(1) REMOVE LOCK files from current directory and ARUN.APL_PHONON_ directories e,g. rm -rf LOCK* ARUN.APL_PHONON_*/LOCK* 
      	    [**WARNING DO NOT REMOVE LOCK from anyother directories except current and ARUN.APL_PHONON_ directories**]

  	(2) aflow -run=multi -rsort -D ./ > out &

  	(3) If you are using a differently named aflow.in and LOCK files (e.g. apl_aflow.in and apl_LOCK) then use the following commands

    	    (a) rm *.bz2 apl_LOCK  ARUN.APL_PHONON_*/apl_LOCK*
    	    (b) aflow --use_aflow.in=apl_aflow.in --use_LOCK=apl_LOCK -run=multi -rsort -D ./ > out & 

   	Output file names:

        [aflow.apl.FVT.out] => free energies, volume, data for user defined temperatures
        [aflow.apl.eos.out] => EOS for all user defined temperatures
        [aflow.apl.err_fit.out] => errors associated with fitting methods. 

   7.5) Data fitting can be performed either using lmsder [this is a robust and efficient version of the 
         Levenberg-Marquardt algorithm] or lmder [this is an unscaled version of the Levenberg-Marquardt algorithm].
         There are 3 fitting functions available
         [AFLOW_QHA]FITTING_TYPE=BM1   // Murnaghan EOS (Murnaghan F D, Am. J. Math. 49, p235 (1937))
         [AFLOW_QHA]FITTING_TYPE=BM2   // Birch-Murnaghan 3rd-order EOS (Birch F, Phys. Rev. 71, p809 (1947))
         [AFLOW_QHA]FITTING_TYPE=BM3   // Birch-Murnaghan 4th-order EOS, (Journal of Computational Chemistry 2012, 33, 1554–1563)
         default option is BM1 

   7.6) Example qsub file file

          #!/bin/bash
          #
          # These commands set up the Grid Environment for your job:
          #PBS -q debug
          #PBS -l nodes=1:ppn=32
          #PBS -l walltime=2:00:00
          #PBS -N Si_test2
          #PBS -o out.$PBS_JOBID
          #PBS -e err.$PBS_JOBID

          AFLOW="AFLOW_PATH/"
          source /opt/intel/bin/compilervars.sh intel64
          export PATH=$PATH:/MAIN/bin/VASP:$AFLOW

          cd $PBS_O_WORKDIR

          # choose random directory
          touch log.$PBS_JOBID
          $AFLOW/aflow  --monitor >> log.$PBS_JOBID &
          $AFLOW/aflow --force -multi  --machine=qrats  -D ./ >> log.$PBS_JOBID

          #post processing
          # remove all bz2 files except qsub1.sh
          if [ -f qsub1.sh.bz2 ]; then bzip2 -d qsub1.sh.bz2; fi
          if [ -f aflow.fgroup.out.bz2 ]; then rm *.bz2; fi

          #(1) rename all phonon LOCK to LOCK.old
          for i in $( ls LOCK* ARUN.APL_PHONON_*/LOCK* ); do mv $i $i.old; done

          #(2) run aflow-multi with rsort option
          $AFLOW/aflow --force -multi  -rsort --machine=qrats  -D ./ >> log.$PBS_JOBID

	
   7.7) Example aflow.in for quasi-harmonic calculations:


[AFLOW]SYSTEM=Si2
#[AFLOW] single element calculation
[AFLOW] **************************************************************************************************************************
[AFLOW] input file for aflow
[AFLOW_MODE=VASP]
[AFLOW] **************************************************************************************************************************
[AFLOW_MODE_ZIP=bzip2]
[AFLOW_MODE_BINARY=vasp53s]
[AFLOW] **************************************************************************************************************************
[AFLOW] **************************************************************************************************************************
#[AFLOW_MODE_MPI]
[AFLOW_MODE_MPI_MODE]NCPUS=MAX
[AFLOW_MODE_MPI_MODE]COMMAND ="mpirun -np"
[AFLOW_MODE_MPI_MODE]AUTOTUNE
[AFLOW_MODE_MPI_MODE]BINARY="mpivasp53s"
[AFLOW] **************************************************************************************************************************
[AFLOW_SYMMETRY]CALC
#[AFLOW_SYMMETRY]SGROUP_WRITE
#[AFLOW_SYMMETRY]SGROUP_RADIUS=7.77
[AFLOW] **************************************************************************************************************************
#[AFLOW_NEIGHBOURS]CALC
[AFLOW_NEIGHBOURS]RADIUS=7.7
[AFLOW_NEIGHBOURS]DRADIUS=0.1
[AFLOW] **************************************************************************************************************************
[AFLOW_QHA]CALC
[AFLOW_QHA]ENGINE=DM
[AFLOW_QHA]DMAG=0.015
[AFLOW_QHA]MINATOMS=100
#[AFLOW_QHA]SUPERCELL=2x2x2
[AFLOW_QHA]POLAR=y
[AFLOW_QHA]DISPLACEMENTS=y
[AFLOW_QHA]PROJECTION_DIR=1:1:1                   
[AFLOW_QHA]GP_VOL_DISTORTION_PERCENTAGE=0.03      
[AFLOW_QHA]EOS=n
[AFLOW_QHA]EOS_VOLRANGE_DIST=-2:4:0.5
[AFLOW_QHA]EOS_KPOINTS_MODE=32768:10000:20:100000
[AFLOW_QHA]DC=y
[AFLOW_QHA]DPM=y
[AFLOW_QHA]ZEROSTATE=y
#[AFLOW_QHA]DCINITSG=227
#[AFLOW_QHA]DCUSERPATH=Gamma-X|X-U|K-Gamma|Gamma-L
[AFLOW_QHA]DOS=y
[AFLOW_QHA]TP=y
[AFLOW_QHA]TPT=0:2000:10
AFLOW] **************************************************************************************************************************
[VASP_FORCE_OPTION]CHGCAR=OFF                         // ON | OFF (default ON)
[VASP_FORCE_OPTION]SYM=ON                             // ON | OFF  (default ON)
[VASP_FORCE_OPTION]PREC=PHONONS 
[VASP_FORCE_OPTION]AUTO_PSEUDOPOTENTIALS=potpaw_PBE   // pot_LDA | pot_GGA | potpaw_LDA | potpaw_GGA | potpaw_PBE
[VASP_FORCE_OPTION]ALGO=NORMAL                        // (NORMAL | VERYFAST | FAST | ALL | DAMPED), PRESERVED (default=NORMAL)
[VASP_FORCE_OPTION]TYPE=DEFAULT                       // METAL | INSULATOR | SEMICONDUCTOR | DEFAULT (default DEFAULT)
[VASP_FORCE_OPTION]CONVERT_UNIT_CELL=SPRIM            // SPRIM, SCONV, NIGGLI, MINK, INCELL, COMPACT, WS, CART, FRAC, PRES
[AFLOW] **************************************************************************************************************************
[VASP_INCAR_MODE_EXPLICIT]START
[AFLOW]SYSTEM=Si2
NELM = 120
NELMIN=2
LPLANE=.TRUE.
LREAL=.FALSE.
LSCALU=.FALSE.
#NBANDS=XX                                            # for hand modification
#IALGO=48                                             # for hand modification
[VASP_INCAR_MODE_EXPLICIT]STOP
[AFLOW] **************************************************************************************************************************
[VASP_KPOINTS_MODE_IMPLICIT]
[VASP_KPOINTS_FILE]KSCHEME=M
[VASP_KPOINTS_FILE]KPPRA=400
[VASP_KPOINTS_FILE]STATIC_KSCHEME=M
[VASP_KPOINTS_FILE]STATIC_KPPRA=400
[VASP_KPOINTS_FILE]BANDS_LATTICE=FCC
[VASP_KPOINTS_FILE]BANDS_GRID=20
[AFLOW] **************************************************************************************************************************
[VASP_POTCAR_MODE_IMPLICIT]
[VASP_POTCAR_FILE]Si
[AFLOW] potpaw_PBE: Si
[AFLOW] COMPOSITION_PP=|Si2|
[AFLOW] COMPOSITION=|Si2|
[AFLOW] VOLUME(A^3)=|0|
[AFLOW] MASS(amu)=|28.0855|
[AFLOW] **************************************************************************************************************************
[AFLOW] Aflow automatically generated (aflow_avasp.cpp)
[AFLOW] **************************************************************************************************************************
[AFLOW] aflow/aconvasp/aqe/apennsy in
[AFLOW] ************************************************************************************************************************
[VASP_POSCAR_MODE_EXPLICIT]START
Si1   [FCC,FCC,cF8] (STD_PRIM doi:10.101
   1.22474500000000
     0.0000000000000000    2.2320315365242354    2.2320315365242354
     2.2320315365242354    0.0000000000000000    2.2320315365242354
     2.2320315365242354    2.2320315365242354    0.0000000000000000
   Si
   2
Direct
  0.0000000000000000  0.0000000000000000  0.0000000000000000
  0.7500000000000000  0.7500000000000000  0.7500000000000000
[VASP_POSCAR_MODE_EXPLICIT]STOP
[AFLOW] **************************************************************************************************************************

8) AAPL Automatic Anharmonic Phonon Library (Lattice thermal conductivity calculations).
   AAPL predicts lattice thermal conductivity commputing the anharmoic force constants and solving the Boltzmann Transport Equation (BTE)
   This implementation in AFLOW was developed by Jose J. Plata
   and is described in "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX".

For lattice thermal conductivity calculation, generally, you have to setup
3 things: 1) Tell to aflow you want a lattice thermal conductivity calculation. 2) Define the parameters for computing the harmonic and anharmonic force constants 
3) Choose how to solve the BTE.

NOTE: Lattice thermal conductivity prediction requires accurate forces (harmonic and anharmonic force constants). This means that the wave function
must be very well converged. A new accuracy flag has been included for this purpose:
	
	[VASP_FORCE_OPTION]PREC=PHONONS

   8.1) Lattice thermal conductivity calculation ([AFLOW_AAPL]CALC).
        Calculate lattice thermal conductivity solving the BTE. Harmonic force constants are also required.

        [AFLOW_AAPL]CALC

   8.2) Anharmonic force constants calculations.
        Different parameters should be include to setup the anharmonic force constants calculations
         
      	8.2.1)	Supercell. The most important setting is the supercell dimension.
         	Supercell can be set in three ways as it was explained previously.
         	EXAMPLES:
         	[AFLOW_AAPL]SUPERCELL=4x4x4
        
   	8.2.2)	Distorsion magnitude. It controls de magnitude of the distorsion of the atoms 
		with respect to their equilibrium positions [Angstrom units].
         	EXAMPLES:
         	[AFLOW_AAPL]TDMAG=0.015 //Default: 0.015

	8.2.3)	Cutoff. A cutoff for the anharmonic force constant should be applied 
		because of the high computational cost. This radious can be define as a distance,
		CUT_RAD [Angstrom], or as coordination shells, CUT_SHELL [number of neighnor shells]. 
         	EXAMPLES:
         	[AFLOW_AAPL]CUT_RAD=5.3 //Default: 5.0
         	[AFLOW_AAPL]CUT_SHELL=3   //Default: 3

	8.2.4)	Sumrule. Sumrules are critical to obtain accurate
		lattice thermal conductivity values. The variable SUMRULE defines the 
		thershold [ev/Angstrom/Angstrom/Angstrom].
         	EXAMPLES:
		[AFLOW_AAPL]SUMRULE=1E-7 //Default: 1E-5  
		 
   8.3) Solving the BTE

	8.3.1)	Solver. There are two pissible oprtions: the relaxation time approximation (RTA)
     		and the full solution using a self consistent field (FULL) using as initial iteration the 
		RTA solution.	
         	EXAMPLES:
		[AFLOW_AAPL]BTE=FULL //Default: RTA 

	8.3.2)	Grid. Phonon scattering is evaluated in a discrete mesh of q-points. The higher the density
 		the more accurate the final values will be. However, a very dense mesh can create 
		memory issues.
         	EXAMPLES:
		[AFLOW_AAPL]THERMALGRID=21x21x21 // Default: 16x16x16 

       	8.3.3)	Elastic/isotopic effects. Elastic scattering is produce by the abundance of the 
		isotopes of each element. The kind of scattering can be included in the calculation of the 
		relaxation time using the ISOTOPE flag 
         	EXAMPLES:
		[AFLOW_AAPL]ISOTOPE=y //Default: no

        8.3.4)	Scattering at boundaries. Two different methods are implemented to predict the lattice 
                thermal conductivity for polycrystalline samples with an average size of the crystals.
                Both methods cannot be used at the same time. Computing the cummulative kappa, (CUMULATIVE)
                aflow compute the lattice thermal conductivity for different size based on the RTA solution.
         	      Using the BOUNDARY flag you will include an extra term for the scattering time. However,
                a specific average size [nm] must be specified.  
         	EXAMPLES:
		[AFLOW_AAPL]CUMULATIVE=y //Default: no
                 OR
		[AFLOW_AAPL]BOUNDARY=y   //Default: no
		[AFLOW_AAPL]NANO_SIZE=y  //Default: no
           
	8.3.5) 	Temperature range Generally, AAPL is calculating the 
               	lattice thermal conductivity for a given range of temperature. This range can be set as
 
	    	[AFLOW_AAPL]TCT=0:2000:10 // Default 100:700:100
 
         	what means it will calculate lattice thermal conductivity from 0 to 2000 K by step of 10 K. This
         	is also a default setting.
 
   8.4) Example aflow.in for lattice thermal conductivity calculation:


[AFLOW]SYSTEM=Si2
#[AFLOW] single element calculation
[AFLOW] **************************************************************************************************************************
[AFLOW] input file for aflow
[AFLOW_MODE=VASP]
[AFLOW] **************************************************************************************************************************
[AFLOW_MODE_ZIP=bzip2]
[AFLOW_MODE_BINARY=vasp53s]
[AFLOW] **************************************************************************************************************************
[AFLOW] **************************************************************************************************************************
#[AFLOW_MODE_MPI]
[AFLOW_MODE_MPI_MODE]NCPUS=MAX
[AFLOW_MODE_MPI_MODE]COMMAND ="mpirun -np"
[AFLOW_MODE_MPI_MODE]AUTOTUNE
[AFLOW_MODE_MPI_MODE]BINARY="mpivasp53s"
[AFLOW] **************************************************************************************************************************
[AFLOW_SYMMETRY]CALC
#[AFLOW_SYMMETRY]SGROUP_WRITE
#[AFLOW_SYMMETRY]SGROUP_RADIUS=7.77
[AFLOW] **************************************************************************************************************************
#[AFLOW_NEIGHBOURS]CALC
[AFLOW_NEIGHBOURS]RADIUS=7.7
[AFLOW_NEIGHBOURS]DRADIUS=0.1
[AFLOW] **************************************************************************************************************************
[AFLOW_AAPL]CALC
#[AFLOW_AAPL] HARMONIC
[AFLOW_AAPL]ENGINE=DM
[AFLOW_AAPL]DMAG=0.015
[AFLOW_AAPL]MINATOMS=100
#[AFLOW_AAPL]SUPERCELL=4x4x4
[AFLOW_AAPL]POLAR=y
[AFLOW_AAPL]DC=y
[AFLOW_AAPL]DPM=y
[AFLOW_AAPL]ZEROSTATE=y
[AFLOW_AAPL]DOS=y
[AFLOW_AAPL]TP=y
[AFLOW_AAPL]TPT=400:2000:10
#[AFLOW_AAPL] ANHARMONIC
[AFLOW_AAPL]TDMAG=0.015
[AFLOW_AAPL]CUT_SHELL=4
#[AFLOW_AAPL]CUT_RAD=4.5
[AFLOW_AAPL]SUMRULE=1E-5
[AFLOW_AAPL]BTE=FULL
[AFLOW_AAPL]THERMALGRID=21x21x21
[AFLOW_AAPL]ISOTOPE=y
#[AFLOW_AAPL]CUMULATIVEK=y
#[AFLOW_AAPL]BOUNDARY=n
#[AFLOW_AAPL]NANO_SIZE=100
[AFLOW_AAPL]TCT=200:700:20
AFLOW] **************************************************************************************************************************
[VASP_FORCE_OPTION]CHGCAR=OFF                         // ON | OFF (default ON)
[VASP_FORCE_OPTION]SYM=ON                             // ON | OFF  (default ON)
[VASP_FORCE_OPTION]PREC=PHONONS 
[VASP_FORCE_OPTION]AUTO_PSEUDOPOTENTIALS=potpaw_PBE   // pot_LDA | pot_GGA | potpaw_LDA | potpaw_GGA | potpaw_PBE
[VASP_FORCE_OPTION]ALGO=NORMAL                        // (NORMAL | VERYFAST | FAST | ALL | DAMPED), PRESERVED (default=NORMAL)
[VASP_FORCE_OPTION]TYPE=DEFAULT                       // METAL | INSULATOR | SEMICONDUCTOR | DEFAULT (default DEFAULT)
[VASP_FORCE_OPTION]CONVERT_UNIT_CELL=SPRIM            // SPRIM, SCONV, NIGGLI, MINK, INCELL, COMPACT, WS, CART, FRAC, PRES
[AFLOW] **************************************************************************************************************************
[VASP_INCAR_MODE_EXPLICIT]START
[AFLOW]SYSTEM=Si2
NELM = 120
NELMIN=2
LPLANE=.TRUE.
LREAL=.FALSE.
LSCALU=.FALSE.
#NBANDS=XX                                            # for hand modification
#IALGO=48                                             # for hand modification
[VASP_INCAR_MODE_EXPLICIT]STOP
[AFLOW] **************************************************************************************************************************
[VASP_KPOINTS_MODE_IMPLICIT]
[VASP_KPOINTS_FILE]KSCHEME=M
[VASP_KPOINTS_FILE]KPPRA=400
[VASP_KPOINTS_FILE]STATIC_KSCHEME=M
[VASP_KPOINTS_FILE]STATIC_KPPRA=400
[VASP_KPOINTS_FILE]BANDS_LATTICE=FCC
[VASP_KPOINTS_FILE]BANDS_GRID=20
[AFLOW] **************************************************************************************************************************
[VASP_POTCAR_MODE_IMPLICIT]
[VASP_POTCAR_FILE]Si
[AFLOW] potpaw_PBE: Si
[AFLOW] COMPOSITION_PP=|Si2|
[AFLOW] COMPOSITION=|Si2|
[AFLOW] VOLUME(A^3)=|0|
[AFLOW] MASS(amu)=|28.0855|
[AFLOW] **************************************************************************************************************************
[AFLOW] Aflow automatically generated (aflow_avasp.cpp)
[AFLOW] **************************************************************************************************************************
[AFLOW] aflow/aconvasp/aqe/apennsy in
[AFLOW] ************************************************************************************************************************
[VASP_POSCAR_MODE_EXPLICIT]START
Si1   [FCC,FCC,cF8] (STD_PRIM doi:10.101
   1.22474500000000
     0.0000000000000000    2.2320315365242354    2.2320315365242354
     2.2320315365242354    0.0000000000000000    2.2320315365242354
     2.2320315365242354    2.2320315365242354    0.0000000000000000
   Si
   2
Direct
  0.0000000000000000  0.0000000000000000  0.0000000000000000
  0.7500000000000000  0.7500000000000000  0.7500000000000000
[VASP_POSCAR_MODE_EXPLICIT]STOP
[AFLOW] **************************************************************************************************************************

****************************************************************************************************
*                                                                                                  *
*                        aflow - STEFANO CURTAROLO Duke University 2003-2018                       *
*                            High-Throughput ab-initio Computing Project                           *
*                                                                                                  *
****************************************************************************************************
