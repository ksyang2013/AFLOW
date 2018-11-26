****************************************************************************************************
*                                                                                                  *
*                        aflow - STEFANO CURTAROLO Duke University 2003-2018                       *
*                            High-Throughput ab-initio Computing Project                           *
*                                                                                                  *
****************************************************************************************************
LATEST VERSION OF THE FILE:                            materials.duke.edu/AFLOW/README_AFLOW_APL.TXT
****************************************************************************************************

Phonon calculations with deformations written by Michal Jahnatek, Pinku Nath, Jose J. Plata, 
Corey Oses, Marco Esters, and Stefano Curtarolo.

The current version of the Automatic Phonon Library (APL) can calculate the phonon dispersion
curves, phonon density of states, and thermal properties. All are calculated in the harmonic
approximation. Two methods of calculation can be selected: the direct (supercell) method (DM) and
the linear response method (LR). The resulting output files can be transformed into agr, ps, eps,
or png by the tool apl2agr (located in SOURCES/APL2AGR).

The Automatic Anharmonic Phonon Library (AAPL) is an expansion of APL and can calculate the lattice
thermal conductivity of a material via the Boltzmann Transport Equation using anharmonic
interatomic force constants. By default, it uses only three-phonon processes but four-phonon
processes can be captured as well.

The default values of APL and AAPL can be changed by editing the appropriate values in aflow.rc.

----------------------------------------------------------------------------------------------------

Table of Contents:

1.    Running APL

2.    Choosing a Phonon Calculation Method
2.1.    Common Parameters Used by Both Methods
2.1.1.    Supercell parameters
2.1.2.    k-Point Parameters
2.1.3.    Polar Corrections
2.2.    Direct Method Parameters
2.3.    Linear Response Method Parameters

3.    Calculating Phonon Properties
3.1.    Phonon Dispersion Curves
3.2.    Phonon Density of States
3.3.    Thermodynamic Properties

4.    Saving and Loading APL Results (Hibernation)

5.    Inheritance from AFLOW

6.    Example APL Input File

7.    The Automatic Anharmonic Phonon Library
7.1     Running AAPL Calculations
7.2.    Anharmonic Interatomic Force Constants
7.2.1.    Cluster Parameters
7.2.2.    Sum Rule Parameters
7.2.3.    Four-Phonon Processes
7.3.    Lattice Thermal Conductivity Calculations
7.4.    Example AAPL Input File

----------------------------------------------------------------------------------------------------

1. Running APL

APL expects a well relaxed structure on input since there is no additional relaxation inside APL.
Such a structure needs to be prepared in a separate AFLOW run.

To perform a phonon calculation, the following line has to be present in the aflow.in file:

    [AFLOW_APL]CALC

If it does not exist, the job is not interpreted as a phonon calculation and no other phonon
settings are read.


2. Choosing a Phonon Calculation Method

The user has the possibility to choose two methods of calculation using the keyword ENGINE: the
direct method (DM), also known as supercell method, and the linear response method (LR). To 
choose the DM, one has to write:

    [AFLOW_APL]ENGINE=DM

This is also a default setting. For the linear response method, it is:

    [AFLOW_APL]ENGINE=LR

Each method has its own parameters, but there are also some common ones that can be used to tune
the whole calculation.

2.1. Common Parameters Used by Both Methods

2.1.1. Supercell parameters

The most important setting is the supercell dimensions (a supercell is required by both methods).
The supercell can be set by the following commands:

    [AFLOW_APL]SUPERCELL=numberXnumberXnumber
    [AFLOW_APL]MINSHELL=number
    [AFLOW_APL]MINATOMS=number
    [AFLOW_APL]MINATOMS_RESTRICTED=number

Examples:

    * Direct specification of the supercell dimensions, e.g.:
            
        [AFLOW_APL]SUPERCELL=4x4x4

    * Specifying the minimal number of nearest neighbor shells, e.g.:

        [AFLOW_APL]MINSHELL=6

      The default settings is MINSHELL=6 and APL will construct the appropriate supercell. HOWEVER,
      THIS IS NOT THE RECOMMENDED WAY BECAUSE THE INCLUSION OF 6 SHELLS IS SUFFICIENT ONLY FOR VERY
      SIMPLE MATERIALS. There is unfortunately no rule on how to find the ideal number of shells
      for the calculation.

    * Require a minimum number of atoms, e.g.:

        [AFLOW_APL]MINATOMS=100

      The code will homogeneously increase the volume of the supercell until it finds a
      SUPERCELL(i,j,k) containing at least the number of atoms specified.

      If it is desired that the supercell is expanded equally in all directions, i.e. to use
      SUPERCELL(i,i,i), use the flag MINATOMS_RESTRICTED as in:

        [AFLOW_APL]MINATOMS_RESTRICTED=100

For best results, it is recommended to specify the supercell dimensions directly. We found that
4x4x4 is usually enough for all fcc systems; 3x3x2 for hcp; bcc cells are complicated, but
generally, it should be more that 3x3x3. Typically, the supercell with about 100 atoms should be
fine for the majority of crystal systems. The MINATOMS option might be useful for large unit cells.

2.1.2. k-Point Parameters

The k-point mesh can be set separately for APL calculations. The settings are similar to the
general AFLOW k-point parameters:

    * The number of k-points per reciprocal atoms can be set via:

        [AFLOW_APL]KPPRA=number

      By default, AFLOW takes the number of KPOINTS (KPPRA) from the "[VASP_KPOINTS_FILE]KPPRA"
      entry. The "[AFLOW_APL]KPPRA" flag overrides this for APL calculations. See README_AFLOW.TXT
      for more information on KPPRA.

    * The k-point scheme can be set using:

        [AFLOW_APL]KSCHEME=scheme

      By default, AFLOW takes the KPOINTS scheme from the "[VASP_KPOINTS_FILE]KSCHEME" entry. The
      "[AFLOW_APL]KSCHEME" overrides this setting. See README_AFLOW.TXT for more information on
      KSCHEME.

    * The k-point parity can be set using:

        [AFLOW_APL]KPOINTS=EVEN | ODD

      By default, AFLOW takes the KPOINTS parity from "[VASP_FORCE_OPTION]KPOINTS" entry, if
      specified. "[AFLOW_APL]KPOINTS" overrides this setting. See README_AFLOW.TXT for more
      information on KPOINTS.

2.1.3. Polar Corrections

If the studied system is polar, dipole-dipole interactions need to be taken into account. This
can be specified by the keyword POLAR. This keyword can be set in the aflow.in file as follows:

    [AFLOW_APL]POLAR=y

This instructs APL to calculate the Born effective charge tensor and the dielectric constant
tensor using density functional perturbation theory (DFPT). The default setting is POLAR=n, i.e.
no polar calculations are performed.

DFPT calculations require VASP 5.2 or higher. APL will check if the appropriate binary file is
present upon start-up, and will refuse to continue if the appropriate VASP version is not found.
They also require large stack sizes, so it is critical that the stack size is set to unlimited.
Otherwise, VASP will likely crash at some point during the calculation. To set the stack size to
unlimited, the following line needs to be written either in the .bashrc file or the job script:

    ulimit -s unlimited

The dielectric constant tensor and the Born effective charge tensor is by default calculated using
the tag LCALCEPS in the VASP INCAR file, but can be calculated using LEPSILON as well. To switch to
LEPSILON, the following lines in aflow_apl.h need to be changed from:

    #undef   AFLOW_APL_VASP_USE_LEPSILON
    #define  AFLOW_APL_VASP_USE_LCALCEPS

To:

    #define  AFLOW_APL_VASP_USE_LEPSILON
    #undef   AFLOW_APL_VASP_USE_LCALCEPS

2.2. Direct Method Parameters

The direct method (DM) requires the calculation of the forces for each unique atom in the primitive
cell along three independent directions. However, some of these directions are related by symmetry
and it is thus unnecessary to perform DFT calculations for all three distortions. For example, if
one distortion direction can by symmetry produce three linearly independent vectors, the atom only
needs to be distorted in one direction. The remaining two can be recovered by the symmetry of the
system as is the case for an fcc lattice. APL will try to determine these distortions to reduce
the number of DFT calculations that need to be performed.

Generally, APL will distort an atom along the [100], [010], [001], [110], [101], [011], and the
[111] direction. It will apply all site symmetry operations, and the direction that produces the
highest number of new independent vectors will by used for the DFT calculations.

To filter out potential anharmonic effects from the forces, APL distorts the atom into the positive
and negative direction of the chosen distortions. The final force acting on each is then calculated
as 0.5 * (f(+) - f(-)).

There are multiple parameters that control how the distortions for the direct method are generated.

    * The distortion magnitude in Angstrom is set via the DMAG entry, e.g.:

        [AFLOW_APL]DMAG=0.015

      The default is 0.015 Angstrom, which is generally accepted to be sufficient for all systems.

    * If the atoms shall only be distorted along the cell axes ([100], [010], [001]), the DXYZONLY
      flag needs to be set to true, i.e.:

        [AFLOW_APL]DXYZONLY=y

      The default is DXYZONLY=n, which will also search along the face and body diagonals. DXYZONLY
      is particularly useful for cells where one primitive cell vector is much longer than the
      others.

    * To prevent the automatic generation of negative distortions, the DPM entry needs to be unset:

        [AFLOW_APL]DPM=n

      In this case, APL will determine whether the negative distortions are required or if they
      can be recovered by symmetry. The default is DPM=y, i.e. negative distortions are always
      created.

If the primitive structure is not well relaxed or cannot be relaxed within a given precision, the
unstrained state has some non-zero forces. In this case, it is possible to subtract these forces
from all distorted states. This can be done using the keyword:

    [AFLOW_APL]ZEROSTATE=y

If activated, APL will generate a calculation of the undistorted supercell and the obtained forces
will be subtracted from all other strained calculations. The default is ZEROSTATE=n.
WARNING: This method is not well tested yet and should be used with caution.

2.3.  Linear Response Method Parameters

The use of the linear response method as implemented in VASP requires VASP version 5.2 or higher.
APL will check if the appropriate binary file is present upon start-up, and will refuse to continue
if the appropriate VASP version is not found.

LR calculations require large stack sizes, so it is critical that the stack size is set to
unlimited. Otherwise, VASP will likely crash at some point during the calculation. To set the stack
size to unlimited, the following line needs to be written either in the .bashrc or the job script:

    ulimit -s unlimited

In contrast to the DM method, APL is not responsible for finding the appropriate set of distortion
vectors for LR calculations as it is all done by VASP. Hence, long run DFT calculations are to be
expected. Please note, the VASP is using the mode XYZONLY so that more distortions are generated
than in the DM method.

APL will use the file DYNMAT, which is generated by VASP, as the basis for post-processing.


3. Calculating Phonon Properties

APL can calculate phonon dispersion curves, phonon density of states, and various thermodynamic
properties. They are activated by the following keywords:

    [AFLOW_APL]DC=y
    [AFLOW_APL]DOS=y
    [AFLOW_APL]TP=y

The default is no for all. Each calculation mode has its own set of parameters, which are described
in the following sections.

3.1. Phonon Dispersion Curves

The calculated dispersion curves are saved in the file aflow.apl.phonon_dispersion.out, where the
1st column is the ID of points in the sub-path, the 2nd is the position in the whole path, and the
remaining columns are the energies/frequencies of the phonon branches.

The phonon dispersion curve calculations have three adjustable parameters: the dispersion path,
the point density, and the frequency unit.

    * There are multiple ways to set the path of the phonon dispersion curve. By default, APL 
      calculates the symmetry of the cell via AFLOW-SYM and uses the default path as introduced
      in [doi=10.1016/j.commatsci.2010.05.010].

      There are two ways to specify a different path. The first is to provide the lattice or the
      spacegroup number of the desired path. The second is to specify the high symmetry points of
      the path explicitly.

      1. Specify the path via the lattice is done using the keyword DCINITLATTICE:

          [AFLOW_APL]DCINITLATTICE=lattice

        Where lattice is the keyword for the Bravais lattice (e.g. RHL). The default paths for each
        lattice is published in [doi=10.1016/j.commatsci.2010.05.010].
        Alternatively, the space group number can also be used using DCINITSG, e.g.:

          [AFLOW_APL]DCINITSG = 166

      2. To specify the path explicitly, the special points need to be defined first. This is done
         by defining the coordinates of each path end point and by assigning labels to these
         points. To define the coordinates, fractional or Cartesian coordinates can be used. They
         are defined using the DCINITCOORDSFRAC and DCINITCOORDSCART keywords, respectively:

             [AFLOW_APL]DCINITCOORDSFRAC=0,0,0;0.5,0.5,0.5;...
             [AFLOW_APL]DCINITCOORDSCART=0,0,0;0.5,0.5,0.5;...

         Each coordinate triplet has to be separated by a semicolon. The labels are defined using
         DCINITCOORDSLABELS as in:

             [AFLOW_APL]DCINITCOORDSLABELS=G,L,M,X...

         Where the first label corresponds to the first coordinate triplet in DCINITCOORDSFRAC or
         DCINITCOORDSCART. The labels have to be separated by a comma, not a semicolon. Finally,
         the path can be specified using the aforementioned labels using the following format:

             [AFLOW_APL]DCUSERPATH = G-X|X-U|K-G|G-L

    * The number of points for each sub-path is set with the DCPOINTS keyword, e.g.:

        [AFLOW_APL]DCPOINTS=100

    * For the frequency, a variety of options are available. They can be set as:

        [AFLOW_APL]FREQFORMAT="format"

      Where "format" needs to be replaced with the appropriate format specifier. These are "HERTZ",
      "THZ", "MEV", and "RECIPROCAL_CM" or "CM-1".

      These frequencies do not include the factor 2 pi. To include it, "OMEGA" needs to be added
      to the format specifier, e.g.:

        [AFLOW_APL]FREQFORMAT="THZ|OMEGA"

      A colon (:), a semicolon(;), or a comma (,) may also be used as delimiters.

      At last, it can be specified whether imaginary frequencies should be set to 0.0 or set to a
      negative value. If they are set to negative values, "ALLOW_NEGATIVE" needs to be added to the
      format specifier, e.g.:

        [AFLOW_APL]FREQFORMAT="THZ|ALLOW_NEGATIVE"

      If imaginary frequencies should be set to 0.0, "ALLOW_NEGATIVE" must not be included, e.g.:

        [AFLOW_APL]FREQFORMAT="THZ"

      Sets all imaginary frequencies to 0.0 while all real frequencies are given in THz.
      The "ALLOW_NEGATIVE" specifier may be combined with "OMEGA", e.g. in:

        [AFLOW_APL]FREQFORMAT="THZ|OMEGA|ALLOW_NEGATIVE"

      The default specifier is "THZ|ALLOW_NEGATIVE".

3.2  Phonon Density of States

The phonon density of states (pDOS) are saved into the file aflow.apl.phonon_dos.out. Besides the
pDOS, this file also contains the frequencies in THz, the wavenumbers in inverse cm, and the
energies in meV.

There are multiple parameters for the pDOS: the Brillouin zone integration method, the k-point
mesh, the number of energy points, and the smearing width.

    * There are two implementations for the integration of the Brillouin zone: the Linear Tetrahedron
      method (LT) and the Root Sampling method (RS). They can be set as:

        [AFLOW_APL]DOSMETHOD=LT(RS)

      Default is DOSMETHOD=LT.

    * The mesh of q-points used for numerical integration should be set as:

        [AFLOW_APL]DOSMESH=21x21x21

      This is also the default setting. APL will automatically find the irreducible part of the
      of BZ (the list of irreducible q-points with appropriate summation weights). Only for these
      q-points is a dynamical matrix constructed and diagonalized. For the LT method, the list of
      irreducible tetrahedrons with their weights is constructed and used in summation.

    * The number of points in the pDOS can be set via the keyword DOSPOINTS, e.g.:

        [AFLOW_APL]DOSPOINTS=2000

      This is also the default setting.

    * The pDOS can be smeared using a Gaussian function. To specify the width of this function, the
      DOSSMEAR tag can be set, e.g.:

        [AFLOW_APL]DOSSMEAR=0.05

      This is especially important for the RS method. The default is 0.05 for the RS method and
      0.0 for the LT method (i.e. no smearing). 

3.3. Thermodynamic Properties

APL can calculate a variety of thermodynamic properties. These are the zero point vibrational
energy (U0), the internal energy (U), the vibrational free energy (Fv), the vibrational
entropy (Sv), and the isochoric specific heat (cv) for a given temperature range. The results are
written into the file aflow.apl.thermodynamic_properties.out.

The temperature range can be set as:

        [AFLOW_APL]TPT=TSTART:TEND:TSTEP

Where TSTART and TEND are the initial and final temperature, respectively, and TSTEP is the
temperature step size. The default setting is 0:2000:10, i.e. it calculates the properties from 0 K
to 2000 K in 10 K steps.

4. Saving and Loading APL Results APL Calculations (Hibernation)

In order to avoid recalculating the force constant matrix with every run, APL has implemented a
hibernation system which saves the forces, the force constants, the Born effective charge tensor,
and the dielectric constant tensor into an output file (aflow.apl.harmonicIFCs.xml).

This behavior can be switched on(off) using the HIBERNATE keyword: 

    [AFLOW_APL]HIBERNATE=y(n)

It is turned on by default.


5. Inheritance from AFLOW

The generated ARUN.APL***/aflow.in input files inherit some options from the master aflow.in
If present, the following options are passed to the local aflow.in:

    [VASP_FORCE_OPTION]AUTO_PSEUDOPOTENTIALS=...
    [VASP_FORCE_OPTION]SPIN=...
    [VASP_FORCE_OPTION]TYPE=...
    [VASP_FORCE_OPTION]AUTO_MAGMOM=...
    [VASP_FORCE_OPTION]BADER=...
    [VASP_FORCE_OPTION]ELF=...
    [VASP_FORCE_OPTION]ABMIX=...
    [VASP_FORCE_OPTION]KPOINTS=...
    [VASP_FORCE_OPTION]IGNORE_AFIX=...

Other options are added by default:

    [VASP_FORCE_OPTION]WAVECAR=OFF
    [VASP_FORCE_OPTION]CHGCAR=OFF
    [VASP_FORCE_OPTION]PREC=ACCURATE
    [VASP_FORCE_OPTION]ALGO=NORMAL

To avoid NPAR warnings from VASP 5.2+, the following option is added:

    [VASP_FORCE_OPTION]IGNORE_AFIX=NPARC

6)  Example Input File

The following is an example APL input file for PbS using the linear response method. For the
direct method, use [AFLOW_APL]ENGINE=DM and uncomment #[AFLOW_APL]DMAG=0.015.

    [AFLOW] **********************************************************************************
    [AFLOW]SYSTEM=PbS
    [AFLOW] **********************************************************************************
    [AFLOW_MODE=VASP]
    [AFLOW_MODE_ZIP=xz]
    [AFLOW_MODE_BINARY=vasp54s]
    [AFLOW] **********************************************************************************
    #[AFLOW_MODE_MPI]
    [AFLOW_MODE_MPI_MODE]NCPUS=MAX
    [AFLOW_MODE_MPI_MODE]COMMAND ="mpirun -np"
    [AFLOW_MODE_MPI_MODE]AUTOTUNE
    [AFLOW_MODE_MPI_MODE]BINARY="mpivasp54s"
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

***************************************************************************************************
NB:
Any shared tags between APL/QHA/AAPL can be equivalently specified with [AFLOW_APL], [AFLOW_QHA],
[AFLOW_AAPL] in the aflow.in. This allows users the flexibility to first run an APL calculation and
then proceed with a subsequent analysis (QHA/AAPL) without the need to change the tags.
***************************************************************************************************


7. The Automatic Anharmonic Phonon Library 

7.1. Running AAPL Calculations

The Automatic Anharmonic Phonon Library (AAPL) is the anharmonic extension to APL and calculates
anharmonic interatomic force constants (IFCs). These can then be used to calculate the lattice
thermal conductivity by solving the Boltzmann Transport Equation (BTE). AAPL can use three-phonon
and four-phonon processes to calculate thermal conductivities.

To use AAPL, the following needs to be added to the aflow.in file:

    [AFLOW_AAPL]CALC

This overrides all other APL CALC options. Lattice thermal conductivity calculations require very
accurate forces to get accurate harmonic and anharmonic IFCs. This requires the primitive cell to
be relaxed with very tight convergence criteria. To do this, the following accuracy needs to be
chosen for the relaxation calculations:

    [VASP_FORCE_OPTION]PREC=PHONONS

An AAPL calculation requires three separate AFLOW runs:

    1. In the first run, the distortions for the anharmonic IFCs (and harmonic IFCs, if not
       already present) are created.

    2. The DFT calculations are performed. However, due to the large number of DFT calculations,
       it is recommended to run each calculation in its own separate AFLOW run.

    3. Post-processing: the harmonic and anharmonic IFCs and the thermal conductivity tensor are
       calculated.

AAPL shares the supercell parameters and the HIBERNATE keyword with APL. AAPL hibernate file store
the clusters and the anharmonic IFCs. IT IS STRONGLY RECOMMENDED TO USE THE HIBERNATE FUNCTION! The
number of calculations for AAPL is very large, so processing the DFT calculations can take a long
time whereas reading the hibernate files takes only a couple of seconds. Moreover, if LR is chosen
for the APL calculations, the DMAG keyword still needs to be specified for AAPL.

7.2. Anharmonic Interatomic Force Constants

Anharmonic IFCs are calculated using the direct method. It is necessary to define a cut-off radius
beyond which the interactions between atoms are considered negligible. Within this shell, clusters
of atoms (three atoms for third order IFCs, four atoms for fourth order IFCs, etc.) are created and
distorted. The number of clusters and distortions is then reduced by symmetry to save computational
resources. The resulting distorted structures are then output into their own aflow.in files inside
a subfolder. These folder all have the name:

    ARUN.AAPL.CN_XX_A1_at1_D1_d1_A2_at2_D2_d2_...

Where the N in CN is replaced with the phonon process order; at1, at2, etc. are the atomic indices
of the distorted atoms; and d1, d2, etc. are the distortion indices corresponding to the atom. They
range from 0 to 5 and correspond to Cartesian +x, +y, +z, -x, -y, -z. Note that when one distortion
index appears twice on the same atom, the atom is not distorted twice into the same direction.

Once the DFT calculations are finished, the forces are read and the IFCs are determined using a
self-consistent field. 

7.2.1. Cluster Parameters

The supercell and the distortion magnitude are defined the same way as in APL. See Section 2.1.1.
for the supercell parameters and Section 2.2. for the distortion magnitude.

The cut-off radius described above can be defined in two ways:

  * Define a cut-off radius directly:

        [AFLOW_AAPL]CUT_RAD=5.0

    This is also the default value.

  * Define the number of coordination shells that should be included in the interaction volume.
    This can be achieved via the CUT_SHELL parameter.

        [AFLOW_AAPL]CUT_SHELL=4

    This is also the default value. The cut-off radius is then determined automatically.

It is possible to use both parameters. In this case, AAPL determines the cut-off radius based on
the CUT_SHELL parameter and chooses whichever value is greater.

7.2.2. Sum Rule Parameters

The iterative procedure behind the calculation of the anharmonic IFCs has a few adjustable
parameters.

    * The most important parameter is the convergence criterion for the sum rule, which can
      be set using the SUMRULE parameter, e.g.:

          [AFLOW_AAPL]SUMRULE=1E-7

      The default value is 1E-5. It causes quick convergence, but may lead to underestimated
      thermal conductivity values.

    * The maximum number of iterations can be set using the parameter SUMRULE_MAX_ITER:

          [AFLOW_AAPL]SUMRULE_MAX_ITER=2000

      This is also the default value, which should be more than enough.

    * The last adjustable parameter is the mixing coefficient, which controls how much of the
      previous iteration step is mixed into the new IFC values:

          [AFLOW_AAPL]MIXING_COEFFICIENT=0.0

      It defaults to 0.0, i.e. there is no mixing between two iteration steps. However, if
      convergence cannot be reached because the corrections are too large, this parameter should
      be adjusted.

7.2.3. Four-Phonon Processes

Four-phonon processes play a critical role at high temperatures. Without including these processes,
thermal conductivities are often overestimated. However, they are also very expensive to compute.

To include four-phonon processes, the FOURTH_ORDER keyword needs to be set:

    [AFLOW_AAPL]FOURTH_ORDER=y

The default is FOURTH_ORDER=n, i.e. only three-phonon processes are considered.

CUT_SHELL and CUT_RAD can be set separately for third and fourth order clusters. To do this, the
flags need to be written as:

    [AFLOW_AAPL]CUT_RAD=5,3

If no second value is provided, AAPL will set the fourth order parameters to the third order
parameters and display a warning.

7.3. Lattice Thermal Conductivity Calculations

The lattice thermal conductivity tensor is calculated by solving the Boltzmann Transport
Equation (BTE) for a specific temperature range. To set the temperature range, the parameter TCT
needs to be set similar to TPT in APL (see Section 3.3.):

    [AFLOW_AAPL]TCT=50:550:10

The default setting is 50:550:50. The starting temperature cannot be set to 0 K because the thermal
conductivity is infinite. If it is set to 0 K, AAPL will simply skip it.

The BTE can be solved using either the Relaxation Time Approximation (RTA) or an iterative
procedure, which is done by the BTE keyword:

    [AFLOW_AAPL]BTE=RTA|FULL

The default is BTE=RTA.

The BTE gets solved along a q-point grid which can be defined using THERMALGRID:

    [AFLOW_AAPL]THERMALGRID=q1Xq2Xq3

Where q1, q2, and q3 are the number of q-points along each dimension. The default is 16X16X16.

There are various factors which can influence the phonon scattering rates besides anharmonic phonon
processes. Two mechanisms that are implemented in AAPL: isotope scattering and grain boundary
scattering.

    * Isotope scattering can be switched on by setting the parameter ISOTOPE:

          [AFLOW_AAPL]ISOTOPE=y

      The default is ISOTOPE=n.

    * Grain boundary scattering can be included in two ways:

        1. By introducing a grain boundary scattering term into the scattering rate. This is done
           by setting BOUNDARY to true:

              [AFLOW_AAPL]BOUNDARY=y

           The default is BOUNDARY=n.

        2. By only including phonon modes with a mean free path below a threshold into the thermal
           conductivity calculation. This is done using the keyword CUMULATIVEK:

               [AFLOW_AAPL]CUMULATIVEK=y

           The default is CUMULATIVEK=n.

        If both are set to true, BOUNDARY will be used, i.e. CUMULATIVEK will be set to false.

        Both methods require an additional parameter, NANO_SIZE:

            [AFLOW_AAPL]NANO_SIZE=300

        Which corresponds to the grain size/mean free path threshold in nm. The default is 100 nm.

7.4. Example AAPL Input File

The following is an example AAPL input file for NaF. Four-phonon processes are included.

[AFLOW] *******************************************************************************************
[AFLOW]System=NaF
[AFLOW] *******************************************************************************************
[AFLOW_MODE=VASP]
[AFLOW_MODE_ZIP=bzip2]
[AFLOW_MODE_BINARY=vasp46s]
[AFLOW_MODE_MPI_MODE]NCPUS=MAX
[AFLOW_MODE_MPI_MODE]COMMAND="mpirun -np"
[AFLOW_MODE_MPI_MODE]AUTOTUNE
[AFLOW_MODE_MPI_MODE]BINARY="mpivasp46s"
[AFLOW] *******************************************************************************************
[AFLOW_AAPL]CALC
[AFLOW_APPL]FOURTH_ORDER=y
[AFLOW_AAPL]CUT_SHELL=4,2
#[AFLOW_AAPL]CUT_RAD=5.0,4.0
[AFLOW_AAPL]SUMRULE=1E-7
[AFLOW_AAPL]SUMRULE_MAX_ITER=1E-7
[AFLOW_AAPL]MIXING_COEFFICIENT=1.0
[AFLOW_AAPL]THERMALGRID=21x21x21
[AFLOW_AAPL]ISOTOPE=Y
[AFLOW_AAPL]BTE=FULL
[AFLOW_AAPL]TCT=50:550:10
#[AFLOW_AAPL]CUMULATIVEK=y
#[AFLOW_AAPL]BOUNDARY=n
#[AFLOW_AAPL]NANO_SIZE=100
[AFLOW] *******************************************************************************************
[AFLOW_APL]ENGINE=DM
[AFLOW_APL]DMAG=0.015
#[AFLOW_APL]POLAR=Y
[AFLOW_APL]MINATOMS=100
[AFLOW_APL]HIBERNATE=Y
[AFLOW_APL]KPPRA=8000
[AFLOW_APL]KSCHEME=G
[AFLOW_APL]KPOINTS=EVEN
[AFLOW] *******************************************************************************************
[VASP_FORCE_OPTION]SYM=ON
[VASP_FORCE_OPTION]AUTO_PSEUDOPOTENTIALS=potpaw_PBE
[VASP_FORCE_OPTION]NBANDS
[VASP_FORCE_OPTION]SPIN=OFF
[VASP_FORCE_OPTION]TYPE=DEFAULT
[VASP_FORCE_OPTION]PREC=PHONON
[VASP_FORCE_OPTION]ALGO=NORMAL
[AFLOW] *******************************************************************************************
[VASP_POTCAR_MODE_IMPLICIT]
[VASP_POTCAR_FILE]F
[VASP_POTCAR_FILE]Na
[AFLOW] *******************************************************************************************
[VASP_INCAR_MODE_EXPLICIT]START
System = F1Na1_ICSD_262837
NELM = 120
NELMIN = 2
LPLANE= .TRUE.
LREAL = .FALSE.
LSCALU = .FALSE.
[VASP_INCAR_MODE_EXPLICIT]STOP
[AFLOW] *******************************************************************************************
[AFLOW_SYMMETRY]CALC
[VASP_POSCAR_MODE_EXPLICIT]START
F1Na1  [FCC,FCC,cF8] (STD_PRIM doi:10.1
1.224745
   0.00000000000000   1.91118771995223   1.91118771995223
   1.91118771995223   0.00000000000000   1.91118771995223
   1.91118771995223   1.91118771995223   0.00000000000000
F Na 
1 1 
Direct(2) [A1B1] 
   0.00000000000000   0.00000000000000   0.00000000000000  F     
   0.50000000000000   0.50000000000000   0.50000000000000  Na    
[VASP_POSCAR_MODE_EXPLICIT]STOP
[AFLOW] *******************************************************************************************

****************************************************************************************************
*                                                                                                  *
*                        aflow - STEFANO CURTAROLO Duke University 2003-2018                       *
*                            High-Throughput ab-initio Computing Project                           *
*                                                                                                  *
****************************************************************************************************
