****************************************************************************************************
*                                                                                                  *
*                        aflow - STEFANO CURTAROLO Duke University 2003-2018                       *
*                            High-Throughput ab-initio Computing Project                           *
*                                                                                                  *
****************************************************************************************************
LATEST VERSION OF THE FILE:                       materials.duke.edu/AFLOW/README_AFLOW_SYM.TXT
****************************************************************************************************

AFLOW-SYM README:

README written by: David Hicks (d.hicks@duke.edu)

Citation info: D. Hicks, C. Oses, E. Gossett, G. Gomez, R.H. Taylor 
               C. Toher, M.J. Mehl, O. Levy, and S. Curtarolo, 
               AFLOW-SYM: platform for the complete, automatic and self-consistent 
               symmetry analysis of crystals, Acta Cryst. A74, 184-203 (2018).


****************************************************************************************************
GENERAL OVERVIEW

AFLOW-SYM is a robust symmetry analysis tool catered to - but not limited to - 
high-throughput frameworks.  With built-in self-consistency symmetry validation and 
an adaptive tolerance scheme, AFLOW-SYM resolves a globally consistent symmetry 
profile.

The AFLOW-SYM functions fully characterize the symmetry of a crystal, including 
the space group, lattice types, Pearson symbol, Wyckoff positions, point group, 
factor group, symmetrically equivalent atoms, etc.  The AFLOW-SYM commands, along 
with important usage information, are listed below.

README Sections:
  - LIST OF AFLOW-SYM COMMANDS
  - SUPPORTED GEOMETRY FILE FORMATS
  - SELF-CONSISTENT SYMMETRY
  - TOLERANCE
  - MAGNETIC ANALYSIS
  - CELL/SETTING CHOICES
  - OUTPUT FORMAT
  - AFLOW COMMANDS USAGE


****************************************************************************************************
LIST OF AFLOW-SYM COMMANDS:

  aflow --aflow-sym            : Returns symmetry operations of the point groups (lattice, crystal, 
                                 reciprocal lattice, and dual of crystal), factor group, space group, 
                                 and site symmetry; along with symmetrically equivalent/unique atoms.
  aflow --aflowSG              : Returns space group in ITC convention
  aflow --cif                  : Returns structure in CIF file format
  aflow --edata                : Returns extended crystallographic symmetry data (edata) of the 
                                 crystal, lattice, reciprocal lattice, and superlattice
  aflow --equivalent           : Returns symmetrically equivalent/unique atoms 
  aflow --factorgroup          : Returns the coset representatives of the factor group (unit cell 
                                 symmetry operations) 
  aflow --lattice_type         : Returns the lattice type corresponding to the crystal (lattice + 
                                 atoms) symmetry
  aflow --lattice_lattice_type : Returns lattice type corresponding to the lattice symmetry
  aflow --pearson_symbol       : Returns Pearson symbol 
  aflow --pointgroup           : Returns point group operations of the lattice
  aflow --pointgroup_crystal   : Returns point group operations of the crystal
  aflow --pointgroupkcrystal   : Returns dual of crystal point group operations
  aflow --pointgroupklattice   : Returns point group operations of the reciprocal lattice 
  aflow --sc                   : Returns standard conventional cell in the AFLOW standard 
  aflow --sgdata               : Returns extended space group information (number, labels, classes, 
                                 settings, Wyckoff positions, etc.)
  aflow --sitepointgroup       : Returns site symmetry operations for atoms in the unit cell 
  aflow --sp                   : Returns standard primitive cell in the AFLOW standard 
  aflow --spacegroup           : Returns space group operations out to a certain radius 
  aflow --wyccar               : Returns Wyckoff positions in a POSCAR file format

For details on usage, see the "AFLOW COMMANDS USAGE" section below.


****************************************************************************************************
SUPPORTED GEOMETRY FILE FORMATS:

AFLOW-SYM reads the following geometry file types:

  - Crystallographic Input File (CIF)
  - VASP (POSCAR)
  - QuantumESPRESSO
  - ABINIT
  - FHI-AIMS

More file types will be added in later versions.


****************************************************************************************************
SELF-CONSISTENT SYMMETRY:

AFLOW-SYM includes self-consistency symmetry checks to identify a globally 
commensurate symmetry description.

The following commands provide the most extensive self-consistency checks and are 
recommended for finding the symmetry of a material:

  aflow --aflow-sym : Validates symmetry operations between point groups, factor group, 
                      space group, site symmetry, and equivalent/unique atoms 
  aflow --edata     : Validates symmetry between the space group, lattice types, Pearson 
                      symbol, etc.

The remaining commands are performed independently and are not validated with the full symmetry 
profile.


****************************************************************************************************
TOLERANCE:

AFLOW-SYM offers an adaptive tolerance scheme to combat issues arising from 
numerical instabilities.  If any of the symmetry analyses yield incommensurate 
symmetry, a radial tolerance scan is performed around the inital input tolerance. 
All symmetry properties are then recalculated from the new tolerance.

The initial tolerance can be specified by appending it to the command with an "=".

  Ex: aflow --aflow-sym=<tolerance>

  Preset AFLOW tolerances:
    tight : minimum_interatomic_distance/100.0 (default)
    loose : minimum_interatomic_distance/10.0

    Ex: aflow --aflow-sym=tight  OR  aflow --aflow-sym=loose

  User defined tolerances:
    Any value between 1e-10 and the minimum interatomic distance
    
    Ex: aflow --aflow-sym=1e-2

No tolerance scan option:
  The tolerance scan functionality can be suppresed by appending the option flag 
  [--no_scan] to the command.


****************************************************************************************************
MAGNETIC ANALYSIS:

The symmetry of magnetic systems is calculated by adding the magnetic moment 
information (symmetry of lattice, atoms, and spin). Works for both collinear 
and non-collinear spin.

The magnetic symmetry analysis is performed by specifying the magnetic moment(s)
via the optional [--magmom] flag.
 
Magnetic moment input formats:
   1) Explicitly via comma separated magnetic moments (m1, m2, ..., mn) in the same 
      order as the input file
        - Ex: aflow --aflow-sym --magmom=1.0,-1.0,1.0
   2) Read from the VASP INCAR by specifying path to file
        - Ex: aflow --aflow-sym --magmom=path_to_incar
   3) Read from the VASP OUTCAR by specifying path to file
        - Ex: aflow --aflow-sym --magmom=path_to_outcar
   4) Magnetic moment readers for other ab-initio codes will be added in later versions.

   AFLOW autodetects collinear or non-collinear based on the input.

Relevant for the following commands:

  aflow --aflow-sym
  aflow --aflowSG
  aflow --edata
  aflow --equivalent
  aflow --factorgroup
  aflow --pointgroup_crystal
  aflow --pointgroupkcrystal   
  aflow --sgdata              
  aflow --sitepointgroup     
  aflow --spacegroup        
  aflow --wyccar           

NOTE: This analysis does not account for time reversal symmetry (spin flips), which is 
      necessary to find the magnetic space group.  The magnetic symmetry analysis 
      here considers the symmetry-breaking feature of the spin-degree of freedom in 
      magnetic systems.
      The magnetic space groups will be added in a later version.

****************************************************************************************************
CELL/SETTING CHOICES:

Upon performing the symmetry analysis, AFLOW-SYM can represent the structure with respect to 
common cell/setting choices listed in the ITC with the option flag [--setting=1 | =2].
Different cell/setting choice exist for:
  - centrosymmetric space groups with two origin choices (24 space groups, e.g., sg=227): 
      --setting=1 : origin centered on high site symmetry (not inversion site)
      --setting=2 : origin centered on inversion site
  - rhombohedral space groups (e.g., space group 166): 
      --setting=1 : rhombohedral setting 
      --setting=2 : hexagonal setting
  - monoclinic space groups (e.g., space group 5): 
      --setting=1 : unique axis-b 
      --setting=2 : unique axis-c

More information on these setting choices is located in ITC Vol. A on pg. 17 (5th Edition).

Relevant for the following commands (i.e., when the conventional cell/Wyckoff positions is shown):
  aflow --cif
  aflow --edata
  aflow --sgdata
  aflow --wyccar


****************************************************************************************************
OUTPUT FORMAT:

AFLOW-SYM contains functionality to print the results in a standard text output 
and JavaScript Object Notation (JSON).

The output format can be specified by appending the optional flag [--print=<format>].

  Text output (default):
    Ex: aflow --edata --print=txt

  JSON output: 
    Ex: aflow --edata --print=json 


****************************************************************************************************
AFLOW COMMANDS USAGE:
  
  aflow --aflow-sym|--AFLOW-SYM|--AFLOWSYM|--aflowSYM|--aflowsym|--full_symmetry|--full_sym|--fullsym[=tolerance| =tight| =loose] [--no_scan] [--print=txt| =json] [--screen_only] [--mag|--magnetic|--magmom=[m1,m2,...|INCAR|OUTCAR]] < POSCAR
        Calculates the full symmetry of a crystal: point group lattice, point group crystal, factor group,
        space group, site point group, inequivalent/equivalent atoms, and the point group of the klattice.
        This routine calculates the full symmetry suite at a certain tolerance. If symmetry rules are broken,
        the code automatically changes tolerance and recalculates from the beginning.
        Tolerance:
           The tolerance is given in Angstroms.  There are two preset tolerances:
             tight: minimum_interatomic_distance/100.0
             loose: minimum_interatomic_distance/10.0
           The default tolerance is the "tight" tolerance value.
           A value can also be specified. (Note, it must not be below the minimum 
           interatomic distance, otherwise an error will be thrown.)
        Magnetic analysis:
          The magnetic tag specifies the magnetic moment on each of the atoms.  The
          two input formats are a comma separated string of the magnetic moments -
          in the same order as the input file, the location of the INCAR, or the location
          of the OUTCAR, which is then read and the magnetic moment for each atom is found
          (must be a static calculation, otherwise this information isn't specified).
          Works for both collinear and non-collinear systems (AFLOW autodetects
          collinear/non-collinear based on the input).
        Options:
          [--no_scan] : Will not perform tolerance scan.
          [--print=txt| =json] : Specify output format.
          [--screen_only] : Prints symmetry elements to string and does not write to a file.
          [--mag|--magnetic|--magmom=[m1,m2,...|INCAR|OUTCAR]]: Perform symmetry analysis of spin system, where spin degree of freedom breaks symmetry 

  aflow --aflowSG[_label,_number][=tolerance | =tight | =loose] [--no_scan] [--mag|--magnetic|--magmom=[m1,m2,...|INCAR|OUTCAR]] < POSCAR
        Calculates the space group of the crystal in the ITC convention.
        No assumptions are made on the input POSCAR-- may be supercell, primitive, etc.
        Gives the Bravais lattice, crystal system, point group, space group number and wyckoff set.
        The origin choice is always the first choice in ITC 5th edition. This is the highest symmetry
        choice. (However, for Monoclinic systems, different unique axes and cell choices are examined to
        ensure a spacegroup is identified.) This routine consists of an adaptive tolerance; if any symmetry 
        rules are broken, the tolerance is changed and everything is recalculated.
        Tolerance:
           The tolerance is given in Angstroms.  There are two preset tolerances:
             tight: minimum_interatomic_distance/100.0
             loose: minimum_interatomic_distance/10.0
           The default tolerance is the "tight" tolerance value.
           A value can also be specified. (Note, it must be below the minimum 
           interatomic distance, otherwise an error will be thrown.)
        Can be called inside aflow via:
              xstructure.SpaceGroup_ITC()        (default tolerance)
              xstructure.SpaceGroup_ITC(tol)     (user defined tolerance)
        AFLOW counterpart to platonSG and findsymSG.
        Magnetic analysis:
          The magnetic tag specifies the magnetic moment on each of the atoms.  The
          two input formats are a comma separated string of the magnetic moments -
          in the same order as the input file, the location of the INCAR, or the location 
          of the OUTCAR, which is then read and the magnetic moment for each atom is found 
          (must be a static calculation, otherwise this information isn't specified).
          Works for both collinear and non-collinear systems (AFLOW autodetects 
          collinear/non-collinear based on the input).
        Options:
          [--no_scan] : Will not perform tolerance scan.
          [--mag|--magnetic|--magmom=[m1,m2,...|INCAR|OUTCAR]]: Perform symmetry analysis of spin system, where spin degree of freedom breaks symmetry 
  
  aflow --cif[=tolerance | =tight | =loose] [--setting=1| =2] [--no_symmetry] < POSCAR
        Outputs to standard out a Crystallographic Information File
        (cif) format file based on the POSCAR input file.
        This can be used as input for Mercury, and is the preferred
        format in Gerd's structure collection.
        Tolerance:
           The tolerance is given in Angstroms.  There are two preset tolerances:
             tight: minimum_interatomic_distance/100.0
             loose: minimum_interatomic_distance/10.0
           The default tolerance is the "tight" tolerance value.
           A value can also be specified. (Note, it must not be below the minimum
           interatomic distance, otherwise an error will be thrown.)
        Space group setting:
          The space group setting can be specified with the optional flag [--setting=1| =2] (the default is setting 1).
          For rhombohedral systems, setting 1=rhl and 2=hex.
          For monoclinic systems, setting 1=unique axis-b and 2=unique axis-c.
        Options:
          [--no_symmetry] : Will not perform symmetry analysis and return the CIF with space group 1.
 
  aflow --edata[=tolerance| =tight| =loose] [--setting=1| =2] [--no_scan] [--print=txt| =json] < POSCAR
        Outputs to standard the information given by --data in addition
        to lattice-, superlattice- and reciprocal-lattice types.
        This routine consists of symmetry consistency checks. If symmetry rules are broken,
        the code automatically changes tolerance and recalculates from the beginning.
        Tolerance:
           The tolerance is given in Angstroms.  There are two preset tolerances:
             tight: minimum_interatomic_distance/100.0
             loose: minimum_interatomic_distance/10.0
           The default tolerance is the "tight" tolerance value.
           A value can also be specified. (Note, it must not be below the minimum
           interatomic distance, otherwise an error will be thrown.)
        Space group setting:
          The space group setting can be specified with the optional flag [--setting=1| =2] (the default is setting 1).
          For rhombohedral systems, setting 1=rhl and 2=hex.
          For monoclinic systems, setting 1=unique axis-b and 2=unique axis-c.
        Magnetic analysis:
          The magnetic tag specifies the magnetic moment on each of the atoms.  The
          two input formats are a comma separated string of the magnetic moments -
          in the same order as the input file, the location of the INCAR, or the location
          of the OUTCAR, which is then read and the magnetic moment for each atom is found
          (must be a static calculation, otherwise this information isn't specified).
          Works for both collinear and non-collinear systems (AFLOW autodetects
          collinear/non-collinear based on the input).
        Options:
          [--no_scan] : Will not perform tolerance scan.
          [--print=txt| =json] : Specify output format.
          [--mag|--magnetic|--magmom=[m1,m2,...|INCAR|OUTCAR]]: Perform symmetry analysis of spin system, where spin degree of freedom breaks symmetry 
  
  aflow --equivalent | --equiv | --iatoms[=tolerance| =tight| =loose] [--no_scan] [--print=txt| =json] [--mag|--magnetic|--magmom=[m1,m2,...|INCAR|OUTCAR]] < POSCAR
        Calculate point/factor/space group and use them to label equivalent
        and inequivalent atoms. On the output, each structure atom has
        number_label_of_atom[equivalent_to_label]* (* if inequivalent).
        (SC1107).
        This routine consists of symmetry consistency checks. If symmetry rules are broken,
        the code automatically changes tolerance and recalculates from the beginning.
        Tolerance:
           The tolerance is given in Angstroms.  There are two preset tolerances:
             tight: minimum_interatomic_distance/100.0
             loose: minimum_interatomic_distance/10.0
           The default tolerance is the "tight" tolerance value.
           A value can also be specified. (Note, it must not be below the minimum
           interatomic distance, otherwise an error will be thrown.)
        Magnetic analysis:
          The magnetic tag specifies the magnetic moment on each of the atoms.  The
          two input formats are a comma separated string of the magnetic moments -
          in the same order as the input file, the location of the INCAR, or the location
          of the OUTCAR, which is then read and the magnetic moment for each atom is found
          (must be a static calculation, otherwise this information isn't specified).
          Works for both collinear and non-collinear systems (AFLOW autodetects
          collinear/non-collinear based on the input).
        Options:
          [--no_scan] : Will not perform tolerance scan.
          [--print=txt| =json] : Specify output format.
          [--mag|--magnetic|--magmom=[m1,m2,...|INCAR|OUTCAR]]: Perform symmetry analysis of spin system, where spin degree of freedom breaks symmetry 
  
  aflow --factorgroup[=tolerance| =tight| =loose] [--no_scan] [--print=txt| =json] [--screen_only] [--mag|--magnetic|--magmom=[m1,m2,...|INCAR|OUTCAR]] < POSCAR
        Calculates factor group symmetry of the cell {R|t} and writes
        in the aflow.fgroup.out file. See documentation of aflow.
        The point group is required for the factor group, therefore
        the aflow.pgroup.out file will be generated as well.
        This routine consists of symmetry consistency checks. If symmetry rules are broken,
        the code automatically changes tolerance and recalculates from the beginning.
        Tolerance:
           The tolerance is given in Angstroms.  There are two preset tolerances:
             tight: minimum_interatomic_distance/100.0
             loose: minimum_interatomic_distance/10.0
           The default tolerance is the "tight" tolerance value.
           A value can also be specified. (Note, it must not be below the minimum
           interatomic distance, otherwise an error will be thrown.)
        Magnetic analysis:
          The magnetic tag specifies the magnetic moment on each of the atoms.  The
          two input formats are a comma separated string of the magnetic moments -
          in the same order as the input file, the location of the INCAR, or the location
          of the OUTCAR, which is then read and the magnetic moment for each atom is found
          (must be a static calculation, otherwise this information isn't specified).
          Works for both collinear and non-collinear systems (AFLOW autodetects
          collinear/non-collinear based on the input).
        Options:
          [--no_scan] : Will not perform tolerance scan.
          [--print=txt| =json] : Specify output format.
          [--screen_only] : Prints symmetry elements to string and does not write to a file.
          [--mag|--magnetic|--magmom=[m1,m2,...|INCAR|OUTCAR]]: Perform symmetry analysis of spin system, where spin degree of freedom breaks symmetry 
  
  aflow --lattice_type | --lattice | --lattice_crystal < POSCAR
        Returns the lattice type and the conventional lattice type of the CRYSTAL following the tables
        of Setyawayn-Curtarolo [https://doi.org/10.1016/j.commatsci.2010.05.010].

  aflow --lattice_lattice_type | --lattice_lattice < POSCAR
        Returns the lattice type and the conventional lattice type of the LATTICE following the tables
        of Setyawayn-Curtarolo [https://doi.org/10.1016/j.commatsci.2010.05.010].
  
  aflow --pearson_symbol | --pearson < POSCAR
        Returns the Pearson symbol of the structure.

  aflow --pointgroup | --pgroup[=tolerance| =tight| =loose] [--no_scan] [--print=txt| =json] [--screen_only] < POSCAR
        Calculates the point group symmetry of the lattice {R} and writes
        it in the aflow.pgroup.out file. See documentation of aflow.
        This routine consists of symmetry consistency checks. If symmetry rules are broken,
        the code automatically changes tolerance and recalculates from the beginning.
        Tolerance:
           The tolerance is given in Angstroms.  There are two preset tolerances:
             tight: minimum_interatomic_distance/100.0
             loose: minimum_interatomic_distance/10.0
           The default tolerance is the "tight" tolerance value.
           A value can also be specified. (Note, it must not be below the minimum
           interatomic distance, otherwise an error will be thrown.)
        Options:
          [--no_scan] : Will not perform tolerance scan.
          [--print=txt| =json] : Specify output format.
          [--screen_only] : Prints symmetry elements to string and does not write to a file.

  aflow --pointgroup_crystal | --pgroup_crystal | --pgroup_xtal[=tolerance| =tight| =loose] [--no_scan] [--print=txt| =json] [--screen_only] [--mag|--magnetic|--magmom=[m1,m2,...|INCAR|OUTCAR]] < POSCAR
        Calculate the point group symmetry of the crystal {R+atoms} and
        writes it the aflow.pgroup_xtal.file.  
        This routine consists of symmetry consistency checks. If symmetry rules are broken,
        the code automatically changes tolerance and recalculates from the beginning.
        Tolerance:
           The tolerance is given in Angstroms.  There are two preset tolerances:
             tight: minimum_interatomic_distance/100.0
             loose: minimum_interatomic_distance/10.0
           The default tolerance is the "tight" tolerance value.
           A value can also be specified. (Note, it must not be below the minimum
           interatomic distance, otherwise an error will be thrown.)
        Magnetic analysis:
          The magnetic tag specifies the magnetic moment on each of the atoms.  The
          two input formats are a comma separated string of the magnetic moments -
          in the same order as the input file, the location of the INCAR, or the location
          of the OUTCAR, which is then read and the magnetic moment for each atom is found
          (must be a static calculation, otherwise this information isn't specified).
          Works for both collinear and non-collinear systems (AFLOW autodetects
          collinear/non-collinear based on the input).
        Options:
          [--no_scan] : Will not perform tolerance scan.
          [--print=txt| =json] : Specify output format.
          [--screen_only] : Prints symmetry elements to string and does not write to a file.
          [--mag|--magnetic|--magmom=[m1,m2,...|INCAR|OUTCAR]]: Perform symmetry analysis of spin system, where spin degree of freedom breaks symmetry 

  aflow --pointgroupkcrystal | --pointgroupk_xtal[=tolerance| =tight| =loose] [--no_scan] [--print=txt| =json] [--screen_only] [--mag|--magnetic|--magmom=[m1,m2,...|INCAR|OUTCAR] < POSCAR
        Calculates the dual of the crystal point group symmetry operations and writes
        it in the aflow.pgroupk_xtal.out file. See documentation of aflow.
        This routine consists of symmetry consistency checks. If symmetry rules are broken,
        the code automatically changes tolerance and recalculates from the beginning.
        Tolerance:
           The tolerance is given in Angstroms.  There are two preset tolerances:
             tight: minimum_interatomic_distance/100.0
             loose: minimum_interatomic_distance/10.0
           The default tolerance is the "tight" tolerance value.
           A value can also be specified. (Note, it must not be below the minimum
           interatomic distance, otherwise an error will be thrown.)
        Magnetic analysis:
          The magnetic tag specifies the magnetic moment on each of the atoms.  The
          two input formats are a comma separated string of the magnetic moments -
          in the same order as the input file, the location of the INCAR, or the location
          of the OUTCAR, which is then read and the magnetic moment for each atom is found
          (must be a static calculation, otherwise this information isn't specified).
          Works for both collinear and non-collinear systems (AFLOW autodetects
          collinear/non-collinear based on the input).
        Options:
          [--no_scan] : Will not perform tolerance scan.
          [--print=txt| =json] : Specify output format.
          [--screen_only] : Prints symmetry elements to string and does not write to a file.
          [--mag|--magnetic|--magmom=[m1,m2,...|INCAR|OUTCAR]]: Perform symmetry analysis of spin system, where spin degree of freedom breaks symmetry 

  aflow --pointgroupklattice | --pgroupk[=tolerance| =tight| =loose] [--no_scan] [--print=txt| =json] [--screen_only] < POSCAR
        Calculates the point group symmetry of the klattice {K} and writes
        it in the aflow.pgroupk.out file. See documentation of aflow.
        This routine consists of symmetry consistency checks. If symmetry rules are broken,
        the code automatically changes tolerance and recalculates from the beginning.
        Tolerance:
           The tolerance is given in Angstroms.  There are two preset tolerances:
             tight: minimum_interatomic_distance/100.0
             loose: minimum_interatomic_distance/10.0
           The default tolerance is the "tight" tolerance value.
           A value can also be specified. (Note, it must not be below the minimum
           interatomic distance, otherwise an error will be thrown.)
        Options:
          [--no_scan] : Will not perform tolerance scan.
          [--print=txt| =json] : Specify output format.
          [--screen_only] : Prints symmetry elements to string and does not write to a file.
 
  aflow --sc | --std_conv | --standard_conventional | --sconv < POSCAR
    	  Output POSCAR in a standard conventional lattice
	  (use --sp | --std_prim | --standard_primitive | --sprim	
	  to output POSCAR in a standard primitive lattice) (WS & SC Nov09)
        REF: Setyawan Curtarolo, DOI: 10.1016/j.commatsci.2010.05.010

  aflow --sgdata | --space_group_data[=tolerance| =tight| =loose] [--setting=1| =2] [--no_scan] [--print=txt| =json] [--screen_only] [--mag|--magnetic|--magmom=[m1,m2,...|INCAR|OUTCAR]] < POSCAR
        Outputs the following space group information: space group number, space group labels
        (Hermann Mauguin, Hall, Schoenflies), Laue class, crystal class, ITC setting, ITC origin shift,
        general Wyckoff position, and the WYCCAR (similar to --edata).
        This routine consists of symmetry consistency checks. If symmetry rules are broken,
        the code automatically changes tolerance and recalculates from the beginning.
        Tolerance:
           The tolerance is given in Angstroms.  There are two preset tolerances:
             tight: minimum_interatomic_distance/100.0
             loose: minimum_interatomic_distance/10.0
           The default tolerance is the "tight" tolerance value.
           A value can also be specified. (Note, it must not be below the minimum
           interatomic distance, otherwise an error will be thrown.)
        Space group setting:
          The space group setting can be specified with the optional flag [--setting=1| =2] (the default is setting 1).
          For rhombohedral systems, setting 1=rhl and 2=hex.
          For monoclinic systems, setting 1=unique axis-b and 2=unique axis-c.
        Magnetic analysis:
          The magnetic tag specifies the magnetic moment on each of the atoms.  The
          two input formats are a comma separated string of the magnetic moments -
          in the same order as the input file, the location of the INCAR, or the location
          of the OUTCAR, which is then read and the magnetic moment for each atom is found
          (must be a static calculation, otherwise this information isn't specified).
          Works for both collinear and non-collinear systems (AFLOW autodetects
          collinear/non-collinear based on the input).
        Options:
          [--no_scan] : Will not perform tolerance scan.
          [--print=txt| =json] : Specify output format.

  aflow --sitepointgroup | --agroup[=tolerance| =tight| =loose] [--no_scan] [--print=txt| =json] [--screen_only] [--mag|--magnetic|--magmom=[m1,m2,...|INCAR|OUTCAR]] < POSCAR
        Calculates the site point group symmetry for every atom in the
        unit cell and writes it in the aflow.agroup.out file.
        This routine consists of symmetry consistency checks. If symmetry rules are broken,
        the code automatically changes tolerance and recalculates from the beginning.
        Tolerance:
           The tolerance is given in Angstroms.  There are two preset tolerances:
             tight: minimum_interatomic_distance/100.0
             loose: minimum_interatomic_distance/10.0
           The default tolerance is the "tight" tolerance value.
           A value can also be specified. (Note, it must not be below the minimum
           interatomic distance, otherwise an error will be thrown.)
        Magnetic analysis:
          The magnetic tag specifies the magnetic moment on each of the atoms.  The
          two input formats are a comma separated string of the magnetic moments -
          in the same order as the input file, the location of the INCAR, or the location
          of the OUTCAR, which is then read and the magnetic moment for each atom is found
          (must be a static calculation, otherwise this information isn't specified).
          Works for both collinear and non-collinear systems (AFLOW autodetects
          collinear/non-collinear based on the input).
        Options:
          [--no_scan] : Will not perform tolerance scan.
          [--print=txt| =json] : Specify output format.
          [--screen_only] : Prints symmetry elements to string and does not write to a file.
          [--mag|--magnetic|--magmom=[m1,m2,...|INCAR|OUTCAR]]: Perform symmetry analysis of spin system, where spin degree of freedom breaks symmetry 
	   See documentation of aflow.
  
  aflow --sp | --std_prim | --standard_primitive | --sprim < POSCAR
   	  Output POSCAR in a standard primitive lattice
	  (use --sc | --std_conv | --standard_conventional | --sconv	
	  to output POSCAR in a standard conventional lattice) (WS & SC Nov09)
        REF: Setyawan Curtarolo, DOI: 10.1016/j.commatsci.2010.05.010
  
  aflow --spacegroup [=tolerance| =tight| =loose] [--no_scan] [--print=txt| =json] [--screen_only] [--radius] [--mag|--magnetic|--magmom=[m1,m2,...|INCAR|OUTCAR]]< POSCAR
        Calculates space group symmetry of the cell {R|t+T} with
        translations as big as |T| and writes it in the aflow.sgroup.out
        file. See documentation of aflow.
        Be careful because the size of the space group increases
        as the radius^3 times the size of the factor roup.
        The point and factopr groups are required for the space group,
        therefore the aflow.pgroup.out and aflow.fgroup.out files will
        be generated as well.
        This routine consists of symmetry consistency checks. If symmetry rules are broken,
        the code automatically changes tolerance and recalculates from the beginning.
        Tolerance:
           The tolerance is given in Angstroms.  There are two preset tolerances:
             tight: minimum_interatomic_distance/100.0
             loose: minimum_interatomic_distance/10.0
           The default tolerance is the "tight" tolerance value.
           A value can also be specified. (Note, it must not be below the minimum
           interatomic distance, otherwise an error will be thrown.)
        Magnetic analysis:
          The magnetic tag specifies the magnetic moment on each of the atoms.  The
          two input formats are a comma separated string of the magnetic moments -
          in the same order as the input file, the location of the INCAR, or the location
          of the OUTCAR, which is then read and the magnetic moment for each atom is found
          (must be a static calculation, otherwise this information isn't specified).
          Works for both collinear and non-collinear systems (AFLOW autodetects
          collinear/non-collinear based on the input).
        Options:
          [--no_scan] : Will not perform tolerance scan.
          [--print=txt| =json] : Specify output format.
          [--screen_only] : Prints symmetry elements to string and does not write to a file.
          [--mag|--magnetic|--magmom=[m1,m2,...|INCAR|OUTCAR]]: Perform symmetry analysis of spin system, where spin degree of freedom breaks symmetry 

  aflow --wyccar[=tolerance | =tight | =loose] [--no_scan] [--setting=1| =2] [--mag|--magnetic|--magmom=[m1,m2,...|INCAR|OUTCAR]] < POSCAR
        Prints the wyccar file (the poscar file with wyckoff positions) for the input POSCAR
        The wyccar format presents the lattice in standard conventional form and the wyckoff
        positions with their associated labels and site symmetry. If variability in the
        Wyckoff positions exists, the algorithm finds the Wykoff positions with the "smallest"
        lettering scheme ("a" being the smallest). (R. Taylor/ D. Hicks)
        Tolerance:
           The tolerance is given in Angstroms.  There are two preset tolerances:
             tight: minimum_interatomic_distance/100.0
             loose: minimum_interatomic_distance/10.0
           The default tolerance is the "tight" tolerance value.
           A value can also be specified. (Note, it must not be below the minimum
           interatomic distance, otherwise an error will be thrown.)
        Space group setting:
          The space group setting can be specified with the optional flag [--setting=1| =2] (the default is setting 1).
          For rhombohedral systems, setting 1=rhl and 2=hex.
          For monoclinic systems, setting 1=unique axis-b and 2=unique axis-c.
        Magnetic analysis:
          The magnetic tag specifies the magnetic moment on each of the atoms.  The
          two input formats are a comma separated string of the magnetic moments -
          in the same order as the input file, the location of the INCAR, or the location
          of the OUTCAR, which is then read and the magnetic moment for each atom is found
          (must be a static calculation, otherwise this information isn't specified).
          Works for both collinear and non-collinear systems (AFLOW autodetects
          collinear/non-collinear based on the input).
        Options:
          [--no_scan] : Will not perform tolerance scan.
          [--mag|--magnetic|--magmom=[m1,m2,...|INCAR|OUTCAR]]: Perform symmetry analysis of spin system, where spin degree of freedom breaks symmetry
        Check aflow --readme=symmetry


For additional information contact: David Hicks (d.hicks@duke.edu)

****************************************************************************************************
*                                                                                                  *
*                        aflow - STEFANO CURTAROLO Duke University 2003-2018                       *
*                            High-Throughput ab-initio Computing Project                           *
*                                                                                                  *
****************************************************************************************************
