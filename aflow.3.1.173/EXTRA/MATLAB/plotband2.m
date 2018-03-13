%script to plot band structures
%the following files must reside in the directory
%  DOSCAR.static
%  POSCAR.bands
%  EIGENVAL.bands
%  KPOINTS.bands
%set the ddir variable below to set the working directory
%
%written by wahyu@alumni.duke.edu
clear;
%-----TO SET-------------------------
ddir='MATLAB_FUNCS_DIRECTORY'; %working directory
ytop=12;   %max energy in band structure plot
ybot=-10;  %min energy in band structure plot
mode='c';  %'c'=cartesian (unit of kpoint segment is absolute)
           %'d'=direct (unit of kpoint segment is fractional of k-lattices)
%---------------------------------------
klabel_tri=['\Gamma-X-V-Y-\Gamma-Z-U-R-T-_{Z|X}-_{U|Y}-_{T|V}-R'];
klabel_mcl=['\Gamma-B-A-Y-\Gamma-Z-D-E-C-_{Z|B}-_{D|Y}-_{C|A}-E'];
klabel_mclc=['\Gamma-A-L-V-\Gamma-Y-M-_{A|V}-_{Y|L}-M'];
klabel_orc=['\Gamma-X-S-Y-\Gamma-Z-U-R-T-_{Z|X}-_{U|Y}-_{T|S}-R'];
klabel_orcc=['\Gamma-Y-T-Z-\Gamma-S-R-Z'];
klabel_orcf=['\Gamma-Y-T-Z-\Gamma-L-_{Y|Z}-L-T'];
klabel_orci=['\Gamma-X-T-\Gamma-R-W-S-\Gamma-W-T'];
klabel_tet=['\Gamma-X-M-\Gamma-Z-R-A-_{Z|X}-_{R|M}-A'];
klabel_bct=['\Gamma-X-M-\Gamma-N-P-_{X|\Gamma}-P'];
klabel_rhl=['\Gamma-M-K-\Gamma-A-L-H-_{A|M}-_{L|K}-H'];
klabel_hex=['\Gamma-M-K-\Gamma-A-L-H-_{A|M}-_{L|K}-H'];
klabel_cub=['\Gamma-X-M-\Gamma-R-_{X|M}-R'];
klabel_fcc=['\Gamma-X-W-K-\Gamma-L-U-_{X|U}-W-L-_{K|W}-\Gamma'];
klabel_bcc=['\Gamma-H-N-\Gamma-P-_{H|P}-N'];
%---------------------------------------
% POSCAR.bands
%---------------------------------------
%calculating recip. lattice vectors b1 b2 b3 from POSCAR
pin=fopen([ddir 'POSCAR.bands'],'r');tmp=fgetl(pin);
LattConst=fscanf(pin,'%f',1); tmp=fgetl(pin);
a1(1:3)=fscanf(pin,'%f',3); tmp=fgetl(pin);
a2(1:3)=fscanf(pin,'%f',3); tmp=fgetl(pin);
a3(1:3)=fscanf(pin,'%f',3); tmp=fgetl(pin);
a1=a1*LattConst; a2=a2*LattConst; a3=a3*LattConst;
VCell=dot(cross(a2,a3),a1);
b1=2*pi*cross(a2,a3)/VCell; b2=2*pi*cross(a3,a1)/VCell; b3=2*pi*cross(a1,a2)/VCell;
%getting species and the number of species
tmp=fgetl(pin);
[species,Nspec]=sscanf(tmp,'%d');
fclose(pin);
%---------------------------------------
% DOSCAR.static
%---------------------------------------
%getting compound name
din=fopen([ddir 'DOSCAR.static'],'r');
tmp=fgetl(din); tmp=fgetl(din); tmp=fgetl(din); tmp=fgetl(din);
name2=fscanf(din,'%s',1); tmp=fgetl(din);
ii=0;name=name2;
funderscore=false;
for i=1:length(name2)
    if(name2(i)>'0' && name2(i)<'9' && ~funderscore)
        ii=ii+1; name(ii)='_';
        ii=ii+1; name(ii)=name2(i);
    elseif(name2(i)=='_')
        funderscore=true;
        ii=ii+1; name(ii)=' ';
    else
        ii=ii+1; name(ii)=name2(i);
    end    
end
%adjust EF to be at midway in the gap
Emax=fscanf(din,'%f',1);Emin=fscanf(din,'%f',1);
NgridDOS=fscanf(din,'%d',1);EF=fscanf(din,'%f',1); tmp=fgetl(din);
Egrid=(Emax-Emin)/NgridDOS;
tol=5e-3;
delE=1e6;delEE=delE;
E=zeros(NgridDOS,1);
D=zeros(NgridDOS,1);
tmp=fgetl(din);
[tmpV,count]=sscanf(tmp,'%f'); %getting ispin
ispin=(count-1)/2;             %getting ispin
isMetal=false;
iclosest=1;
E(1)=tmpV(1);  D(1)=tmpV(2);
for i=2:NgridDOS
    E(i)=fscanf(din,'%f',1); D(i)=fscanf(din,'%f',1); tmp=fgetl(din);
end
tmp=fgetl(din);tmp=fgetl(din); %getting Norb
[tmpV,count]=sscanf(tmp,'%f'); %getting Norb
Norb=(count-1)/ispin;          %getting Norb
switch Norb
    case 1
        orbital='s';
    case 2
        orbital='s p';
    case 3
        orbital='s p d';
    case 3
        orbital='s p d f';
end
for i=1:NgridDOS
    if(abs(EF-E(i))<Egrid)
        if(D(i)>tol && D(i+1)>tol) isMetal=true; break; %metal
        end
    end
end
for i=1:NgridDOS
    Etmp=E(i); Dtmp=D(i);
    if(Dtmp<tol) delEE=abs(EF-Etmp);
    end
    if(delEE<delE) iclosest=i; delE=delEE;
    end
end
for i=iclosest-1:-1:1
  if(D(i)>tol) break;
  end
end
Ev=E(i);
for i=iclosest+1:NgridDOS
  if(D(i)>tol) break;
  end
end
Ec=E(i);
if(isMetal)
    Ec=EF; Ev=EF;
end
EF=(Ev+Ec)/2;
fclose(din);
%---------------------------------------
% KPOINTS.bands
%---------------------------------------
kin=fopen([ddir 'KPOINTS.bands'],'r');
latt_type=fscanf(kin,'%s',1); tmp=fgetl(kin);
%getting Ngrid per segment in band structure
Ngrid=fscanf(kin,'%d',1); tmp=fgetl(kin);
%-----block to make KBAND similar to KBAND.out------
  NSEGMAX=50;
  ksegdirect=zeros(NSEGMAX,4);
  tmpstr=fgetl(kin);
  tmpstr=fgetl(kin);
  i=0;
  while(~feof(kin))
    tmpstr=fgetl(kin);tmpstr2=tmpstr;
    [ftmp,count]=sscanf(tmpstr,'%s');
    if(count>3)
      i=i+1;
      ftmp=sscanf(tmpstr2,'%f',3);
      ksegdirect(i,1)=ftmp(1);
      ksegdirect(i,2)=ftmp(2);
      ksegdirect(i,3)=ftmp(3);
    end
  end
  Nseg=i/2;
  %constructing kline from the norm of kpoints cascaded for continuous plot
  kline=zeros(Nseg*Ngrid,1);
  koffset=0; ind=1;
  for iseg=1:Nseg
    i = 2*iseg-1;
    kx1=ksegdirect(i,1);    ky1=ksegdirect(i,2);    kz1=ksegdirect(i,3);
    kx2=ksegdirect(i+1,1);  ky2=ksegdirect(i+1,2);  kz2=ksegdirect(i+1,3);
    dkx = kx2-kx1;
    dky = ky2-ky1;
    dkz = kz2-kz1;
    dk = (sqrt(dkx*dkx + dky*dky + dkz*dkz))/(Ngrid-1);
    kline(ind)=koffset; ind=ind+1;
    for j=2:Ngrid
      kline(ind) = koffset+(j-1)*dk; ind=ind+1;
    end
    koffset = koffset+dk*(Ngrid-1);
  end
  fclose(kin);
%  //creating overall kband data
%  //format:
%  //             spinup                 spindown
%  //kline[1] band1 band2 band3 ...  band1 band2 band3 ...
%  //kline[2] band1 band2 band3 ...  band1 band2 band3 ...
%  //...
%  //this format is nice for plotting with kpoints as the x-axis
%-------------------------------
% EIGENVAL.bands
%-------------------------------
  ein=fopen('EIGENVAL.bands','r');
  tmpstr=fgetl(ein);
  tmpstr=fgetl(ein);
  tmpstr=fgetl(ein);
  tmpstr=fgetl(ein);
  tmpstr=fgetl(ein);
  itmp=fscanf(ein,'%d',1); itmp=fscanf(ein,'%d',1); 
  %if(Nseg~=itmp/Ngrid) {cerr<<"error Nseg inconsistent between KPOINTS and EIGENVAL"<<endl;; return;}
  Nband=fscanf(ein,'%d',1);  tmpstr=fgetl(ein);
  KBAND=zeros(Nseg*Ngrid,ispin*Nband+1);
  for i=1:Nseg*Ngrid
    ftmp=fscanf(ein,'%f',1); ftmp=fscanf(ein,'%f',1);
    ftmp=fscanf(ein,'%f',1); ftmp=fscanf(ein,'%f',1);
    KBAND(i,1)=kline(i);
    for j=1:Nband
      ftmp=fscanf(ein,'%f',1);  
      KBAND(i,j+1)=fscanf(ein,'%f',1);
      if(ispin==2) KBAND(i,Nband+j+1)=fscanf(ein,'%f',1);
      end
    end
  end
  fclose(ein);
%------------------------------
%finding max valence band
klinedirect=KBAND(:,1);
nk=size(KBAND,1);
maxEv=EF-1000;
for i=2:Nband
  testband=KBAND(:,i);
  maxtest=max(testband);
  if(maxtest>maxEv && maxtest<EF) maxEv=maxtest; end
end
%---------------------------------------
%calculating the boundary k(s) for each segment in cartesian coordinate
for iseg=1:Nseg*2
  dummy=ksegdirect(iseg,1:3);
  dummy=dummy(1)*b1+dummy(2)*b2+dummy(3)*b3;
  ksegcart(iseg,1:3)=dummy;
end
klinecart=zeros(Nseg*Ngrid,1);
kposxdirect=zeros(Nseg,1);
kposxcart=zeros(Nseg,1);
tnormcart=norm(ksegcart(1,1:3));
for iseg=1:Nseg
  ka=ksegdirect(2*iseg-1,1:3);     kb=ksegdirect(2*iseg,1:3);
  kposxdirect(iseg+1)=kposxdirect(iseg)+norm(kb-ka);
  ka=ksegcart(2*iseg-1,1:3);     kb=ksegcart(2*iseg,1:3);
  kposxcart(iseg+1)=kposxcart(iseg)+norm(kb-ka);
  klinecart((iseg-1)*Ngrid+1:iseg*Ngrid)=kposxcart(iseg)+linspace(0,norm(kb-ka),Ngrid);
end
k=klinecart;
kposx=kposxcart;
if(mode=='d')
  k=klinedirect;
  kposx=kposxdirect;
end
%---------------------------------------
% PLOT   PLOT   PLOT
%---------------------------------------
figure(1);close(1);figure(1);
axorig=[0.13 0.11 0.775 0.815];
ax1=axes('position',[0.1 0.1 0.59 0.8]);
set(gca,'fontsize',16);
set(gca,'linewidth',0.8);
for i=2:Nband
	hold on; plot(k,KBAND(:,i)-maxEv,'k','linewidth',0.8);
end
drawnow;
axis([0 max(k) ybot ytop]);
box on;
ylabel('E (eV)','fontsize',18);
switch(latt_type)
    case 'TRI'
        klabel=klabel_tri;
    case 'MCL'
        klabel=klabel_mcl;
    case 'MCLC'
        klabel=klabel_mclc;
    case 'ORC'
        klabel=klabel_orc;
    case 'ORCC'
        klabel=klabel_orcc;
    case 'ORCF'
        klabel=klabel_orcf;
    case 'ORCI'
        klabel=klabel_orci;
    case 'TET'
        klabel=klabel_tet;
    case 'BCT'
        klabel=klabel_bct;
    case 'RHL'
        klabel=klabel_rhl;
    case 'HEX'
        klabel=klabel_hex;
    case 'CUB'
        klabel=klabel_cub;
    case 'FCC'
        klabel=klabel_fcc;
    case 'BCC'
        klabel=klabel_bcc;
end
iLbot=1; iLtop=1;
for i=1:Nseg+1
    if i==1
        iLbot=0;
    end
    fLtop=0;
    for iL=iLbot+1:length(klabel)
        if(klabel(iL)=='-')
            fLtop=1;
            iLtop=iL;
            break;
        end
    end
    if i==Nseg+1
        iLtop=length(klabel)+1;
    end
    kklabel=klabel(iLbot+1:iLtop-1);
    iLbot=iLtop;
    text(kposx(i)-(max(k)-min(k))*0.015,ybot-(ytop-ybot)*0.035,kklabel,'fontsize',18);
end
title([name ' (' latt_type ')'],'fontsize',18);
set(gca,'xtick',[-1 10]);
for i=2:length(kposx)-1
    hold on; plot([kposx(i) kposx(i)],[ybot ytop],'k--','linewidth',0.8);
end
%------------------------------------
% EDOS.out
%------------------------------------
DO_S=false;
DO_P=false;
DO_D=false;
DO_F=false;
if(Norb>=1) DO_S=true; end
if(Norb>=2) DO_P=true; end
if(Norb>=3) DO_D=true; end
if(Norb>=4) DO_F=true; end
%//processing DOSCAR
din=fopen('DOSCAR.static','r');
tmpstr=fgetl(din); tmpstr=fgetl(din); tmpstr=fgetl(din);
tmpstr=fgetl(din); tmpstr=fgetl(din); 
maxE=fscanf(din,'%f',1);
minE=fscanf(din,'%f',1);
nE=fscanf(din,'%d',1); tmpstr=fgetl(din);  
%//energy loop for DOS
%//only DOS will be considered and written out because the Integrated DOS
%//can be calculated from DOS and energy bins
Ncol=1;
if(DO_S) Ncol=1+1*Nspec; end
if(DO_P) Ncol=1+2*Nspec; end
if(DO_D) Ncol=1+3*Nspec; end
if(DO_F) Ncol=1+4*Nspec; end
if(ispin) Ncol=Ncol*2; end
Ncol=Ncol+1;
EDOS=zeros(nE,Ncol);
for i=1:nE
   EDOS(i,1)=fscanf(din,'%f',1);
   EDOS(i,2)=fscanf(din,'%f',1);
   if(ispin==2) EDOS(i,3)=-fscanf(din,'%f',1); end
   tmpstr=fgetl(din);
end
%//energy loop for DOS for s,p,d,f
%//sum over all atoms for the same SPECIES! 
%//note that we read up to the highest between s,p,d,f
%//even though not all of them will be outputed
%//We will output only the orbitals that are requested
%//at the prompt input. 
ioffset=0;
for i=1:Nspec
  ioffset=(i-1)*Norb*ispin;
  for j=1:species(i)
    tmpstr=fgetl(din);
    for k=1:nE
      ftmp=fscanf(din,'%f',1);%discard energy
      if(ispin==1)
        for iorb=1:Norb
          ftmp=fscanf(din,'%f',1); 
          EDOS(k,ioffset+iorb+2)=EDOS(k,ioffset+iorb+2)+ftmp;%//s
        end
        tmpstr=fgetl(din);
      end
      if(ispin==2)
        for iorb=1:Norb
          ftmp=fscanf(din,'%f',1);
          EDOS(k,ioffset+(iorb-1)*2+4)=EDOS(k,ioffset+(iorb-1)*2+4)+ftmp;%//up
          ftmp=fscanf(din,'%f',1);
          EDOS(k,ioffset+(iorb-1)*2+5)=EDOS(k,ioffset+(iorb-1)*2+5)-ftmp;%//down
        end
        tmpstr=fgetl(din);
      end
    end
  end
end
fclose(din);
%------------------------------------
%finding the appropriate scaling for DOS plotting
istop=size(EDOS,2);
xdosmax=0;
for i=1:nE
    if(EDOS(i,1)>(ybot-EF) && EDOS(i,1)<(ytop-EF))
        istart=3;
        if(ispin==2) istart=istart+1; end
        linetmp=abs(EDOS(i,istart:istop));
        maxtmp=max(linetmp);
        if(maxtmp>xdosmax) xdosmax=maxtmp; end
    end
end
if(ispin==2) xdosmax=xdosmax*2;
end
%for i=1:nE
%    if(EDOS(i,1)>maxEv)
%        istart=3;
%        if(ispin==2) istart=istart+1; end
%        linetmp=abs(EDOS(i,istart:istop));
%        maxtmp=max(linetmp);
%        if(maxtmp>xdosmax) xdosmax=maxtmp; end
%    end
%end
%xdosmax=xdosmax*3;
%------------------------------------
DOSup=EDOS(:,2);
DOSdw=zeros(nE,1);
ic=2;
if(ispin==2) ic=ic+1; DOSdw=EDOS(:,ic); end
%----------S-------
ic=ic+1;
Sup=0;Sdown=0;
if(ispin==1)
  for ispec=1:Nspec 
    Sup=Sup+EDOS(:,ic+(ispec-1)); 
  end
end
if(ispin==2)
  for ispec=1:Nspec
    Sup=Sup+EDOS(:,ic+(ispec-1)*2);
    Sdown=Sdown+EDOS(:,ic+(ispec-1)*2+1);
  end
end
ic=ic+Nspec*(ispin);
%---------P--------
Pup=0; Pdown=0;
if(Norb>1)
  if(ispin==1)
    for ispec=1:Nspec
      Pup=Pup+EDOS(:,ic+(ispec-1));
    end
  end
  if(ispin==2)
    for ispec=1:Nspec
      Pup=Pup+EDOS(:,ic+(ispec-1)*2);
      Pdown=Pdown+EDOS(:,ic+(ispec-1)*2+1);
    end
  end
  ic=ic+Nspec*(ispin);
end
%----------D---------
Dup=0; Ddown=0;
if(Norb>2)
  if(ispin==1)
    for ispec=1:Nspec
      Dup=Dup+EDOS(:,ic+(ispec-1));
    end
  end
  if(ispin==2)
    for ispec=1:Nspec
      Dup=Dup+EDOS(:,ic+(ispec-1)*2);
      Ddown=Ddown+EDOS(:,ic+(ispec-1)*2+1);
    end
  end
  ic=ic+Nspec*(ispin);
end
%----------F---------
Fup=0; Fdown=0;
if(Norb==4)
  if(ispin==1)
  for ispec=1:Nspec
    Fup=Fup+EDOS(:,ic+(ispec-1));
    end
  end
  if(ispin==2)
    for ispec=1:Nspec
      Fup=Fup+EDOS(:,ic+(ispec-1)*2);
      Fdown=Fdown+EDOS(:,ic+(ispec-1)*2+1);
    end
  end
  ic=ic+Nspec*(ispin);
end

ax2=axes('position',[0.7 0.1 0.2 0.8]);
set(gca,'fontsize',16);
set(gca,'linewidth',0.8);
%%plot(DOSup-DOSdw-EDOS(:,1),EDOS(:,1),'k','linewidth',0.8);
%plot(DOSup-DOSdw,EDOS(:,1)-maxEv,'k','linewidth',0.8);
%axis([0 25 ybot ytop]);
%set(gca,'ytick',[ybot-10 ytop+10]);
%set(gca,'xtick',[-100 100]);
%box on;
%xlabel('N(E)','fontsize',18);
%print -depsc fig.eps
%print -dtiff fig.tiff

plot(Sup-Sdown,EDOS(:,1)-maxEv,'k:','linewidth',1.2);
if(Norb>1)
  hold on;plot(Pup-Pdown,EDOS(:,1)-maxEv,'k-.','linewidth',0.8);
end
if(Norb>2)
  hold on;plot(Dup-Ddown,EDOS(:,1)-maxEv,'k--','linewidth',0.8);
end
if(Norb==4)
  hold on;plot(Fup-Fdown,EDOS(:,1)-maxEv,'k-','linewidth',0.8);
end
drawnow;

if(Norb==2)
  legend('s','p');
end
if(Norb==3);
  legend('s','p','d');%,'location','southeast');
end
if(Norb==4)
  legend('s','p','d','f');
end

axis([0 xdosmax ybot ytop]);
set(gca,'ytick',[ybot-10 ytop+10]);
set(gca,'xtick',[-300 300]);
box on;
xlabel('N(E)','fontsize',18);

eval(['print -deps ' 'figband_' name2 '.eps']);
system(['convert ' 'figband_' name2 '.eps ' ' ' name2 '.jpg']);
system(['convert ' 'figband_' name2 '.eps ' ' ' name2 '.pdf']);

