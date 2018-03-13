%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elems={'Y','Sc','La','Zr','Hf','Ti','Nb','Ta','V','Mo','W','Cr','Re','Mn','Tc','Fe','Os','Ru','Co','Ir','Rh','Ni','Pt','Pd','Au','Ag','Cu','Hg','Cd','Zn'};
%      Y   Sc   La   Zr   Hf   Ti   Nb   Ta    V   Mo    W   Cr   Re   Mn   Tc   Fe   Os   Ru   Co   Ir   Rh   Ni   Pt   Pd   Au   Ag   Cu   Hg   Cd   Zn
%%%%%%%%% T DISORDER %%%%%%%%%%%%%%
T=[   0    0    0    0    0    0    0    0    0    0    0    0  230    0  191   80  331  377  198  874  862  497 1252 1065  943  368  308  525  376  458
      0    0    0   36   15    0    0    0    0   10    0    0  371  154  331  305  435  539  382 1032 1034  572 1294 1055  894  336  284  438  298  403
      0    0    0    0    0    0    0    0    0    0    0    0    0    0   54    0  199  328  138  795  791  443 1197 1014  894  334  249  542  419  486
      0   36    0    0   23    0    0    0    0  150  157   53  389  208  356  315  528  646  341  874  846  570 1270 1007  631  127  197  216  141  351
      0   15    0   23    0   11    0    0    0  185  186  130  429  291  481  385  709  818  406  986  972  671 1355 1085  616  133  204  131   96  253
      0    0    0    0   11    0    0    0    0  198  126  132  431  301  492  418  713  763  385  882  823  604 1065  798  467   65  147  122   77  244
      0    0    0    0    0    0    0    9   60  132   94   53  306  166  364  163  339  274  218  773  675  389  835  558  161    0   31    0    0  197
      0    0    0    0    0    0    9    0  132  195  123  143  421  253  500  234  407  347  294  847  752  453  887  649  103    0    0    0    0   93
      0    0    0    0    0    0   60  132    0  136  112   97  367  286  377  217  445  365  245  622  485  307  604  355   53    0    0    0    0   60
      0   10    0  150  185  198  132  195  136    0    7    0    2   74    0   10   64   69   63  408  305  123  398  173    0    0    0    0    0   51
      0    0    0  157  186  126   94  123  112    7    0    0    0   99    0   36   70   80  103  430  337  152  402  246    0    0    0    0    0    0
      0    0    0   53  130  132   53  143   97    0    0    0    0    0    0    0   28    0    0  281  159   33  322  102    0    0    0    0    0    0
    230  371    0  389  429  431  306  421  367    2    0    0    0  151    5   26   95   99   85  273  201  143  285   70    0    0    0    0    0    0
      0  154    0  208  291  301  166  253  286   74   99    0  151    0  102    0   42   21   23  209  187  142  446  290  120    0    0    0    0   31
    191  331   54  356  481  492  364  500  377    0    0    0    5  102    0    3   88   77   65  287  194  131  329   90    0    0    0    0    0  102
     80  305    0  315  385  418  163  234  217   10   36    0   26    0    3    0    0    0   66   77   68  119  251  163    0    0    0    0    0   28
    331  435  199  528  709  713  339  407  445   64   70   28   95   42   88    0    0   15    0   15    8    0    0    0    0    0    0    0    0    0
    377  539  328  646  818  763  274  347  365   69   80    0   99   21   77    0   15    0    0   64    8    0   32    0    0    0    0    0    0  223
    198  382  138  341  406  385  218  294  245   63  103    0   85   23   65   66    0    0    0    0    0   23  119   13    0    0    0    0    0   71
    874 1032  795  874  986  882  773  847  622  408  430  281  273  209  287   77   15   64    0    0   20   37    0    0    0    0    0    0    0  276
    862 1034  791  846  972  823  675  752  485  305  337  159  201  187  194   68    8    8    0   20    0    0   23    0    0    0    8    5  180  432
    497  572  443  570  671  604  389  453  307  123  152   33  143  142  131  119    0    0   23   37    0    0   99    8    0    0   12    0   20  255
   1252 1294 1197 1270 1355 1065  835  887  604  398  402  322  285  446  329  251    0   32  119    0   23   99    0   37    0   45  177   66  344  569
   1065 1055 1014 1007 1085  798  558  649  355  173  246  102   70  290   90  163    0    0   13    0    0    8   37    0   99   72  132  174  422  571
    943  894  894  631  616  467  161  103   53    0    0    0    0  120    0    0    0    0    0    0    0    0    0   99    0   86   69    0  191  221
    368  336  334  127  133   65    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0   45   72   86    0    0    0   69   68
    308  284  249  197  204  147   31    0    0    0    0    0    0    0    0    0    0    0    0    0    8   12  177  132   69    0    0    0    0   93
    525  438  542  216  131  122    0    0    0    0    0    0    0    0    0    0    0    0    0    0    5    0   66  174    0    0    0    0    3    0
    376  298  419  141   96   77    0    0    0    0    0    0    0    0    0    0    0    0    0    0  180   20  344  422  191   69    0    3    0    0
    458  403  486  351  253  244  197   93   60   51    0    0    0   31  102   28    0  223   71  276  432  255  569  571  221   68   93    0    0    0
];
%%%%%%%%% X EXPERIMENTS %%%%%%%%%%%%%%
X=[   0    0    1    0    0    0    0    0    0    0    0    0    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1
      0    0    0    0    0    0    0    0    0   -1    0    0    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1
      1    0    0    0    0    0    0    0    0    0    0    0    0    0   -1    0    1    1    1    1    1    1    1    1    1    1    1    1    1    1
      0    0    0    0    0    0    0    0    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1
      0    0    0    0    0    0    0    0    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1
      0    0    0    0    0    0    0    0    0    0    0    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1
      0    0    0    0    0    0    0    0    0    0    0    1    1    1    1    1    1    1    1    1    1    1    1    1    1    0    0    0    1    1
      0    0    0    0    0    0    0    0    1    0    0    1    1    1    1    1    1    1    1    1    1    1    1    1    1    0    0    1   -1    1
      0    0    0    1    1    0    0    1    0    0    0    0    1    1    1    1    1    1    1    1    1    1    1    1    1    0    0    0    0    1
      0   -1    0    1    1    0    0    0    0    0    0    1    1    1    1    1    1    1    1    1    1    1    1    1    0    0    0    0    0    1
      0    0    0    1    1    0    0    0    0    0    0    0    1    0    1    1    1    1    1    1    1    1    0    0    0    0    0    0    0    0
      0    0    0    1    1    1    1    1    0    1    0    0    1    1    1    1    1    1    1    1    1    1    1    1    1    0    0    0    0    1
      1    1    0    1    1    1    1    1    1    1    1    1    0    1    0    1    0    0    0    0    0    0    0    0    0    0    0    0    0   -1
      1    1    0    1    1    1    1    1    1    1    0    1    1    0    1    0   -1    0    1    1    1    1    1    1    1    0    1    1    0    1
      1    1   -1    1    1    1    1    1    1    1    1    1    0    1    0    1    0    0    0    0    0    0    0    0    0    0    0   -1   -1    1
      1    1    0    1    1    1    1    1    1    1    1    1    1    0    1    0    0    0    1    0    1    1    1    1    0    0    0    0    0    1
      1    1    1    1    1    1    1    1    1    1    1    1    0   -1    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
      1    1    1    1    1    1    1    1    1    1    1    1    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0   -1    1
      1    1    1    1    1    1    1    1    1    1    1    1    0    1    0    1    0    0    0    0    0    0    1    0    0    0    0    0    0    1
      1    1    1    1    1    1    1    1    1    1    1    1    0    1    0    0    0    0    0    0    0    0    0    0    0    1    0    0   -1    1
      1    1    1    1    1    1    1    1    1    1    1    1    0    1    0    1    0    0    0    0    0    0    0    0    0    0    0    1    0    1
      1    1    1    1    1    1    1    1    1    1    1    1    0    1    0    1    0    0    0    0    0    0    1    0    0    0    0    1    1    1
      1    1    1    1    1    1    1    1    1    1    0    1    0    1    0    1    0    0    1    0    0    1    0    0    0    1    1    1    1    1
      1    1    1    1    1    1    1    1    1    1    0    1    0    1    0    1    0    0    0    0    0    0    0    0    0    0    1    1    1    1
      1    1    1    1    1    1    1    1    1    0    0    1    0    1    0    0    0    0    0    0    0    0    0    0    0    0    1    1    1    1
      1    1    1    1    1    1    0    0    0    0    0    0    0    0    0    0    0    0    0    1    0    0    1    0    0    0    0    1    1    1
      1    1    1    1    1    1    0    0    0    0    0    0    0    1    0    0    0    0    0    0    0    0    1    1    1    0    0    1    1    1
      1    1    1    1    1    1    0    1    0    0    0    0    0    1   -1    0    0    0    0    0    1    1    1    1    1    1    1    0    1    1
      1    1    1    1    1    1    1   -1    0    0    0    0    0    0   -1    0    0   -1    0   -1    0    1    1    1    1    1    1    1    0    0
      1    1    1    1    1    1    1    1    1    1    0    1   -1    1    1    1    0    1    1    1    1    1    1    1    1    1    1    1    0    0
];

[i N]=size(elems);  
horX=+0.05;horY=0.4;
verX=-0.60;verY=0.80;
shiftX=0.1;shiftY=-0.5;


%N=20

FONTDIM=16;
FONTDIMsmall=14;
FONTDIMbig=18;
COLOR=[0 0 0];
LINEWIDTH=1.2;
MARKERSIZE=25;
PAPERPOSITION=0.4*[0 0 N N]; % standard ratio is 16/9;
PAPERSIZE=[8.5 11];
POSITION=[420 560 560 420]; % video
AXISWIDTH=1;
AXISPOSITION=[0.0 1 1 0.0]*N;
VERBOSE=1; 


H_FIGURE=figure(1);clf;
%set(H_FIGURE,'PaperUnits','inches');
%set(H_FIGURE,'PaperOrientation','portrait');
set(H_FIGURE,'PaperPosition',PAPERPOSITION); % left top width height
%set(H_FIGURE,'PaperSize',PAPERSIZE);
%set(H_FIGURE,'PaperType','usletter');
%set(H_FIGURE,'Position',POSITION);




greylight2=[hex2dec('f0') hex2dec('f0') hex2dec('f0')]/255;
greylight=[hex2dec('d9') hex2dec('d9') hex2dec('d9')]/255;
greydark=[hex2dec('7f') hex2dec('7f') hex2dec('7f')]/255;
green=[hex2dec('4f')         hex2dec('93')         hex2dec('3e')]/255;
%green=[hex2dec('00') hex2dec('aa') hex2dec('00')]/255;
red=[hex2dec('ff') hex2dec('00') hex2dec('00')]/255;
bluelight=[hex2dec('8a') hex2dec('96') hex2dec('e2')]/255;
bluedark=[hex2dec('38') hex2dec('48') hex2dec('b5')]/255;
pink=[hex2dec('ff') hex2dec('00') hex2dec('ff')]/255;
yellow=[hex2dec('ff') hex2dec('ff') hex2dec('00')]/255;

hsvM=hsv;boneM=bone;hotM=hot;jetM=jet;
pinkM=pink;coolM=cool;prismM=prism;flagM=flag;
grayM=gray;colormapM=colormap;


cshade=0.97;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cla;
%axis ij;
axis off;

t=sqrt((T/max(max(T))));

jet;

xlim([0,2*N+1]);ylim([0,2*N+1]);
sizeC=1.5;
sizeR=1.4;
for i=1:N
  for j=1:N
    if i<j 
     posX=2*i-1+shiftX;posY=2*j-1+shiftY;
      if t(i,j)>=0.0 & t(i,j)<=0.01
        %circle greylight
        rectangle('Position',[posX,posY,sizeC,sizeC],'Curvature',[1,1],'LineWidth',1,'EdgeColor',cshade*greylight2,'FaceColor',greylight2);
        %daspect([1,1,1]);
      end
      if t(i,j)>0.01 & t(i,j)<=1.00
        color=((interp1((0:63)/63,colormapM,1.0-t(i,j))));%color(1)=color(1)*0.10;color(2)=color(2)*0.85;color(3)=1-0.8*(1-color(3));
        %circle grey light
        rectangle('Position',[posX,posY,sizeC,sizeC],'Curvature',[1,1],'LineWidth',1,'EdgeColor',cshade*color,'FaceColor',color);
        %daspect([1,1,1]);
      end
      
      %  if t(i,j)>0.5 & t(i,j)<=0.8
      %     %circle greydark
      %	   rectangle('Position',[2*i-1+shiftX,2*j-1+shiftY,sizeC,sizeC],'Curvature',[1,1],'LineWidth',1,'EdgeColor',greydark,'FaceColor',greydark);
      %      %daspect([1,1,1]);
      %   end
      %   if t(i,j)>0.8 & t(i,j)<=0.9
      %     %rectangle green
      % 		h=rectangle('Position',[2*i-1+shiftX,2*j-1+shiftY,sizeR,sizeR],'Curvature',[0,0],'LineWidth',1,'EdgeColor',green,'FaceColor',green);
      %      %daspect([1,1,1]);
      %   end
      %   if t(i,j)>0.9 & t(i,j)<=1.0
      %     %rectangle red
      %		h=rectangle('Position',[2*i-1+shiftX,2*j-1+shiftY,sizeR,sizeR],'Curvature',[0,0],'LineWidth',1,'EdgeColor',red,'FaceColor',red);
      %       %daspect([1,1,1]);
      %    end
      
    end
  end
end

%colormap;


for i=1:N
  if i==1 | i==9  | i==11  
    text(2*i-1+horX+shiftX+0.33,horY+shiftY+2*N+1,elems(i),'FontName','Courier','FontSize',FONTDIMsmall,'FontWeight','bold');
    text(verX+shiftX,2*i-1+verY+shiftY,elems(i),'FontName','Courier','FontSize',FONTDIMsmall,'FontWeight','bold');
  else
    text(2*i-1+horX+shiftX,horY+shiftY+2*N+1,elems(i),'FontName','Courier','FontSize',FONTDIMsmall,'FontWeight','bold');
    text(verX+shiftX-0.33,2*i-1+verY+shiftY,elems(i),'FontName','Courier','FontSize',FONTDIMsmall,'FontWeight','bold');
  end    
end


%h=gca;
%get(h);
%set(h,'Visible','off');
%set(h,'Color',[1 1 0]);
%set(h,'XColor',[1 1 1]);

print -depsc orderTR.eps
%unix('convert orderTR.eps orderTR.pdf');
unix('convert orderTR.eps orderTR.png');

%%%%%%%%%%%%%%%%%%%%%%%%%%
hold off

cla;

%axis ij;
axis off;

x=X;


xlim([0,2*N+1]);ylim([0,2*N+1]);
sizeC=1.5;
sizeR=1.4;
daspect([1,1,1]);

for i=1:N
  for j=1:N
    if i<j 
      taU=sizeC*[1 1 0];tbU=sizeC*[0 1 1];
      taB=sizeC*[0 1 0];tbB=sizeC*[0 0 1];
      posX=2*i-1+shiftX;posY=2*j-1+shiftY;
      %   trianglesupX=[0 1E-8 0]+2*i-1+shiftX;trianglesupY=[0 0
      %   1E-8]+2*j-1+shiftY;
        color=[1 1 1];Curvature=[1,1];
   %     rectangle('Position',[posX,posY,sizeC,sizeC],'Curvature',Curvature,'LineWidth',1,'EdgeColor',cshade*color,'FaceColor',color);
        hold on
      if (t(i,j)==0.0) & x(i,j)==0
        %circle greylight
        rectangle('Position',[posX,posY,sizeC,sizeC],'Curvature',[1,1],'LineWidth',1,'EdgeColor',cshade*greylight2,'FaceColor',greylight2);
        %daspect([1,1,1]);
      end      
     if t(i,j)>0.0 & x(i,j)==1  % MIX
        color=green; 
        rectangle('Position',[posX,posY,sizeC,sizeC],'Curvature',Curvature,'LineWidth',1,'EdgeColor',cshade*color,'FaceColor',color);
        hold on
     end
     if t(i,j)>0.0 & x(i,j)~=1   % MIX HT NO-MIX EXPERIMENTS
        color=red; 
        fill(posX+taU,posY+tbU,color,'EdgeColor',cshade*color,'FaceColor',color);%,'LineWidth',1,'EdgeColor',cshade*color);
  %     rectangle('Position',[posX,posY,sizeC,sizeC],'Curvature',Curvature,'LineWidth',1,'EdgeColor',cshade*color,'FaceColor',color);
        hold on
      end
      if (x(i,j)==1 | x(i,j)==-1) & t(i,j)==0.0   % EXPERIMENTS
        if x(i,j)==1 color=bluelight; end
        if x(i,j)==-1 color=yellow; end
        z=fill(posX+taB,posY+tbB,color,'EdgeColor',cshade*color,'FaceColor',color);%,'LineWidth',1,'EdgeColor',cshade*color);
        hold on
        % daspect([1,1,1]);
      end
  %      if t(i,j)>0.01 & t(i,j)<=1.00
%        color=((interp1((0:63)/63,grayM,1.0-t(i,j))));color(2)=color(2)*0.10;color(3)=color(3)*0.15;
%        %circle grey light
%        rectangle('Position',[posX,posY,sizeC,sizeC],'Curvature',[1,1],'LineWidth',1,'EdgeColor',cshade*color,'FaceColor',color);
%        %daspect([1,1,1]);
%      end
      
%      if x(i,j)==0    % nomix
%        color=[0.1 0.1 0.1];Curvature=[1,1];
%        %            rectangle('Position',[posX,posY,sizeC,sizeC],'Curvature',Curvature,'LineWidth',1,'EdgeColor',cshade*color,'FaceColor',color);
%      end
%      if x(i,j)==-1     % unknown
                        %   color=[0 0 0.5];Curvature=[0,0];
%      end
%      
%      daspect([1,1,1]);
      
      
    end
  end
end


for i=1:N
  if i==1 | i==9  | i==11  
    text(2*i-1+horX+shiftX+0.33,horY+shiftY+2*N+1,elems(i),'FontName','Courier','FontSize',FONTDIMsmall,'FontWeight','bold');
    text(verX+shiftX,2*i-1+verY+shiftY,elems(i),'FontName','Courier','FontSize',FONTDIMsmall,'FontWeight','bold');
  else
    text(2*i-1+horX+shiftX,horY+shiftY+2*N+1,elems(i),'FontName','Courier','FontSize',FONTDIMsmall,'FontWeight','bold');
    text(verX+shiftX-0.33,2*i-1+verY+shiftY,elems(i),'FontName','Courier','FontSize',FONTDIMsmall,'FontWeight','bold');
  end    
end


taU=sizeC*[1 1 0];tbU=sizeC*[0 1 1];
taB=sizeC*[0 1 0];tbB=sizeC*[0 0 1];
i=14;j=10;
j=j-1;posX=2*i-1+shiftX;posY=2*j-1+shiftY;
color=green; 
rectangle('Position',[posX,posY,sizeC,sizeC],'Curvature',Curvature,'LineWidth',1,'EdgeColor',cshade*color,'FaceColor',color);
text(posX+2,posY+0.8,'mix {HT=exp.}','FontName','Courier','FontSize',FONTDIMsmall,'FontWeight','bold');
j=j-1;posX=2*i-1+shiftX;posY=2*j-1+shiftY;
color=green; 
rectangle('Position',[posX,posY,sizeC,sizeC],'Curvature',[1,1],'LineWidth',1,'EdgeColor',cshade*greylight2,'FaceColor',greylight2);
text(posX+2,posY+0.8,'no-mix {HT=exp.}','FontName','Courier','FontSize',FONTDIMsmall,'FontWeight','bold');
j=j-1;posX=2*i-1+shiftX;posY=2*j-1+shiftY;
color=red; 
fill(posX+taU,posY+tbU,color,'EdgeColor',cshade*color,'FaceColor',color);
text(posX+2,posY+0.8,'mix {HT \neq no-mix exp.}','FontName','Courier','FontSize',FONTDIMsmall,'FontWeight','bold');
j=j-1;posX=2*i-1+shiftX;posY=2*j-1+shiftY;
color=bluelight; 
fill(posX+taB,posY+tbB,color,'EdgeColor',cshade*color,'FaceColor',color);
text(posX+2,posY+0.8,'no-mix HT \neq mix exp.','FontName','Courier','FontSize',FONTDIMsmall,'FontWeight','bold');
j=j-1;posX=2*i-1+shiftX;posY=2*j-1+shiftY;
color=yellow; 
fill(posX+taB,posY+tbB,color,'EdgeColor',cshade*color,'FaceColor',color);
text(posX+2,posY+0.8,'no-mix HT \neq unk. exp.','FontName','Courier','FontSize',FONTDIMsmall,'FontWeight','bold');


%h=gca;
%get(h);
%set(h,'Visible','off');
%set(h,'Color',[1 1 0]);
%set(h,'XColor',[1 1 1]);

print -depsc experimentsTL.eps
unix('convert experimentsTL.eps experimentsTL.png');

%%%%%%%%%%%%%%%%%%%%%%%%%%
hold off

cla;

%axis ij;
axis off;

x=X;


xlim([0,2*N+1]);ylim([0,2*N+1]);
sizeC=1.5;
sizeR=1.4;
daspect([1,1,1]);

for i=1:N
  for j=1:N
    if i>j 
      taU=sizeC*[1 1 0];tbU=sizeC*[0 1 1];
      taB=sizeC*[0 1 0];tbB=sizeC*[0 0 1];
      posX=2*i-1+shiftX;posY=2*j-1+shiftY;
      %   trianglesupX=[0 1E-8 0]+2*i-1+shiftX;trianglesupY=[0 0
      %   1E-8]+2*j-1+shiftY;
        color=[1 1 1];Curvature=[1,1];
   %     rectangle('Position',[posX,posY,sizeC,sizeC],'Curvature',Curvature,'LineWidth',1,'EdgeColor',cshade*color,'FaceColor',color);
        hold on
      if (t(i,j)==0.0) & x(i,j)==0
        %circle greylight
        rectangle('Position',[posX,posY,sizeC,sizeC],'Curvature',[1,1],'LineWidth',1,'EdgeColor',cshade*greylight2,'FaceColor',greylight2);
        %daspect([1,1,1]);
      end      
     if t(i,j)>0.0 & x(i,j)==1  % MIX
        color=green; 
        rectangle('Position',[posX,posY,sizeC,sizeC],'Curvature',Curvature,'LineWidth',1,'EdgeColor',cshade*color,'FaceColor',color);
        hold on
     end
     if t(i,j)>0.0 & x(i,j)~=1   % MIX HT NO-MIX EXPERIMENTS
        color=red; 
        fill(posX+taU,posY+tbU,color,'EdgeColor',cshade*color,'FaceColor',color);%,'LineWidth',1,'EdgeColor',cshade*color);
  %     rectangle('Position',[posX,posY,sizeC,sizeC],'Curvature',Curvature,'LineWidth',1,'EdgeColor',cshade*color,'FaceColor',color);
        hold on
      end
      if (x(i,j)==1 | x(i,j)==-1) & t(i,j)==0.0   % EXPERIMENTS
        if x(i,j)==1 color=bluelight; end
        if x(i,j)==-1 color=yellow; end
        z=fill(posX+taB,posY+tbB,color,'EdgeColor',cshade*color,'FaceColor',color);%,'LineWidth',1,'EdgeColor',cshade*color);
        hold on
        % daspect([1,1,1]);
      end
  %      if t(i,j)>0.01 & t(i,j)<=1.00
%        color=((interp1((0:63)/63,grayM,1.0-t(i,j))));color(2)=color(2)*0.10;color(3)=color(3)*0.15;
%        %circle grey light
%        rectangle('Position',[posX,posY,sizeC,sizeC],'Curvature',[1,1],'LineWidth',1,'EdgeColor',cshade*color,'FaceColor',color);
%        %daspect([1,1,1]);
%      end
      
%      if x(i,j)==0    % nomix
%        color=[0.1 0.1 0.1];Curvature=[1,1];
%        %            rectangle('Position',[posX,posY,sizeC,sizeC],'Curvature',Curvature,'LineWidth',1,'EdgeColor',cshade*color,'FaceColor',color);
%      end
%      if x(i,j)==-1     % unknown
                        %   color=[0 0 0.5];Curvature=[0,0];
%      end
%      
%      daspect([1,1,1]);
      
      
    end
  end
end

for i=1:N
  if i==1 | i==9  | i==11  
    text(2*i-1+horX+shiftX+0.33,horY+shiftY+0.7+0*2*N-1,elems(i),'FontName','Courier','FontSize',FONTDIMsmall,'FontWeight','bold');
    text(verX+shiftX+2*N+2,2*i-1+verY+shiftY,elems(i),'FontName','Courier','FontSize',FONTDIMsmall,'FontWeight','bold');
  else
    text(2*i-1+horX+shiftX,horY+shiftY+0.7-1,elems(i),'FontName','Courier','FontSize',FONTDIMsmall,'FontWeight','bold');
    text(verX+shiftX+2*N+2-0.33,2*i-1+verY+shiftY,elems(i),'FontName','Courier','FontSize',FONTDIMsmall,'FontWeight','bold');
  end    
end


%h=gca;
%get(h);
%set(h,'Visible','off');
%set(h,'Color',[1 1 0]);
%set(h,'XColor',[1 1 1]);

print -depsc experimentsBR.eps
unix('convert experimentsBR.eps experimentsBR.png');



%%%%%%%%%%%%%%%%%%%%%%%%%%
hold off

cla;

%axis ij;
axis off;

x=X;


xlim([0,2*N+1]);ylim([0,2*N+1]);
sizeC=1.5;
sizeR=1.4;
daspect([1,1,1]);

for i=1:N
  for j=1:N
    posX=2*i-1+shiftX;posY=2*j-1+shiftY;
    if i<j 
      if t(i,j)>=0.0 & t(i,j)<=0.01
        %circle greylight
        rectangle('Position',[posX,posY,sizeC,sizeC],'Curvature',[1,1],'LineWidth',1,'EdgeColor',cshade*greylight2,'FaceColor',greylight2);
        %daspect([1,1,1]);
      end
      if t(i,j)>0.01 & t(i,j)<=1.00
        color=((interp1((0:63)/63,colormapM,1.0-t(i,j))));%color(1)=color(1)*0.10;color(2)=color(2)*0.85;color(3)=1-0.8*(1-color(3));
        %circle grey light
        rectangle('Position',[posX,posY,sizeC,sizeC],'Curvature',[1,1],'LineWidth',1,'EdgeColor',cshade*color,'FaceColor',color);
        %daspect([1,1,1]);
      end
      
%      value=3*i/N;
%      if value <=1.0 
%	color=((interp1((0:63)/63,colormapM,value)));
%        rectangle('Position',[posX,posY,sizeC,sizeC],'Curvature',[1,1],'LineWidth',1,'EdgeColor',cshade*color,'FaceColor',color);     
%      else
%	color=greylight2;
%        rectangle('Position',[posX,posY,sizeC,sizeC],'Curvature',[1,1],'LineWidth',1,'EdgeColor',cshade*color,'FaceColor',color);     
%     end
 
      %  if t(i,j)>0.5 & t(i,j)<=0.8
      %     %circle greydark
      %	   rectangle('Position',[2*i-1+shiftX,2*j-1+shiftY,sizeC,sizeC],'Curvature',[1,1],'LineWidth',1,'EdgeColor',greydark,'FaceColor',greydark);
      %      %daspect([1,1,1]);
      %   end
      %   if t(i,j)>0.8 & t(i,j)<=0.9
      %     %rectangle green
      % 		h=rectangle('Position',[2*i-1+shiftX,2*j-1+shiftY,sizeR,sizeR],'Curvature',[0,0],'LineWidth',1,'EdgeColor',green,'FaceColor',green);
      %      %daspect([1,1,1]);
      %   end
      %   if t(i,j)>0.9 & t(i,j)<=1.0
      %     %rectangle red
      %		h=rectangle('Position',[2*i-1+shiftX,2*j-1+shiftY,sizeR,sizeR],'Curvature',[0,0],'LineWidth',1,'EdgeColor',red,'FaceColor',red);
      %       %daspect([1,1,1]);
      %    end
    end
       
    if i>j 
      taU=sizeC*[1 1 0];tbU=sizeC*[0 1 1];
      taB=sizeC*[0 1 0];tbB=sizeC*[0 0 1];
      posX=2*i-1+shiftX;posY=2*j-1+shiftY;
      %   trianglesupX=[0 1E-8 0]+2*i-1+shiftX;trianglesupY=[0 0
      %   1E-8]+2*j-1+shiftY;
        color=[1 1 1];Curvature=[1,1];
   %     rectangle('Position',[posX,posY,sizeC,sizeC],'Curvature',Curvature,'LineWidth',1,'EdgeColor',cshade*color,'FaceColor',color);
        hold on
      if (t(i,j)==0.0) & x(i,j)==0
        %circle greylight
        rectangle('Position',[posX,posY,sizeC,sizeC],'Curvature',[1,1],'LineWidth',1,'EdgeColor',cshade*greylight2,'FaceColor',greylight2);
        %daspect([1,1,1]);
      end      
     if t(i,j)>0.0 & x(i,j)==1  % MIX
        color=green; 
        rectangle('Position',[posX,posY,sizeC,sizeC],'Curvature',Curvature,'LineWidth',1,'EdgeColor',cshade*color,'FaceColor',color);
        hold on
     end
     if t(i,j)>0.0 & x(i,j)~=1   % MIX HT NO-MIX EXPERIMENTS
        color=red; 
        fill(posX+taU,posY+tbU,color,'EdgeColor',cshade*color,'FaceColor',color);%,'LineWidth',1,'EdgeColor',cshade*color);
  %     rectangle('Position',[posX,posY,sizeC,sizeC],'Curvature',Curvature,'LineWidth',1,'EdgeColor',cshade*color,'FaceColor',color);
        hold on
      end
      if (x(i,j)==1 | x(i,j)==-1) & t(i,j)==0.0   % EXPERIMENTS
        if x(i,j)==1 color=bluelight; end
        if x(i,j)==-1 color=yellow; end
        z=fill(posX+taB,posY+tbB,color,'EdgeColor',cshade*color,'FaceColor',color);%,'LineWidth',1,'EdgeColor',cshade*color);
        hold on
        % daspect([1,1,1]);
      end
  %      if t(i,j)>0.01 & t(i,j)<=1.00
%        color=((interp1((0:63)/63,grayM,1.0-t(i,j))));color(2)=color(2)*0.10;color(3)=color(3)*0.15;
%        %circle grey light
%        rectangle('Position',[posX,posY,sizeC,sizeC],'Curvature',[1,1],'LineWidth',1,'EdgeColor',cshade*color,'FaceColor',color);
%        %daspect([1,1,1]);
%      end
      
%      if x(i,j)==0    % nomix
%        color=[0.1 0.1 0.1];Curvature=[1,1];
%        %            rectangle('Position',[posX,posY,sizeC,sizeC],'Curvature',Curvature,'LineWidth',1,'EdgeColor',cshade*color,'FaceColor',color);
%      end
%      if x(i,j)==-1     % unknown
                        %   color=[0 0 0.5];Curvature=[0,0];
%      end
%      
%      daspect([1,1,1]);
      
      
    end
  end
end

for i=1:N
  if i==1 | i==9  | i==11  
    text(2*i-1+horX+shiftX+0.33,horY+shiftY+2*N+1,elems(i),'FontName','Courier','FontSize',FONTDIMsmall,'FontWeight','bold');
    text(verX+shiftX,2*i-1+verY+shiftY,elems(i),'FontName','Courier','FontSize',FONTDIMsmall,'FontWeight','bold');
    text(2*i-1+horX+shiftX+0.33,horY+shiftY+0.7+0*2*N-1,elems(i),'FontName','Courier','FontSize',FONTDIMsmall,'FontWeight','bold');
    text(verX+shiftX+2*N+2,2*i-1+verY+shiftY,elems(i),'FontName','Courier','FontSize',FONTDIMsmall,'FontWeight','bold');
  else
    text(2*i-1+horX+shiftX,horY+shiftY+2*N+1,elems(i),'FontName','Courier','FontSize',FONTDIMsmall,'FontWeight','bold');
    text(verX+shiftX-0.33,2*i-1+verY+shiftY,elems(i),'FontName','Courier','FontSize',FONTDIMsmall,'FontWeight','bold');
    text(2*i-1+horX+shiftX,horY+shiftY+0.7-1,elems(i),'FontName','Courier','FontSize',FONTDIMsmall,'FontWeight','bold');
    text(verX+shiftX+2*N+2-0.33,2*i-1+verY+shiftY,elems(i),'FontName','Courier','FontSize',FONTDIMsmall,'FontWeight','bold');
  end    
end


%h=gca;
%get(h);
%set(h,'Visible','off');
%set(h,'Color',[1 1 0]);
%set(h,'XColor',[1 1 1]);

print -depsc order_experiments.eps
unix('convert order_experiments.eps order_experiments.png');



