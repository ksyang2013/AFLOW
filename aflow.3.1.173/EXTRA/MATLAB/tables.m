elements={'Ag','Au','Cd','Co','Cr','Cu','Fe','Hf','Hg','Ir','La','Mn','Mo','Nb','Ni','Os','Pd','Pt','Re','Rh','Ru','Sc','Ta','Tc','Ti','V','W','Y','Zn','Zr'};   
%elements={'Ag','Au','Cd','Co','Cr','Cu','Fe','Hf','Hg','Ir','La','Mn','Mo','Nb','Ni'};   

[i N]=size(elements);  
horX=+0.15;horY=0.2;
verX=-0.45;verY=0.70;
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


x=rand(N,N);
x=x.*x';

greylight=[hex2dec('d9') hex2dec('d9') hex2dec('d9')]/255;
greydark=[hex2dec('7f') hex2dec('7f') hex2dec('7f')]/255;
green=[hex2dec('4f') hex2dec('93') hex2dec('3e')]/255;
%green=[hex2dec('00') hex2dec('aa') hex2dec('00')]/255;
red=[hex2dec('ff') hex2dec('00') hex2dec('00')]/255;

hsvM=hsv;boneM=bone;hotM=hot;jetM=jet;
pinkM=pink;coolM=cool;prismM=prism;flagM=flag;
grayM=gray;




cla;
axis ij;
axis off;

xlim([0,2*N+1]);ylim([0,2*N+1]);
sizeC=1.5;
sizeR=1.4;
for i=1:N
   for j=1:N
      if i~=j 
      if x(i,j)>0.0 & x(i,j)<=0.1
         %circle grey light
   	   rectangle('Position',[2*i-1+shiftX,2*j-1+shiftY,sizeC,sizeC],'Curvature',[1,1],'LineWidth',1,'EdgeColor',greylight,'FaceColor',greylight);
         %daspect([1,1,1]);
      end
      if x(i,j)>0.1 & x(i,j)<=1.00
         color=((interp1((0:63)/63,grayM,x(i,j))));color(2)=color(2)*0.66;color(3)=color(3)*0.5;

         %circle grey light
   	   rectangle('Position',[2*i-1+shiftX,2*j-1+shiftY,sizeC,sizeC],'Curvature',[1,1],'LineWidth',1,'EdgeColor',color,'FaceColor',color);
         %daspect([1,1,1]);
      end
    
    %  if x(i,j)>0.5 & x(i,j)<=0.8
    %     %circle greydark
   %	   rectangle('Position',[2*i-1+shiftX,2*j-1+shiftY,sizeC,sizeC],'Curvature',[1,1],'LineWidth',1,'EdgeColor',greydark,'FaceColor',greydark);
   %      %daspect([1,1,1]);
   %   end
   %   if x(i,j)>0.8 & x(i,j)<=0.9
 	%     %rectangle green
  % 		h=rectangle('Position',[2*i-1+shiftX,2*j-1+shiftY,sizeR,sizeR],'Curvature',[0,0],'LineWidth',1,'EdgeColor',green,'FaceColor',green);
   %      %daspect([1,1,1]);
   %   end
   %   if x(i,j)>0.9 & x(i,j)<=1.0
 	%     %rectangle red
   %		h=rectangle('Position',[2*i-1+shiftX,2*j-1+shiftY,sizeR,sizeR],'Curvature',[0,0],'LineWidth',1,'EdgeColor',red,'FaceColor',red);
  %       %daspect([1,1,1]);
  %    end
       
      end
    end
end


for i=1:N
   text(2*i-1+horX+shiftX,horY+shiftY,elements(i),'FontName','Courier','FontSize',11);
   text(verX+shiftX,2*i-1+verY+shiftY,elements(i),'FontName','Courier','FontSize',11);
end

%h=gca;
%get(h);
%set(h,'Visible','off');
%set(h,'Color',[1 1 0]);
%set(h,'XColor',[1 1 1]);

print -depsc test.eps

