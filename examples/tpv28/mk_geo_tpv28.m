clc
clear
close all

faultsize=200;
boundsize=10e3;

x = -18e3:50:18e3;
nx = length(x);
y = zeros(size(x));

for i = 1:nx
     x1 = x(i);
    if (x1>7.5e3 && x1<13.5e3)
        y(i) = 300*(1+cos(pi*abs(x1-10.5e3)/3e3));
    end
    if (x1<-7.5e3 && x1>-13.5e3)
        y(i) = 300*(1+cos(pi*abs(x1+10.5e3)/3e3));
    end
end

nx = length(x);
plot(x,y,'o-')
axis equal


fnm = 'tpv28_2d.geo'

fid = fopen(fnm,'w');
fprintf(fid,'SetFactory("OpenCASCADE");\n');
for i = 1:nx
   fprintf(fid,'Point(%d)={%g,%g,0,1.0};\n',i,x(i),y(i));
end
fprintf(fid,'Spline(1) = {%d:%d};\n',1,nx);

% add four corners
LenX = 50e3; LenY = 40e3;
c = [-LenX, -LenY;
      LenX, -LenY;
      LenX,  LenY;
     -LenX,  LenY];
hold on
plot(c(:,1),c(:,2),'ro')

fprintf(fid,'Rectangle(1) = {-50e3, -40e3, 0, 100e3, 80e3, 0};\n');

fprintf(fid,'BooleanFragments{ Surface{1}; Delete; }{ Curve{1}; Delete; }\n',1,i);
fprintf(fid,'MeshSize{%d:%d} = %g;\n',1,nx,faultsize*2);
fprintf(fid,'MeshSize{%d:%d} = %g;\n',nx+1,nx+4,boundsize*2);
fprintf(fid,'Mesh.Algorithm = 6;\n');
%fprintf(fid,'Mesh.RecombinationAlgorithm = 2;\n');
%fprintf(fid,'Recombine Surface{:};\n');
fprintf(fid,'Mesh.SubdivisionAlgorithm = 1;\n');
fclose(fid);
