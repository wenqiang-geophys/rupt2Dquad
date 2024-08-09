clc
clear
%close all

nproc = 4;

u = cell(nproc,1);
x = cell(nproc,1);
y = cell(nproc,1);
tri = cell(nproc,1);

data_dir = 'data';


for iproc = 1:nproc
 
fid = fopen([data_dir,'/vertices',num2str(iproc-1,'%06d')],'r');
a = fscanf(fid,'%d %d %d', 3); % NGLL NGLL nelem
NGLL = a(1);
nelem = a(3);
v = fscanf(fid,'%f',2*NGLL*NGLL*nelem);
v = reshape(v,2,NGLL*NGLL*nelem);
fclose(fid);

x{iproc} = reshape(v(1,:),NGLL,NGLL,nelem);
y{iproc} = reshape(v(2,:),NGLL,NGLL,nelem);

if 1
tic
% fast plot using trisurf
tri2 = [];
for i = 1:nelem
    tri1 = delaunay(x{iproc}(:,:,i),y{iproc}(:,:,i)) + (i-1)*NGLL*NGLL;
    tri2 = [tri2;tri1];
end

tri{iproc} = tri2;
toc
end
%return

outdir = data_dir;
%addpath('~/works/colormap');

%%xr = 0e3;
%%yr = 7e3;
%r = sqrt((x-xr).^2+(y-yr).^2);
%[~,idx] = min(r(:));
%recv = [];

fnm = [outdir, '/pdy', num2str(iproc-1,'%06d')];
usize = NGLL*NGLL*nelem*1;

%usize = 1*1*nelem*1;

fid = fopen(fnm,'rb');
tmp = fread(fid,usize,'float32');
u{iproc} = reshape(tmp,NGLL,NGLL,nelem);
%u{iproc} = repmat(reshape(tmp,1,1,nelem),NGLL,NGLL,1);
fclose(fid);
end


figure
for iproc = 1:nproc
    u1 = u{iproc}(:,:,:);
    plotSolutionFast(tri{iproc},x{iproc},y{iproc},u1)
    hold on
end
hold off
