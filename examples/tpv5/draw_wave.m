clc
clear

addpath(genpath('../../mscripts'))

nproc = 4;

data_dir = 'data';

%it = 1;
%x = []; y = []; v = [];
%for i = 1:nproc
%    fnm = [data_dir,'/wave_mpi',num2str(i-1,'%06d'),'.nc'];
%    if (exist(fnm,'file'))
%        disp(fnm)
%        x1 = ncread(fnm, 'x');
%        y1 = ncread(fnm, 'y');
%        v1 = ncread(fnm, 'Vy',[1,1,1,it],[Inf Inf Inf 1]);
%        t = ncread(fnm, 'time');
%        x = cat(3,x,x1);
%        y = cat(3,y,y1);
%        v = cat(3,v,v1);
%    end
%end
[x,y,tvec] = gather_wave_coord(data_dir, nproc);

[NGLL,~,nelem] = size(x);

tic
% fast plot using trisurf
tri2 = [];
for i = 1:nelem
    tri1 = delaunay(x(:,:,i),y(:,:,i)) + (i-1)*NGLL*NGLL;
    tri2 = [tri2;tri1];
end
tri = tri2;
toc


%figure('Position',[100 300 800 400])
figure
for it =   1:1:1000


    [v1,t] = gather_wave_snap(data_dir, nproc, 'Vx', it);
    %v1 = v(:,:,:,it);
    vmax = max(abs(v1(:)));
    plotSolutionFast(tri,x*1e-3,y*1e-3,v1,1   )
  

    colormap rdbu

    %hold off
    xlabel('Horizontal Distance (km)')
    %ylabel('Depth (km)')
    ylabel('Normal Distance (km)')
    %vmax = 60;
    caxis([-1 1]*vmax/2)
    title(['T = ',num2str(t),' sec'])

    
    
    %ylim([-1 0]*50)
    %xlim([-1 1]*80)


    pause(0.005)

end


