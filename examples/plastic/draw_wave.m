clc
clear

addpath(genpath('../../mscripts'))

par = ReadYaml('parameters.yaml');
nproc = par.nproc;
data_dir = par.data_dir;
varnm = 'damage';
%varnm = 'Vx';
%data_dir = 'data_Psi9'
data_dir = 'data_Psi20'
%data_dir = 'data_Psi50'
data_dir = 'data'

[x,y,tvec] = gather_wave_coord(data_dir, nproc);
nt = length(tvec);
[NGLL,~,nelem] = size(x);

% fast plot using trisurf
tri = [];
for i = 1:nelem
    tri1 = delaunay(x(:,:,i),y(:,:,i)) + (i-1)*NGLL*NGLL;
    tri = cat(1,tri,tri1);
end

figure
for it = 1:1:nt

    [v1,t] = gather_wave_snap(data_dir, nproc, varnm, it);

    if strcmp(varnm, 'damage')
        %v1(v1<1e-16)=1e-16;
        v1 = log10(v1);
    end
    vmax = max(abs(v1(:)));
    disp(['it = ',num2str(it),' vmax = ',num2str(vmax)]);
    plotSolutionFast(tri,x*1e-3,y*1e-3,v1 )
    %colormap rdbu

    %hold off
    xlabel('X (km)')
    %ylabel('Depth (km)')
    ylabel('Y (km)')
    %vmax = 60;
    %vmax=1;
    vmax = max(vmax,1e-16);
    %caxis([-1 1]*vmax/2)
    title(['T = ',num2str(t),' sec'])
    cb = colorbar;
    title(cb,'log_{10}\eta')

    %axis([-41 41  -10 10])
    if strcmp(varnm, 'damage')
       caxis([-8 -6])
    end
    hold on
    %plot3([-15 15],[0 0],[1 1]*1e30,'k','LineWidth',1.5)

    hold off

    pause(0.005)

end
