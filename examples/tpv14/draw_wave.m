clc
clear

addpath(genpath('../../mscripts'))

par = ReadYaml('parameters.yaml');
nproc = par.nproc;
data_dir = par.data_dir;

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

    [v1,t] = gather_wave_snap(data_dir, nproc, 'Vx', it);
    vmax = max(abs(v1(:))); vmax = max(vmax,1e-16);
    disp(['it = ',num2str(it),' vmax = ',num2str(vmax)]);
    plotSolutionFast(tri,x*1e-3,y*1e-3,v1 )
    colormap rdbu
    xlabel('X (km)')
    ylabel('Y (km)')
    %vmax=1;
    caxis([-1 1]*vmax/2)
    title(['T (',num2str(it),') = ',num2str(t),' sec'])

    if 0
    axis([-2 2 -1 1]*10)
    hold on
    plot3([-16 12],[0 0],[1 1]*1e30,'k','LineWidth',1.5)
    plot3([0 12*cosd(30)],[0 -12*sind(30)],[1 1]*1e30,'k','LineWidth',1.5)
    hold off
    end

    pause(0.005)

end
