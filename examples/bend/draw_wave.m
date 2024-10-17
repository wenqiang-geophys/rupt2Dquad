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
for it = 20%1:1:nt

    [v1,t] = gather_wave_snap(data_dir, nproc, 'Vx', it);
    vmax = max(abs(v1(:)));
    disp(['it = ',num2str(it),' vmax = ',num2str(vmax)]);
    plotSolutionFast(tri,x*1e-3,y*1e-3,v1,1)
    colormap rdbu

    %hold off
    xlabel('X (km)')
    %ylabel('Depth (km)')
    ylabel('Y (km)')
    %vmax = 60;
    %vmax=1;
    vmax = max(vmax,1e-16);
    caxis([-1 1]*vmax/5)
    title(['T = ',num2str(t),' sec'])

    %axis([-2 2 -1 1]*40)
    hold on
    plot3([-1 1/2]*40,[0 0],[1 1]*1e30,'k','LineWidth',1.)
    plot3([0 20*cosd(10)]+20,[0 20*sind(10)],[1 1]*1e30,'k','LineWidth',1.)


    hold off

    pause(0.005)

end
