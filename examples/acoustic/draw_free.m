clc
%close all
clear


nproc = 4;
varnm = 'Vy';
data_dir = 'data';

x = []; y = []; v = [];
nx = []; ny = [];
Nt = 35000;
for i = 1:nproc

    fnm = [data_dir,'/free_mpi',num2str(i-1,'%06d'),'.nc'];
    disp(fnm)
    if (exist(fnm,'file'))
    x1 = ncread(fnm, 'x');
    y1 = ncread(fnm, 'y');
    nx1 = ncread(fnm, 'nx');
    ny1 = ncread(fnm, 'ny');
    v1 = ncread(fnm, varnm, [1 1 1],[Inf,Inf,Nt]);
    t = ncread(fnm, 'time', [1],[Nt]);
    end

    %x = [x,x1];
    %y = [y,y1];
    %v = [v,v1];

    x = cat(2,x,x1);
    y = cat(2,y,y1);
    nx = cat(2,nx,nx1);
    ny = cat(2,ny,ny1);
    v = cat(2,v,v1);

end

[Nfp,Nface,Nt ] = size(v);

v = reshape(v,[Nfp*Nface,Nt]);
%return;

x = x(:)*1e-3;
y = y(:)*1e-3;

[x,idx] = sort(x);
v = v(idx,:);

dt = t(2)-t(1);
u = cumtrapz(v,2)*dt;

figure
pcolor(x,t, u')
shading interp
c=colorbar;

colormap RdBu

vm = max(abs(v(:)));
%caxis([-1 1]*vm/50)
xlabel('X (km)')
ylabel('T (sec)')
set(gca,'FontSize',12)
%fnm = ['free_x_t_mesh_asymm_mixed_',varnm];
%fnm = ['free_x_t_mesh_',key1,'_',varnm];
%print( '-r300', '-dpng', fnm)
