clc
close all
clear

par = ReadYaml('parameters.yaml');
nproc = par.nproc;
data_dir = par.data_dir;
varnm = 'rate';

x = []; y = []; v = [];
for i = 1:nproc

    fnm = [data_dir,'/fault_mpi',num2str(i-1,'%06d'),'.nc'];
    if (exist(fnm,'file'))
        disp(fnm)
    x1 = ncread(fnm, 'x');
    y1 = ncread(fnm, 'y');
    v1 = ncread(fnm, varnm);
    t = ncread(fnm, 'time');

    x = cat(2,x,x1);
    y = cat(2,y,y1);
    v = cat(2,v,v1);
    end

end

nt = length(t);
nflt = length(x(:));
v = reshape(v,[nflt,nt]);
[x,idx] = sort(x(:));
v = v(idx,:);

figure
pcolor(x,t,v')
shading interp
c=colorbar;
colormap jet
xlabel('X (m)')
ylabel('T (sec)')
set(gca,'FontSize',12)
