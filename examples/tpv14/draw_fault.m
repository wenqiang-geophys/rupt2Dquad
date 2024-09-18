clc
close all
clear

par = ReadYaml('parameters.yaml');
nproc = par.nproc;
data_dir = par.data_dir;

%varnm = 'stress';
varnm = 'rate';
fault_segment_id = 1;

x = []; y = []; v = [];
fault_id = [];
for i = 1:nproc

    fnm = [data_dir,'/fault_mpi',num2str(i-1,'%06d'),'.nc'];
    if (exist(fnm,'file'))
        disp(fnm)
    x1 = ncread(fnm, 'x');
    y1 = ncread(fnm, 'y');
    v1 = ncread(fnm, varnm);
    t = ncread(fnm, 'time');
    fault_id1 = ncread(fnm, 'fault_id');
    fault_id = cat(1,fault_id,fault_id1);

    x = cat(2,x,x1);
    y = cat(2,y,y1);
    v = cat(2,v,v1);
    end

end

idx = find(fault_id<3);
x = x(:,idx);
y = y(:,idx);
v = v(:,idx,:);

%return

nt = length(t);
nflt = length(x(:));
v = reshape(v,[nflt,nt]);
x = x(:)*1e-3;
y = y(:)*1e-3;
[x,idx] = sort(x);
v = v(idx,:);

if strcmp(varnm, 'stress')
    v = -v*1e-6;
end

figure
pcolor(x,t,v')
shading interp
%caxis([0 3])
c=colorbar;

title(c,'m/s')
if strcmp(varnm, 'stress')
    title(c,'MPa')
    %caxis([63 81.24])
end
if strcmp(varnm, 'rate')
    title(c,'m/s')
    %caxis([0 3])
end
%colormap jet
xlabel('X (km)')
ylabel('T (sec)')
set(gca,'FontSize',12)
%fnm = ['fault_x_t_mesh_asymm_mixed_',varnm];
%fnm = ['fault_x_t_mesh_',key1,'_',varnm];
%print( '-r300', '-dpng', fnm)
