clc
close all
clear

nproc = 4;
varnm = 'stress';
%data_dir = 'data_500_20000_v1'
data_dir = 'data'
%data_dir = 'data_500_20000_llf_pOrder8'

key1 = 'symm_upwind'
data_dir = ['data_',key1];

x = []; y = []; v = [];
for i = 1:nproc

    fnm = [data_dir,'/fault_mpi',num2str(i-1,'%06d'),'.nc'];
    if (exist(fnm))
    x1 = ncread(fnm, 'x');
    y1 = ncread(fnm, 'y');
    v1 = ncread(fnm, varnm);
    t = ncread(fnm, 'time');
    end

    x = [x,x1];
    y = [y,y1];
    v = [v,v1];

end

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
%plot(x(:), y(:), 'x')
pcolor(x,t, v')
shading interp
%caxis([0 3])
c=colorbar;

title(c,'m/s')
if strcmp(varnm, 'stress')
    title(c,'MPa')
    caxis([63 81.24])
end
if strcmp(varnm, 'rate')
    title(c,'m/s')
    caxis([0 3])
end
colormap jet
xlabel('X (km)')
ylabel('T (sec)')
set(gca,'FontSize',12)
%fnm = ['fault_x_t_mesh_asymm_mixed_',varnm];
fnm = ['fault_x_t_mesh_',key1,'_',varnm];
print( '-r300', '-dpng', fnm)