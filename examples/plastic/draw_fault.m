clc
close all
clear

addpath(genpath('../../mscripts'))

par = ReadYaml('parameters.yaml');
nproc = par.nproc;
data_dir = par.data_dir;
data_dir = 'data_psi20'

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

    nt = length(t);
    %nt = 700;
    t = t(1:nt);
    v1 = v1(:,:,1:nt);

    x = cat(2,x,x1);
    y = cat(2,y,y1);
    v = cat(2,v,v1);
    end

end

nt = length(t);
nflt = length(x(:));
v = reshape(v,[nflt,nt]);
x = x(:)*1e-3;
y = y(:)*1e-3;
[x,idx] = sort(x);
v = v(idx,:);

if strcmp(varnm, 'tau')
    v = -v*1e-6;
end


figure

% for bending fault
faultBendAngle = 10;
% x(x>20) = 20+(x(x>20)-20)/cosd(faultBendAngle);

pcolor(x,t, (v'))
shading interp
c=colorbar;

title(varnm)

title(c,'m/s')
if strcmp(varnm, 'tau')
    title(c,'MPa')
    %caxis([63 81.24])
end
if strcmp(varnm, 'rate')
    %title(c,'log_{10} V [m/s]')
    title(c,'m/s')
    %ylabel(c,'slip rate')
    %caxis([0 10])
end
colormap  whitejet
xlabel('X (km)')
ylabel('Time (sec)')
