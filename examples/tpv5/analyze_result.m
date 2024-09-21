clc
clear

fnm = 'data_tpv5.nc';

%% draw snapshot of wavefield
x = ncread(fnm, '/wave/x');
y = ncread(fnm, '/wave/y');
z = x * 0;
tri = ncread(fnm, '/wave/tri');
t = ncread(fnm, '/wave/t');

t1 = 2;

[~,it] = min(abs(t-t1));

v = squeeze(ncread(fnm, '/wave/Vx', [1,1,1,it], [Inf,Inf,Inf,1]));

figure
trisurf(tri,x,y,z,v);view(2)
axis image
shading interp
colorbar
%colormap jet
colormap seismic

%% draw fault
x = ncread(fnm, '/fault/x');
y = ncread(fnm, '/fault/y');
t = ncread(fnm, '/fault/t');
nt = length(t);

x = x(:);
[x,idx] = sort(x);

v = ncread(fnm, '/fault/rate');
v = reshape(v,[],nt);
v = v(idx,:);

figure
pcolor(x(:),t(:),v')
shading interp
colorbar
