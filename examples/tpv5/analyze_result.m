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
% Vx or Vy
v = squeeze(ncread(fnm, '/wave/Vx', [1,1,1,it], [Inf,Inf,Inf,1]));

figure
trisurf(tri,x,y,z,v);view(2);shading interp
axis image;colorbar;colormap jet
title(['Vx snapshot at ',num2str(t1),' sec'])
xlabel('X (m)'); ylabel('Y (m)'); axis([-1 1 -1 1]*15e3)

%% draw fault
x = ncread(fnm, '/fault/x');
y = ncread(fnm, '/fault/y');
t = ncread(fnm, '/fault/t');
nt = length(t);

x = x(:);
[x,idx] = sort(x);

% rate, slip, tau, sigma
v = ncread(fnm, '/fault/rate');
v = reshape(v,[],nt);
v = v(idx,:);

figure
pcolor(x(:),t(:),v')
shading interp;colorbar;colormap jet
xlabel('X (m)')

figure
X1 = 3e3;
[~,ix] = min(abs(x(:)-X1));
plot(t,v(ix,:),'LineWidth',1)
xlabel('Time (sec)')
title(['at X = ',num2str(X1), ' m'])

figure
t1 = 3;
[~,it] = min(abs(t(:)-t1));
plot(x,v(:,it),'LineWidth',1)
xlabel('Time (sec)')
title(['at T = ',num2str(t1), ' sec'])