clc
clear

fnm = 'data_tpv5.nc';

%% draw snapshot of wavefield
x = ncread(fnm, '/wave/x')*1e-3; % m to km
y = ncread(fnm, '/wave/y')*1e-3;
z = x * 0;
tri = ncread(fnm, '/wave/tri');
t = ncread(fnm, '/wave/t');

t1 = 3;
[~,it] = min(abs(t-t1));
Vx = squeeze(ncread(fnm, '/wave/Vx', [1,1,1,it], [Inf,Inf,Inf,1]));
Vy = squeeze(ncread(fnm, '/wave/Vy', [1,1,1,it], [Inf,Inf,Inf,1]));

figure
v = Vx;
trisurf(tri,x,y,z,v);view(2);shading interp
axis image;colorbar;colormap jet;%colormap seismic
title(['Vx snapshot at ',num2str(t(it)),' sec'])
xlabel('X (km)'); ylabel('Y (km)'); %axis([-1 1 -1 1]*50e3)

%% draw fault
x = ncread(fnm, '/fault/x')*1e-3; % convert to km
y = ncread(fnm, '/fault/y')*1e-3;
t = ncread(fnm, '/fault/t');
nt = length(t);

mu_s = ncread(fnm, '/fault/mu_s');
mu_d = ncread(fnm, '/fault/mu_d');
Dc = ncread(fnm, '/fault/Dc');
C0 = ncread(fnm, '/fault/C0');

figure
plot(x(:),mu_s(:),'-','LineWidth',1)
hold on
plot(x(:),mu_d(:),'--','LineWidth',1)
hold off
xlabel('X (km)')
ylabel('Friction')
legend('\mu_f','\mu_d','FontSize',20,'Location','Best')

% sliprate, slip, tau, sigma
sliprate = ncread(fnm, '/fault/sliprate');
slip = ncread(fnm, '/fault/slip');
tau = ncread(fnm, '/fault/tau');
sigma = ncread(fnm, '/fault/sigma');

% plot initial stresses
figure
v = -tau(:,:,1)*1e-6;
plot(x(:),v(:),'-','LineWidth',1)
hold on
v = -sigma(:,:,1)*1e-6;
plot(x(:),v(:),'--','LineWidth',1)
hold off
xlabel('X (km)')
ylabel('Initial Stresses (Pa)')
legend('\tau','-\sigma','FontSize',20,'Location','Best')

v = reshape(sliprate,[],nt);

figure
pcolor(x(:),t(:),v')
shading interp;colorbar;colormap jet
xlabel('X (km)')

figure
X1 = 6;
[~,ix] = min(abs(x(:)-X1));
plot(t,v(ix,:),'LineWidth',1)
xlabel('Time (sec)')
title(['at X = ',num2str(X1), ' km'])

figure
t1 = 4;
[~,it] = min(abs(t(:)-t1));
plot(x(:),v(:,it),'LineWidth',1)
xlabel('X (km)')
title(['at T = ',num2str(t1), ' sec'])