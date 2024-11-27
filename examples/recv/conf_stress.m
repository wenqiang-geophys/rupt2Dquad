%clc
clear
%close all

%%%%% Change this:
tau_2 = 63e6;   % Pa
%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('../../mscripts'));
myconstants;

par = ReadYaml('parameters.yaml');
nproc = par.nproc;
mesh_dir = par.mesh_dir;
ForcedRup = par.ForcedRup;


tau = []; sig = []; xf = []; mu_1 = []; mu_2 = [];

for iproc = 0:nproc-1

    fnm_out = [mesh_dir,'/meshVar',num2str(iproc,'%06d'),'.nc'];

    node = ncread(fnm_out,'node');
    elem = ncread(fnm_out,'elem');
    Sxx  = ncread(fnm_out,'Sxx0');
    Syy  = ncread(fnm_out,'Syy0');
    Sxy  = ncread(fnm_out,'Sxy0');
    mu_s = ncread(fnm_out,'mu_s');
    mu_d = ncread(fnm_out,'mu_d');
    Dc   = ncread(fnm_out,'Dc'  );
    C0   = ncread(fnm_out,'C0'  );
    bctype = ncread(fnm_out,'bctype');
    fault2wave = ncread(fnm_out,'fault2wave');

    Nelem = size(elem,2);
    Nnode = size(node,2);
    Nelem_fault = numel(fault2wave);
    fprintf('iproc=%d,Nelem=%d,Nnode=%d,Nelem_fault=%d\n', ...
        iproc,Nelem,Nnode,Nelem_fault);

    X = Sxx * 0;
    for ief = 1:Nelem_fault
    %for ie = 1:Nelem
        ie = fault2wave(ief);
        for is = 1:Nfaces
            if (bctype(is,ie)>=BC_FAULT)
                %ief = ief + 1;
                xc = mean(node(1,elem(FtoV(is,:),ie)));
                yc = mean(node(2,elem(FtoV(is,:),ie)));

                for i = 1:2
                    j = FtoV(is,i);
                    x = node(1,elem(j,ie));
                    y = node(2,elem(j,ie));

                    sxx = 0e6;
                    syy = -120e6;
                    sxy = 70e6;
                    if (ForcedRup == 0)
                    if (abs(x)<1.5e3)
                        sxy = 81.6e6;
                    end
                    end

                    %if (abs(x-7.5e3)<1.5e3)
                    %    sxy = 60e6;
                    %end
                    %if (abs(x+7.5e3)<1.5e3)
                    %    sxy = 80e6;
                    %end

                    if x > 10e3
                        sxy = tau_2;
                    end
                    if x < 0e3
                        sxy = tau_2;
                    end

                    X(i,is,ief) = x;

                    Sxx(i,is,ief) = sxx;
                    Syy(i,is,ief) = syy;
                    Sxy(i,is,ief) = sxy;

                    mu_s1 = 0.677;
                    mu_d1 = 0.525;

                    if x < 0
                        mu_s1 = mu_d1;
                    end

                    mu_s(i,is,ief) = mu_s1;   % static friction
                    mu_d(i,is,ief) = mu_d1;   % dynamic friciton
                    Dc(i,is,ief) = 0.4;       % slip-weakening distance
                    C0(i,is,ief) = 0;

                    tau = [tau,sxy];
                    sig = [sig,syy];
                    xf = [xf,x];
                    mu_1 = [mu_1,mu_s1];
                    mu_2 = [mu_2,mu_d1];
                end

            end
        end
    end

    


    ncid = netcdf.open(fnm_out,'WRITE');

    var1 = netcdf.inqVarID(ncid,'Sxx0');
    var2 = netcdf.inqVarID(ncid,'Syy0');
    var3 = netcdf.inqVarID(ncid,'Sxy0');
    var4 = netcdf.inqVarID(ncid,'mu_s');
    var5 = netcdf.inqVarID(ncid,'mu_d');
    var6 = netcdf.inqVarID(ncid,'Dc');
    var7 = netcdf.inqVarID(ncid,'C0');

    if (Nelem_fault > 0)
    netcdf.putVar(ncid,var1,Sxx);
    netcdf.putVar(ncid,var2,Syy);
    netcdf.putVar(ncid,var3,Sxy);
    netcdf.putVar(ncid,var4,mu_s);
    netcdf.putVar(ncid,var5,mu_d);
    netcdf.putVar(ncid,var6,Dc);
    netcdf.putVar(ncid,var7,C0);
    end

    netcdf.close(ncid);


end

% check initial settings for stresses and frictions
[xf,idx] = sort(xf);
tau = tau(idx);
mu_1 = mu_1(idx);
mu_2 = mu_2(idx);

figure
plot(xf,tau,'x-','LineWidth',1)
hold on
plot(xf,-sig,'o-','LineWidth',1)
xlabel('X (m)')
legend('\tau','-\sigma','Location','best','FontSize',20)
ylabel('Initial Stress (Pa)')
%print -dpng -r150 stress

figure
plot(xf,mu_1,'x-','LineWidth',1)
hold on
plot(xf,mu_2,'o-','LineWidth',1)
xlabel('X (m)')
legend('\mu_s','\mu_d','Location','best','FontSize',20)
ylabel('Friction')
%print -dpng -r150 friction
