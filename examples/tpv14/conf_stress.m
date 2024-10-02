%clc
clear
%close all

addpath(genpath('../../mscripts'));
myconstants;

par = ReadYaml('parameters.yaml');
nproc = par.nproc;
mesh_dir = par.mesh_dir;
%data_dir = par.data_dir;
ForcedRup = par.ForcedRup;

for iproc = 0:nproc-1

    fnm_out = [mesh_dir,'/meshVar',num2str(iproc,'%06d'),'.nc'];

    node = ncread(fnm_out,'node');
    elem = ncread(fnm_out,'elem');
    sigma0 = ncread(fnm_out,'sigma0');
    tau0 = ncread(fnm_out,'tau0');
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


    for ief = 1:Nelem_fault
        ie = fault2wave(ief);
        % calculate normal vectors
        for is = 1:Nfaces
            j1=is;j2=is+1;
            if (j1==Nfaces)
                j2=1;
            end
            va = node(:,elem(j1,ie));
            vb = node(:,elem(j2,ie));
            vec_m = vb-va;
            vec_m = vec_m/norm(vec_m);
            vec_n = [vec_m(2),-vec_m(1)];

            if (bctype(is,ie)>=BC_FAULT)
                xc = mean(node(1,elem(FtoV(is,:),ie)));
                yc = mean(node(2,elem(FtoV(is,:),ie)));

                for i = 1:2
                    j = FtoV(is,i);
                    x = node(1,elem(j,ie));
                    y = node(2,elem(j,ie));

                    sigma = -120e6;
                    tau = -70e6;
                    if (ForcedRup == 0)
                    if (abs(x+8e3)<1.5e3)
                        tau = -81.6e6;
                    end
                    end

                    %if (vec_n(2)<0)
                        %tau = -tau;
                        %sigma = -sigma;
                    %end

                    sigma0(i,is,ief) = sigma;
                    tau0(i,is,ief) = tau;

                    mu_s(i,is,ief) = 0.677;
                    mu_d(i,is,ief) = 0.525;
                    Dc(i,is,ief) = 0.4;
                    C0(i,is,ief) = 0;
                end

            end
        end
    end


    ncid = netcdf.open(fnm_out,'WRITE');

    var1 = netcdf.inqVarID(ncid,'sigma0');
    var2 = netcdf.inqVarID(ncid,'tau0');
    var4 = netcdf.inqVarID(ncid,'mu_s');
    var5 = netcdf.inqVarID(ncid,'mu_d');
    var6 = netcdf.inqVarID(ncid,'Dc');
    var7 = netcdf.inqVarID(ncid,'C0');

    if Nelem_fault > 0
    netcdf.putVar(ncid,var1,sigma0);
    netcdf.putVar(ncid,var2,tau0);
    netcdf.putVar(ncid,var4,mu_s);
    netcdf.putVar(ncid,var5,mu_d);
    netcdf.putVar(ncid,var6,Dc);
    netcdf.putVar(ncid,var7,C0);
    end

    netcdf.close(ncid);


end
