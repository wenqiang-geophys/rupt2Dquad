clear
close all
clc

par = ReadYaml('parameters.yaml');
nproc = par.nproc;
data_dir = par.data_dir;

[x,y,tvec] = gather_wave_coord(data_dir, nproc);
nt = length(tvec);
[NGLL,~,Nelem] = size(x);

% fast plot using trisurf
tri = [];
for i = 1:Nelem
    tri1 = delaunay(x(:,:,i),y(:,:,i)) + (i-1)*NGLL*NGLL;
    tri = cat(1,tri,tri1);
end

Vx = [];Vy = [];
for iproc = 0:nproc-1
    fnm = sprintf('%s/wave_mpi%06d.nc',data_dir,iproc);
    if (exist(fnm,'file'))
        v1 = ncread(fnm, 'Vx', [1,1,1,1], [Inf,Inf,Inf,Inf], [1,1,1,1]);
        Vx = cat(3,Vx,v1);
        v1 = ncread(fnm, 'Vy', [1,1,1,1], [Inf,Inf,Inf,Inf], [1,1,1,1]);
        Vy = cat(3,Vy,v1);
        t = ncread(fnm, 'time');
    end
end

Vx = [];Vy = [];
rate = [];stress=[];
xf = []; yf = [];
for iproc = 0:nproc-1
    fnm = sprintf('%s/wave_mpi%06d.nc',data_dir,iproc);
    if (exist(fnm,'file'))
        v1 = ncread(fnm, 'Vx'); Vx = cat(3,Vx,v1);
        v1 = ncread(fnm, 'Vy'); Vy = cat(3,Vy,v1);
        t = ncread(fnm, 'time');
    end
    fnm = sprintf('%s/fault_mpi%06d.nc',data_dir,iproc);
    if (exist(fnm,'file'))
        tf = ncread(fnm, 'time');
        v1 = ncread(fnm, 'rate'); rate = cat(2,rate,v1);
        v1 = ncread(fnm, 'stress'); stress = cat(2,stress,v1);
        v1 = ncread(fnm, 'x'); xf = cat(2,xf,v1);
        v1 = ncread(fnm, 'y'); yf = cat(2,yf,v1);
    end
end
[~,Nelem_fault,Nt_fault] = size(rate);

fnm = 'data_tpv5.nc';
system(['rm -r ',fnm]);
ncid = netcdf.create(fnm,'netcdf4');%'CLOBBER'

waveGrpId = netcdf.defGrp(ncid,"wave");
faultGrpId = netcdf.defGrp(ncid,"fault");

dimid2 = netcdf.defDim(ncid,'two',2);
dimid3 = netcdf.defDim(ncid,'three',3);
dimid4 = netcdf.defDim(ncid,'four',4);
dimid_N = netcdf.defDim(ncid,'Nfp',NGLL);

dimid_elem = netcdf.defDim(waveGrpId,'Nelem',Nelem);
dimid_tris = netcdf.defDim(waveGrpId,'Ntri',size(tri,1));
dimid_t = netcdf.defDim(waveGrpId,'Nt',length(t));

dimid_elem2 = netcdf.defDim(faultGrpId,'Nelem',Nelem_fault);
dimid_t2 = netcdf.defDim(faultGrpId,'Nt',length(tf));

varid1 = netcdf.defVar(waveGrpId,'x','NC_DOUBLE',[dimid_N,dimid_N,dimid_elem]);
varid2 = netcdf.defVar(waveGrpId,'y','NC_DOUBLE',[dimid_N,dimid_N,dimid_elem]);
varid3 = netcdf.defVar(waveGrpId,'tri','NC_INT',[dimid_tris,dimid3]);
varid4 = netcdf.defVar(waveGrpId,'t','NC_DOUBLE',[dimid_t]);
varid5 = netcdf.defVar(waveGrpId,'Vx','NC_FLOAT',[dimid_N,dimid_N,dimid_elem,dimid_t]);
varid6 = netcdf.defVar(waveGrpId,'Vy','NC_FLOAT',[dimid_N,dimid_N,dimid_elem,dimid_t]);
netcdf.endDef(waveGrpId);

varid7 = netcdf.defVar(faultGrpId,'x','NC_DOUBLE',[dimid_N,dimid_elem2]);
varid8 = netcdf.defVar(faultGrpId,'y','NC_DOUBLE',[dimid_N,dimid_elem2]);
varid9 = netcdf.defVar(faultGrpId,'t','NC_DOUBLE',[dimid_t2]);
varid10 = netcdf.defVar(faultGrpId,'rate','NC_FLOAT',[dimid_N,dimid_elem2,dimid_t2]);
varid11 = netcdf.defVar(faultGrpId,'stress','NC_FLOAT',[dimid_N,dimid_elem2,dimid_t2]);
netcdf.endDef(faultGrpId);

netcdf.putVar(waveGrpId,varid1,x);
netcdf.putVar(waveGrpId,varid2,y);
netcdf.putVar(waveGrpId,varid3,tri);
netcdf.putVar(waveGrpId,varid4,t);
netcdf.putVar(waveGrpId,varid5,Vx);
netcdf.putVar(waveGrpId,varid6,Vy);

netcdf.putVar(faultGrpId,varid7,xf);
netcdf.putVar(faultGrpId,varid8,yf);
netcdf.putVar(faultGrpId,varid9,tf);
netcdf.putVar(faultGrpId,varid10,rate);
netcdf.putVar(faultGrpId,varid11,stress);

netcdf.close(ncid);
