function [ v,t ] = gather_wave_snap( data_dir, nproc, varnm, itstep )

v = [];
for iproc = 0:nproc-1
    fnm = sprintf('%s/wave_mpi%06d.nc',data_dir,iproc);
    if (exist(fnm,'file'))
        v1 = ncread(fnm, varnm, [1,1,1,itstep], [Inf,Inf,Inf,1], [1,1,1,1]);
        v = cat(3,v,v1);
        t = ncread(fnm, 'time', [itstep], [1], [1]);
    end
end

end
