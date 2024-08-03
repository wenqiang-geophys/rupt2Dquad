function [ x, y, t ] = gather_wave_coord( data_dir, nproc)

x = []; y = [];
for iproc = 0:nproc-1
    fnm = sprintf('%s/wave_mpi%06d.nc',data_dir,iproc);
    if (exist(fnm,'file'))
        t = ncread(fnm, 'time');
        x1 = ncread(fnm, 'x');
        y1 = ncread(fnm, 'y');
        x = cat(3,x,x1);
        y = cat(3,y,y1);
    end
end

end
