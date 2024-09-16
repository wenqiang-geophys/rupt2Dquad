clc
clear
close all

addmypath;

if_plot_mesh = 1;

fnm_in = 'tpv5.exo';
fnm_out = 'mesh.nc';

%% check the tet_id, surf_id
% use: check_surface_id(fnm_in,surf_id);
tet_id = 24;
faultsurf_id = [17];
freesurf_id = 20;
%tet_id = 39;
%faultsurf_id = [29,30,31,32];
%freesurf_id = 35;
msh = build_mesh(fnm_in,tet_id,faultsurf_id,freesurf_id);
msh.node = msh.node * 1e3; % convert to m

%% write mesh.nc
disp('writing mesh ...')
write_mesh_nc(fnm_out,msh);

if if_plot_mesh
%% checking
% check fault connectivity
fault_tri = get_fault_connectivity(msh.elem',msh.fault2wave,msh.bctype);
figure
trisurf(fault_tri,...
    msh.node(1,:),...
    msh.node(2,:),...
    msh.node(3,:))
end
