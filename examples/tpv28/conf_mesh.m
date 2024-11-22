clc
clear
%close all

addpath(genpath('../../mscripts'))
myconstants;

flag_plot_mesh = 1;

par = ReadYaml('parameters.yaml');
nproc = par.nproc;
data_dir = par.data_dir;

fnm = 'tpv28_2d.exo';
coord = ncread(fnm,'coord');
node = coord(:,1:2);
Nnode = size(node,1);

quad_list = [6];
elem = [];
for i = 1:length(quad_list)
    elem1 = ncread(fnm,['connect',num2str(quad_list(i))]);
    elem = cat(1, elem, elem1);
end
[~,Nelem] = size(elem);

elem = ccw_sort(node,elem'); elem = elem';

fault_bdr_list = [1];
%elem_fault = [];
bc1 = zeros(4,Nelem);
for i = 1:length(fault_bdr_list)
    elem1 = ncread(fnm,['connect',num2str(fault_bdr_list(i))]);
    bc2 = set_bctype_from_curve(elem,elem1,BC_FAULT+(i-1));
    bc1 = bc1 + bc2;
end
bctype = bc1;

%elem = elem'; % Nelem x 4

mesh.node = node;
mesh.elem = elem';
mesh = init_mesh(mesh,1);

elem = mesh.elem;
node = mesh.node;

fluxtype = zeros(4,Nelem);
elemtype = zeros(1,Nelem);
elemtype(:) = ELEM_SOLID;
vp = zeros(1,Nelem);
vs = zeros(1,Nelem);
rho = zeros(1,Nelem);
vp(:) = 6000;
vs(:) = 3464;
rho(:) = 2670;

mesh.bctype = bctype;
mesh.fluxtype = fluxtype;
mesh.elemtype = elemtype;
mesh.rho = rho;
mesh.vp = vp;
mesh.vs = vs;
%mesh.zp = rho.*vp;
%mesh.zs = rho.*vs;

system('mkdir -p data');
db = gen_mesh_mpi_v6(mesh, nproc);

%% check local node and elem ...

check_mesh(mesh)
