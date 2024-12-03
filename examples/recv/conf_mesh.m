clc
clear
%close all

addpath(genpath('../../mscripts'))
myconstants;

flag_plot_mesh = 1;

pOrder = 3; % should be consistent with Makefile

par = ReadYaml('parameters.yaml');
nproc = par.nproc;
data_dir = par.data_dir;

fnm = 'tpv5_2d.exo';
try
    coord = ncread(fnm,'coord');
catch
    coordx = ncread(fnm,'coordx');
    coordy = ncread(fnm,'coordy');
    coordz = ncread(fnm,'coordz');
    coord = [coordx,coordy,coordz];
end
node = coord(:,1:2);
Nnode = size(node,1);

quad_list = [1];
 
%quad_list = [10];
elem = [];
for i = 1:length(quad_list)
    elem1 = ncread(fnm,['connect',num2str(quad_list(i))]);
    elem = cat(1, elem, elem1);
end
[~,Nelem] = size(elem);

elem = ccw_sort(node,elem'); elem = elem';

if 0
    % use gmsh
    fault_bdr_list = [1];
    elem_fault = [];
    bc1 = zeros(4,Nelem);
    for i = 1:length(fault_bdr_list)
        elem1 = ncread(fnm,['connect',num2str(fault_bdr_list(i))]);
        elem_fault = cat(2,elem_fault,elem1);
        bc2 = set_bctype_from_curve(elem,elem1,BC_FAULT+(i-1));
        bc1 = bc1 + bc2;
    end
    fnodes = elem_fault;
else
    % use cubit
    fnodes = ncread(fnm,['node_ns1']);
    bc1 = set_bctype_from_nodes(elem,fnodes,BC_FAULT+(i-1));
end
bctype = bc1;

%elem = elem'; % Nelem x 4

mesh.node = node;
mesh.elem = elem';
mesh = init_mesh(mesh,1);

elem = mesh.elem;
node = mesh.node;

fluxtype = zeros(4,Nelem);
fluxtype = set_fluxtype_quad(elem, fnodes);
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

% body recvs

recv_x = linspace(0,40e3,21);
recv_y = linspace(1e3,1e3,21);
nrecv = length(recv_x);

if 0
nrecv = 0;
recv_x = [];
recv_y = [];
end

recv_i = zeros(nrecv,1);
recv_j = zeros(nrecv,1);
recv_ie = zeros(nrecv,1);
recv_refx = zeros(nrecv,1);
recv_refy = zeros(nrecv,1);
NGLL = pOrder + 1;
for i = 1:nrecv
    x = recv_x(i);
    y = recv_y(i);
    r = sqrt((mesh.x(:)-x).^2+(mesh.y(:)-y).^2);
    [~,idx]=min(r);
    %[recv_i(i),recv_j(i),recv_ie(i)] = ind2sub([NGLL,NGLL,nelem],idx);
    for ie = 1:mesh.Nelem
        VX = mesh.node(1,mesh.elem(:,ie));
        VY = mesh.node(2,mesh.elem(:,ie));
        in = inpolygon(x,y,VX,VY);
        if in
            recv_ie(i) = ie;
            [recv_refx(i),recv_refy(i)]=xy2rs(recv_x(i),recv_y(i),VX,VY);
            break;
        end
    end
    %[recv_refx(i),recv_refy(i)]=xy2rs(recv(i,1),recv(i,2),VX,VY);
end

mesh.body_nrecv = nrecv;
mesh.body_recv_x = recv_x;
mesh.body_recv_y = recv_y;
mesh.body_recv_i = recv_i;
mesh.body_recv_j = recv_j;
mesh.body_recv_ie = recv_ie;
mesh.body_recv_refx = recv_refx;
mesh.body_recv_refy = recv_refy;

system('mkdir -p data');
db = gen_mesh_mpi_v6(mesh, nproc);

%% check local node and elem ...

check_mesh(mesh)
