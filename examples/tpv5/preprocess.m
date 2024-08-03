clc
clear
%close all

addpath(genpath('../../mscripts'))

myconstants;

flag_plot_mesh = 1;

nproc = 4;

TPV14 = 0;
TPV10 = 0;
TPV5  = 1;

order = 1;
%[node,elem] = read_mesh2(basename);

fnm = 'tpv5_2d.exo';
elem = ncread(fnm,'connect6');
elem_fault = ncread(fnm,'connect1');
fnm = 'tpv5_2d_symm.exo';
elem = ncread(fnm,'connect10');
elem_fault = ncread(fnm,'connect1');
coord = ncread(fnm,'coord');
node = coord(:,1:2);
elem = elem';

%return

L = max(node(:,1))-min(node(:,2));
for i = 1:size(node,1)
x = node(i,1);
y = node(i,1);
end

mesh.node = node*1e3;
mesh.elem = elem;
mesh = init_mesh(mesh,order);

mesh.node = mesh.node * 1e-0;
mesh.x = mesh.x * 1e-0;
mesh.y = mesh.y * 1e-0;
elem = mesh.elem;
node = mesh.node;
EtoE = mesh.EtoE;
EtoF = mesh.EtoF;

nelem = mesh.nelem;
nnode = mesh.nnode;

bctype = zeros(4,mesh.nelem);
elemtype = zeros(1,mesh.nelem);
elemtype(:) = ELEM_SOLID;
vp = zeros(1,mesh.nelem);
vs = zeros(1,mesh.nelem);
rho = zeros(1,mesh.nelem);
vp(:) = 6000;
vs(:) = 3464;
rho(:) = 2670;
bctype = set_bctype_from_curve(elem,elem_fault);
%for ie = 1:mesh.nelem
%    for is = 1:4
%        is1 = is + 1;
%        if (is1 > 4)
%            is1 = is1-4;
%        end
%        x = node(1,elem([is,is1],ie));
%        y = node(2,elem([is,is1],ie));
%
%        if (norm(y) < 1e1 && min(x) >= -15e3 && max(x)<= 15e3)
%            bctype(is,ie) = BC_FAULT;
%        end
%    end
%end

fnodes = unique(sort(elem_fault(:)));
fluxtype = set_fluxtype_quad(elem,fnodes);



mesh.bctype = bctype;
mesh.fluxtype = fluxtype;
mesh.elemtype = elemtype;
mesh.rho = rho;
mesh.vp = vp;
mesh.vs = vs;
mesh.zp = rho.*vp;
mesh.zs = rho.*vs;

a = [cosd(60),-sind(60)];
recv = [...
    a*0;
    a*1.5;
    a*3;
    a*4.5;
    a*7.5;
    a*12]*1e3;
if TPV14
a = [cosd(30),-sind(30)];
recv = [...
    -2,0;
    2,0;
    5,0;
    9,0;
    a*2;
    a*5;
    a*9;]*1e3;
end
if TPV5
   recv = [0 0;
       -4.5 0;
       -7.5 0;
       -12 0;
       -12.5 0;
       4.5 0;
       7.5 0;
       12 0;
       12.5 0;
       ]*1e3;

    recv = [cosd(30),sind(30);
    cosd(60),sind(60);
    cosd(90),sind(90);
    cosd(120),sind(120);
    cosd(150),sind(150);
    ]*15e3;
    %recv=[...
    %0,0;
    %5,0;
    %7.5,2.5;
    %12.0,0];
    recv = [0,5]*1e3;
%     recv = [...
%     cosd(60),sind(60);
%     cosd(90),sind(90);
%     cosd(120),sind(120);
%     ]*25e3-[0,20e3];
    recv = [-2e3,0;
    0e3,0;
    0,5e3];
    recv = [0,0;
    3e3 600];
    recv = [0 0;
    3e3 0
    4.5e3 0;
    7.5e3 0];
    recv = [0,0];
    %recv=[...
    %-10,0;
    %-5,0;
    %0,0;
    %2, .8*sin(2/L*2*pi*2);
    %4, .8*sin(4/L*2*pi*2);
    %6, .8*sin(6/L*2*pi*2);
    %8, .8*sin(8/L*2*pi*2);
    %10,.8*sin(10/L*2*pi*2);
    %12,.8*sin(12/L*2*pi*2);
    %14,.8*sin(14/L*2*pi*2);
    %16,.8*sin(16/L*2*pi*2);
    %13,0]*1e3;
    %7.5,5-sqrt(1.5^2+5^2);
    %recv = [0 15]*1e3;
end

nrecv = size(recv,1)
recv_idx = zeros(nrecv,1);
recv_ie = zeros(nrecv,1);
recv_refx = zeros(nrecv,1);
NGLL = order + 1;
for i = 1:nrecv
    x = recv(i,1);
    y = recv(i,2);
    r = sqrt((mesh.x(:)-x).^2+(mesh.y(:)-y).^2);
    [~,idx]=min(r);
    [recv_i,recv_j,recv_ie(i)] = ind2sub([NGLL,NGLL,nelem],idx);
    %%%%recv_is = find(bctype(:,recv_ie(i))==-3,1);
    recv_is = find(bctype(:,recv_ie(i))>=BC_FAULT,1);
    if ((recv_is==1) | (recv_is == 3))
        recv_idx(i) = recv_i;
    else
        recv_idx(i) = recv_j;
    end
    if (recv_is==1)
        x0=mesh.x(1,  1,recv_ie(i));
        y0=mesh.y(1,  1,recv_ie(i));
        x1=mesh.x(end,1,recv_ie(i));
        y1=mesh.y(end,1,recv_ie(i));
    elseif (recv_is==2)
        x0=mesh.x(end,  1,recv_ie(i));
        y0=mesh.y(end,  1,recv_ie(i));
        x1=mesh.x(end,end,recv_ie(i));
        y1=mesh.y(end,end,recv_ie(i));
    elseif (recv_is==3)
        x0=mesh.x(1,  end,recv_ie(i));
        y0=mesh.y(1,  end,recv_ie(i));
        x1=mesh.x(end,end,recv_ie(i));
        y1=mesh.y(end,end,recv_ie(i));
    elseif (recv_is==4)
        x0=mesh.x(1,  1,recv_ie(i));
        y0=mesh.y(1,  1,recv_ie(i));
        x1=mesh.x(1,end,recv_ie(i));
        y1=mesh.y(1,end,recv_ie(i));
    end

    r1 = sqrt((x1-x0)^2+(y1-y0)^2);
    rr = sqrt((x -x0)^2+(y -y0)^2);
    recv_refx(i) = rr/r1;

end

mesh.nrecv = nrecv;
mesh.recv_i = recv_idx;
mesh.recv_ie = recv_ie;
mesh.recv_refx = recv_refx;

%% body recvs
%recv = [
%    -1,0;
%    +1,0;
%    -2,0;
%    +2,0;
%    -3,0;
%    +3,0;
%    -0.5,-0.3;
%    13,0]*1e3;
%    %7.5,5-sqrt(1.5^2+5^2);
%    %recv = [0 15]*1e3;
%%end
%
%nrecv = size(recv,1);
%recv_idx = zeros(nrecv,1);
%recv_ie = zeros(nrecv,1);
%recv_refx = zeros(nrecv,1);
%NGLL = order + 1;
%for i = 1:nrecv
%    x = recv(i,1);
%    y = recv(i,2);
%    r = sqrt((mesh.x(:)-x).^2+(mesh.y(:)-y).^2);
%    [~,idx]=min(r);
%    [recv_i,recv_j,recv_ie(i)] = ind2sub([NGLL,NGLL,nelem],idx);
%    %%%%recv_is = find(bctype(:,recv_ie(i))==-3,1);
%    recv_is = find(bctype(:,recv_ie(i))==-3,1);
%    if ((recv_is==1) | (recv_is == 3))
%        recv_idx(i) = recv_i;
%    else
%        recv_idx(i) = recv_j;
%    end
%    if (recv_is==1)
%        x0=mesh.x(1,  1,recv_ie(i));
%        y0=mesh.y(1,  1,recv_ie(i));
%        x1=mesh.x(end,1,recv_ie(i));
%        y1=mesh.y(end,1,recv_ie(i));
%    elseif (recv_is==2)
%        x0=mesh.x(end,  1,recv_ie(i));
%        y0=mesh.y(end,  1,recv_ie(i));
%        x1=mesh.x(end,end,recv_ie(i));
%        y1=mesh.y(end,end,recv_ie(i));
%    elseif (recv_is==3)
%        x0=mesh.x(1,  end,recv_ie(i));
%        y0=mesh.y(1,  end,recv_ie(i));
%        x1=mesh.x(end,end,recv_ie(i));
%        y1=mesh.y(end,end,recv_ie(i));
%    elseif (recv_is==4)
%        x0=mesh.x(1,  1,recv_ie(i));
%        y0=mesh.y(1,  1,recv_ie(i));
%        x1=mesh.x(1,end,recv_ie(i));
%        y1=mesh.y(1,end,recv_ie(i));
%    end
%
%    r1 = sqrt((x1-x0)^2+(y1-y0)^2);
%    rr = sqrt((x -x0)^2+(y -y0)^2);
%    recv_refx(i) = rr/r1;
%
%end
%
%mesh.nrecv = nrecv;
%mesh.recv_i = recv_idx;
%mesh.recv_ie = recv_ie;
%mesh.recv_refx = recv_refx;

% body recvs
recv = [
    -1,0;
    +1,0;
    -2,0;
    +2,0;
    -3,0;
    +3,0;
    -0.5,-0.3;
    +0.5,-0.3;
    -1,-0.3;
    +1,-0.3;
    ]*1e3;

recv = [0,0];

if TPV14
recv = [...
    2,-0.6;
    5,-1.4 ;
    8,-2.3 ;
    -2, -3;
    2, -4.2;
    5, -5.9;
    8, -7.6;
    -2, 3;
    2, 3;
    5, 3;
    8, 3;
    ]*1e3;
end

if TPV5
recv = [...
-12,-3;
+12,-3;
]*1e3;
recv = [0,0];
end

nrecv = size(recv,1)
recv_i = zeros(nrecv,1);
recv_j = zeros(nrecv,1);
recv_ie = zeros(nrecv,1);
recv_refx = zeros(nrecv,1);
recv_refy = zeros(nrecv,1);
NGLL = order + 1;
for i = 1:nrecv
    x = recv(i,1);
    y = recv(i,2);
    r = sqrt((mesh.x(:)-x).^2+(mesh.y(:)-y).^2);
    [~,idx]=min(r);
    %[recv_i(i),recv_j(i),recv_ie(i)] = ind2sub([NGLL,NGLL,nelem],idx);
    for ie = 1:mesh.nelem
        VX = mesh.node(1,mesh.elem(:,ie));
        VY = mesh.node(2,mesh.elem(:,ie));
        in = inpolygon(x,y,VX,VY);
        if in
            recv_ie(i) = ie;
            [recv_refx(i),recv_refy(i)]=xy2rs(recv(i,1),recv(i,2),VX,VY);
            break;
        end
    end
    %[recv_refx(i),recv_refy(i)]=xy2rs(recv(i,1),recv(i,2),VX,VY);
end

mesh.body_nrecv = nrecv;
mesh.body_recv_i = recv_i;
mesh.body_recv_j = recv_j;
mesh.body_recv_ie = recv_ie;
mesh.body_recv_refx = recv_refx;
mesh.body_recv_refy = recv_refy;

if 0
   % check searched reference coordinates
   for i = 1:nrecv
       x = recv(i,1);
       y = recv(i,2);
       VX = mesh.node(1,mesh.elem(:,mesh.body_recv_ie(i)));
       VY = mesh.node(2,mesh.elem(:,mesh.body_recv_ie(i)));
       figure
       plot(VX([1:4,1]),VY([1:4,1]),'o-')
       hold on
       plot(x,y,'rv')
       [x1,y1]=rs2xy(mesh.body_recv_refx(i),mesh.body_recv_refy(i),VX,VY);
       plot(x1,y1,'b*')
       hold off
   end
end


db = gen_mesh_mpi_v5(mesh,nproc);


%% check local node and elem ...
if flag_plot_mesh
figure


%subplot(121)
%for iproc = 1:nproc
%    idx = find(part==iproc);
elem1 = elem(1:4,:);
node1 = node;
for ie = 1:size(elem1,2)
    patch(node1(1,elem1(:,ie)),node1(2,elem1(:,ie)),[1 1 1 1]*-1);%,'FaceColor','none');
    for is = 1:4
        if (bctype(is,ie) >= BC_FAULT)
            is1 = is + 1;
            if (is1 > 4)
                is1 = is1-4;
            end
            xx = node(1,elem([is,is1],ie));
            yy = node(2,elem([is,is1],ie));
            hold on;
            plot(xx,yy,'r','linewidth',2);
        end
        if (bctype(is,ie) == BC_FREE)
            is1 = is + 1;
            if (is1 > 4)
                is1 = is1-4;
            end
            xx = node(1,elem([is,is1],ie));
            yy = node(2,elem([is,is1],ie));
            hold on;
            plot(xx,yy,'g','linewidth',2);
        end
    end
    for is = 1:4
        if (fluxtype(is,ie) == 1)
            c='b';
            is1 = is + 1;
            if (is1 > 4)
                is1 = is1-4;
            end
            xx = node(1,elem([is,is1],ie));
            yy = node(2,elem([is,is1],ie));
            hold on;
            plot(xx,yy,c,'linewidth',1.5);
        end
    end
end
%end
axis image;
hold off
colormap default
end

if 0
subplot(122)

colors=meshcolormap(nproc);
for iproc = 1:nproc
    elem1 = db(iproc).elem;
    node1 = db(iproc).node*1e-3;
    for i = 1:size(elem1,2)
        %patch(node1(1,elem1(:,i)),node1(2,elem1(:,i)),[1 1 1 1]+iproc);
        plot(node1(1,elem1([1:4,1],i)),...
            node1(2,elem1([1:4,1],i)),...
            'color',colors(iproc,:));
        hold on
    end
end
axis image;
xlabel('Horizontal Distance (km)')
ylabel('Depth (km)')
set(gcf,'renderer','Painters')
%print('-depsc','mesh_partition_8')
%hold on
colormap(meshcolormap(nproc))
end

%%
%if 0
%for iproc = 1
%    xx = zeros(order+1,db(iproc).mpi_ne,db(iproc).mpi_nn);
%    yy = zeros(order+1,db(iproc).mpi_ne,db(iproc).mpi_nn);
%    for i = 1:db(iproc).mpi_nn
%        for ie = 1:db(iproc).mpi_ne
%            je = db(iproc).mpi_connection(i,ie,1);
%            fc = db(iproc).mpi_connection(i,ie,2);
%            if (je > 0)
%                xx(:,ie,i) = get_face(db(iproc).x(:,:,je),fc);
%                yy(:,ie,i) = get_face(db(iproc).y(:,:,je),fc);
%                plot(xx(1,ie,i),yy(1,ie,i),'o','linewidth',2); hold on
%            end  
%        end
%    end
%end
%
%
%axis equal
%hold on
%%return
% 
%for iproc = 1
%    for ie = 1:db(iproc).nelem
%        for is = 1:4
%            neigh = db(iproc).neighbor(is,ie);
%            if (neigh == -1)
%                mpi_e = db(iproc).mpi_ibool(is,ie);
%                mpi_n = db(iproc).mpi_interface(4,is,ie);
%                if mpi_n == 1
%                xx1 = xx(1,mpi_e,mpi_n);
%                yy1 = yy(1,mpi_e,mpi_n);
%                plot(xx1,yy1,'rx','linewidth',2)
%                hold on
%                end
%            end
%        end 
%    end  
%end
%end
