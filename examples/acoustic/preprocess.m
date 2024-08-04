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
%basename = 'mesh/2layer_2'
%[node,elem] = read_mesh2(basename);

fnm = 'twolayer.exo';
elem1 = ncread(fnm,'connect8');
elem2 = ncread(fnm,'connect9');
elem = [elem1,elem2];

coord = ncread(fnm,'coord');
coord = coord * 1e3; % km to m
node = coord(:,1:2);
elem = elem';

%return;

%L = max(node(:,1))-min(node(:,2));
%for i = 1:size(node,1)
%x = node(i,1);
%y = node(i,1);
%if (x > 7.5e3 && x < 13.5e3)
%node(i,2) = y+800*sin(x/L*pi*2*2);
%end
%end

mesh.node = node;
mesh.elem = elem;
mesh = init_mesh(mesh,order);

mesh.node = mesh.node;
mesh.x = mesh.x * 1e-0;
mesh.y = mesh.y * 1e-0;
elem = mesh.elem;
node = mesh.node;
EtoE = mesh.EtoE;
EtoF = mesh.EtoF;

nelem = mesh.nelem;
nnode = mesh.nnode;

bctype = zeros(4,mesh.nelem);
fluxtype = zeros(4,mesh.nelem);
elemtype = zeros(1,mesh.nelem);
elemtype(:) = -999;
bctype(:,:) = -999;
vp = zeros(1,mesh.nelem);
vs = zeros(1,mesh.nelem);
rho = zeros(1,mesh.nelem);

seabed_depth = -6e3;
for ie = 1:mesh.nelem
%     x4 = node(1,elem(:,ie));
     y4 = node(2,elem(:,ie));
%
%     if (mean(y4)>0)
%         elemtype(ie) = ELEM_FLUID;
%     else
%         elemtype(ie) = ELEM_SOLID;
%     end

    for is = 1:4
        is1 = is + 1;
        if (is1 > 4)
            is1 = is1-4;
        end
        x = node(1,elem([is,is1],ie));
        y = node(2,elem([is,is1],ie));

        y2 = 1000*( cos(pi*x/5e3));
        y2(:) = seabed_depth;

        %if (mean(y)>0)
        if (mean(y-y2)>=0)
            bctype(is,ie) = BC_FLUID_FLUID;
            elemtype(ie) = ELEM_FLUID;
        else
            bctype(is,ie) = BC_SOLID_SOLID;
            elemtype(ie) = ELEM_SOLID;
        end

        % fluid solid interface
        if (mean(abs(y-y2))<0.1)
            bctype(is,ie) = BC_SOLID_FLUID;
        end

        if 0
        % reset to solid for testing
        bctype(is,ie) = BC_SOLID_SOLID;
        elemtype(ie) = ELEM_SOLID;

        % reset to fluid for testing
        bctype(is,ie) = BC_FLUID_FLUID;
        elemtype(ie) = ELEM_FLUID;
        %% set all to solid
        %bctype(is,ie) = BC_SOLID_SOLID;
        %elemtype(ie) = ELEM_SOLID;
        end

        %if (mean(abs(y-40e3))<1e-1)
        if (mean(abs(y-0.0))<1e-1)
            %bctype(is,ie) = BC_FREE_G;
            bctype(is,ie) = BC_FREE;
        end


        if (elemtype(ie) == ELEM_SOLID)
        %if ((mean(y4-0))<=0)
            %vp(ie) = 6000;
            %vs(ie) = 3464;
            %rho(ie) = 2670;
            vp(ie) = 1500*1.5 ;
            vs(ie) = 900*1.5 ;
            rho(ie) = 1000*1.5 ;
                        
            vp(ie) = 8000 ;
            vs(ie) = 4600 ;
            rho(ie) = 3200 ;
            
        %elseif (elemtype(ie) == ELEM_FLUID)
        else
            %vp(ie) = 6000;
            %vs(ie) = 3464;
            %rho(ie) = 2670;
            vp(ie) = 1500;
            vs(ie) = 900*0;
            rho(ie) = 1000;
        end
    end
end
if 0
elemtype(:) = ELEM_FLUID;
bctype(:) = BC_FLUID_FLUID;
end
if 0
elemtype(:) = ELEM_SOLID;
bctype(:) = BC_SOLID_SOLID;
end

%bctype(:) = -3;
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
    recv = [];
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
    recv = [];
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
    recv_is = find(bctype(:,recv_ie(i))==-3,1);
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


system('mkdir -p data');
db = gen_mesh_mpi_v5(mesh,nproc);

disp('gen_mesh_mpi done!')

%% check local node and elem ...
if flag_plot_mesh
figure


%subplot(121)
%for iproc = 1:nproc
%    idx = find(part==iproc);
elem1 = elem(1:4,:);
node1 = node*1e-3;
for ie = 1:size(elem1,2)
    if elemtype(ie) == ELEM_FLUID
        c = [127,205,255]/255;
    elseif elemtype(ie) == ELEM_SOLID
        c = [1 1 1]*0.9;
    end
    patch(node1(1,elem1(:,ie)),node1(2,elem1(:,ie)),...
        [1 1 1 1]*vp(ie),'FaceColor',c);
        %[1 1 1 1]*elemtype(ie));
    for is = 1:4
        if (bctype(is,ie) == BC_IN)
            c = 'r'; lw = 1;
        elseif (bctype(is,ie) == BC_FLUID_FLUID)
            c = 'c'; lw = 1;
        elseif (bctype(is,ie) == BC_FAULT)
            c = 'y'; lw = 2;
        elseif (bctype(is,ie) == BC_SOLID_FLUID)
            c = [0,100,0]/255; lw = 3;
        elseif (bctype(is,ie) == BC_FREE)
            c = 'b'; lw = 1;
        elseif (bctype(is,ie) == BC_FREE_G)
            c = 'k'; lw = 2;
        elseif (bctype(is,ie) == -999)
            c = 'w'; lw = 10;
        end
        %if (bctype(is,ie) == -3)
            is1 = is + 1;
            if (is1 > 4)
                is1 = is1-4;
            end
            xx = node1(1,elem([is,is1],ie));
            yy = node1(2,elem([is,is1],ie));
            hold on;
            plot(xx,yy,'Color',c,'linewidth',lw);
        %end
        %if (bctype(is,ie) == -1)
        %    is1 = is + 1;
        %    if (is1 > 4)
        %        is1 = is1-4;
        %    end
        %    xx = node(1,elem([is,is1],ie));
        %    yy = node(2,elem([is,is1],ie));
        %    hold on;
        %    plot(xx,yy,'g','linewidth',2);
        %
        %end
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
