function db = gen_mesh_mpi_v6(mesh,nproc)
% add elemtype
% add media and impedance mpi
% add fluxtype
% add media mpi for rho,vp,vs
% add netcdf output (v6)

myconstants

%%[EtoE,EtoF] = tiConnect2D_quad(elem);
mesh.node = mesh.node;
elem = mesh.elem;
node = mesh.node;
EtoE = mesh.EtoE;
EtoF = mesh.EtoF;
Nelem = mesh.Nelem;
Nnode = mesh.Nnode;

%nproc = 8;
if nproc > 1
    part = metis_part(elem,nproc) + 1;
else
    part = zeros(Nelem,1)+1;
end

% clean neighbor array
for ie = 1:size(EtoE,2)
    for is = 1:4
        if (EtoE(is,ie) == ie)
            EtoE(is,ie) = 0;
        end
    end
end
neighbor = EtoE;
face = EtoF;

%tic

max_neighbor = 10;
nsize = 100;

elmnts = elem;
elmnts = elmnts(:);

% local numbering of elements
glob2loc_elmnts = zeros(Nelem,1);
num_loc = ones(Nelem,1);
for i = 1:Nelem
    ipart = part(i);
    glob2loc_elmnts(i) = num_loc(ipart);
    num_loc(ipart) = num_loc(ipart) + 1;
end

% local node numbering
Nnodes_elmnts = zeros(Nnode,1);
nodes_elmnts = zeros(Nnode*nsize,1);

for i= 1 : 4*Nelem
    % nodes is like a matrix with notes as rows and nsize elements as colums
    nodes_elmnts((elmnts(i)-1)*nsize + Nnodes_elmnts(elmnts(i))+1) = 1+idivide(int32(i-1),int32(4));
    Nnodes_elmnts(elmnts(i)) = Nnodes_elmnts(elmnts(i)) + 1;
end
% now create the local node numbering
glob2loc_nodes_nparts = zeros(Nnode+1,1);
part_nodes = zeros(nproc,1);
num_parts = zeros(nproc,1);

size_glob2loc_nodes = 1;
part_nodes(:) = 0;

% go through all coordinates
for in = 1:Nnode
    glob2loc_nodes_nparts(in) = size_glob2loc_nodes;
    for ie = 1:Nnodes_elmnts(in)
        % nodes_elmnts((ie)+nsize*(in-1)) gibt eins der maximal _nsize_ elemente zum knoten _in_
        % gibt an in welchen partitionen der knoten alles liegt
        part_nodes(part(nodes_elmnts((ie)+nsize*(in-1))))=1;
    end
    for num_part = 1:nproc % gehe durch die partitionen
        if (part_nodes(num_part) == 1)
            % get number of nodes in the global array, there might be some nodes at the interfaces doubble
            size_glob2loc_nodes = size_glob2loc_nodes + 1;
            part_nodes(num_part) = 0;
        end
    end
end
glob2loc_nodes_nparts(Nnode+1) = size_glob2loc_nodes;
glob2loc_nodes_parts = zeros(size_glob2loc_nodes,1);
glob2loc_nodes       = zeros(size_glob2loc_nodes,1);

glob2loc_nodes(:) = 1;
part_nodes(:) = 0;
num_parts(:) = 1;
size_glob2loc_nodes = 1;


for in = 1:Nnode
    for ie = 1:Nnodes_elmnts(in)
        % nodes_elmnts((ie)+nsize*(in-1)) gibt eins der maximal _nsize_ elemente zum knoten _in_
        % gibt an in welchen partitionen der knoten alles liegt
        part_nodes(part(nodes_elmnts((ie)+nsize*(in-1))))=1;
    end
    for num_part = 1:nproc
        if (part_nodes(num_part) == 1)
            % build arrays with local nodes ordering
            glob2loc_nodes_parts(size_glob2loc_nodes) = num_part;
            glob2loc_nodes(size_glob2loc_nodes) = num_parts(num_part);
            size_glob2loc_nodes = size_glob2loc_nodes + 1;
            num_parts(num_part) = num_parts(num_part) + 1;
            part_nodes(num_part) = 0;
        end
    end
end

%db = {[],[],[],[]}


for iproc = nproc:-1:1

    db(iproc).Nelem = num_loc(iproc)-1;
    db(iproc).Nnode = num_parts(iproc)-1;
    db(iproc).elem = zeros(4,db(iproc).Nelem);
    db(iproc).node = zeros(2,db(iproc).Nnode);
    % build local coordinates
    for i = 1:Nnode
        for j = glob2loc_nodes_nparts(i) : glob2loc_nodes_nparts(i+1)-1
            if (glob2loc_nodes_parts(j) == iproc)
                k=glob2loc_nodes(j);
                db(iproc).node(1,k)=node(1,i);
                db(iproc).node(2,k)=node(2,i);
                db(iproc).loc2glob_nodes(k)=i;
            end
        end
    end

    loc_nodes = zeros(4,1);
    % build local elements / bulid local node indizes
    for i = 1:Nelem
        if (part(i) == iproc)
            for j = 1:4
                l=elmnts((i-1)*4+j);
                m=elmnts((i-1)*4+j);
                for k=glob2loc_nodes_nparts(l) : glob2loc_nodes_nparts(m+1)-1
                    if (glob2loc_nodes_parts(k) == iproc)
                        loc_nodes(j) = glob2loc_nodes(k);
                    end
                end
            end
            k = glob2loc_elmnts(i);
            db(iproc).elem(1,k) = loc_nodes(1);
            db(iproc).elem(2,k) = loc_nodes(2);
            db(iproc).elem(3,k) = loc_nodes(3);
            db(iproc).elem(4,k) = loc_nodes(4);
        end
    end

    % find mpi neighbors and the neigbor element
    db(iproc).mpi_interface = zeros(4,4,db(iproc).Nelem);
    db(iproc).neighbor = zeros(4,db(iproc).Nelem);
    db(iproc).face     = zeros(4,db(iproc).Nelem);
    db(iproc).direction= zeros(4,db(iproc).Nelem);
    db(iproc).bctype   = zeros(4,db(iproc).Nelem);
    db(iproc).fluxtype = zeros(4,db(iproc).Nelem);
    db(iproc).elemtype = zeros(1,db(iproc).Nelem);
    db(iproc).rho      = zeros(1,db(iproc).Nelem);
    db(iproc).vp       = zeros(1,db(iproc).Nelem);
    db(iproc).vs       = zeros(1,db(iproc).Nelem);
    db(iproc).zp       = zeros(1,db(iproc).Nelem);
    db(iproc).zs       = zeros(1,db(iproc).Nelem);
    %db(iproc).x     = zeros(order+1,order+1,db(iproc).Nelem);
    %db(iproc).y     = zeros(order+1,order+1,db(iproc).Nelem);
    %db(iproc).rx    = zeros(order+1,order+1,db(iproc).Nelem);
    %db(iproc).ry    = zeros(order+1,order+1,db(iproc).Nelem);
    %db(iproc).sx    = zeros(order+1,order+1,db(iproc).Nelem);
    %db(iproc).sy    = zeros(order+1,order+1,db(iproc).Nelem);
    %db(iproc).detJ  = zeros(order+1,order+1,db(iproc).Nelem);
    %db(iproc).nx    = zeros(4,order+1,db(iproc).Nelem);
    %db(iproc).ny    = zeros(4,order+1,db(iproc).Nelem);
    %db(iproc).norm_n= zeros(4,order+1,db(iproc).Nelem);

    db(iproc).neighbor(:,:) = 0;
    db(iproc).mpi_interface(:,:,:) = 0;
    db(iproc).pinterfaces = 0;

    for ie=1:Nelem
        for j=1:4
            k=neighbor(j,ie);
            if (k>0)
                if (part(ie) == iproc)
                    if (part(k) == iproc)
                        l = glob2loc_elmnts(ie);
                        m = glob2loc_elmnts(k);
                        db(iproc).neighbor(j,l) = m;
                    else
                        l = glob2loc_elmnts(ie);
                        m = glob2loc_elmnts(k);
                        db(iproc).neighbor(j,l) = -1;% ! means mpi neighbor
                        db(iproc).mpi_interface(1,j,l) = part(k);% ! neighbor partition
                        db(iproc).mpi_interface(2,j,l) = m;% ! element in partition
                        db(iproc).pinterfaces=db(iproc).pinterfaces+1;
                    end
                end
            end
        end
    end

    % set arrays to local element numbering
    for ie=1:Nelem
        if (part(ie)== iproc)
            l = glob2loc_elmnts(ie);
            db(iproc).face(:,l)=face(:,ie);
            db(iproc).direction(:,l)=mesh.direction(:,ie);
            db(iproc).bctype(:,l)=mesh.bctype(:,ie);
            db(iproc).fluxtype(:,l)=mesh.fluxtype(:,ie);
            db(iproc).elemtype(l)=mesh.elemtype(ie);
            db(iproc).rho(l) = mesh.rho(ie);
            db(iproc).vp(l) = mesh.vp(ie);
            db(iproc).vs(l) = mesh.vs(ie);
            db(iproc).zp(l) = mesh.zp(ie);
            db(iproc).zs(l) = mesh.zs(ie);
            % db(iproc).x(:,:,l)=mesh.x(:,:,ie);
            % db(iproc).y(:,:,l)=mesh.y(:,:,ie);
            % db(iproc).rx(:,:,l)=mesh.rx(:,:,ie);
            % db(iproc).ry(:,:,l)=mesh.ry(:,:,ie);
            % db(iproc).sx(:,:,l)=mesh.sx(:,:,ie);
            % db(iproc).sy(:,:,l)=mesh.sy(:,:,ie);
            % db(iproc).detJ(:,:,l)=mesh.detJ(:,:,ie);
            % db(iproc).nx(:,:,l)=mesh.nx(:,:,ie);
            % db(iproc).ny(:,:,l)=mesh.ny(:,:,ie);
            % db(iproc).norm_n(:,:,l)=mesh.norm_n(:,:,ie);
        end
    end
end


icom = zeros(nproc,nproc);
icom(:,:) = 0;
for iproc=1:nproc
    for ie=1:db(iproc).Nelem
        for i=1:4
            if (db(iproc).mpi_interface(1,i,ie) > 0)
                icom(iproc,db(iproc).mpi_interface(1,i,ie))=1;
            end
        end
    end
    db(iproc).mpi_nn=sum(icom(iproc,:));
end
% set up mpi_neighbor array
for iproc=1:nproc
    db(iproc).mpi_neighbor = zeros(db(iproc).mpi_nn,1);
    c=1;
    for i=1:nproc
        if (icom(iproc,i) == 1)
            db(iproc).mpi_neighbor(c)=i;
            c=c+1;
        end
    end
end

% how many mpi interfaces do we have?
for iproc=1:nproc
    tempv = zeros(db(iproc).mpi_nn,1);
    db(iproc).mpi_ninterface = zeros(db(iproc).mpi_nn,1);
    tempv(:)=0;
    c=1;
    for i=1:db(iproc).mpi_nn
        in=db(iproc).mpi_neighbor(i);
        for ie=1:db(iproc).Nelem
            for is=1:4
                if ( (db(iproc).neighbor(is,ie) == -1) && (db(iproc).mpi_interface(1,is,ie) == in))
                    tempv(i)=tempv(i)+1;
                end
            end
        end
        c=c+1;
        db(iproc).mpi_ninterface(i)=tempv(i);
    end
end

temp2=0;
for iproc=1:nproc
    temp1=db(iproc).Nelem;
    if (temp2<temp1); temp2=temp1; end;
end

% set up mpi_ibool with max mpi neighbors
for iproc=1:nproc
    db(iproc).mpi_nemax=temp2;
    %allocate(db(iproc).mpi_ibool(4,db(iproc).mpi_nemax));
    db(iproc).mpi_ibool = zeros(4,db(iproc).mpi_nemax);
    db(iproc).mpi_ibool(:,:)=0;
end

temp2=0;
for iproc=1:nproc
    temp1=max(db(iproc).mpi_ninterface(:));
    if (temp2<temp1); temp2=temp1; end;
end

% build interface array
%"--------------------------------------------------------------------------------------"
% create main mpi reconnection arrays.
%"--------------------------------------------------------------------------------------"
mpi_nn_vec = zeros(nproc,1);
for iproc=1:nproc
    mpi_nn_vec(iproc) = db(iproc).mpi_nn;
end
for iproc=1:nproc
    db(iproc).mpi_nnmax=max(mpi_nn_vec);
    %allocate(db(iproc).mpi_icon(db(iproc).mpi_nnmax));
    db(iproc).mpi_icon = zeros(db(iproc).mpi_nnmax,1);
    db(iproc).mpi_icon(:) = 0;
    db(iproc).mpi_ne=temp2;
    % arbeite mit max values
    %allocate(db(iproc)%mpi_connection(db(iproc)%mpi_nn,db(iproc)%mpi_ne,2))
    db(iproc).mpi_connection = zeros(db(iproc).mpi_nn,db(iproc).mpi_ne,2);
    db(iproc).mpi_connection(:,:,:)=0;
    for i=1:db(iproc).mpi_nn
        in=db(iproc).mpi_neighbor(i);
        c=1;
        for ie=1:db(iproc).Nelem
            for is=1:4
                if ( (db(iproc).neighbor(is,ie) == -1) && (db(iproc).mpi_interface(1,is,ie) == in))
                    db(iproc).mpi_connection(i,c,1)=ie; % mpi interface to global element
                    db(iproc).mpi_connection(i,c,2)=is; % which side of the element is the interface
                    db(iproc).mpi_interface(3,is,ie)=c; % global element to local mpi interface
                    db(iproc).mpi_interface(4,is,ie)=i; % which local interface?
                    c=c+1;
                end
            end
        end
    end
end

% create impedance mpi connection array
for iproc = 1:nproc
    %db(iproc).mpi_zp = zeros(db(iproc).pinterfaces,3);
    %db(iproc).mpi_zs = zeros(db(iproc).pinterfaces,3);
    db(iproc).mpi_rho = zeros(db(iproc).pinterfaces,3);
    db(iproc).mpi_vp = zeros(db(iproc).pinterfaces,3);
    db(iproc).mpi_vs = zeros(db(iproc).pinterfaces,3);
    i = 1;
    for ie = 1:db(iproc).Nelem
        for is = 1:4
            if (db(iproc).neighbor(is,ie) == -1)
                idb = db(iproc).mpi_interface(1,is,ie);
                iz  = db(iproc).mpi_interface(2,is,ie);
                %db(iproc).mpi_zp(i,1) = db(idb).zp(iz);
                %db(iproc).mpi_zp(i,2) = idb;
                %db(iproc).mpi_zp(i,3) = iz;
                %db(iproc).mpi_zs(i,1) = db(idb).zs(iz);
                %db(iproc).mpi_zs(i,2) = idb;
                %db(iproc).mpi_zs(i,3) = iz;
                db(iproc).mpi_rho(i,1) = db(idb).rho(iz);
                db(iproc).mpi_rho(i,2) = idb;
                db(iproc).mpi_rho(i,3) = iz;
                db(iproc).mpi_vp(i,1) = db(idb).vp(iz);
                db(iproc).mpi_vp(i,2) = idb;
                db(iproc).mpi_vp(i,3) = iz;
                db(iproc).mpi_vs(i,1) = db(idb).vs(iz);
                db(iproc).mpi_vs(i,2) = idb;
                db(iproc).mpi_vs(i,3) = iz;
                i = i + 1;
            end
        end
    end
end

for iproc = 1:nproc
    for i = 1:db(iproc).mpi_nn
        l = db(iproc).mpi_neighbor(i);
        for ie = 1:db(iproc).mpi_ne
            if ( db(iproc).mpi_connection(i,ie,1) > 0)
                is = db(iproc).mpi_connection(i,ie,2);
                ee = db(iproc).mpi_connection(i,ie,1);
                k = db(iproc).mpi_interface(2,is,ee);
                c = db(l).mpi_interface(3,db(iproc).face(is,ee),k);
                db(iproc).mpi_icon(i)=db(l).mpi_interface(4,db(iproc).face(is,ee),k);
                %very important for correct mpi boundarys
                db(iproc).mpi_ibool(is,ee) = c;
                %for j=1,NGLL
                %    mpi_ti = db(iproc)%ibt( j,is,ee)
                %    mpi_ni = db(l)%ibt(j,db(iproc)%face(is,ee),k)
                %    ! LL LL 24.4.14
                %    if ((abs(db(iproc)%vx(mpi_ti) - db(l)%vx(mpi_ni) ) < epsi)&
                %    .and. (abs(db(iproc)%vz(mpi_ti) - db(l)%vz(mpi_ni))  < epsi)) then
                %        db(l)%mpi_ibt(j,db(iproc)%face(is,ee),k)=j
                %    else
                %        db(l)%mpi_ibt(j,db(iproc)%face(is,ee),k)=NGLL+1-j
                %    end if
                %end do
                %do j=1,NGLL
                %    mpi_ti = db(iproc)%ibt( j,is,ee)
                %    mpi_ni = db(l)%ibt(db(l)%mpi_ibt(j,db(iproc)%face(is,ee),k),db(iproc)%face(is,ee),k)
                %end do
            end
        end
    end
end

% fault recvs
nrecv = mesh.nrecv;
db(iproc).recv_fid = [];
db(iproc).recv_i = [];
db(iproc).recv_refx = [];
db(iproc).recv_ie = [];

for iproc = 1:nproc
    db(iproc).recv_number = 0;
    k = 0;
    for i = 1:nrecv
        ie = mesh.recv_ie(i);
        if (part(ie)== iproc)
            k = k + 1;
            db(iproc).recv_fid(k)=i;
            db(iproc).recv_i(k)=mesh.recv_i(i);
            db(iproc).recv_refx(k)=mesh.recv_refx(i);
            %db(iproc).recv_j(k)=mesh.recv_j(i);
            db(iproc).recv_ie(k)= glob2loc_elmnts(ie);
        end
        db(iproc).recv_number = k;
    end
end
% body recvs
nrecv = mesh.body_nrecv;
for iproc = 1:nproc
    db(iproc).body_recv_number = 0;
    k = 0;
    for i = 1:nrecv
        ie = mesh.body_recv_ie(i);
        if (part(ie)== iproc)
            k = k + 1;
            db(iproc).body_recv_fid(k)=i;
            db(iproc).body_recv_i(k)=mesh.body_recv_i(i);
            db(iproc).body_recv_j(k)=mesh.body_recv_j(i);
            db(iproc).body_recv_refx(k)=mesh.body_recv_refx(i);
            db(iproc).body_recv_refy(k)=mesh.body_recv_refy(i);
            db(iproc).body_recv_ie(k)= glob2loc_elmnts(ie);
        end
        db(iproc).body_recv_number = k;
    end
end

% counting fault elems
for iproc = 1:nproc
    db(iproc).Nelem_fault = 0;
    for ie = 1:db(iproc).Nelem
        isfault = 0;
        for is = 1:Nfaces
            if(db(iproc).bctype(is,ie) >= BC_FAULT)
                isfault=1;
            end
        end
        if (isfault==1)
            db(iproc).Nelem_fault = db(iproc).Nelem_fault+1;
        end
    end

    % fault2wave, wave2fault

    disp(['iproc=',num2str(iproc),';Nelem_fault=',num2str(db(iproc).Nelem_fault)]);
    db(iproc).fault2wave = zeros(db(iproc).Nelem_fault,1);
    db(iproc).wave2fault = zeros(db(iproc).Nelem,1);
    k = 0;
    for ie = 1:db(iproc).Nelem
        isfault=0;
        for is = 1:Nfaces
            if (db(iproc).bctype(is,ie) >= BC_FAULT)
                isfault=1;
            end
        end
        if (isfault)
            k = k + 1;
            db(iproc).fault2wave(k) = ie;
            db(iproc).wave2fault(ie) = k;
        end
    end
end

if 0
    for iproc = 1:nproc
        fnm = sprintf('data/meshVar%06d',iproc-1);
        fid = fopen(fnm, 'w');

        fprintf(fid, '%d\n',db(iproc).Nnode);
        fprintf(fid, '%g ',db(iproc).node); fprintf(fid,'\n');
        fprintf(fid, '%d\n',db(iproc).Nelem);
        fprintf(fid, '%d ',db(iproc).elem); fprintf(fid, '\n');
        fprintf(fid, '%d ',db(iproc).neighbor); fprintf(fid, '\n');
        fprintf(fid, '%d ',db(iproc).face); fprintf(fid, '\n');
        fprintf(fid, '%d ',db(iproc).direction); fprintf(fid, '\n');
        fprintf(fid, '%d ',db(iproc).bctype); fprintf(fid, '\n');
        fprintf(fid, '%d ',db(iproc).fluxtype); fprintf(fid, '\n');
        fprintf(fid, '%d ',db(iproc).elemtype); fprintf(fid, '\n');
        fprintf(fid, '%d ',db(iproc).rho); fprintf(fid, '\n');
        fprintf(fid, '%d ',db(iproc).vp); fprintf(fid, '\n');
        fprintf(fid, '%d ',db(iproc).vs); fprintf(fid, '\n');

        fprintf(fid, '%d\n',db(iproc).mpi_nn);
        fprintf(fid, '%d\n',db(iproc).mpi_ne);
        fprintf(fid, '%d\n',db(iproc).mpi_nemax);
        fprintf(fid, '%d\n',db(iproc).pinterfaces);
        fprintf(fid, '%d ',db(iproc).mpi_neighbor); fprintf(fid, '\n');
        fprintf(fid, '%d ',db(iproc).mpi_connection); fprintf(fid, '\n');
        %fprintf(fid, '%d ',db(iproc).loc2glob_nodes); fprintf(fid, '\n');
        fprintf(fid, '%d ',db(iproc).mpi_ibool); fprintf(fid, '\n');
        fprintf(fid, '%d ',db(iproc).mpi_interface); fprintf(fid, '\n');
        if (nproc > 1)
            fprintf(fid, '%d ',db(iproc).mpi_rho); fprintf(fid, '\n');
            fprintf(fid, '%d ',db(iproc).mpi_vp); fprintf(fid, '\n');
            fprintf(fid, '%d ',db(iproc).mpi_vs); fprintf(fid, '\n');
        end

        fprintf(fid, '%d ',db(iproc).recv_number); fprintf(fid, '\n');
        if (db(iproc).recv_number > 0)
            fprintf(fid, '%d ',db(iproc).recv_fid); fprintf(fid, '\n');
            fprintf(fid, '%d ',db(iproc).recv_i); fprintf(fid, '\n');
            %fprintf(fid, '%d ',db(iproc).recv_j); fprintf(fid, '\n');
            fprintf(fid, '%d ',db(iproc).recv_ie); fprintf(fid, '\n');
            fprintf(fid, '%d ',db(iproc).recv_refx); fprintf(fid, '\n');
        end

        fprintf(fid, '%d ',db(iproc).body_recv_number); fprintf(fid, '\n');
        if (db(iproc).body_recv_number > 0)
            fprintf(fid, '%d ',db(iproc).body_recv_fid); fprintf(fid, '\n');
            fprintf(fid, '%d ',db(iproc).body_recv_i); fprintf(fid, '\n');
            fprintf(fid, '%d ',db(iproc).body_recv_j); fprintf(fid, '\n');
            fprintf(fid, '%d ',db(iproc).body_recv_ie); fprintf(fid, '\n');
            fprintf(fid, '%g ',db(iproc).body_recv_refx); fprintf(fid, '\n');
            fprintf(fid, '%g ',db(iproc).body_recv_refy); fprintf(fid, '\n');
        end
        fclose(fid);
    end
end

for iproc = 1:nproc
    fnm = sprintf('data/meshVar%06d.nc',iproc-1);
    ncid = netcdf.create(fnm,'CLOBBER');
    dimid_node = netcdf.defDim(ncid,'Nnode',db(iproc).Nnode);
    dimid_elem = netcdf.defDim(ncid,'Nelem',db(iproc).Nelem);
    dimid_fault_elem = netcdf.defDim(ncid,'Nelem_fault',db(iproc).Nelem_fault);

    dimid2 = netcdf.defDim(ncid,'two',2);
    dimid3 = netcdf.defDim(ncid,'three',3);
    dimid4 = netcdf.defDim(ncid,'four',4);
    dimid_face = netcdf.defDim(ncid,'Nfaces',4);
    dimid_vertex = netcdf.defDim(ncid,'Nvertex',4);
    dimid_mpi_nn = netcdf.defDim(ncid,'mpi_nn',db(iproc).mpi_nn);
    dimid_mpi_ne = netcdf.defDim(ncid,'mpi_ne',db(iproc).mpi_ne);
    dimid_mpi_mx = netcdf.defDim(ncid,'mpi_nemax',db(iproc).mpi_nemax);
    dimid_mpi_pf = netcdf.defDim(ncid,'pinterfaces',db(iproc).pinterfaces);

    varid1 = netcdf.defVar(ncid,'node','NC_DOUBLE',[dimid2,dimid_node]);
    varid2 = netcdf.defVar(ncid,'elem','NC_INT',[dimid_vertex,dimid_elem]);
    varid3 = netcdf.defVar(ncid,'neighbor','NC_INT',[dimid_face,dimid_elem]);
    varid4 = netcdf.defVar(ncid,'face','NC_INT',[dimid_face,dimid_elem]);
    varid5 = netcdf.defVar(ncid,'direction','NC_INT',[dimid_face,dimid_elem]);
    varid6 = netcdf.defVar(ncid,'bctype','NC_INT',[dimid_face,dimid_elem]);
    varid7 = netcdf.defVar(ncid,'fluxtype','NC_INT',[dimid_face,dimid_elem]);
    varid8 = netcdf.defVar(ncid,'elemtype','NC_INT',[dimid_elem]);
    varid9 = netcdf.defVar(ncid,'rho','NC_DOUBLE',[dimid_elem]);
    varid10 = netcdf.defVar(ncid,'vp','NC_DOUBLE',[dimid_elem]);
    varid11 = netcdf.defVar(ncid,'vs','NC_DOUBLE',[dimid_elem]);
    varid32 = netcdf.defVar(ncid,'fault2wave','NC_INT',[dimid_fault_elem]);
    varid33 = netcdf.defVar(ncid,'wave2fault','NC_INT',[dimid_elem]);
    % fault congfiguration
    varid12 = netcdf.defVar(ncid,'Sxx0'  ,'NC_DOUBLE',[dimid2,dimid_face,dimid_fault_elem]);
    varid13 = netcdf.defVar(ncid,'Syy0'  ,'NC_DOUBLE',[dimid2,dimid_face,dimid_fault_elem]);
    varid14 = netcdf.defVar(ncid,'Sxy0'  ,'NC_DOUBLE',[dimid2,dimid_face,dimid_fault_elem]);
    varid30 = netcdf.defVar(ncid,'sigma0','NC_DOUBLE',[dimid2,dimid_face,dimid_fault_elem]);
    varid31 = netcdf.defVar(ncid,'tau0'  ,'NC_DOUBLE',[dimid2,dimid_face,dimid_fault_elem]);
    varid15 = netcdf.defVar(ncid,'mu_s'  ,'NC_DOUBLE',[dimid2,dimid_face,dimid_fault_elem]);
    varid16 = netcdf.defVar(ncid,'mu_d'  ,'NC_DOUBLE',[dimid2,dimid_face,dimid_fault_elem]);
    varid17 = netcdf.defVar(ncid,'Dc'    ,'NC_DOUBLE',[dimid2,dimid_face,dimid_fault_elem]);
    varid18 = netcdf.defVar(ncid,'C0'    ,'NC_DOUBLE',[dimid2,dimid_face,dimid_fault_elem]);
    varid19 = netcdf.defVar(ncid,'a'     ,'NC_DOUBLE',[dimid2,dimid_face,dimid_fault_elem]);
    varid20 = netcdf.defVar(ncid,'b'     ,'NC_DOUBLE',[dimid2,dimid_face,dimid_fault_elem]);
    varid21 = netcdf.defVar(ncid,'Vw'    ,'NC_DOUBLE',[dimid2,dimid_face,dimid_fault_elem]);
    varid22 = netcdf.defVar(ncid,'state' ,'NC_DOUBLE',[dimid2,dimid_face,dimid_fault_elem]);
    % mpi
    varid23 = netcdf.defVar(ncid,'mpi_rho','NC_DOUBLE',[dimid_mpi_pf,dimid3]);
    varid24 = netcdf.defVar(ncid,'mpi_vp' ,'NC_DOUBLE',[dimid_mpi_pf,dimid3]);
    varid25 = netcdf.defVar(ncid,'mpi_vs' ,'NC_DOUBLE',[dimid_mpi_pf,dimid3]);
    varid26 = netcdf.defVar(ncid,'mpi_neighbor','NC_DOUBLE',[dimid_mpi_nn]);
    varid27 = netcdf.defVar(ncid,'mpi_connection','NC_DOUBLE',[dimid_mpi_nn,dimid_mpi_ne,dimid2]);
    varid28 = netcdf.defVar(ncid,'mpi_ibool','NC_DOUBLE',[dimid4,dimid_mpi_mx]);
    varid29 = netcdf.defVar(ncid,'mpi_interface','NC_DOUBLE',[dimid4,dimid4,dimid_elem]);
    netcdf.endDef(ncid);

    netcdf.putVar(ncid,varid1,db(iproc).node);
    netcdf.putVar(ncid,varid2,db(iproc).elem);
    netcdf.putVar(ncid,varid3,db(iproc).neighbor);
    netcdf.putVar(ncid,varid4,db(iproc).face);
    netcdf.putVar(ncid,varid5,db(iproc).direction);
    netcdf.putVar(ncid,varid6,db(iproc).bctype);
    netcdf.putVar(ncid,varid7,db(iproc).fluxtype);
    netcdf.putVar(ncid,varid8,db(iproc).elemtype);
    netcdf.putVar(ncid,varid9,db(iproc).rho);
    netcdf.putVar(ncid,varid10,db(iproc).vp);
    netcdf.putVar(ncid,varid11,db(iproc).vs);
    if (db(iproc).Nelem_fault > 0)
        netcdf.putVar(ncid,varid32,db(iproc).fault2wave);
    end
    netcdf.putVar(ncid,varid33,db(iproc).wave2fault);
    % mpi
    netcdf.putVar(ncid,varid23,db(iproc).mpi_rho);
    netcdf.putVar(ncid,varid24,db(iproc).mpi_vp);
    netcdf.putVar(ncid,varid25,db(iproc).mpi_vs);
    netcdf.putVar(ncid,varid26,db(iproc).mpi_neighbor);
    netcdf.putVar(ncid,varid27,db(iproc).mpi_connection);
    netcdf.putVar(ncid,varid28,db(iproc).mpi_ibool);
    netcdf.putVar(ncid,varid29,db(iproc).mpi_interface);

    netcdf.close(ncid);


end
