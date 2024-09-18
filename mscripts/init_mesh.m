function mesh = init_mesh(mesh, order)

% Nodes and quadrature points
[xnode,wts] = lglnodes(order);
%[xint,wts2] = lglnodes(order+0);
%W     = diag(wts2);
%[G,D] = lagint(xnode, xint);
%[~,D1] = lagint(xnode,xnode);
%G = eye(order+1);
%D = derivative_GL(order,xnode,wts);
%D = D';
%D1 = D;

%mesh.G = G;
%mesh.D = D;
%mesh.W = W;

%%

%basename = 'mesh/circle_dh4';
%[node,elem] = read_mesh2(basename);
node = mesh.node;
elem = mesh.elem;
Nelem = size(elem,1);
Nnode = size(node,1);

%elem = ccw_sort(node,elem);

[EtoE,EtoF] = tiConnect2D_quad(elem);
% clean neighbor array
for ie = 1:size(EtoE,1)
    for is = 1:4
        if (EtoE(ie,is) == ie)
            EtoE(ie,is) = 0;
        end
    end
end

direction = get_face_direction(node,elem,EtoE,EtoF);
[x,y] = nodes2d(order,node,elem,xnode);

mesh.Nelem = Nelem;
mesh.Nnode = Nnode;

mesh.node = node';
mesh.elem = elem';
mesh.EtoE = EtoE';
mesh.EtoF = EtoF';
mesh.direction = direction';

mesh.x = x;
mesh.y = y;

%disp('calculating metrics...')
%xr = zeros(order+1,order+1,Nelem);
%xs = zeros(order+1,order+1,Nelem);
%yr = zeros(order+1,order+1,Nelem);
%ys = zeros(order+1,order+1,Nelem);
%detJ = zeros(order+1,order+1,Nelem);
%J = zeros(2,2,order+1,order+1,Nelem);
%invJ = zeros(2,2,order+1,order+1,Nelem);
%
%[r,s]=meshgrid(xnode,xnode);
%r = r'; s = s';
%for i = 1:Nelem
%%    xr(:,:,i) = D1 * squeeze(x(:,:,i));
%%    xs(:,:,i) = transpose(D1*transpose(squeeze(x(:,:,i))));
%%    yr(:,:,i) = D1 * squeeze(y(:,:,i));
%%    ys(:,:,i) = transpose(D1*transpose(squeeze(y(:,:,i))));
%    vx = node(elem(i,:),1);
%    vy = node(elem(i,:),2);
%    [c1,c2,c3,c4] = cal_covariant(vx,vy,r,s);
%    xr(:,:,i) = c1;
%    xs(:,:,i) = c2;
%    yr(:,:,i) = c3;
%    ys(:,:,i) = c4;  
%end
%
%
%% (xr yr)
%% (xs ys)
%%
%% (rx sx)
%% (ry sy)
%J(1,1,:,:,:) = xr;
%J(2,1,:,:,:) = xs;
%J(1,2,:,:,:) = yr;
%J(2,2,:,:,:) = ys;
%
%for ex = 1:Nelem
%    for i = 1:order+1
%        for j = 1:order+1
%            detJ(i,j,ex) = det(squeeze(J(:,:,i,j,ex)));
%            invJ(:,:,i,j,ex) = inv(squeeze(J(:,:,i,j,ex)));
%        end
%    end
%end
%
%rx = squeeze(invJ(1,1,:,:,:));
%ry = squeeze(invJ(2,1,:,:,:));
%sx = squeeze(invJ(1,2,:,:,:));
%sy = squeeze(invJ(2,2,:,:,:));
%
%nx = zeros(4,order+1,Nelem);
%ny = zeros(4,order+1,Nelem);
%norm_n = zeros(4,order+1,Nelem);
%
%%       3
%%    4------3
%%    |      | 2
%%  4 |      | 
%%    1------2
%%       1
%nx(4,:,:) = -detJ(1,:,:)    .*rx(1,:,:)  ;
%ny(4,:,:) = -detJ(1,:,:)    .*ry(1,:,:)  ;
%nx(2,:,:) =  detJ(end,:,:)  .*rx(end,:,:);
%ny(2,:,:) =  detJ(end,:,:)  .*ry(end,:,:);
%nx(1,:,:) = -detJ(:,1,:)    .*sx(:,1,:)  ;
%ny(1,:,:) = -detJ(:,1,:)    .*sy(:,1,:)  ;
%nx(3,:,:) =  detJ(:,end,:)  .*sx(:,end,:);
%ny(3,:,:) =  detJ(:,end,:)  .*sy(:,end,:);
%
%for ie = 1:Nelem
%    for i = 1:order+1
%        for is = 1:4
%            if (is == 1) 
%                rx1   = sx  (i,:,ie);
%                ry1   = sy  (i,:,ie);
%                detJ1 = detJ(i,:,ie);
%                nx1 = -detJ1.*rx1;
%                ny1 = -detJ1.*ry1;
%                nx2 = lagrange_interpol(order+1,0,xnode,nx1);
%                ny2 = lagrange_interpol(order+1,0,xnode,ny1);
%            elseif (is == 2)
%                rx1   = rx  (:,i,ie);
%                ry1   = ry  (:,i,ie);
%                detJ1 = detJ(:,i,ie);
%                nx1 = detJ1.*rx1;
%                ny1 = detJ1.*ry1;
%                nx2 = lagrange_interpol(order+1,1,xnode,nx1);
%                ny2 = lagrange_interpol(order+1,1,xnode,ny1);
%            elseif (is == 3)
%                rx1   = sx  (i,:,ie);
%                ry1   = sy  (i,:,ie);
%                detJ1 = detJ(i,:,ie);
%                nx1 = detJ1.*rx1;
%                ny1 = detJ1.*ry1;
%                nx2 = lagrange_interpol(order+1,1,xnode,nx1);
%                ny2 = lagrange_interpol(order+1,1,xnode,ny1);
%            elseif (is == 4)
%                rx1   = rx  (:,i,ie);
%                ry1   = ry  (:,i,ie);
%                detJ1 = detJ(:,i,ie);
%                nx1 = -detJ1.*rx1;
%                ny1 = -detJ1.*ry1;
%                nx2 = lagrange_interpol(order+1,0,xnode,nx1);
%                ny2 = lagrange_interpol(order+1,0,xnode,ny1);
%            end
%            nx(is,i,ie) = nx2;
%            ny(is,i,ie) = ny2;
%            n = [nx(is,i,ie),ny(is,i,ie)];
%            norm_n(is,i,ie) = norm(n);
%        end
%    end
%end
%nx = nx./norm_n;
%ny = ny./norm_n;
%
%mesh.rx = rx;
%mesh.ry = ry;
%mesh.sx = sx;
%mesh.sy = sy;
%mesh.detJ = detJ;
%%mesh.invJ = invJ;
%mesh.nx = nx;
%mesh.ny = ny;
%mesh.norm_n = norm_n;

end
