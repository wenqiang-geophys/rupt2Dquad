function [bctype] = set_bctype_from_curve(elem,tri,id)

% elem(4,Nelem) quad
% tri(2,Nface) line

BC_FAULT = 100;

if (nargin < 3 )
    id = BC_FAULT;
end

tri = sort(tri);

fnode1 = tri(1,:);
fnode2 = tri(2,:);
%fnode3 = tri(3,:);

Nelem = size(elem,2); % elem(4,nelem);
Nnode = max(elem(:));

bctype = zeros(4,Nelem);

if(isempty(bctype))
    return
end

FToV = [ ...
    1,2,3;
    1,2,4;
    2,3,4;
    1,3,4];
FToV = [ ...
    1,2;
    2,3;
    3,4;
    4,1; ];

node_flag1 = zeros(Nnode,1);
node_flag2 = zeros(Nnode,1);
%node_flag3 = zeros(Nnode,1);
node_flag1(fnode1)=1;
node_flag2(fnode2)=1;
%node_flag3(fnode3)=1;

%elem4 = reshape(elem(FToV',:),[3,4*Nelem]);
elem4 = reshape(elem(FToV',:),[2,4*Nelem]);
elem4 = sort(elem4);

fn1 = node_flag1(elem4(1,:));
fn2 = node_flag2(elem4(2,:));
%fn3 = node_flag3(elem4(3,:));
%fn = fn1 & fn2 & fn3;
fn = fn1 & fn2;
bctype(:,:) = reshape(fn,[4,Nelem]);

bctype1 = zeros(size(bctype));
bctype1(bctype==1) = id;

bctype = bctype1;

end
