function [x,y] = nodes2d(N,node,elem,xnode)
% node(nnode,2)
% elem(nelem,4)
% x or y (order+1,order+1,nelem)

nelem = size(elem,1);

%[r,w]=get_gll(N);
%[r,w]=lglnodes(N);
%[r,w]=glnodes(N);
%r = linspace(0,1,N+1);
r = xnode;
[r,s]=meshgrid(r,r);
r = r'; s = s'; % (nx by ny)
%r = r(:); s = s(:);
Np=(N+1)*(N+1);
x = zeros(N+1,N+1,nelem);
y = zeros(N+1,N+1,nelem);
%r = (r+1)/2; % [-1,1] to [0,1]
%s = (s+1)/2;
for i = 1:nelem
VX = node(elem(i,:),1);
VY = node(elem(i,:),2);
x(:,:,i) = (1-r).*(1-s)*VX(1)+r.*(1-s)*VX(2)+r.*s*(VX(3))+(1-r).*s*VX(4);
y(:,:,i) = (1-r).*(1-s)*VY(1)+r.*(1-s)*VY(2)+r.*s*(VY(3))+(1-r).*s*VY(4);
end
end
