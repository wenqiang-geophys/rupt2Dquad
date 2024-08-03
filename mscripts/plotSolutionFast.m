function plotSolutionFast(tri,x,y,u,varargin)
% u(o+1,o+1,nelem)
% x(o+1,o+1,nelem)
% y(o+1,o+1,nelem)
% nelem = size(x,3);
% tri is assembled element by element
%tri = [];
%for i = 1:nelem
%    tri1 = delaunay(x(:,:,i),y(:,:,i)) + (i-1)*NGLL*NGLL;
%    tri = [tri;tri1];
%end
%toc

%u = squeeze(u);
trisurf(tri,x(:),y(:),u(:))

shading interp

flag_line = 0;
if nargin == 5
    flag_line = 1;
end


if flag_line
hold on
xx = [x(1,1,:);x(end,1,:);x(end,end,:);x(1,end,:);x(1,1,:)];
yy = [y(1,1,:);y(end,1,:);y(end,end,:);y(1,end,:);y(1,1,:)];
xx = squeeze(xx);
yy = squeeze(yy);
plot3(xx,yy,1e30*ones(size(xx)),'color',[1,1,1]*.5,'LineWidth',0.25)
hold off % dont forget this
end

view(2)
daspect([1 1 1])
axis image
%xmin = min(x(:));
%xmax = max(x(:));
%ymin = min(y(:));
%ymax = max(y(:));
umax = max(u(:));
%axis([xmin xmax ymin ymax])
%caxis([-1 1]*umax/2);
colormap('jet')
colorbar
grid off

end
