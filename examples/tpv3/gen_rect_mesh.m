function [node,elem,bctype] = gen_rect_mesh()

myconstants;

dx = 200;
dy1 = 200;
dy2 = 1*dy1;

nx = 100;
ny1 = 60;
ny2 = ny1/1;

Nx = 2*nx+1;
Ny = ny1+ny2+1;
X = zeros(Nx,Ny);
Y = zeros(Nx,Ny);
for i = -nx:nx
    for j = -ny1:ny2
        if (j > 0)
            dy = dy2;
        else
            dy = dy1;
        end
        X(i+nx+1,j+ny1+1) = dx*i;
        Y(i+nx+1,j+ny1+1) = dy*j;
    end
end

for i = -nx:nx
    for j = 1:ny2
        dy = dy2;
        x = dx*i+dy*j*0;
        y0 = dy*j;
        X(i+nx+1,j+ny1+1) = x;
        dymin = dy2;
        dymax = 2*dy2;
        amp = (1+cos(x*2*pi/5e3))/2; % 0~1
        damp = exp(-(y0/2e3).^2);
        dy = dymin+damp*amp*(dymax-dymin);
        Y(i+nx+1,j+ny1+1) = dy+Y(i+nx+1,j+ny1);
    end
end


if 0
    figure
    pcolor(X*1e-3,Y*1e-3,zeros(size(X)))
    axis image
    caxis([-1e30 0])
    colormap gray
    xlabel('X (km)')
    ylabel('Y (km)')
end

Nnode = Nx*Ny;
Nelem = (Nx-1)*(Ny-1);
node = zeros(2,Nnode);
elem = zeros(4,Nelem);

for i = 1:Nx
    for j = 1:Ny
        in = i+(j-1)*Nx;
        node(1,in) = X(i,j);
        node(2,in) = Y(i,j);
    end
end

for i = 1:Nx-1
    for j = 1:Ny-1
        ie = i+(j-1)*(Nx-1);
        in = i+(j-1)*Nx;
        elem(:,ie) = [in,in+1,in+Nx+1,in+Nx];
    end
end

bctype = zeros(4,Nelem);

for ie = 1:Nelem
    for is = 1:4
        is1 = mod(is,4)+1;
        x2 = node(1,elem([is,is1],ie));
        y2 = node(2,elem([is,is1],ie));
        xc = mean(x2);
        yc = mean(y2);
        if (abs(yc)<1e-16 && abs(xc)<10e3)
            bctype(is,ie) = BC_FAULT;
        end

    end
end
