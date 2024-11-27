function check_mesh(m)
elem = m.elem;
node = m.node;
bctype = m.bctype;
c = m.vs; % color

rx = m.body_recv_x;
ry = m.body_recv_y;

myconstants

figure
x4 = reshape(node(1,elem),size(elem));
y4 = reshape(node(2,elem),size(elem));
% plot Vp
c4 = repmat(c,4,1);
patch(x4,y4,c4)

hold on;
for ie = 1:size(elem,2)
    for is = 1:4
        is1 = mod(is,4)+1;
        x2 = node(1,elem([is,is1],ie));
        y2 = node(2,elem([is,is1],ie));
        if (bctype(is,ie) >= BC_FAULT) % fault surface
            plot(x2,y2,'r','linewidth',1.5);
        end
        if (bctype(is,ie) == BC_FREE) % free surface
            plot(x2,y2,'b','linewidth',1.5);
        end
    end
end
plot(rx, ry, 'bx', 'linewidth', 1)

hold off
axis image;
%colormap sky
cb = colorbar;
ylabel(cb, 'Vs (m/s)')
xlabel('X (m)'); ylabel('Y (m)')
end
