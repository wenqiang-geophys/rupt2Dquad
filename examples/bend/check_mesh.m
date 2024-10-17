



% check local node and elem ...
if flag_plot_mesh
    figure

    elem1 = elem(1:4,:);
    node1 = node;
    for ie = 1:size(elem1,2)
        %patch(node1(1,elem1(:,ie)),node1(2,elem1(:,ie)),[1 1 1 1]*-1);
        plot(node1(1,elem1([1:4,1],ie)),node1(2,elem1([1:4,1],ie)),'k','LineWidth',0.5)
        hold on
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
    end
    %end
    axis image;
    hold off
    colormap default
end
