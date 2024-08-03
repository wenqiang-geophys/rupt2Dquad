function direction = get_face_direction(node,elem,EtoE,EtoF)
nelem = size(elem,1);
direction = zeros(nelem,4);
% check neigh ...
for ie = 1:nelem
    for is = 1:4
        
        
        %if (is == 1)
        %    x1 = x(:,1,ie);
        %    y1 = y(:,1,ie);
        %elseif (is == 2)
        %    x1 = x(end,:,ie);
        %    y1 = y(end,:,ie);
        %elseif (is == 3);
        %    x1 = x(:,end,ie);
        %    y1 = y(:,end,ie);
        %elseif (is == 4);
        %    x1 = x(1,:,ie);
        %    y1 = y(1,:,ie);
        %end
        
        xx = node(elem(ie,:),1);
        yy = node(elem(ie,:),2);
        
        if (is == 1)
            x1 = [xx(1),xx(2)];
            y1 = [yy(1),yy(2)];
        elseif (is == 2)
            x1 = [xx(2),xx(3)];
            y1 = [yy(2),yy(3)];
        elseif (is == 3)
            x1 = [xx(4),xx(3)];
            y1 = [yy(4),yy(3)];
        elseif (is == 4)
            x1 = [xx(1),xx(4)];
            y1 = [yy(1),yy(4)];
        end
        
        j = EtoE(ie,is); % neigh
        face = EtoF(ie,is);
        
        if (j > 0)
        xx = node(elem(j,:),1);
        yy = node(elem(j,:),2);
        else
            xx = [1 1 1 1]*1e100;
            yy = xx;
        end
        
        if (face == 1)
            x2 = [xx(1),xx(2)];
            y2 = [yy(1),yy(2)];
        elseif (face == 2)
            x2 = [xx(2),xx(3)];
            y2 = [yy(2),yy(3)];
        elseif (face == 3)
            x2 = [xx(4),xx(3)];
            y2 = [yy(4),yy(3)];
        elseif (face == 4)
            x2 = [xx(1),xx(4)];
            y2 = [yy(1),yy(4)];
        end
        
        %j = EtoE(ie,is); % neigh
        %face = EtoF(ie,is);
        %if (face == 1)
        %    x2 = x(:,1,j);
        %    y2 = y(:,1,j);
        %elseif (face == 2)
        %    x2 = x(end,:,j);
        %    y2 = y(end,:,j);
        %elseif (face == 3)
        %    x2 = x(:,end,j);
        %    y2 = y(:,end,j);
        %elseif (face == 4)
        %    x2 = x(1,:,j);
        %    y2 = y(1,:,j);
        %end
        
        x1 = reshape(x1,[],1);
        x2 = reshape(x2,[],1);
        y1 = reshape(y1,[],1);
        y2 = reshape(y2,[],1);
        
        x2r = x2(end:-1:1);
        y2r = y2(end:-1:1);
        err = sqrt((x1-x2).^2+(y1-y2).^2);
        rerr = sqrt((x1-x2r).^2+(y1-y2r).^2);
        
        if norm(err) > 1e-6
            direction(ie,is) = 1;
            %fprintf('ie = %d neigh = %d is = %d face = %d error = %g r = %g\n',...
            %    ie,j,is,face,norm(err),norm(rerr));
        end
    end
end

end