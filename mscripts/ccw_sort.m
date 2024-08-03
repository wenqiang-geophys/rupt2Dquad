function [elem1] = ccw_sort(node,elem)

% elem (nelem,4)
% node (nnode,2)

nelem = size(elem,1);
%nnode = size(node,1);

elem1 =  elem;

for i = 1:nelem
    x = node(elem(i,:),1);
    y = node(elem(i,:),2);
    
    x0 = mean(x);
    y0 = mean(y);
    
    ang = atan2(y-y0,x-x0);
    
    [~,idx]=sort(ang);
    
    elem1(i,:) = elem(i,idx);
    
end

end