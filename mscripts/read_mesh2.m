function [node,elem] = read_mesh2(fnm)
%node = load([fnm,'.node']);
%elem = load([fnm,'.elem']);
%
%node = node(:,2:3);
%if elem(1,2) == 2
%    elem = elem(:,6:8);
%    nface = 3;
%elseif elem(1,2) == 3
%    elem = elem(:,6:9);
%    nface = 4;
%end

fid = fopen([fnm,'.msh']);

for i = 1:4
fgetl(fid);
end

nnode=fscanf(fid,'%d\n',1);
node = zeros(nnode,2);
for i = 1:nnode
a=fscanf(fid,'%d %g %g %g\n',4);
node(i,:) = a(2:3);
end

for i = 1:2
fgetl(fid);
end

nelem=fscanf(fid,'%d\n',1);
elem = zeros(nelem,4);

k = 0;
for i = 1:nelem
   a = fgetl(fid);
   b = sscanf(a,'%d',2);
   if(b(2) == 3)
       k = k+1;
       c = sscanf(a,'%d',9);
       elem(k,:) = c(6:9);
   end
end

nelem = k;
elem = elem(1:nelem,:);

fclose(fid);

end
