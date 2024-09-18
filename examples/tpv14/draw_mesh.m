clc
clear
close all

%% determine quad lists (QUAD4) and fault boundary curves (BAR2)
% first, ncdisp('your_mesh_file.exo') to list all "BAR2" and "QUAD4"
% then, use the following scripts to determine which connectX is the fault
% boundary curve
% alternatively, you can use paraview to do this, but I think this script
% is enough for users

fnm = 'tpv14_2d.exo';

coord = ncread(fnm,'coord');
node = coord(:,1:2)';
Nnode = size(node,1);

quad_list = [5];
elem = [];
for i = 1:length(quad_list)
    elem1 = ncread(fnm,['connect',num2str(quad_list(i))]);
    elem = cat(2, elem, elem1);
end

bdr_list = [1,2,3];

elem2 = [];
for i = 1:length(bdr_list)
    elem1 = ncread(fnm,['connect',num2str(bdr_list(i))]);
    elem2 = cat(2, elem2, elem1);
end
elem5 = elem([1:4,1],:);
x = reshape(node(1,elem5),size(elem5));
y = reshape(node(2,elem5),size(elem5));
plot(x,y,'k')
hold on
axis image
x = reshape(node(1,elem2),size(elem2));
y = reshape(node(2,elem2),size(elem2));
plot(x,y,'r','LineWidth',1)