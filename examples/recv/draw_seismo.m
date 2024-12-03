clear
close all

id_list = [1:21];

for fid = id_list
    
    fnm = ['data/body_recv_',num2str(fid,'%06d'),'.txt'];

    a = load(fnm);
    t = a(:,1);
    vx = a(:,2);
    vy = a(:,3);

    plot(t,vy,'LineWidth',1)
    hold on

end

xlabel('Time (sec)')
ylabel('Vy (m/s)')