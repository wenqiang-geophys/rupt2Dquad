function [r1,s1] = xy2rs(x,y,VX,VY)

% syms r s;
% 
% f = (1-r)*(1-s)*VX(1)+r*(1-s)*VX(2)+r*s*(VX(3))+(1-r)*s*VX(4) - x;
% g = (1-r)*(1-s)*VY(1)+r*(1-s)*VY(2)+r*s*(VY(3))+(1-r)*s*VY(4) - y;
% 
% sol = solve(f,g,r,s);
% 
% r1 = double(sol.r);
% s1 = double(sol.s);


F = @(r)[...
    (1-r(1))*(1-r(2))*VX(1)+r(1)*(1-r(2))*VX(2)+r(1)*r(2)*(VX(3))+(1-r(1))*r(2)*VX(4)-x;
    (1-r(1))*(1-r(2))*VY(1)+r(1)*(1-r(2))*VY(2)+r(1)*r(2)*(VY(3))+(1-r(1))*r(2)*VY(4)-y;
    ];

rs=fsolve(F,[0.5;0.5]);
%rs=bisection(F,[0;0],[1;1],[0;0]);
%rs=bisection(F,[0;0],[1;1]);
r1=rs(1);
s1=rs(2);

end
