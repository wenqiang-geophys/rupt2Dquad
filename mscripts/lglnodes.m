function [xi,wts] = lglnodes(N)
% Legendre-Gauss-Lobatto nodes and weights
N1 = N+1;
x = cos( pi*(0:N)/(N) )';
P = zeros(N1,N1);

xold = 2;

while max(abs(x-xold)) > 1e-8
    xold = x;
    P(:,1) = 1;
    P(:,2) = x;

    for k = 2:N
        P(:,k+1) = ((2*k-1)*x.*P(:,k) - (k-1)*P(:,k-1))/k;
    end

    x = xold - (x.*P(:,N1) - P(:,N)) ./ (N1*P(:,N1));
end

w = 2.0./(N*N1*P(:,N1).^2);
xi = 0.5*(x(end:-1:1)+1); % convert [-1,1] to [0,1]
wts = 0.5*w(end:-1:1);
end