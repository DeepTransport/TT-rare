% Log-likelihood of SIR
function [lF] = sir_ll(betagamma, data, W, sigma_n, x0, tobs)
d = size(W,1);
beta = betagamma(:,1:2:2*d-1)';
beta = reshape(beta, 1, d, []);
gamma = betagamma(:,2:2:2*d)';
gamma = reshape(gamma, 1, d, []);
x0 = repmat(x0, 1, size(beta,3));
[t,x] = ode45(@(t,x)sir_rhs(x,W,beta,gamma), tobs, x0, odeset('AbsTol', 1e-6, 'RelTol', 1e-6));
x = x(2:end, :);
x = reshape(x, [], 3*d, size(beta,3));
x = x(:, 2+(0:d-1)*3, :);
x = reshape(x, [], size(beta,3));

lF = -sum((x - data).^2, 1)/(2*sigma_n^2);
lF = lF(:);

end
