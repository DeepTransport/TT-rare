% Quantity of interest function (Infected) in full space for Cross Entropy
function [x] = InfectedFun(u, W, x0, tfine)
betagamma = 1 + erf(u/sqrt(2)); 

d = size(W,1);
beta = betagamma(:,1:2:2*d-1)';
beta = reshape(beta, 1, d, []);
gamma = betagamma(:,2:2:2*d)';
gamma = reshape(gamma, 1, d, []);
x0 = repmat(x0, 1, size(beta,3));
[t,x] = ode45(@(t,x)sir_rhs(x,W,beta,gamma), tfine, x0, odeset('AbsTol', 1e-6, 'RelTol', 1e-6));
x = reshape(x, [], 3*d, size(beta,3));
x = x(:, 2+(d-1)*3, :);
x = max(x,[],1); % Check for the maximum number of infected in last compart
x = x(:);
end
