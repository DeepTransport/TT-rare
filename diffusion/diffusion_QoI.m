function [Q, C, U] = diffusion_QoI(y, phil, bound, W2, Mass_summed)
% Sample QoI on the random vectors y (size I x d)

C = y.';
if (~isempty(phil))
    C = phil * C;
end
C = exp(C);
C = sparse(C);

I = size(C, 2);
n = round(sqrt(size(C, 1)));

% Boundary condition
bound_l = bound(1:numel(bound)/2);
bound_r = bound((numel(bound)/2+1):end);
y = (0:(1/(n-1)):1)';
left_bc = 1+0.5*y;
right_bc = -sin(2*pi*y)-1;
% left_bc = 1+0*y;
% right_bc = 0*y;
% Prolong matrices
Pb = speye(n^2);
Pb(:,bound) = [];
ud = zeros(n^2,1);
ud(bound_l) = left_bc;
ud(bound_r) = right_bc;

if (nargout==1)
    Q = zeros(size(Mass_summed,1), I);
else
    U = zeros(n^2, I);
end

% Prevent OOM
blocksize = 8192;
for nblock = 0:floor(I/blocksize)
    ind = (nblock*blocksize+1) : min((nblock+1)*blocksize, I);
    
    Ball = W2*C(:,ind);
    
    for i=ind
        % Assemble matrix from W2g and coeff
        B = reshape(Ball(:,i-ind(1)+1), n^2, n^2);
        % Eliminate the BC to the RHS
        g = B(:,bound_l)*left_bc + B(:,bound_r)*right_bc;
        g = -g;
        g(bound) = [];
        B(bound,:) = [];
        B(:,bound) = [];
        
        % Solve
        u = B\g;
        u = Pb*u;
        u = u+ud;
        
        if (nargout==1)
            Q(:,i) = Mass_summed * u;
        else
            U(:,i) = u;
        end
    end
    
end

if (nargout>1)
    Q = Mass_summed * U;
end

Q = Q.'; % we need I x R output

if (nargout>1)
    C = full(C);    
end
end
