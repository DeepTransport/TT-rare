function [F] = flux_interpolate_vec(X, a1, u1)
% Interpolates flux at the given 2D point X.
% a1, u1 should be M x n^2
% X should be M x 2

n = round(sqrt(size(a1,2)));
M = size(a1,1);
xsf = linspace(0,1,n)';

% Binary search for the closest to X points
yk = X(:); % size M*2
i0 = ones(M*2,1);
i2 = n*ones(M*2,1);
active_ind = (1:(2*M))';
y_active = yk;
ia = i0;
while (~isempty(active_ind))  % (any((i2-i0)>1))
    i1 = floor((ia+i2)/2);
    ind_left = y_active > xsf(i1); % at these indices i0 becomes i1
    ia(ind_left) = i1(ind_left);
    i2(~ind_left) = i1(~ind_left);
    i0(active_ind) = ia;
    not_active = (i2-ia)<2;
    active_ind(not_active) = [];
    i2(not_active) = [];
    ia(not_active) = [];
    y_active(not_active) = [];
end

x1 = xsf(i0);  % grid points at i0, i0+1
x2 = xsf(i0+1);
h3 = x2 - x1; % grid intervals
% Linear interpolation coefficients into xk
Aq = (x2-yk)./h3;
Bq = (yk-x1)./h3;
% Gradient coefficients into xk
gAq = -1./h3;
gBq =  1./h3;

i0 = reshape(i0, M, 2);
Aq = reshape(Aq, M, 2);
Bq = reshape(Bq, M, 2);
gAq = reshape(gAq, M, 2);
gBq = reshape(gBq, M, 2);

i0lb = i0(:,1) + (i0(:,2)-1)*n;
i0rb = i0(:,1)+1 + (i0(:,2)-1)*n;
i0lt = i0(:,1) + (i0(:,2)+1-1)*n;
i0rt = i0(:,1)+1 + (i0(:,2)+1-1)*n;

aint = tracemult(a1, i0lb);  % aint = a1(i0(1,:), i0(2,:), :);
ulb = tracemult(u1, i0lb); % u1(i0(1,:),i0(2,:),:)
urb = tracemult(u1, i0rb); % u1(i0(1,:)+1,i0(2,:),:)
ult = tracemult(u1, i0lt); % u1(i0(1,:),i0(2,:)+1,:)
urt = tracemult(u1, i0rt); % u1(i0(1,:)+1,i0(2,:)+1,:)

F = zeros(M,2);
% du / dx
F(:,1) = ulb.*gAq(:,1).*Aq(:,2) + urb.*gBq(:,1).*Aq(:,2) ...
       + ult.*gAq(:,1).*Bq(:,2) + urt.*gBq(:,1).*Bq(:,2);
% du / dy
F(:,2) = ulb.*Aq(:,1).*gAq(:,2) + urb.*Bq(:,1).*gAq(:,2) ...
       + ult.*Aq(:,1).*gBq(:,2) + urt.*Bq(:,1).*gBq(:,2);
% Multiply with a - use piecewise constant now

F = - F .* aint;
end
