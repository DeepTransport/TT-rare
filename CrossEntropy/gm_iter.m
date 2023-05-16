function gm = gm_iter(gm, samples, weights)
% given a set of samples (nxs matrix), and the gment mixture
% apply one iteration of EM to work out the next mixture
%
% next.alpha: non-negative weights (Kx1)
% next.mean:  mean vectors (Kxs)
% next.cov:   covariance (cell array of Kx1)
%
% We have a defensive term (K+1) which is not updated during EM
%

a = gm.alpha(end);

s = size(gm.mean, 1);
K = length(gm.alpha)-1;
n = size(samples, 2);

lf = zeros(K+1, n);
for i = 1:(K+1)
    v = gm.L{i}\(samples-gm.mean(:,i));
    lf(i,:) = - 0.5*sum(v.^2, 1) - sum(log(gm.S{i})) + log(gm.alpha(i));
end
P = exp(lf);
P = P./sum(P,1);

if nargin <= 2
    weights = ones(1,n);
end
weights = reshape(weights/sum(weights), 1, []);

gm.alpha = sum(P.*weights,2)';
gm.alpha(1:K) = gm.alpha(1:K)/sum(gm.alpha(1:K))*(1-a);
gm.alpha(end) = a;


for i = 1:K
    if gm.alpha(i) < max(a, gm.tol)
        gm.mean(:,i) = zeros(s,1);
        gm.cov{i} = eye(s);
        gm.L{i}   = eye(s);
        gm.S{i}   = ones(s,1);
    else
        tmp = (P(i,:).*weights)/sum(P(i,:).*weights);
        gm.mean(:,i) = samples*tmp(:);
        gm.cov{i} = (samples.*tmp)*samples' - gm.mean(:,i)*gm.mean(:,i)';
        [U,S,~] = svd(gm.cov{i});
        gm.S{i} = (diag(S) + gm.tol).^(0.5);
        gm.L{i} = U*diag(gm.S{i});
    end
end

end


