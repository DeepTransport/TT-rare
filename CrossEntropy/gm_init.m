function gm = gm_init(s, K, def)

gm.tol = 1E-4;

gm.alpha = ones(1,K+1);
gm.alpha(1:K) = gm.alpha(1:K)/sum(gm.alpha(1:K))*(1-def);
gm.alpha(end) = max(def, gm.tol);

for i = 1:K
    gm.mean(:,i) = randn(s,1);
    gm.cov{i} = eye(s);
    gm.L{i}   = eye(s);
    gm.S{i}   = ones(s,1);
end

gm.mean(:,K+1) = zeros(1,s);
gm.cov{K+1} = eye(s);
gm.L{K+1}   = eye(s);
gm.S{K+1}   = ones(s,1);
    
end