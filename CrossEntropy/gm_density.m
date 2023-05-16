function [f, fi] = gm_density(gm, samples)

[s,n] = size(samples);
K = length(gm.alpha);

lf = zeros(K, n);
for i = 1:K
    %L = chol(gm.cov{i} + eye(s)*1E-8)';
    %v = L\samples';
    v = gm.L{i}\(samples-gm.mean(:,i));
    lf(i,:) = -0.5*sum(v.^2, 1) - sum(log(gm.S{i})) + log(gm.alpha(i)) - (s/2)*log(2*pi);
end
%ref = max(lf, [], 2);
%logf = log(sum(exp(lf - ref),2)) + ref;

fi = exp(lf);
f  = sum(fi,1);
fi = fi./f;

end