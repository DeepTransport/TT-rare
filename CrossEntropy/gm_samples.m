function samples = gm_samples(gm, N)

s  = size(gm.mean,1);
NN = mnrnd(N, gm.alpha);
samples = zeros(s, sum(NN));

n = 0;
for i = 1:length(gm.alpha)
    ind = (1:NN(i)) + n;
    samples(:,ind) = gm.mean(:,i) + (gm.L{i}*randn(s,NN(i)));
    n = n+NN(i);
end

end