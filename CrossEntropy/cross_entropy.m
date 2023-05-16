function [gm,prob,ess,k] = cross_entropy(qoi, target, a, gm, em_iter, n, rho, delta)

ess_prev = 0;
for k=1:em_iter
    % sample the current gm
    samples = gm_samples(gm, n);
    
    % compute weights
    gmf = gm_density(gm, samples);
    tf = target(samples);
    log_weights = log(tf) - log(gmf);
    %
    mlw = max(log_weights);
    weights = exp(log_weights - mlw);
        
    % compute the performance function
    Qs = qoi(samples);
    
    % ess for the QoI
    weights_ess = weights.*(reshape(Qs,1,[])>a);
    ess = mean(weights_ess).^2/mean(weights_ess.^2);
    
    % 1-rho quantile
    a_k = quantile(Qs, 1-rho);
    a_k = min(a_k, a);
    
    fprintf('\n\nk=%2d, a_k=%3.3e, ess/N=%3.3e \n', ...
        k, a_k, ess);
    
    weights = weights.*(reshape(Qs,1,[])>a_k);
    
    % EM iterations
    for i = 1:em_iter
        gm = gm_iter(gm, samples, weights);
    end
    
    if (ess<ess_prev)
        weights = weights_prev;
        mlw = mlw_prev;
        break;
    end
    ess_prev = ess;
    weights_prev = weights;
    mlw_prev = mlw;
    
    
%     if a_k >= (a-1E-10)
%         break;
%     end
        
end

%{
% sample the current gm
samples = gm_samples(gm, n);

% compute weights
gmf = gm_density(gm, samples);
tf = target(samples);
log_weights = log(tf) - log(gmf);
weights = exp(log_weights - max(log_weights));

% compute the performance function
Qs = qoi(samples);

% 1-rho quantile
a_k = quantile(Qs, 1-rho);

% ess for the QoI
weights_ess = weights.*(reshape(Qs,1,[])>a);
ess = mean(weights_ess).^2/mean(weights_ess.^2);
%}

ess = mean(weights).^2/mean(weights.^2);

prob = mean(weights)*exp(mlw);

fprintf('\n\nfinal ess/N=%3.3e, (unnormalised) failure prob=%3.3e \n', ess, prob);


end
