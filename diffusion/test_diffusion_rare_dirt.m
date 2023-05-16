% DIRT Diffusion rare event sampling test
function test_diffusion_rare_dirt(varargin)
% Check for and download TT-Toolbox
mydir = fileparts(mfilename('fullpath'));
try
    check_ttirt;
catch
    cd(mydir); cd('..'); cd('..'); cd('..'); cd('utils'); check_ttirt;
end
check_tt;
cd(mydir);

% Parse parameters or ask a user for them
params = parse_diffusion_inputs(varargin{:});
% Extra parameters (only for DIRT)
if (~isfield(params, 'npi'))
    params.npi = input('Uniform grid size for the posterior npi = ? (default 17): ');
    if (isempty(params.npi))
        params.npi = 17;
    end
end
if (~isfield(params, 'rpi'))
    params.rpi = input('TT rank of ratios R = ? (default 5): ');
    if (isempty(params.rpi))
        params.rpi = 5;
    end
end
if (~isfield(params, 'beta'))
    params.beta = input('Tempering powers (default 10.^(-4:0.5:0)): ');
    if (isempty(params.beta))
        params.beta = 10.^(-4:0.5:0);
    end
end
if (~isfield(params, 'pp'))
    params.pp = input('Power of prior tempering (default 0.2): ');
    if (isempty(params.pp))
        params.pp = 0.2;
    end
end


% Build the discretization and KLE
[~,pm,bound,W1g,W1m,spind,Pua,phi,lambda,Mass, W2] = build_kle_eig(2, 'DN', params.nu, params.corr_length, params.tol_kle, params.m0);

Mass_summed = cellfun(@(M)sum(M,1), Mass, 'uni', 0);
Mass_summed = reshape(Mass_summed, [], 1);
Mass_summed = cell2mat(Mass_summed);

% weighted KLE components
L = numel(lambda);
phil = full(phi*spdiags(sqrt(lambda), 0, L, L));

% Simulate some observations
if (exist(sprintf('Q_obs_nu%g_kle%g_sigman%g_m0%d_ytrue%s.mat', params.nu, -log10(params.tol_kle), params.sigma_n, params.m0, string(params.y0)), 'file')>0)
    % Load the same observations for all experiments
    fprintf('Found Q_obs file for nu=%g, tol_kle=%g, sn=%g, m0=%d, ytrue=%s\nRemove it to regenerate the observations\n', params.nu, params.tol_kle, params.sigma_n, params.m0, string(params.y0));
    load(sprintf('Q_obs_nu%g_kle%g_sigman%g_m0%d_ytrue%s.mat', params.nu, -log10(params.tol_kle), params.sigma_n, params.m0, string(params.y0)));
else
    fprintf('Generating Q_obs from the true parameters\n');
    if (params.y0=='b')
        % Synthetic semi-annulus coefficient
        C = ((pm(1,:)-0.7).^2 + (pm(2,:)-0.3).^2>0.4^2) ...
            & ((pm(1,:)-0.7).^2 + (pm(2,:)-0.3).^2<=0.5^2) ...
            & (pm(1,:)<=0.7) & (pm(2,:)>0.3);  % Semi-annulus inside the domain
        C = 0.01 + 100*double(C);
        Q_obs = diffusion_QoI(log(C), [], bound, W2, Mass_summed);        
    elseif (params.y0=='c')
        C = ((pm(1,:)-0.7).^2 + (pm(2,:)-0.3).^2>0.4^2) ...
            & ((pm(1,:)-0.7).^2 + (pm(2,:)-0.3).^2<=0.5^2) ...
            & (pm(1,:)<=0.7) & (pm(2,:)>0.3);
        C = C | ( (pm(1,:)<=0.35) & pm(2,:)>0.45 & pm(2,:) < 0.55 ) ...
            | ( (pm(1,:)>0.7) & pm(2,:)>0.72 & pm(2,:) < 0.78 );   % Channels to boundaries
        C = 0.01 + 100*double(C);
        Q_obs = diffusion_QoI(log(C), [], bound, W2, Mass_summed);        
    else
        Q_obs = diffusion_QoI(params.y0*ones(1,L), phil, bound, W2, Mass_summed);
        if (~isinf(params.sigma_n))
            Q_obs = Q_obs + randn(1,params.m0^2)*sqrt(params.sigma_n);
        end
    end
    save(sprintf('Q_obs_nu%g_kle%g_sigman%g_m0%d_ytrue%s.mat', params.nu, -log10(params.tol_kle), params.sigma_n, params.m0, string(params.y0)), 'Q_obs');
end

% Grid in parameters
ys = 20*sqrt(params.sigma)/(params.npi-1);
ys = ((-10*sqrt(params.sigma)) : ys : (10*sqrt(params.sigma)))';
ys = repmat({ys}, L, 1);


% Log posterior function
lpfun = @(theta,beta0,beta)-sum((diffusion_QoI(theta, phil, bound, W2, Mass_summed) - Q_obs).^2, 2)*(beta-beta0)/(2*params.sigma_n) - sum(theta.^2,2)*(beta^params.pp-beta0^params.pp)/(2*params.sigma);

% Rare event Log biasing function
lpfun_rare = @(theta,beta0,beta)- sum(theta.^2,2)*(beta^params.pp-beta0^params.pp)/(2*params.sigma) ...
                                + loglike_time(theta, beta0, beta, params.thres, phil, bound, W2, Mass_summed, Q_obs, params.sigma_n, params.gamma);
                            

Z_post = zeros(params.runs, 1);
ttimes_post = zeros(params.runs, 1);
evalcnt_post = zeros(params.runs, 1);
ess_post = zeros(params.runs, 1);
tau_post = zeros(params.runs, 1);
hell_post = zeros(params.runs, 1);

ttimes_event = zeros(params.runs, 1);
evalcnt_event = zeros(params.runs, 1);
ess_event = zeros(params.runs, 1);
tau_event = zeros(params.runs, 1);
hell_event = zeros(params.runs, 1);

P_event = zeros(params.runs, 1);


% Fire up Bayesian runs
for irun=1:params.runs
    %% Sample Posterior
    tic;
    IRTstruct = tt_dirt_approx(ys, lpfun, params.beta, 'testsamples', 1e3, ...
                               'nswp', 1, 'y0', ~isinf(params.sigma_n)*params.rpi + isinf(params.sigma_n), ...
                               'kickrank', 0, 'boundary', true, 'reference', 'n3', ...
                               'trunctol', 0, 'interpolation', 's', 'IRTdenom', false);
                           
    q = randref(IRTstruct.reference, 2^params.log2N, L);
                           
    [z,lFapp,lFex] = tt_dirt_sample(IRTstruct, q, @(x)lpfun(x,0,1));
    ttimes_post(irun) = toc;
    
    Z_post(irun) = mean(exp(lFex-lFapp))
    
    % Reject (a.k.a. check error)
    [z2,lFex2,lFapp2] = mcmc_prune(z,lFex,lFapp);
    
    tau_post(irun) = statsiact(z2)
    ess_post(irun) = essinv(lFex,lFapp)
    hell_post(irun) = hellinger(lFex,lFapp)
    evalcnt_post(irun) = sum(IRTstruct.evalcnt - 1e3);
    
    
    
    %% Rare event approximation
    tic;
    [IRTstruct2] = tt_dirt_approx(ys, lpfun_rare, params.beta, 'testsamples', 1e3, ...
                                   'nswp', 1, 'y0', params.rpi, 'kickrank', 0, 'boundary', true, ...
                                   'reference', 'n3', 'interpolation', 's');
                                                              
    [z2,lFapp2,lFex2] = tt_dirt_sample(IRTstruct2, q, @(x)lpfun_rare(x,0,1));
    ttimes_event(irun) = toc;
    
    [z3,lFex3,lFapp3] = mcmc_prune(z2,lFex2,lFapp2);
    tau_event_approx(irun) = statsiact(z3)
    ess_event_approx(irun) = essinv(lFex2,lFapp2)
    hell_event_approx(irun) = hellinger(lFex2,lFapp2)
    
    [lFex2, tz2] = loglike_time(z2, 0, 1, params.thres, phil, bound, W2, Mass_summed, Q_obs, params.sigma_n, 0);
    lFex2 = lFex2 - sum(z2.^2,2)/(2*params.sigma);   
            
    [z3,lFex3,lFapp3] = mcmc_prune(z2,lFex2,lFapp2);
    tau_event(irun) = statsiact(z3)
    ess_event(irun) = essinv(lFex2,lFapp2)
    hell_event(irun) = hellinger(lFex2,lFapp2)
    evalcnt_event(irun) = sum(IRTstruct2.evalcnt - 1e3);
    
    P_event(irun) = mean(exp(lFex2 - lFapp2)) / Z_post(irun)    
end % irun



% Print some statsy information
fprintf('DIRT Diffusion completed. Some average values:\n');
fprintf('\tCPU time of posterior DIRT: %g\n', mean(ttimes_post));
fprintf('\tNumber of evaluations in posterior DIRT: %g\n', mean(evalcnt_post));
fprintf('\tIACT for posterior sampling: %g\n', mean(tau_post));
fprintf('\tN/ESS for posterior sampling: %g\n', mean(ess_post));
fprintf('\tD_H error for posterior sampling: %g\n', mean(hell_post));
fprintf('\n');
fprintf('\tCPU time of event DIRT: %g\n', mean(ttimes_event));
fprintf('\tNumber of evaluations in event DIRT: %g\n', mean(evalcnt_event));
fprintf('\tIACT for event sampling: %g\n', mean(tau_event));
fprintf('\tN/ESS for event sampling: %g\n', mean(ess_event));
fprintf('\tD_H error for event sampling: %g\n', mean(hell_event));
fprintf('\n');
fprintf('P[t<thres]: %g\n', mean(P_event));

% Copy vars to main space
vars = whos;
for i=1:numel(vars)
    if (exist(vars(i).name, 'var'))
        assignin('base', vars(i).name, eval(vars(i).name));
    end
end
end
