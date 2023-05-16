% Cross Entropy Diffusion rare event sampling test
function test_diffusion_rare_CE(varargin)
% Check for and download TT-IRT
mydir = fileparts(mfilename('fullpath'));
try
    check_ttirt;
catch
    warning('TT-IRT is not found. Running install from the root folder TT-rare...');
    cd('..');
    install;
end
cd(mydir);

% Parse parameters or ask a user for them
params = parse_diffusion_inputs(varargin{:});
% Extra parameters (only for CE)
if (~isfield(params, 'K'))
    params.K = input('Number of Gaussians in the mixture K = ? (default 1): ');
    if (isempty(params.K))
        params.K = 1;
    end
end
if (~isfield(params, 'em_iter'))
    params.em_iter = input('Maximal number of iterations em_iter = ? (default 30): ');
    if (isempty(params.em_iter))
        params.em_iter = 30;
    end
end
if (~isfield(params, 'rho'))
    params.rho = input('1-quantile rho = ? (default 0.1): ');
    if (isempty(params.rho))
        params.rho = 0.1;
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


% Log posterior function
postfun = @(theta)exp(-sum((diffusion_QoI(theta', phil, bound, W2, Mass_summed) - Q_obs).^2, 2)/(2*params.sigma_n) - sum((theta').^2,2)/(2*params.sigma))';

% Quantity of interest function for CE: negative escape time
qoi = @(theta)m_escape_time(theta', phil, bound, W2, Mass_summed)';

Z_post = zeros(params.runs, 1);
ttimes_post = zeros(params.runs, 1);
n_iter_post = zeros(params.runs, 1);
ess_post = zeros(params.runs, 1);

ttimes_event = zeros(params.runs, 1);
n_iter_event = zeros(params.runs, 1);
ess_event = zeros(params.runs, 1);

P_event = zeros(params.runs, 1);


% Fire up Bayesian runs
for irun=1:params.runs
    tic;
    % Posterior
    gm0 = gm_init(L, params.K, 1E-2);
    qoi0 = @(theta)ones(1,size(theta,2));
    [gm0,Z_post(irun),ess_post(irun),n_iter_post(irun)] = cross_entropy(qoi0, postfun, 0, gm0, params.em_iter, params.Nsamples, params.rho, []);
    ttimes_post(irun) = toc;
 
    % CE Event
    tic;
    gm1 = gm_init(L, params.K, 1E-2);
    [gm1,P_event(irun),ess_event(irun),n_iter_event(irun)] = cross_entropy(qoi, postfun, -params.thres, gm1, params.em_iter, params.Nsamples, params.rho, []);
    P_event(irun) = P_event(irun) / Z_post(irun)
    ttimes_event(irun) = toc;
end % irun



% Print some statsy information
fprintf('Cross Entropy Diffusion completed. Some average values:\n');
fprintf('\tCPU time of posterior CE: %g\n', mean(ttimes_post));
fprintf('\tNumber of iterations in posterior CE: %g\n', mean(n_iter_post));
fprintf('\tESS/N for posterior sampling: %g\n', mean(ess_post));
fprintf('\n');
fprintf('\tCPU time of event CE: %g\n', mean(ttimes_event));
fprintf('\tNumber of iterations in event CE: %g\n', mean(n_iter_event));
fprintf('\tESS/N for event sampling: %g\n', mean(ess_event));
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


% Negative escape time for Cross Entropy
function [mt] = m_escape_time(theta, phil, bound, W2, Mass_summed)
[~, a, u] = diffusion_QoI(theta, phil, bound, W2, Mass_summed);
a = a.';
u = u.';  % size N x n^2

mt = escape_time_vec(a, u);
mt = -mt;
end

