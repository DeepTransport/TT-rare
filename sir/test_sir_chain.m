% DIRT SIR rare event sampling test for chain topology
function test_sir_chain(varargin)
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
params = parse_sir_inputs(varargin{:});
% Extra parameters (only for DIRT)
if (~isfield(params, 'npi'))
    params.npi = input('Uniform grid size for the posterior npi = ? (default 17): ');
    if (isempty(params.npi))
        params.npi = 17;
    end
end
if (~isfield(params, 'rpi'))
    params.rpi = input('TT rank of ratios R = ? (default 7): ');
    if (isempty(params.rpi))
        params.rpi = 7;
    end
end
if (~isfield(params, 'beta'))
    params.beta = input('Tempering powers (default 10.^(-4:1/3:0)): ');
    if (isempty(params.beta))
        params.beta = 10.^(-4:1/3:0);
    end
end


% number of compartments
d = params.K; 
% Adjacency matrix of the diffusion
W = spdiags(ones(d,1)*[1 -2 1], -1:1, d, d);
W(1,d) = 1;
W(d,1) = 1;
W = W*0.5;
if (d==2)
    W = [-1 1; 1 -1];
end
if (d==1)
    W = 0;
end
d = size(W,1);

% Simulate the exact trajectory from some known parameters
beta = 0.1*ones(1,d);
gamma = 1*ones(1,d);
% Try different initial states in different compartments
x0 = [99; 0; 0] + [0;1;0]*(d:-1:1);
x0(1,:) = 100-x0(2,:);
x0 = x0(:);
% Observation time points
tobs = linspace(0,5,7);
[t,x] = ode45(@(t,x)sir_rhs(x,W,beta,gamma), tobs, x0, odeset('AbsTol', 1e-6, 'RelTol', 1e-6));

% fine discretization we need for events
tfine = linspace(0,5,25*6+1);
ind_obs = 25+1:25:25*6+1;
assert(norm(tfine(ind_obs)-tobs(2:end))==0)

data = x(2:end, 2+(0:d-1)*3);
data = data(:);

% Grid in parameters
xsf = repmat({linspace(0,2,params.npi)'}, 2*d, 1);

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

for irun=1:params.runs
    %% Posterior density
    tic;
    IRT = tt_dirt_approx(xsf, @(bg,b1,b2)sir_ll(bg, data, W, params.sigma_n, x0, tobs)*(b2-b1), ...
                         params.beta, 'nswp', 1, 'kickrank', 0, 'y0', params.rpi, 'testsamples', 1e2, ...
                         'reference', 'n3', 'boundary', true, 'interpolation', 's');

    q = randref(IRT.reference, params.Nsamples, 2*d);
    
    [z,lFapp] = tt_dirt_sample(IRT, q); 
    ttimes_post(irun) = toc;

    lFex = sir_ll(z, data, W, params.sigma_n, x0, tobs);
    % Normalising constant
    Z_post(irun) = mean(exp(lFex-lFapp))
    % Error estimate of the posterior
    z2 = mcmc_prune(z, lFex, lFapp);
    tau_post(irun) = statsiact(z2)
    ess_post(irun) = essinv(lFex,lFapp)
    hell_post(irun) = hellinger(lFex,lFapp)
    evalcnt_post(irun) = sum(IRT.evalcnt - 1e2);


    % Simulate posterior samples
    [t,x] = ode45(@(t,x)sir_rhs(x,W, reshape(z(:,1:2:2*d-1)',1,d,[]) , reshape(z(:,2:2:2*d)',1,d,[]) ), tfine, repmat(x0, 1, size(z,1)), odeset('AbsTol', 1e-6, 'RelTol', 1e-6));
    x = reshape(x, [], 3*d, size(z,1));
    x = x(:, 2+(d-1)*3, :);
    x = reshape(x, [], size(z,1));
    x = max(x, [], 1);
    x = x(:);

    % Direct importance sampling probability estimation
    P_event_direct(irun) = mean(exp(lFex-lFapp).*double(x>params.Imax))/Z_post(irun)
    
    % Direct MCMC probability estimation
    [zmcmc,lFexmcmc] = mcmc_prune(z,[lFex x],lFapp);
    P_event_mcmc(irun) = mean(double(lFexmcmc(:,2)>params.Imax))


    %% Event biasing
    tic;
    IRT2 = tt_dirt_approx(xsf, @(bg,b1,b2)ll_logsigmoid(bg, data, W, params.sigma_n, x0, tfine, ind_obs, params.Imax, params.gammaev, b1, b2), ...
                          params.beta, 'nswp', 1, 'kickrank', 0, 'y0', params.rpi, 'testsamples', 1e2, ...
                          'reference', IRT.reference, 'boundary', true, 'interpolation', 's');
    
    [zR,lFappR] = tt_dirt_sample(IRT2, q);
    ttimes_event(irun) = toc;
    
    lFexR = sir_ll(zR, data, W, params.sigma_n, x0, tobs);
    lFsmooth = ll_logsigmoid(zR, data, W, params.sigma_n, x0, tfine, ind_obs, params.Imax, params.gammaev, 0, 1);
    tau_event_approx(irun) = statsiact(mcmc_prune(zR, lFsmooth, lFappR))
    ess_event_approx(irun) = essinv(lFsmooth, lFappR)
    hell_event_approx(irun) = hellinger(lFsmooth, lFappR)
    evalcnt_event(irun) = sum(IRT2.evalcnt - 1e2);

    % Simulate infected trajectories on zR
    [t,x] = ode45(@(t,x)sir_rhs(x,W, reshape(zR(:,1:2:2*d-1)',1,d,[]) , reshape(zR(:,2:2:2*d)',1,d,[]) ), tfine, repmat(x0, 1, size(zR,1)), odeset('AbsTol', 1e-6, 'RelTol', 1e-6));
    x = reshape(x, [], 3*d, size(zR,1));
    x = x(:, 2+(d-1)*3, :);
    x = reshape(x, [], size(zR,1));
    x = max(x, [], 1);
    x = x(:);

    zmcmc = mcmc_prune(zR, lFexR + log(double(x>params.Imax)), lFappR);
    tau_event(irun) = statsiact(zmcmc)
    ess_event(irun) = essinv(lFexR + log(double(x>params.Imax)), lFappR)
    hell_event(irun) = hellinger(lFexR + log(double(x>params.Imax)), lFappR)

    P_event(irun) = mean(exp(lFexR-lFappR).*double(x>params.Imax))/Z_post(irun)
end

% Print some statsy information
fprintf('DIRT chain SIR completed. Some average values:\n');
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
fprintf('P[I>Imax]: %g\n', mean(P_event));

% Copy vars to main space
vars = whos;
for i=1:numel(vars)
    if (exist(vars(i).name, 'var'))
        assignin('base', vars(i).name, eval(vars(i).name));
    end
end
end
