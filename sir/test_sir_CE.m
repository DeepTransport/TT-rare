% Cross Entropy SIR rare event sampling test for chain topology
function test_sir_CE(varargin)
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
% Extra parameters (only for CE)
if (~isfield(params, 'Ng'))
    params.Ng = input('Number of Gaussians in the mixture Ng = ? (default 4): ');
    if (isempty(params.Ng))
        params.Ng = 4;
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


% Posterior function for CE - it can only sample well in full space -> transform
postfun = @(u)exp(sir_ll(1 + erf(u/sqrt(2))', data(:), W, params.sigma_n, x0, tobs) - 0.5*sum(u.^2, 1)')';

% Quantity of interest function for CE
qoi = @(bg)InfectedFun(bg', W, x0, tfine)';

Z_post = zeros(params.runs, 1);
ttimes_post = zeros(params.runs, 1);
n_iter_post = zeros(params.runs, 1);
ess_post = zeros(params.runs, 1);

ttimes_event = zeros(params.runs, 1);
n_iter_event = zeros(params.runs, 1);
ess_event = zeros(params.runs, 1);

P_event = zeros(params.runs, 1);

for irun=1:params.runs
    % Step 1: EM approximation for the normalising constant
    tic;
    gm0 = gm_init(2*d, params.Ng, 1E-2);
    [gm0,Z_post(irun),ess_post(irun),n_iter_post(irun)] = cross_entropy(@(u)ones(1,size(u,2)), postfun, 0, gm0, params.em_iter, params.Nsamples, params.rho);
    ttimes_post(irun) = toc;

    % Step 1: EM approximation for the event probability
    tic;
    gm1 = gm_init(2*d, params.Ng, 1E-2);
    [gm1,P_event(irun),ess_event(irun),n_iter_event(irun)] = cross_entropy(qoi, postfun, params.Imax, gm1, params.em_iter, params.Nsamples, params.rho);
    P_event(irun) = P_event(irun) / Z_post(irun)
    ttimes_event(irun) = toc;
end

% Print some statsy information
fprintf('Cross Entropy SIR completed. Some average values:\n');
fprintf('\tCPU time of posterior CE: %g\n', mean(ttimes_post));
fprintf('\tNumber of iterations in posterior CE: %g\n', mean(n_iter_post));
fprintf('\tESS/N for posterior sampling: %g\n', mean(ess_post));
fprintf('\n');
fprintf('\tCPU time of event CE: %g\n', mean(ttimes_event));
fprintf('\tNumber of iterations in event CE: %g\n', mean(n_iter_event));
fprintf('\tESS/N for event sampling: %g\n', mean(ess_event));
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



% Log-likelihood transformed into full space
function [lF] = sir_ll_2(u, data, W, sigma_n, x0, tobs)
betagamma = 1 + erf(u/sqrt(2)); 

d = size(W,1);
beta = betagamma(:,1:2:2*d-1)';
beta = reshape(beta, 1, d, []);
gamma = betagamma(:,2:2:2*d)';
gamma = reshape(gamma, 1, d, []);
x0 = repmat(x0, 1, size(beta,3));
[t,x] = ode45(@(t,x)sir_rhs(x,W,beta,gamma), tobs, x0, odeset('AbsTol', 1e-6, 'RelTol', 1e-6));
x = x(2:end, :);
x = reshape(x, [], 3*d, size(beta,3));
x = x(:, 2+(0:d-1)*3, :);
x = reshape(x, [], size(beta,3));

lF = -sum((x - data).^2, 1)/(2*sigma_n^2) - 0.5*sum(u'.^2, 1);
lF = lF(:);

end






