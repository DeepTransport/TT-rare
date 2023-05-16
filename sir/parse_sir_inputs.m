function [params] = parse_sir_inputs(varargin)
% Parse model parameters
params = struct;
for i=1:2:numel(varargin)
    params.(varargin{i}) = varargin{i+1};
end

if (~isfield(params, 'K'))
    params.K = input('Number of compartments K = ? (default 1): ');
    if (isempty(params.K))
        params.K = 1;
    end
end

if (~isfield(params, 'sigma_n'))
    params.sigma_n = input('Standard deviation of the observation noise sigma_n = ? (default 1): ');
    if (isempty(params.sigma_n))
        params.sigma_n = 1; % Noise std
    end
end

params.Imaxdef = 80;
if (isfield(params, 'model')) && (strncmpi(params.model, 'a', 1))
    params.Imaxdef = 69; % 80 is too large for Austrian topology
end

if (~isfield(params, 'Imax'))
    params.Imax = input(sprintf('Exceedance threshold Imax = ? (default %g): ', params.Imaxdef));
    if (isempty(params.Imax))
        params.Imax = params.Imaxdef;
    end
end

if (~isfield(params, 'gammaev'))
    params.gammaev = input('Sigmoid width coefficient gammaev = ? (default 3e3/Imax): ');
    if (isempty(params.gammaev))
        params.gammaev = 3e3/params.Imax;
    end
end

if (~isfield(params, 'Nsamples'))
    params.Nsamples = input('Number of samples Nsamples = ? (default 2^14): ');
    if (isempty(params.Nsamples))
        params.Nsamples = 2^14;
    end
end

if (~isfield(params, 'runs'))
    params.runs = input('Number of test runs = ? (default 1): ');
    if (isempty(params.runs))
        params.runs = 1;
    end
end

end
