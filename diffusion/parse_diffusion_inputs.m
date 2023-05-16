function [params] = parse_diffusion_inputs(varargin)
params = struct;
for i=1:2:numel(varargin)
    params.(varargin{i}) = varargin{i+1};
end
if (~isfield(params, 'sigma'))
    params.sigma = input('Prior model variance sigma = ? (default 1): ');
    if (isempty(params.sigma))
        params.sigma = 1;
    end
end
if (~isfield(params, 'corr_length'))
    params.corr_length = input('Prior model corr_length = ? (default 1): ');
    if (isempty(params.corr_length))
        params.corr_length = 1;
    end
end
if (~isfield(params, 'nu'))
    params.nu = input('Decay rate nu = ? (default 2): ');
    if (isempty(params.nu))
        params.nu = 2;
    end
end


if (~isfield(params, 'tol_kle'))
    params.tol_kle = input('KLE truncation tolerance (default 1e-2): ');
    if (isempty(params.tol_kle))
        params.tol_kle = 1e-2;
    end
end

% Here the inverse problem starts
if (~isfield(params, 'sigma_n'))
    params.sigma_n = input('Noise variance sigma_n = ? (default 1e-2, inf for prior inference): ');
    if (isempty(params.sigma_n))
        params.sigma_n = 1e-2;
    end
end
if (~isfield(params, 'm0'))
    params.m0 = input('Number of measurements in each direction m0 = ? (default 15): ');
    if (isempty(params.m0))
        params.m0 = 15;
    end
end
if (~isfield(params, 'y0'))
    params.y0 = input('Truth value of parameters for synthetic data y0 = ? (default 1.5, ''c'' for channel, ''b'' for barrier): ');
    if (isempty(params.y0))
        params.y0 = 1.5;
    end
end

if (~isfield(params, 'Nsamples'))
    params.Nsamples = input('Number of samples Nsamples = ? (default 10000): ');
    if (isempty(params.Nsamples))
        params.Nsamples = 10000;
    end
end

if (~isfield(params, 'thres'))
    params.thres = input('Threshold time thres defining the event t<thres = ? (default 0.15): ');
    if (isempty(params.thres))
        params.thres = 0.15;
    end
end

if (~isfield(params, 'gamma'))
    params.gamma = input('Sigmoid width coefficient gamma = ? (default 100/thres): ');
    if (isempty(params.gamma))
        params.gamma = 100/params.thres;
    end
end


if (~isfield(params, 'runs'))
    params.runs = input('Number of test runs = ? (default 1): ');
    if (isempty(params.runs))
        params.runs = 1;
    end
end

end
