% Installation script for TT-rare

% Remember our directory
mydir = fileparts(mfilename('fullpath'));

% Check TT-IRT with prerequisites
try
    z = tt_irt_sqr({[0;1;2]; [0;1;2]}, {randn(1,3,2); randn(2,3,1)}, rand(10,2));
catch ME
    if (exist('TT-IRT-master', 'dir')==0)
        if (exist('TT-IRT-master.zip', 'file')==0)
            try
                fprintf('TT-IRT is not found, downloading...\n');
                opts = weboptions; opts.CertificateFilename=('');
                websave('TT-IRT-master.zip', 'https://github.com/dolgov/TT-IRT/archive/master.zip', opts);
            catch ME2
                error('%s. Automatic download failed. Please download TT-IRT from https://github.com/dolgov/TT-IRT', ME2.message);
            end
        end
        try
            unzip('TT-IRT-master.zip');
        catch ME2
            error('%s. Automatic unzipping failed. Please extract TT-IRT-master.zip here', ME2.message);
        end
    end
    cd('TT-IRT-master');
    cd('matlab');
    install;
end
check_tt;
check_mcmc(false); % We need UWerr but not DRAM

cd(mydir);
% Add our own paths
addpath('CrossEntropy');
addpath('diffusion');
addpath('sir');

