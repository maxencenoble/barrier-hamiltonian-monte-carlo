%% This file tests the rhmc package on benchmark instances
% It generates files "result_(modelname)_(dim)" with a struct with the following 
% properties in the folder "rhmc_test".
% 1. dim
% 2. ess
% 3. preparation time
% 4. (total) sample time
% 5. (total) # steps


function testpaper_bias_cluster_GL_cube_10(idx)
filename0=strcat('../PolytopeSamplerMatlab-develop');
addpath(genpath(filename0));
addpath(genpath('../CHRR-test'));
maxNumCompThreads(1);
assert(maxNumCompThreads == 1);

curFolder = fileparts(mfilename('fullpath'));
instances = {'box', 'simplex', 'birkhoff'};


dimBoxSimplex = [10];
dimBirkhoff = [9 16];
numSamples = 1e6;
maxIter = 5e5;
nb_repet = 10;
thin = 1;
thinOutput=false;

%%%
% generalized leapfrog
for c = [idx]
    inst = instances{c};
    if strcmp('box', inst) || strcmp('simplex', inst)
        dimInput = dimBoxSimplex;
    elseif strcmp('birkhoff', inst)
        dimInput = dimBirkhoff;
    else
        error('Not supported')
    end

    for dim = dimInput
        fprintf('%s: %d\n', inst, dim);
        disp(datetime('now'));
        P = loadProblem(['basic/' inst '@' int2str(dim)]);
        unit_vec = 1:dim;
        unit_vec(1) = 0;
        unit_vec = unit_vec'/norm(unit_vec);
        unit_vec(2) = 1;
        unit_vec = 10*unit_vec;
        f = @(x) sum((x-unit_vec).^2)/2;
        df = @(x) (x-unit_vec); 
        
        
        %P.f = @(x) sum(unit_vec .*x);
        %P.df = @(x) unit_vec;
        P.f = f;
        P.df = df;
        opts = default_options_GL();
        opts.module = {'MixingTimeEstimator', 'MemoryStorage', 'DynamicRegularizer', 'DynamicStepSize', 'DynamicWeight', 'ProgressBar', 'DebugLogger'};
        opts.MemoryStorage.memoryLimit = 8*1024*1024*1024;
        opts.simdLen = 1;
        opts.ProgressBar.refreshInterval = 600;
        opts.maxTime = 3600*24;
        opts.maxStep = maxIter;
        opts.initialStepsize = 1e-1;
        opts.thin = thin;
        opts.implicitTol_RC = 10;
        opts.MemoryStorage.thinOutput = thinOutput;
        opts.logging = 'testpaper_bias_GL_benchmark.log'; % Output the debug log to testpaper_benchmark.log
        
        for i=1:nb_repet
            rng(1+i);
            opts.seed = randi(1024)+1+i;
            o = sample_GL(P, numSamples, opts);

            exps = struct; 
            exps.dim = size(P.Aeq, 2); 
            exps.ess = size(o.samples, 2);
            exps.numelA = numel(P.Aeq);
            exps.nnzA = nnz(P.Aeq);
            exps.roundTime = o.prepareTime;
            exps.sampleTime = o.sampleTime;
            exps.step = o.totalStep * o.opts.simdLen;
            exps.samples = o.samples;
            exps.total_samples = o.total_samples;
            fprintf('ESS: %f\n', exps.ess);

            save(fullfile(curFolder, strcat('bias_analysis_local/', 'GL_', inst, int2str(dim),'_rep_',num2str(i)))        , 'exps');
        end
end
end

