% 1 Abiotrophia_defectiva_ATCC_49176.mat: (952, 1069)
% 2 Acidaminococcus_fermentans_DSM_20731.mat: (1009, 1090)
% 3 Acidaminococcus_intestini_RyC_MR95.mat: (917, 994)
% 4 Acidaminococcus_sp_D21.mat: (856, 851)
% 5 Acinetobacter_calcoaceticus_PHEA_2.mat: (1319, 1561)
% 6 Recon1.mat: (2766, 3742)
% 7 Recon2.mat: (5063, 7440)
% 8 Recon3.mat: (8399, 13543)
% 9 cardiac_mit_glcuptake_atpmax.mat: (230, 220)
% 10 ecoli_core_model.mat: (72, 95)
% 11 iAF1260.mat: (1668, 2382)
% 12 iJO1366.mat: (1805, 2583)
% 13 lp_25fv47.mat: (821, 1876) -- unbounded
% 14 lp_80bau3b.mat: (2262, 12061) -- unbounded
% 15 lp_cre_a.mat: (3516, 7248)
% 16 lp_gfrd_pnc.mat: (616, 1160) -- unbounded 
% 17 lp_israel.mat: (174, 316) -- unbounded 
% 18 lp_ken_18.mat: (105127, 154699)
% 19 lp_pilot_ja.mat: (940, 2267)
% 20 lp_sctap2.mat: (1090, 2500) -- unbounded 
% 21 lp_ship08l.mat: (778, 4363) -- unbounded 
% 22 lp_woodw.mat: (1098, 8418)

%% This file tests the rhmc package
% It generates files "result_(modelname)" with a struct with the following 
% properties in the folder "rhmc_test".
% 1. dim
% 2. ess
% 3. sampleTime
% 4. steps


filename0=strcat('../PolytopeSamplerMatlab-develop')
addpath(genpath(filename0));
addpath(genpath('../CHRR-test'));

maxNumCompThreads(1);
assert(maxNumCompThreads == 1);
curFolder = fileparts(mfilename('fullpath'));
datapath = fullfile(fileparts(curFolder), '/Instances/0raw', '*.mat');
matfiles = dir(datapath);

rhmc_result = struct;
idx = {[1:19]}; 


gamma0_list = [1e-1];
eta_list = [10];


for c = idx{1}
	idx_pb = c;
        fprintf('%d: %s\n', c, matfiles(c).name);
        disp(datetime('now'));
        model = load(fullfile(fileparts(curFolder), 'Instances/0raw', matfiles(c).name));
        model = cellfun(@(x)(model.(x)),fieldnames(model));
        P = struct;
        if matfiles(c).name(5:6) == "lp"
            if (~isfield(model,'A') || ~isfield(model, 'b'))
                error('You need to define both model.A and model.b');
            else
                P.Aeq = model.A;
                P.beq = model.b;
            end
            if isfield(model.aux, 'lo')
                P.lb = max(-1e7, model.aux.lo);
            end
            if isfield(model.aux, 'hi')
                P.ub = min(model.aux.hi, 1e7);
            end
        else
            if (~isfield(model,'S') || ~isfield(model,'b'))
                error('You need to define both model.S and model.b');
            else
                P.Aeq = model.S;
                P.beq = model.b;
            end
            if isfield(model,'lb')
                P.lb = model.lb;
            end
            if isfield(model,'ub')
                P.ub = model.ub;
            end
            if isfield(model,'dsense')
                I = (model.dsense == 'E');
                P.Aeq = [P.Aeq; model.C(I,:)];
                P.beq = [P.beq; model.d(I)];
                P.Aineq = model.C(~I,:);
                P.bineq = model.d(~I,:);
                flip = 1-2*(model.dsense(~I) == 'G');
                P.Aineq = flip.*P.Aineq;
                P.bineq = flip.*P.bineq;
            end
        end
        
        fprintf('%d %s: (%d, %d)\n', c, matfiles(c).name, size(P.Aeq, 1), size(P.Aeq, 2));
for eta = eta_list
	for gamma0 = gamma0_list
     try
        opts = default_options_GL();
        opts.initalStepSize = gamma0;
	opts.implicitTol_RC = eta;
        opts.maxStep = +Inf;
        opts.thin = +Inf;
        opts.module = {'MixingTimeEstimator', 'MemoryStorage', 'DynamicRegularizer', 'DynamicStepSize', 'DynamicWeight', 'ProgressBar', 'DebugLogger'};
        opts.MemoryStorage.memoryLimit = 6*1024*1024*1024;
        opts.ProgressBar.refreshInterval = 300;
        if c == 22
            opts.maxTime = 3600*24*3;
        else
            opts.maxTime = 24*3600;
        end
        opts.seed = 1;
        opts.logging = strcat('testpaper_bhmc_v2_', num2str(idx_pb), '.log'); % Output the debug log to testpaper_rhmc.log
        
        o = sample_GL_v2(P, 1000, opts);

        exps = struct; exps.dim = size(P.Aeq, 2); 
        exps.ess = min(effective_sample_size(o.samples));
        exps.roundTime = o.prepareTime;
        exps.sampleTime = o.sampleTime;
        exps.step = o.totalStep * o.opts.simdLen;
        fprintf('ESS: %f\n', exps.ess);
        filename = strcat('gamma0', num2str(gamma0,5));
        save(fullfile(curFolder, strcat('/bhmc_test/result_rhmc_v2_', filename, matfiles(c).name)), 'exps')
        clear('exps');
        clear('o');
    catch e
        fprintf('Error occurs. Skip model %d\n', c);
        fprintf(1,'The identifier was:\n%s\n', e.identifier);
        fprintf(1,'There was an error! The message was:\n%s\n', e.message);
        continue;
        
    end
end
end
end
