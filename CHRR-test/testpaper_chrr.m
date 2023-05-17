% 1 Abiotrophia_defectiva_ATCC_49176.mat: (952, 1069) -> 157
% 2 Acidaminococcus_fermentans_DSM_20731.mat: (1009, 1090) -> 164
% 3 Acidaminococcus_intestini_RyC_MR95.mat: (917, 994) -> 123
% 4 Acidaminococcus_sp_D21.mat: (856, 851) -> 103
% 5 Acinetobacter_calcoaceticus_PHEA_2.mat: (1319, 1561) -> 328
% 6 Recon1.mat: (2766, 3742) -> 932
% 7 Recon2.mat: (5063, 7440) -> 2430
% 8 Recon3.mat: (8399, 13543)
% 9 cardiac_mit_glcuptake_atpmax.mat: (230, 220) -> 12
% 10 ecoli_core_model.mat: (72, 95) -> 24
% 11 iAF1260.mat: (1668, 2382) -> 524
% 12 iJO1366.mat: (1805, 2583) -> 582
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

%% This file tests the chrr package
% It generates files "result_(modelname)" with a struct with the following 
% properties in the folder "chrr_test".
% 1. dim
% 2. ess
% 3. sampleTime
% 4. steps
maxNumCompThreads(1);
assert(maxNumCompThreads == 1);
curFolder = fileparts(mfilename('fullpath'));
datapath = fullfile(fileparts(curFolder), '/Instances/1chrr', '*.mat');
matfiles = dir(datapath);

result = struct;
attempts = 1000;

idx = {[1:21]};
for c = idx{1}
     try
        fprintf('%d: %s\n', c, matfiles(c).name);
        disp(datetime('now'));
        load(fullfile(fileparts(curFolder), 'Instances/1chrr', matfiles(c).name));
        skip = size(P.A, 2)^2;
        fprintf("Skip: %d\n", skip)
        
        rng(1, 'simdTwister');
        P = Sample_time(attempts, skip, P, [], 0.1);
        %assert(norm(P.originalb - P.originalA * P.samples, Inf) <= 1e-8)

        exps = struct; exps.dim = size(P.A,2);
        exps.ess = min(effective_sample_size(P.samples));
        exps.roundTime = P.roundTime;
        exps.sampleTime = P.sampleTime;
        exps.step = P.attempts * skip;
        fprintf('ESS: %f\n', exps.ess);

        save(fullfile(curFolder, strcat('/chrr_test/result_', matfiles(c).name)), 'exps')
     catch e
        fprintf('Error occurs. Skip model %d\n', c);
        fprintf(1,'The identifier was:\n%s\n', e.identifier);
        fprintf(1,'There was an error! The message was:\n%s\n', e.message);
        continue;
     end
end
