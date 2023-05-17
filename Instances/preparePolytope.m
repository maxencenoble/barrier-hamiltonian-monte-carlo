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

addpath(genpath('../PolytopeSamplerMatlab-develop'));
addpath(genpath('../CHRR-test'));
curFolder = fileparts(mfilename('fullpath'));
matfiles = dir(fullfile(curFolder, '/0raw/*.mat'));
numPoly = size(matfiles, 1);


polytopes = struct;

%  1~12: biological
% 13~22: Netlib 

for idx = 21:21
     try 
        %% 0) Load instances
        name = matfiles(idx).name(1:end-4);
        fprintf('** Instance %d: %s\n', idx, name);
        polytopes.(name) = struct;
        disp(datetime('now'));

        model = load(fullfile(curFolder, '/0raw', matfiles(idx).name));
        model = cellfun(@(x)(model.(x)), fieldnames(model));
        P = struct;        
        if matfiles(idx).name(5:6) == "lp"
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
        fprintf('Original Size: (%d, %d)\n', size(P.Aeq, 1), size(P.Aeq, 2));
        polytopes.(name).originalSize = size(P.Aeq);

        %% 1) First preproces by RHMC's code
        % Remove redundant rows and scale instances for numerical stability
        %o = sample(P, 1);
        P = rhmcPreprocess(P);
        fprintf('Size after preprocessing: (%d, %d)\n', size(P.A, 1), size(P.A, 2));
        polytopes.(name).processedSize = size(P.A);

        %% 2) Change into full-dimensional form {Ax \leq b} & Run MVE
        % If n_fulldim >= 10000, don't run CHRR and CDHR due to memory limit
        if size(P.A, 2)-size(P.A, 1) >= 10000
            continue
        end
        
        % Make it full-dimensional + Round by MVE
        mve = 1;
        P = fullify(P, mve);
        if matfiles(idx).name(5:6) == "lp"
            P.originalA = model.A;
            P.originalb = model.b;
        else
            P.originalA = model.S; 
            P.originalb = model.b; 
        end
        
        fprintf('Size after fullify: (%d, %d)\n', size(P.A, 1), size(P.A, 2));
        polytopes.(name).fullSize = size(P.A);

        % Save these rounded polytopes which will be fed into CHRR (Matlab)
        if size(P.A, 2) >= 6000
            save(fullfile(curFolder, strcat('/1chrr/chrr_', matfiles(idx).name)), 'P', '-v7.3')
        else
            save(fullfile(curFolder, strcat('/1chrr/chrr_', matfiles(idx).name)), 'P')
        end

        %% 3) Change it into HPolytope for CDHR (C++)
        poly=cell(2,1);
        poly{1} = P.A; poly{2} = P.b;
        if size(P.A, 2) >= 6000
            save(fullfile(curFolder, strcat('/2cdhr/cdhr_', matfiles(idx).name)), 'poly', '-v7.3')
        else
            save(fullfile(curFolder, strcat('/2cdhr/cdhr_', matfiles(idx).name)), 'poly')
        end
        
        %% 4) Save dimension info
        save(strcat(curFolder, '/polytopes.mat'), 'polytopes')
    catch e
        fprintf('Error occurs. Skip model %d\n', idx);
        fprintf(1,'The identifier was:\n%s\n',e.identifier);
        fprintf(1,'There was an error! The message was:\n%s\n',e.message);
        continue;
    end
end