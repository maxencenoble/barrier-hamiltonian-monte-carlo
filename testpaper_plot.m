curFolder = fileparts(mfilename('fullpath'));
load(strcat(curFolder, '/Instances/polytopes.mat'))
close all

%plotter = str2func('semilogx');
plotter = str2func('loglog');

%%%%%%%%%%%%%%%%%%%%%
% %% 1. Mixing Time vs Dim
 figure;
 quantVSdim(curFolder, 'RHMC-test/bhmc_test', plotter, polytopes, 'ro', 'originalSize', "mixing",25)
 quantVSdim(curFolder, 'RHMC-test/rhmc_test', plotter, polytopes, 'bo', 'originalSize', "mixing",32)
 legend('BHMC', 'RHMC', 'Location', 'northwest')
 xlim([50 2*1e5]); ylim([10 1e9]);
% 
% %% 2. Sample Time vs Dim
 figure;
 quantVSdim(curFolder, 'RHMC-test/bhmc_test', plotter, polytopes, 'ro', 'originalSize', "sample",25)
 quantVSdim(curFolder, 'RHMC-test/rhmc_test', plotter, polytopes, 'bo', 'originalSize', "sample",32)
 legend('BHMC', 'RHMC', 'Location', 'northwest')
 xlim([50 2*1e5]); ylim([1e-5 1e5]);
%
% %% 3. Mixing Time vs NNZ
 figure;
 quantVSnnz(curFolder, 'RHMC-test/bhmc_test', plotter, 'ro', "mixing")
 quantVSnnz(curFolder, 'RHMC-test/rhmc_test', plotter, 'bo', "mixing")
 legend('BHMC', 'RHMC', 'Location', 'northwest')
% 
% %% 4. Sample Time vs NNZ
 figure;    
 quantVSnnz(curFolder, 'RHMC-test/bhmc_test', plotter, 'ro', "sample")
 quantVSnnz(curFolder, 'RHMC-test/rhmc_test', plotter, 'bo', "sample")
 legend('BHMC', 'RHMC', 'Location', 'northwest')

%% 5. Preparation Time (RHMC: presolve & CHRR: Fullify and MVE)
% figure;
% roundResult1(curFolder, 'RHMC-test/rhmc_test', plotter, polytopes, 'ro', 'originalSize')
% roundResult2(curFolder, 'CHRR-test/chrr_test', plotter, polytopes, 'bo', 'originalSize')
% legend('CRHMC', 'CHRR', 'Location', 'northwest')
% hold off;
% xlim([50 2*1e5]);



%%%%%%%%%%%%%%%%%%%%%

%% Quantities vs Dimension 
 figure;
 plotResult(curFolder, 'RHMC-test/bhmc_test', plotter, polytopes, 'r--o', 'processedSize',25)
 plotResult(curFolder, 'RHMC-test/rhmc_test', plotter, polytopes, 'b--o', 'processedSize',32)
 subplot(1,2,1)
 legend('BHMC', 'RHMC', 'Location', 'northwest')
 subplot(1,2,2)
 legend('BHMC', 'RHMC', 'Location', 'northwest')
 hold off;
% 
 figure;
 plotResult(curFolder, 'RHMC-test/bhmc_test', plotter, polytopes, 'r--o', 'fullSize',25)
 plotResult(curFolder, 'RHMC-test/rhmc_test', plotter, polytopes, 'b--o', 'fullSize',32)
 subplot(1,2,1)
 legend('BHMC', 'RHMC', 'Location', 'northwest')
 subplot(1,2,2)
 legend('BHMC', 'RHMC', 'Location', 'northwest')
 hold off;
%%%%%%%

%% Quantities vs nnz(A)
figure;
plotResult2(curFolder, 'RHMC-test/rhmc_test', plotter, 'ro')
%plotResult2(curFolder, 'CHRR-test/chrr_test', plotter, 'bo')
%plotResult2(curFolder, 'Volesti-test/volesti_test', plotter, 'go')
%subplot(1,2,1)
%legend('RHMC', 'CHRR (Ben)', 'CDHR (Volesti)', 'Location', 'northwest')
%subplot(1,2,2)
%legend('RHMC', 'CHRR (Ben)', 'CDHR (Volesti)', 'Location', 'northwest')

%%%%%%%
%% Prepare Time (Presolve for RHMC & MVE for CHRR and CDHR)
% figure;
% roundResult1(curFolder, 'RHMC-test/rhmc_test', plotter, polytopes, 'ro', 'fullSize')
% roundResult2(curFolder, 'CHRR-test/chrr_test', plotter, polytopes, 'bo', 'fullSize')
% legend('RHMC', 'CHRR', 'Location', 'northwest')
% hold off;

function quantVSnnz(curFolder, pathdir, plotter, color, quant)
    matfiles = dir(fullfile(fullfile(curFolder, pathdir), '*.mat'));
    
    datapath = fullfile(curFolder, '/Instances/1chrr/', '*.mat');
    matfiles2 = dir(datapath);

    numPoly = length(matfiles);
    vnnz = zeros(numPoly, 1);
    time = zeros(numPoly, 1);
    step = zeros(numPoly, 1);
    for idx = 1:numPoly
        result = load(strcat(matfiles(1).folder,'/',matfiles(idx).name));
        if idx < 22
            load(fullfile(curFolder, 'Instances/1chrr/', matfiles2(idx).name));

            if result.exps.ess >= 10
                vnnz(idx) = nnz(P.A_eq);
                time(idx) = result.exps.sampleTime/result.exps.ess;
                step(idx) = result.exps.step/result.exps.ess;
            else
                fprintf("%s: Ess %d is too small\n", matfiles(idx).name, result.exps.ess) 
            end
        else
            if result.exps.ess >= 10
                vnnz(idx) = 295946;
                time(idx) = result.exps.sampleTime/result.exps.ess;
                step(idx) = result.exps.step/result.exps.ess;
            else
                fprintf("%s: Ess %d is too small\n", matfiles(idx).name, result.exps.ess) 
            end
        end
    end

    [vnnz, seq] = sort(vnnz); time = time(seq); step = step(seq);
    filter = vnnz>0;
    vnnz= vnnz(filter); time = time(filter); step = step(filter);
    
    if quant == "sample"
        h = plotter(vnnz, time, color);
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
        title('Sampling Time', 'FontSize', 15);
        xlabel('NNZ', 'FontSize', 15); ylabel('Time/Sample (s)', 'FontSize', 15);
        hold on;
        fit = polyfit(log(vnnz), log(time), 1);
        vnnz = [100, 1000, 10000, 100000, 1e6];
        z = polyval(fit, log(vnnz));
        plotter(vnnz, exp(z), color(1));
        
        algo = pathdir(1:4);
        fprintf(strcat(algo, ": Time/NNZ = %f\n"), fit(1));
    elseif quant == "mixing"
        h = plotter(vnnz, step, color);
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
        title('Mixing Time', 'FontSize', 15);
        xlabel('NNZ', 'FontSize', 15); ylabel('Step/Sample', 'FontSize', 15);
        hold on;
        fit = polyfit(log(vnnz), log(step), 1);
        vnnz = [100, 1000, 10000, 100000, 1e6];
        z = polyval(fit, log(vnnz));
        plotter(vnnz, exp(z), color(1));
        
        algo = pathdir(1:4);
        fprintf(strcat(algo, ": Step/NNZ= %f\n"), fit(1));
    end
end

function quantVSdim(curFolder, pathdir, plotter, polytopes, color, dimOpt, quant, beginning)
    matfiles = dir(fullfile(fullfile(curFolder, pathdir), '*.mat'));

    numPoly = length(matfiles);
    dim = zeros(numPoly, 1);
    time = zeros(numPoly, 1);
    step = zeros(numPoly, 1);
    for idx = 1:numPoly
        result = load(strcat(matfiles(1).folder,'/',matfiles(idx).name));
        name = matfiles(idx).name(beginning:end-4); %26 for bhmc

        if result.exps.ess >= 10
            dim(idx) = polytopes.(name).(dimOpt)(2);
            time(idx) = result.exps.sampleTime/result.exps.ess;
            step(idx) = result.exps.step/result.exps.ess;
            fprintf(strcat(name, ": Sampling Time= %f\n"), time(idx));
        else
            fprintf("%s: Ess %d is too small\n", matfiles(idx).name, result.exps.ess) 
        end
    end
    
    [dim, seq] = sort(dim); time = time(seq); step = step(seq);
    filter = dim>0;
    dim = dim(filter); time = time(filter); step = step(filter);
    
    if quant == "sample"
        h = plotter(dim, time, color);
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
        title('Sampling Time', 'FontSize', 15);
        xlabel('Dimension', 'FontSize', 15); ylabel('Time/Sample (s)', 'FontSize', 15);
        hold on;
        fit = polyfit(log(dim), log(time), 1);
        newdim = [50, 10, 100, 1000, 10000, 100000, 200000];
        z = polyval(fit, log(newdim));
        plotter(newdim, exp(z), color(1));
        algo = pathdir(1:4);
        fprintf(strcat(algo, ": Time/Dim = %f & Coeff = %f\n"), fit(1), exp(fit(2)));
    elseif quant == "mixing"
        h = plotter(dim, step, color);
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
        title('Mixing Time', 'FontSize', 15);
        xlabel('Dimension', 'FontSize', 15); ylabel('Step/Sample', 'FontSize', 15);
        hold on;
        fit = polyfit(log(dim), log(step), 1);
        newdim = [50, 10, 100, 1000, 10000, 100000, 200000];
        z = polyval(fit, log(newdim));
        plotter(newdim, exp(z), color(1));
        algo = pathdir(1:4);
        fprintf(strcat(algo, ": Step/Dim = %f & Coeff = %f\n"), fit(1), exp(fit(2)));
    end
end

function roundResult2(curFolder, pathdir, plotter, polytopes, color, dimOpt,beginning)
    matfiles = dir(fullfile(fullfile(curFolder, pathdir), '*.mat'));
    datapath = fullfile(curFolder, '/Instances/1chrr/', '*.mat');
    matfiles2 = dir(datapath);
    
    numPoly = length(matfiles2);
    dim = zeros(numPoly, 1);
    time = zeros(numPoly, 1);
    
    for idx = 1:numPoly
        load(fullfile(curFolder, 'Instances/1chrr/', matfiles2(idx).name));
        name = matfiles(idx).name(beginning:end-4); 
        dim(idx) = polytopes.(name).(dimOpt)(2);
        time(idx) = P.roundTime;
    end

    [dim, seq] = sort(dim); time = time(seq);
    filter = dim>0;
    dim = dim(filter); time = time(filter);
    h = plotter(dim, time, color);
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    title('Preparation Time', 'FontSize', 15);
    xlabel('Dimension', 'FontSize', 15); ylabel('Time', 'FontSize', 15);
    hold on;
    fit = polyfit(log(dim), log(time), 1);
    newdim = [60, 100, 1000, 10000, 100000, 1e6];
    z = polyval(fit, log(newdim));
    plotter(newdim, exp(z), color(1));
    algo = pathdir(1:4);
    fprintf(strcat(algo, ": RoundTime/Dim = %f & Coeff = %f\n"), fit(1), exp(fit(2)));
end

function roundResult1(curFolder, pathdir, plotter, polytopes, color, dimOpt,beginning)
    matfiles = dir(fullfile(fullfile(curFolder, pathdir), '*.mat'));

    numPoly = length(matfiles);
    dim = zeros(numPoly, 1);
    time = zeros(numPoly, 1);
    for idx = 1:numPoly
        result = load(strcat(matfiles(1).folder,'/',matfiles(idx).name));
        name = matfiles(idx).name(beginning:end-4); 

        if result.exps.ess >= 10
            dim(idx) = polytopes.(name).(dimOpt)(2);
            time(idx) = result.exps.roundTime;
        else
            fprintf("%s: Ess %d is too small\n", matfiles(idx).name, result.exps.ess) 
        end
    end

    [dim, seq] = sort(dim); time = time(seq);
    filter = dim>0;
    dim = dim(filter); time = time(filter);
    h = plotter(dim, time, color);
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    title('Preparation Time', 'FontSize', 15);
    xlabel('Dimension', 'FontSize', 15); ylabel('Time', 'FontSize', 15);
    hold on;
    fit = polyfit(log(dim), log(time), 1);
    newdim = [60, 100, 1000, 10000, 100000, 1e6];
    z = polyval(fit, log(newdim));
    plotter(newdim, exp(z), color(1));
    algo = pathdir(1:4);
    fprintf(strcat(algo, ": RoundTime/Dim = %f & Coeff = %f\n"), fit(1), exp(fit(2)));
end

function plotResult2(curFolder, pathdir, plotter, color)
    matfiles = dir(fullfile(fullfile(curFolder, pathdir), '*.mat'));
    
    datapath = fullfile(curFolder, '/Instances/1chrr/', '*.mat');
    matfiles2 = dir(datapath);

    numPoly = length(matfiles);
    vnnz = zeros(numPoly, 1);
    time = zeros(numPoly, 1);
    step = zeros(numPoly, 1);
    for idx = 1:numPoly
        result = load(strcat(matfiles(1).folder,'/',matfiles(idx).name));
        if idx < 22
            load(fullfile(curFolder, 'Instances/1chrr/', matfiles2(idx).name));

            if result.exps.ess >= 10
                vnnz(idx) = nnz(P.A_eq);
                time(idx) = result.exps.sampleTime/result.exps.ess;
                step(idx) = result.exps.step/result.exps.ess;
            else
                fprintf("%s: Ess %d is too small\n", matfiles(idx).name, result.exps.ess) 
            end
        else
            if result.exps.ess >= 10
                vnnz(idx) = 295946;
                time(idx) = result.exps.sampleTime/result.exps.ess;
                step(idx) = result.exps.step/result.exps.ess;
            else
                fprintf("%s: Ess %d is too small\n", matfiles(idx).name, result.exps.ess) 
            end
        end
    end

    [vnnz, seq] = sort(vnnz); time = time(seq); step = step(seq);
    filter = vnnz>0;
    vnnz= vnnz(filter); time = time(filter); step = step(filter);
    subplot(1,2,1);
    plotter(vnnz, time, color)
    title('Sampling Time', 'FontSize', 15);
    xlabel('NNZ', 'FontSize', 15); ylabel('Time/Sample (s)', 'FontSize', 15);
    hold on;

    subplot(1,2,2)
    plotter(vnnz, step, color)
    title('Mixing Time', 'FontSize', 15);
    xlabel('NNZ', 'FontSize', 15); ylabel('Step/Sample', 'FontSize', 15);
    hold on;

    fit1 = polyfit(log(vnnz), log(time), 1);
    fit2 = polyfit(log(vnnz), log(step), 1);
    
    algo = pathdir(1:4);
    fprintf(strcat(algo, ": Time/NNZ = %f, Step/NNZ= %f\n"), fit1(1), fit2(1));
end

function plotResult_RHMC(curFolder, pathdir, plotter, polytopes, color,beginning)
    matfiles = dir(fullfile(fullfile(curFolder, pathdir), '*.mat'));

    numPoly = length(matfiles);
    dim = zeros(numPoly, 1);
    time = zeros(numPoly, 1);
    step = zeros(numPoly, 1);
    for idx = 1:numPoly
        result = load(strcat(matfiles(1).folder,'/',matfiles(idx).name));
        name = matfiles(idx).name(beginning:end-4); 

        if result.exps.ess >= 10
            dim(idx) = polytopes.(name).processedSize(2);
            time(idx) = result.exps.sampleTime/result.exps.ess;
            step(idx) = result.exps.step/result.exps.ess;
        else
            fprintf("%s: Ess %d is too small\n", matfiles(idx).name, result.exps.ess) 
        end
    end

    [dim, seq] = sort(dim); time = time(seq); step = step(seq);
    subplot(1,2,1);
    plotter(dim, time, color)
    title('Sampling Time', 'FontSize', 15);
    xlabel('Dimension', 'FontSize', 15); ylabel('Time/Sample (s)', 'FontSize', 15);
    hold on;

    subplot(1,2,2)
    plotter(dim, step, color)
    title('Mixing Time', 'FontSize', 15);
    xlabel('Dimension', 'FontSize', 15); ylabel('Step/Sample', 'FontSize', 15);
    hold on;

    fit1 = polyfit(log(dim), log(time), 1);
    fit2 = polyfit(log(dim), log(step), 1);
    
    algo = pathdir(1:4);
    fprintf(strcat(algo, ": Time/Dim = %f, Step/Dim = %f\n"), fit1(1), fit2(1));
end

function plotResult(curFolder, pathdir, plotter, polytopes, color, dimOpt,beginning)
    matfiles = dir(fullfile(fullfile(curFolder, pathdir), '*.mat'));

    numPoly = length(matfiles);
    dim = zeros(numPoly, 1);
    time = zeros(numPoly, 1);
    step = zeros(numPoly, 1);
    for idx = 1:numPoly
        result = load(strcat(matfiles(1).folder,'/',matfiles(idx).name));
        name = matfiles(idx).name(beginning:end-4); 

        if result.exps.ess >= 10
            dim(idx) = polytopes.(name).(dimOpt)(2);
            time(idx) = result.exps.sampleTime/result.exps.ess;
            step(idx) = result.exps.step/result.exps.ess;
        else
            fprintf("%s: Ess %d is too small\n", matfiles(idx).name, result.exps.ess) 
        end
    end

    [dim, seq] = sort(dim); time = time(seq); step = step(seq);
    filter = dim>0;
    dim = dim(filter); time = time(filter); step = step(filter);
    subplot(1,2,1);
    plotter(dim, time, color)
    title('Sampling Time', 'FontSize', 15);
    xlabel('Dimension', 'FontSize', 15); ylabel('Time/Sample (s)', 'FontSize', 15);
    hold on;

    subplot(1,2,2)
    plotter(dim, step, color)
    title('Mixing Time', 'FontSize', 15);
    xlabel('Dimension', 'FontSize', 15); ylabel('Step/Sample', 'FontSize', 15);
    hold on;

    fit1 = polyfit(log(dim), log(time), 1);
    fit2 = polyfit(log(dim), log(step), 1);
    
    algo = pathdir(1:4);
    fprintf(strcat(algo, ": Time/Dim = %f, Step/Dim = %f\n"), fit1(1), fit2(1));
end

