% problemlist = {'LPnetlib@"lp_agg"'};
% %problemlist = {'Recon@1'};
% result = struct;
% attempts = 5000;
% 
% problem = loadProblem(problemlist{1});
% problem.S = problem.Aeq;
% problem.b = problem.beq;
% [samples, stimes] = chrrSampler(problem, 0, attempts, 0);
% 
% nSamples = size(samples, 2);
% stimes = stimes / nSamples
% steps = attempts / nSamples


%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fullify
%problemlist = {'Recon@1'};
%problemlist = {'LPnetlib@"lp_25fv47"'};%, 'LPnetlib@"lp_80bau3b"', 'LPnetlib@"lp_afiro"', 'LPnetlib@"lp_agg"', 'LPnetlib@"lp_beaconfd"', 'LPnetlib@"lp_blend"', 'LPnetlib@"lp_degen2"', 'LPnetlib@"lp_degen3"', 'LPnetlib@"lp_etamacro"', 'LPnetlib@"lp_scorpion"', 'LPnetlib@"lp_sierra"', 'LPnetlib@"lp_truss"'};
%problemlist = {'LPnetlib@"lp_afiro"'};

problemlist = {'metabolic/Recon1'};
%problemlist = {'metabolic/Recon2', 'netlib/truss'};%, 'netlib/80bau3b'};
%problemlist = {'netlib/truss', 'netlib/80bau3b'};
%problemlist = {'netlib/afiro', 'netlib/blend', 'netlib/beaconfd', 'netlib/scorpion', 'netlib/agg'};
%problemlist = {'netlib/degen2', 'netlib/etamacro'};
%problemlist = {'netlib/25fv47'}
%problemlist = {'netlib/degen3','netlib/sierra','netlib/truss','netlib/80bau3b',};
%problemlist = {'netlib/sierra','netlib/truss'};
%problemlist = {'netlib/etamacro'};
%problemlist = {'tvball/tvball'};
result = struct;
attempts = 10000;
skip = 5000;

for i = 1:length(problemlist)
    problem = loadProblem_rhmc(problemlist{i});
    %problem = problemTVballDifficult(205);
    %problem = loadProblem(problemlist{i});
    P = fullify(problem);%, 1e-25);
    disp(size(P.A));
    
    [o, sampleTime] = Sample_time(attempts, skip, P, [], 0.1);
    nSamples = size(o, 2);
    stimes = sampleTime / nSamples;
    steps = attempts*skip / nSamples;
    
    exps = struct; exps.size = size(P.A); 
    exps.stimes = stimes; exps.steps = steps;
    
    tmp = split(problemlist{i}, '/');
    fieldname = tmp{2};
    result.(fieldname) = exps;
    disp([problemlist{i} ' done']);
end

for i = 1:length(problemlist)
    tmp = split(problemlist{i}, '/');
    fieldname = tmp{2};
    disp(fieldname);
    disp(result.(fieldname));
end

% o_shift = o - mean(o, 2); % now, mean is at the origin
% covmtx = o_shift*o_shift'/(nSamples-1);
% 
% e = eig(covmtx);
% e = e(e>0); % ignore zero eigenvalues
% e = sort(e, 'descend');
% 
% subplot(1, 3, 1);
% plot(e, 'o'); title('eigenvalues in descending order')
% subplot(1, 3, 2);
% plot(e(2:end), 'o'); title('eigenvalues without largest one')
% subplot(1, 3, 3);
% plot(log(e), 'o'); title('log(eigenvalues)')

% num_averages = 100;
% dim = 5000;
% times = zeros(1, dim);
% for i = 10:dim+10
%    time = zeros(1, num_averages);
%    for j = 1:num_averages
%        arr = rand(1, dim);
%        tic;
%        arr = arr(6:dim);
%        arr(3) = arr(3) + 3;
%        time(j) = toc;
%    end
%    times(i-9) = mean(time);
% end
% plot(times)
