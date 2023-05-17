function P = loadProblem(name)
path_size = split(name,'@');
path_folders = split(path_size{1},'/');
curFolder = fileparts(mfilename('fullpath'));
path = fullfile(curFolder, path_folders{:});

% check if the file exists as a mat
if exist([path '.mat'], 'file') && length(path_size) == 1
    model = load([path '.mat']);
    model = cellfun(@(x)(model.(x)),fieldnames(model));
    P = struct;        
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
    %P = problem;
%elseif exist([path '.m'], 'file') && length(path_size) == 2
%    [folder,name] = fileparts(path);
%    prevDir = pwd;
%    cd(folder);
%    h = str2func(['@' name]);
%    cd(prevDir);
%    
%    scurr = rng;
%    rng(123456); % fix the seed for random generator
%    P = h(str2double(path_size{2}));
%    rng(scurr);
else
    error(['Problem ' name ' does not exists']);
end

P = standardize_problem(P);

% netlib LP are not bounded
if contains(name, 'netlib/')
    df = P.df; % the df for netlib problem is a fixed vector
    x = linprog(df, P.Aineq, P.bineq, P.Aeq, P.beq, P.lb, P.ub, struct('Display','none'));
    threshold = df' * x + abs(df)' * abs(x);
    P.f = [];
    P.df = [];
    P.ddf = [];
    P.Aineq = [P.Aineq; df'];
    P.bineq = [P.bineq; threshold];
    P.ub = min(P.ub, 2 * max(abs(x)));
    P.lb = max(P.lb, -2 * max(abs(x)));
end
