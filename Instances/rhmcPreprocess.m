function o = rhmcPreprocess(problem)
%Input: a structure problem with the following fields
%  .Aineq
%  .bineq
%  .Aeq
%  .beq
%  .lb
%  .ub
%
%Output:
%  o - problem structure

%% Initialize parameters and compiling if needed
opts = default_options();
compile_solver(0); compile_solver(opts.simdLen);

%% Presolve
if isempty(opts.seed)
    opts.seed = randi(2^31);
end

if ischar(opts.logging) || isstring(opts.logging) % logging for Polytope
    fid = fopen(opts.logging, 'a');
    opts.presolve.logFunc = @(tag, msg) fprintf(fid, '%s', msg);
elseif ~isempty(opts.logging)
    opts.presolve.logFunc = opts.logging;
else
    opts.presolve.logFunc = @(tag, msg) 0;
end

o = Polytope(problem, opts);

if ischar(opts.logging) || isstring(opts.logging)
    fclose(fid);
end
end