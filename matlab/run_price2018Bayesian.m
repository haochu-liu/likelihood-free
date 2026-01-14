load('matlab/data_poisson.mat');
addpath('matlab');

n = [1; 2; 5; 6; 7; 10; 15; 20];

sz = [8 5]; 
varTypes = ["double", "double", "double", "double", "double"];
varNames = ["n", "acc_rate.fix_var", "ess.fix_var", "acc_rate.est_var" "ess.est_var"];
mcmc_quality1 = table('Size', sz, 'VariableTypes', varTypes, 'VariableNames', varNames);
mcmc_quality1{:, 1} = n;

for i = 1:8
    [mcmc_chain, acc_rate] = bayes_sl_simple(y, 100000, n(i));
    
end

