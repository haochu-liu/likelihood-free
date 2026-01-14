load('matlab/data_poisson.mat');
addpath('matlab');

n = [1; 2; 5; 6; 7; 10; 15; 20];

sz = [8 3]; 
varTypes = ["double", "double", "double"];
varNames = ["n", "acc_rate.fix_var", "acc_rate.est_var"];
mcmc_acc1 = table('Size', sz, 'VariableTypes', varTypes, 'VariableNames', varNames);
mcmc_acc1{:, 1} = n;
chain_est_var = zeros(100000, 8);
chain_fix_var = zeros(100000, 8);

for i = 1:8
    [mcmc_chain, acc_rate] = bayes_sl_simple_fix(y, 100000, n(i));
    mcmc_acc1{i, "acc_rate.fix_var"} = acc_rate;
    chain_fix_var(:, i) = mcmc_chain;

    if i > 1
        [mcmc_chain, acc_rate] = bayes_sl_simple(y, 100000, n(i));
        mcmc_acc1{i, "acc_rate.est_var"} = acc_rate;
        chain_est_var(:, i) = mcmc_chain;
    end

    fprintf('Finish n = %d\n.', n(i));
end

save('price2018Bayesian_results.mat', 'mcmc_acc1', 'chain_fix_var', 'chain_est_var');
