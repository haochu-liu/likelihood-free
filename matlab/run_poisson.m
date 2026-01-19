load('matlab/data_poisson.mat');
addpath('matlab');

n = [1, 2, 5, 6, 7, 10, 15, 20];

acc_est_var = zeros(8, 10);
acc_fix_var = zeros(8, 10);
chain_est_var = zeros(100000, 80);
chain_fix_var = zeros(100000, 80);

for i = 1:8
    fprintf('----------\n');
    for j = 1:10
        [mcmc_chain, acc_rate] = bayes_sl_simple_fix(y, 100000, n(i));
        acc_fix_var(i, j) = acc_rate;
        chain_fix_var(:, j+10*(i-1)) = mcmc_chain;

        if i > 1
            [mcmc_chain, acc_rate] = bayes_sl_simple(y, 100000, n(i));
            acc_est_var(i, j) = acc_rate;
            chain_est_var(:, j+10*(i-1)) = mcmc_chain;
        end
        fprintf('Finish chain %d.\n', j);
    end    
    fprintf('Finish n = %d.\n', n(i));
end

save('price2018Bayesian_results.mat', ...
    'n', 'acc_fix_var', 'acc_est_var', 'chain_fix_var', 'chain_est_var');
