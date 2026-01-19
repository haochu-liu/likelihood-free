load('matlab/data_ricker.mat');
addpath('matlab');

std_rw = [0.15 0.5 0.15];
corr_rw = [1 -0.7 -0.6; -0.7 1 0.4; -0.6 0.4 1];
cov_rw = corr2cov(std_rw,corr_rw);
n = [30, 40, 50, 60, 70];

acc_rate_m = zeros(5, 10);
chain_m = zeros(5000, 3, 80);

for i = 1:5
    fprintf('----------\n');
    for j = 1:10
        [mcmc_chain, acc_rate] = bayes_sl_ricker_wood(y,1,5000,n(i),cov_rw);
        acc_rate_m(i, j) = acc_rate;
        chain_m(:, :, j+10*(i-1)) = mcmc_chain;

        fprintf('Finish chain %d.\n', j);
    end    
    fprintf('Finish n = %d.\n', n(i));
end

save('ricker_results.mat', ...
    'n', 'acc_rate_m', 'chain_m');
