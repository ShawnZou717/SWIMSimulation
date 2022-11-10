function [random_sigma]  = GeneratingInitPhase(M, N)

sigma_et_r = randn(M + 1, N / 2 + 1);
sigma_et1_r = randn(M + 1, N / 2 + 1);
sigma_et_l = rot90(sigma_et_r, 2);
sigma_et1_l = - rot90(sigma_et1_r, 2);
sigma_et = [sigma_et_l, sigma_et_r(1 : M + 1, 2 : N / 2 + 1)];
sigma_et1 = [sigma_et1_l, sigma_et1_r(1 : M + 1, 2 : N / 2 + 1)];
sigma_et(1 : M / 2, N / 2 + 1) = sigma_et(M + 1 : -1 : M / 2 + 2, N / 2 + 1);
sigma_et1(1 : M / 2, N / 2 + 1)= - sigma_et1(M + 1 : -1 : M / 2 + 2, N / 2 + 1);

random_sigma = (sigma_et + 1i * sigma_et1) / sqrt(2);
random_sigma(M / 2 + 1, N / 2 + 1) = real(random_sigma(M / 2 + 1,N / 2 + 1));
random_sigma(M / 2 + 1, N + 1) = real(random_sigma(M / 2 + 1, N + 1));
random_sigma(M + 1, N / 2 + 1) = real(random_sigma(M + 1, N / 2 + 1));
random_sigma(M + 1, N + 1) = real(random_sigma(M + 1, N + 1));  

random_sigma = rot90(random_sigma, 2);
random_sigma_ = random_sigma(1 : M, 1 : N);
random_sigma_(1, end) = random_sigma(1, end);
random_sigma_(end, 1) = random_sigma(end, 1);
random_sigma_(1, 2 : N / 2) = conj(random_sigma_(1, end : -1 : N / 2 + 2));
random_sigma_(2 : M / 2, 1) = conj(random_sigma_(end : -1 : M / 2 + 2, 1));
random_sigma(2 : end, 2 : end)=random_sigma_;
end