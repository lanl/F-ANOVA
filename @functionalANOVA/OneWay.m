function pvalue_matrix = OneWay(self, n_tests, q, eig_gamma_hat, eta_i, eta_grand, params, pair_vec)
%ONEWAY performs homogeneous One-way ANOVA for functional data
% PVALUE_MATRIX = ONEWAY(SELF, N_TESTS, Q, EIG_GAMMA_HAT, ETA_I, ETA_GRAND, PARAMS, PAIR_VEC)
% Internally called through the funtionalANOVA class method: OneWayANOVA.
%
% <strong>Required Input</strong>
%
%       n_tests ([1x1] Numeric)
%               - Number of tests to run. 
%               - For Hypothesis="FAMILY",   n_tests=1 
%               - For Hypothesis="PAIRWISE", n_tests=nchoosek(self.k_groups, 2)
%             q ([1x1] Numeric)
%               - Rank of contrast matrix. q=self.k_group - 1
% eig_gamma_hat ([MxM] Numeric)
%               - Positive eigen values of the compute estimated covariance
%               - M is the number of data points within each function sample
%         eta_i ([MxN] Numeric)
%               - Mean function response for each level within the main group
%     eta_grand ([1x1] Numeric)
%               - Constant used for Hypothesis testing. Typically, c=0 to denote
%        params ([1x8] Cell)
%               - Parameters used in method, see commented code.
%      pair_vec ([Px1] String)
%               - P = n_tests
%               - Strings used to denote the Hypothesis.
%               - For Hypothesis="Family", pair_vec="FAMILY"
%               - For Hypothesis="PAIRWISE", pair_vec are all the pairwise
%                 combinations of the levels within the main group.
%
% <strong>Output</strong>
%
% pvalue_matrix ([PxQ] Numeric)
%               - P = n_tests
%               - Q = Number of methods used in the Analysis
%               - Matrix of p-values from running all the tests and methods
%
% See also FUNCTIONALANOVA

%  History
%  Initial May 26, 2023 Los Alamos National Laboratory, USA

%{
© 2023. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos
National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S.
Department of Energy/National Nuclear Security Administration. All rights in the program are.
reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear
Security Administration. The Government is granted for itself and others acting on its behalf a
nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare.
derivative works, distribute copies to the public, perform publicly and display publicly, and to permit.
others to do so.

Author: Adam Watts (acwatts@lanl.gov)
%}

n_methods = numel(self.ANOVA_Methods_Used);

pvalue_matrix = zeros(n_tests, n_methods);
self.CriticalValues = cell(n_methods, 3);

T_n  = params{1}; % L2 Test Statistic
F_n = params{2};  % F-teype test Statistic
beta_hat = params{3}; % Chi-square mixture Approximation parameters 
kappa_hat = params{4}; % Chi-square mixture Approximation parameters 
beta_hat_unbias = params{5}; % Chi-square mixture Approximation parameters 
kappa_hat_unbias  = params{6} ; % Chi-square mixture Approximation parameters 
C = params{7} ; % Contrast Matrix
D = params{8} ; % Accounts for Sample Size differences

counter = 0;
for method = self.ANOVA_Methods_Used   % Generate Test Statistic
    counter = counter + 1;

    switch method

        case "L2-Simul"
            T_null = functionalANOVA.Chi_sq_mixture(q, eig_gamma_hat, self.N_simul);

            T_NullFitted = fitdist(T_null, "Kernel");

            p_value = zeros(n_tests, 1);
            for cc = 1:n_tests
                %  p_value(cc) = 1 - sum(T_n_Boot(:,cc) < T_n(cc) ) ./ double(N_boot);
                p_value(cc) = 1 - T_NullFitted.cdf(T_n(cc));

                if self.showSimulPlot
                    functionalANOVA.PlotTestStats(p_value(cc), self.alpha, T_null, T_n(cc), method + " test", self.Hypothesis, pair_vec(cc))
                end

            end
            pvalue_matrix(:, counter) = p_value;

            % Critical Values for Plotting
            self.CriticalValues{counter, 1} = method;
            self.CriticalValues{counter, 2} = quantile(T_null, 1 - self.alpha);

            if strcmp(self.Hypothesis, "FAMILY")
                self.updateFamilyTable( method, {T_NullFitted});
            end


        case "L2-Naive"
            p_value = zeros(n_tests, 1);
            for cc = 1:n_tests
                p_value(cc) =  1 - ncx2cdf(T_n(cc) / beta_hat, q * kappa_hat, 0 ); %eq 5.47
            end
            pvalue_matrix(:, counter) = p_value;

            % Critical Values for Plotting
            self.CriticalValues{counter, 1} = method;
            self.CriticalValues{counter, 2} = beta_hat*ncx2inv(1 - self.alpha, q * kappa_hat, 0); % Get 1-alpha value for estiamted chi-square

        case "L2-BiasReduced"
            p_value = zeros(n_tests, 1);
            for cc = 1:n_tests
                p_value(cc) =  1 - ncx2cdf(T_n(cc) / beta_hat_unbias, q * kappa_hat_unbias, 0 );
            end
            pvalue_matrix(:, counter) = p_value;


            % Critical Values for Plotting
            self.CriticalValues{counter, 1} = method;
            self.CriticalValues{counter, 2} = beta_hat_unbias*ncx2inv(1 - self.alpha, q * kappa_hat_unbias, 0); % Get 1-alpha value for estiamted chi-square

        case "L2-Bootstrap"
            T_n_Boot = zeros(self.N_boot, n_tests);
            switch self.Hypothesis
                case "FAMILY"
                    self.hypothesis_LABEL =  pair_vec(1);
                    T = self.setUpTimeBar(method);

                    d_points = self.n_domain_points;
                    k_group =  self.k_groups;
                    n_iii = self.n_i;
                    g_data = self.data;
                    n = self.N;
                    parfor K = 1 : self.N_boot
                        [eta_i_star, eta_grand_star, ~] = functionalANOVA.GroupBooter(g_data, d_points, k_group, n_iii, n);
                        T_n_Boot(K) = functionalANOVA.SSH_n_boot(eta_i, eta_grand, eta_i_star, eta_grand_star, n_iii);
                        T.progress;
                    end
                    T.stop; T.delete;
                    self.updateFamilyTable(method, {self.N_boot});

                case "PAIRWISE"
                    for cc = 1:n_tests
                        self.hypothesis_LABEL =  pair_vec(cc);
                        T = self.setUpTimeBar(method);

                        Ct = C(cc,:);

                        d_points = self.n_domain_points;
                        k_group =  self.k_groups;
                        n_iii = self.n_i;
                        g_data = self.data;
                        n = self.N;
                        parfor K = 1 : self.N_boot
                            [eta_i_star, ~, ~] = functionalANOVA.GroupBooter(g_data, d_points, k_group, n_iii, n);
                            SSH_t = ((Ct*(eta_i_star - eta_i)').^2) .* inv(Ct*D*Ct') ;
                            T_n_Boot(K, cc) = sum(SSH_t);
                            T.progress;
                        end
                        T.stop; T.delete;
                    end
            end

            p_value = zeros(n_tests, 1);
            critVals = zeros(n_tests, 1);

            for cc = 1:n_tests
                p_value(cc) = 1 - sum(T_n_Boot(:,cc) < T_n(cc) ) ./ double(self.N_boot);
                critVals(cc) = quantile(T_n_Boot(:,cc), 1-self.alpha);
            end
            pvalue_matrix(:, counter) = p_value;

            % *********** NEEDS HELP IN PAIRWISE ************
            % Critical Values for Plotting
            self.CriticalValues{counter, 1} = method;
            self.CriticalValues{counter, 2} = critVals(cc);

        case "F-Simul"
            ratio = (self.N - self.k_groups)/ q ;
            T_null = functionalANOVA.Chi_sq_mixture(q, eig_gamma_hat, self.N_simul);
            F_null_denom = functionalANOVA.Chi_sq_mixture(self.N - self.k_groups, eig_gamma_hat, self.N_simul);
            F_null = (T_null ./ F_null_denom) .* ratio;

            F_NullFitted = fitdist(F_null, "Kernel");

            p_value = zeros(n_tests, 1);
            for cc = 1:n_tests
                p_value(cc) = 1 - F_NullFitted.cdf(F_n(cc));
                if self.showSimulPlot
                    functionalANOVA.PlotTestStats(p_value(cc), self.alpha, F_null, F_n(cc), method + " test", self.Hypothesis, pair_vec(cc))
                end
            end
            pvalue_matrix(:, counter) = p_value;

            if strcmp(self.Hypothesis, "FAMILY")
                self.updateFamilyTable(method, {F_NullFitted});
            end

            % Critical Values for Plotting
            self.CriticalValues{counter, 1} = method;
            self.CriticalValues{counter, 3} = quantile(F_null, 1 - self.alpha); % Get 1-alpha value for estiamted chi-square
            self.CriticalValues{counter, 2} = quantile((self.CriticalValues{counter, 3} / ratio) .* F_null_denom, 1 - self.alpha); % Reverse Calculate 0.95 Quantile of T_N

        case "F-Naive"

            p_value = zeros(n_tests, 1);
            for cc = 1:n_tests
                p_value(cc) = 1 - fcdf(F_n(cc), q * kappa_hat, (self.N - self.k_groups) * kappa_hat);
            end
            pvalue_matrix(:, counter) = p_value;

            A = q * kappa_hat;
            B = (self.N - self.k_groups) * kappa_hat;

            self.CriticalValues{counter, 1} = method;
            self.CriticalValues{counter, 3} = finv(1 - self.alpha,A , B); % Get 1-alpha value for estiamted f-dist

            % Unable to back-calculate what the T_n Value is.
            % Since both numerator and denom are reduced

        case "F-BiasReduced"
            p_value = zeros(n_tests, 1);
            for cc = 1:n_tests
                p_value(cc) = 1 - fcdf(F_n(cc), q * kappa_hat_unbias, (self.N - self.k_groups) * kappa_hat_unbias);
            end
            pvalue_matrix(:, counter) = p_value;


            self.CriticalValues{counter, 1} = method;
            self.CriticalValues{counter, 3} = finv(1 - self.alpha, q * kappa_hat_unbias, (self.N - self.k_groups) * kappa_hat_unbias); % Get 1-alpha value for estiamted f-dist


        case "F-Bootstrap"
            F_n_Boot = zeros(self.N_boot, n_tests);
            ratio = (self.N - self.k_groups)/ q ;
            critVals = zeros(n_tests, 1);
            ReversedT_n = zeros(n_tests, 1);

            switch self.Hypothesis
                case "FAMILY"
                    self.hypothesis_LABEL =  pair_vec(1);
                    T = self.setUpTimeBar(method);
                    d_points = self.n_domain_points;
                    k_group =  self.k_groups;
                    n_iii = self.n_i;
                    g_data = self.data;
                    n = self.N;
                    f_n_Denominator_Boot = zeros(self.N_boot, n_tests);
                    parfor K = 1 : self.N_boot
                        [eta_i_star, eta_grand_star, gamma_hat_star] = functionalANOVA.GroupBooter(g_data, d_points, k_group, n_iii, n);
                        f_n_Denominator_Boot(K, 1) = trace(gamma_hat_star) * (n-k_group);
                        T_n_Boot = functionalANOVA.SSH_n_boot(eta_i, eta_grand, eta_i_star, eta_grand_star, n_iii);
                        F_n_Boot(K, :) = (T_n_Boot ./ f_n_Denominator_Boot(K, 1)) .* ratio;
                        T.progress;
                    end
                    T.stop; T.delete;
                    self.updateFamilyTable(method, {self.N_boot});

                case "PAIRWISE"
                    f_n_Denominator_Boot = zeros(self.N_boot, n_tests);

                    for cc = 1:n_tests
                        self.hypothesis_LABEL =  pair_vec(cc);
                        T = self.setUpTimeBar(method);
                        Ct = C(cc,:);

                        d_points = self.n_domain_points;
                        k_group =  self.k_groups;
                        n_iii = self.n_i;
                        g_data = self.data;
                        n = self.N;
                        parfor K = 1 : self.N_boot
                            [eta_i_star, ~, gamma_hat_star] = functionalANOVA.GroupBooter(g_data, d_points, k_group, n_iii, n);
                            f_n_Denominator_Boot(K, cc) = trace(gamma_hat_star) * (n-k_group);
                            SSH_t = ((Ct*(eta_i_star - eta_i)').^2) .* inv(Ct*D*Ct') ;
                            T_n_Boot = sum(SSH_t);
                            F_n_Boot(K, :) = (T_n_Boot ./ f_n_Denominator_Boot(K, cc)) .* ratio;
                            T.progress;
                        end
                        T.stop; T.delete;
                    end
            end

            p_value = zeros(n_tests, 1);
            for cc = 1:n_tests
                p_value(cc) = 1 - sum(F_n_Boot(:,cc) < F_n(cc) ) ./ double(self.N_boot);

                critVals(cc) = quantile(F_n_Boot(:,cc), 1-self.alpha);
                ReversedT_n(cc) = quantile((critVals(cc) / ratio) .* f_n_Denominator_Boot(:, cc), 1-self.alpha);
                % Get Crit Value

            end
            pvalue_matrix(:, counter) = p_value;

            % Critical Values for Plotting
            self.CriticalValues{counter, 1} = method;
            self.CriticalValues{counter, 2} = ReversedT_n(cc);
            self.CriticalValues{counter, 3} = critVals(cc);

    end


end
end

