function pvalue = k_group_cov(self, method, stat, Sigma, V)
%K_GROUP_COV k-Group Equality of covariance function
% PVALUE=K_GROUP_COV(SELF, METHOD, STAT, SIGMA, V) Test of equality of
% k-Group covariance functions for functional data. All tests are versions
% of the L2-Norm based tests. Internally called through the funtionalANOVA
% class method: CovarianceTest
%
% <strong>Required Input</strong>
%
% method ([1x1] String)
%        - Options:  "L2-Simul", "L2-Naive", "L2-BiasReduced",
%          "Bootstrap-Test", "Permutation-Test".
%   stat ([1x1] Numeric)
%        - Test statistic to compare against
%  Sigma ([GxG] Numeric)
%        - N: Number of functional samples/realizations
%        - m: number of data points within a sample/realization
%        - if N > m then G=m, else G=N.
%      V ([Nxm] Numeric)
%        - Concatenated zero meaned data matrix from their 
%          respective level.
%
% <strong>Output</strong>
%
% pvalue ([1x1] Numeric)
%        - p value from running the test
%
% See also FUNCTIONALANOVA, COVARIANCETEST

% History
% Initial   April 27, 2010, National University of Singapore (NUS), Singapore (Jin-Ting Zhang)
% Revised   April 18, 2011, National University of Singapore, Singapore
% Revised   Oct 04,   2011, Princeton University, USA
% Revised   April 30, 2012, Princeton University
% Modified  May 26,   2023, Los Alamos National Laboratory, USA
% Modified  July 31,  2024, Los Alamos National Laboratory, USA

%% L2-norm based test

gsize = self.n_i;
N = self.N;
k = self.k_groups;
p = self.n_domain_points;

switch method
    case "L2-Simul"
        q = self.k_groups - 1;

        v_array = cell(1, self.k_groups);

        for K = 1 : self.k_groups
            v_array{K} = self.data{K} - mean(self.data{K}, 2);  %subject-effect matrix for ith group
        end


        n_array = self.n_i;
        LHS = 0;
        ts = TimedProgressBar(sum(n_array), 35, 'Calculating Simulated Null Distribution: ', '');
        for ii = 1:self.k_groups
            n_i = n_array(ii);
            V = v_array{ii};
            for jj = 1:n_i
                v_ij = V(:, jj);
                LHS = LHS + (v_ij * v_ij') * (v_ij * v_ij');
                ts.progress()
            end

        end
        ts.stop;ts.delete

        % Not Any Faster than doing double loop above but uses more memory

        % LHS_TEST = 0;
        % for ii = 1:self.k_groups
        %     V = v_array{ii};
        %     V = reshape(V, [size(V, 1), 1, size(V, 2)]);
        %     V_outer = pagemtimes(V, 'none',  V, 'ctranspose');
        %     V_outer = pagemtimes(V_outer, 'none', V_outer, 'none');
        %     LHS_TEST = LHS_TEST + V_outer;
        % end
        % LHS_TEST = sum(LHS_TEST, 3);

        LHS = LHS ./ N;

        if size(LHS) == size(Sigma)
            omega_hat = LHS - Sigma*Sigma;  %  matrix multiplication
        else % Recalculate Sigma 
            %% Set up Data Matrix
            vmu=[]; V=[];
            for ii=1:self.k_groups
                yyi=self.data{ii}';
                mui=mean(yyi);
                Vi=yyi-ones(self.n_i(ii),1)*mui;
                vmu=[vmu;mui];
                V=[V;Vi];
            end
            SigmaLarge = V'*V/(self.N-self.k_groups);  %% [m x m] pooled covariance matrix
            omega_hat = LHS - SigmaLarge*SigmaLarge;  %  matrix multiplication
        end

        eig_gamma_hat = real(eig(omega_hat));
        
        eig_gamma_hat = eig_gamma_hat(eig_gamma_hat>0);
        T_null = functionalANOVA.Chi_sq_mixture(q, eig_gamma_hat,  self.N_simul);
        T_NullFitted = fitdist(T_null, "Kernel");
        pvalue = 1 - T_NullFitted.cdf(stat);

    case "L2-Naive"
        A=trace(Sigma);
        B=trace(Sigma^2);

        A2=A^2;
        B2=B;

        A3=B2+A2;
        B3=2*trace(Sigma^4)+2*B2^2;

        df=(k-1)*A3^2/B3;
        alpha=B3/A3;

        pvalue=1-chi2cdf(stat/alpha,df);
%         params=[alpha,df,A,B,A2,B2,A3,B3];
    case "L2-BiasReduced"
        A=trace(Sigma);
        B=trace(Sigma^2);

        A2=(N-k)*(N-k+1)/(N-k-1)/(N-k+2)*(A^2-2*B/(N-k+1));
        B2=(N-k)^2/(N-k-1)/(N-k+2)*(B-A^2/(N-k));

        %A=trace(Sigma^2)+trace(Sigma)^2;
        %B=2*trace(Sigma^4)+2*trace(Sigma^2)^2;
        A3=B2+A2;
        B3=2*trace(Sigma^4)+2*B2^2;

        df=(k-1)*A3^2/B3;alpha=B3/A3;

        pvalue=1-chi2cdf(stat/alpha,df);
%         params=[alpha,df,A,B,A2,B2,A3,B3];
        % Bootstrap  test
    case "Bootstrap-Test"
        aflag  = aflagMaker(gsize);
        aflag0=unique(aflag); %% Levels of Factor A
        ts = TimedProgressBar(self.N_boot, 35, 'Running Bootstrap Test: ', 'Finished Bootstrap: ');
        vstat = nan(self.N_boot, 1);
        parfor ii=1:self.N_boot
            % flag=fix(rand(N,1)*(N-1))+1;
            flag = randsample(N, N, true);
            py=V(flag,:);

            %%Computing pS : Boot-strap Pooled Covariance
            R=[];
            for i=1:k
                iflag=(aflag==aflag0(i));
                yyi=py(iflag,:);
                ni=gsize(i);
                mui=mean(yyi);
                Ri=yyi-ones(ni,1)*mui;
                R=[R;Ri];
            end
            if N>p
                pS=R'*R/(N-k);  %% pxp pooled covariance matrix
            else
                pS=R*R'/(N-k);
            end

            stat0=0;nni=0;
            for i=1:k
                ni=gsize(i);
                flag=(nni+1):(nni+ni);
                Ri=R(flag,:);

                % ith groups covariance
                if N>p
                    pSi=Ri'*Ri/(ni-1);  %% Vi: nixp
                    temp=trace((pSi-pS)^2);
                else
                    pSi=Ri*Ri'/(ni-1);
                    temp=trace(pSi^2)-2*trace(Ri*R'*R*Ri')/(N-k)/(ni-1)+trace(pS^2);
                end
                stat0=stat0+(ni-1)*temp;
                nni=nni+ni;
            end

            vstat(ii) = stat0;
            ts.progress()

        end
        ts.stop;ts.delete
        pvalue= 1 - sum(vstat < stat) / double(self.N_boot);
       
    case "Permutation-Test" % Permutation  test
        aflag  = aflagMaker(gsize);
        aflag0=unique(aflag); %% Levels of Factor A
        ts = TimedProgressBar(self.N_permutations, 35, 'Running Permutation Test: ', 'Finished Permutation: ');
        vstat = nan(self.N_permutations, 1);
        parfor ii=1:self.N_permutations

            flag=randperm(N);
            py=V(flag,:);

            %%Computing pS : Boot-strap Pooled Covariance
            R=[];
            for i=1:k
                iflag=(aflag==aflag0(i));
                yyi=py(iflag,:);
                ni=gsize(i);
                mui=mean(yyi);
                Ri=yyi-ones(ni,1)*mui;
                R=[R;Ri];
            end

            if N>p
                pS=R'*R/(N-k);  %% pxp pooled covariance matrix
            else
                pS=R*R'/(N-k);
            end


            stat0=0;nni=0;
            for i=1:k
                ni=gsize(i);
                flag=(nni+1):(nni+ni);
                Ri=R(flag,:);

                % ith groups covariance
                if N>p
                    pSi=Ri'*Ri/(ni-1);  %% Vi: nixp
                    temp=trace((pSi-pS)^2);
                else
                    pSi=Ri*Ri'/(ni-1);
                    temp=trace(pSi^2)-2*trace(Ri*R'*R*Ri')/(N-k)/(ni-1)+trace(pS^2);
                end
                stat0=stat0+(ni-1)*temp;

                nni=nni+ni;
            end

            vstat(ii) = stat0;
            ts.progress()
        end
        ts.stop;ts.delete

        pvalue= 1 - sum(vstat < stat) / self.N_permutations;
end

end



function aflag  = aflagMaker(n_i)

aflag=[];

for K = 1 :numel(n_i)
    idicator = repmat(K, n_i(K), 1);
    aflag = [aflag; idicator];
end



end



