function pvalue = k_group_cov_pairwise(self, method, y1, y2)
%K_GROUP_COV_PAIRWISE Test of equality of covariance functions from 2 groups
% PVALUE = TWO_GROUP_COV(SELF, METHOD, Y1, Y2) 2-Group Equality of Covariance
% Test of equality of covariance functions from 2 groups for functional data. All
% tests are versions of the L2-Norm based tests. Slightly different than
% two_group_cov for Bootstrap and Permutation methods only as it uses the
% sample space vs feature space.
%
% This method works with a pooled covariance estimate and looks at
% deviations of each group’s covariance from this pooled estimate. The
% two_group_cov compares the difference between the two bootstrapped
% covariance matrices directly (adjusted by the original difference).
%
% Internally called through the funtionalANOVA class method: k_group_cov
%
% <strong>Required Input</strong>
%
% method ([1x1] String)
%        - Options:  "L2-Simul", "L2-Naive", "L2-BiasReduced",
%          "Bootstrap-Test", "Permutation-Test"
% y1 ([Wxm] Numeric)
%         - Sample 1 in short Format
%         - W: Number of functional samples/realizations from sample 1
%         - m: number of data points within a sample/realization
%         - Rows are realizations, columns are the data for a realization
% y2 ([Zxm] Numeric)
%         - Sample 2 in short Format
%         - Z: Number of functional samples/realizations from sample 2
%         - m: number of data points within a sample/realization
%         - Rows are realizations, columns are the data for a realization
%
% <strong>Output</strong>
%
% pvalue ([1x1] Numeric)
%        - p value from running the test
%
% See also FUNCTIONALANOVA, K_GROUP_COV

% History
% Initial   March 13, 2008 National University of Singapore, Singapore (Jin-Ting Zhang)
% Modified  May 26,   2023 Los Alamos National Laboratory, USA
% Modified July 30,   2024 Los Alamos National Laboratory, USA

n1 = size(y1,1);
n2 = size(y2,1);
L = size(y1,2);
gsize = [n1, n2];

Sigma1=cov(y1);
Sigma2=cov(y2);
N=n1+n2;
% Sigma=((n1-1)*Sigma1+(n2-1)*Sigma2)/(N-2);
stat=(n1-1)*(n2-1) / (N-2) * trace((Sigma1-Sigma2)^2);  % Corrected


% TESTING
vmu=[]; V=[];
n_i = [n1, n2];
temp_data = {y1, y2};
for ii=1:2
    yyi=temp_data{ii};
    mui=mean(yyi);
    Vi=yyi-ones(n_i(ii),1)*mui;
    vmu=[vmu;mui];
    V=[V;Vi];
end

m = self.n_domain_points;

if N > m
    Sigma=V'*V/(N-2);  %% [m x m] pooled covariance matrix
else
    Sigma=V*V'/(N-2);  %% [N x N] pooled covariance matrix
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Computing the test statistic  (Same as above)
% stat=0;  % T_n test statistic
% nni=0;
% for ii=1:2
%     ni=n_i(ii);
%     flag=(nni+1):(nni+ni);
%     Vi=V(flag,:);
%     if N > m
%         Si=Vi'*Vi/(ni-1);  %% Vi: nixp
%         temp=trace((Si-Sigma)^2);
%     else
%         Si=Vi*Vi'/(ni-1);
%         temp=trace(Si^2)-2*trace(Vi*V'*V*Vi')/(N-2)/(ni-1)+trace(Sigma^2);
%     end
%
%     stat=stat+(ni-1)*temp;
%     nni=nni+ni;
% end
%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch method
    case "L2-Simul"
        q=1;

        v_1j = y1' - mean(y1, 1)'; %subject-effect matrix for 1st group
        v_2j = y2' - mean(y2, 1)'; %subject-effect matrix for 2nd group

        n_array = [n1, n2];
        v_array = {v_1j, v_2j};
        LHS = 0;
        ts = TimedProgressBar(sum(n_array), 35, 'Calculating Simulated Null Distribution: ', '');
        for ii = 1:2
            n_i = n_array(ii);
            V = v_array{ii};
            for jj = 1:n_i
                v_ij = V(:, jj);
                LHS = LHS + (v_ij * v_ij') * (v_ij * v_ij');
                ts.progress()
            end

        end

        ts.stop;ts.delete

        LHS = LHS ./ N;

        if size(LHS) == size(Sigma)
            omega_hat = LHS - Sigma*Sigma;  %  matrix multiplication
        else
            SigmaLarge=((n1-1)*cov(y1)+(n2-1)*cov(y2))/(N-2);
            omega_hat = LHS - SigmaLarge*SigmaLarge;  %  matrix multiplicatio
        end
        eig_gamma_hat = real(eig(omega_hat));

        eig_gamma_hat = eig_gamma_hat(eig_gamma_hat>0);
        T_null = functionalANOVA.Chi_sq_mixture(q, eig_gamma_hat,  self.N_simul);
        T_NullFitted = fitdist(T_null, "Kernel");
        pvalue = 1 - T_NullFitted.cdf(stat);
    case "L2-BiasReduced"  % Bias Reduced
        A=trace(Sigma^2)+trace(Sigma)^2;
        B=2*trace(Sigma^4)+2*trace(Sigma^2)^2;
        alpha=(N-2)^2/N/(N-3)*(B-A^2/(N-2))/A;
        df=(1+1/(N-2))*(A^2-2*B/(N-1))/(B-A^2/(N-2));

        pvalue=1-chi2cdf(stat/alpha,df);
        %         pstat=[stat,pvalue,alpha,df];

    case "L2-Naive" % Naive Method

        an=trace(Sigma);bn=trace(Sigma^2);
        cn=trace(Sigma^3);dn=trace(Sigma^4);
        An=an;
        Bn=(N-2)^2/(N*(N-3))*(bn-An^2/(N-2));
        Cn=(cn-3/(N-2)*Bn*An)/(1+3/(N-2));
        Dn=(dn-6/(N-2)*Cn*An)/(1+6/(N-2));
        A=Bn+An^2;
        B=2*Dn+2*Bn^2;
        alpha=B/A;df=A^2/B;

        pvalue=1-chi2cdf(stat/alpha,df);
        %         pstat=[stat,pvalue,alpha,df];

    case "Bootstrap-Test" %% Bootstrap test
        ts = TimedProgressBar(self.N_boot, 35, 'Running Bootstrap Test: ', '');
        vstat=nan(self.N_boot, 1);
        k=2;

        parfor ii=1:self.N_boot
            flag1=randsample(n1, n1, true);
            flag2=randsample(n2, n2, true);

            yy1=y1(flag1,:);
            yy2=y2(flag2,:);

            %Computing pS : Boot-strap Pooled Covariance
            mu1=mean(yy1);
            R1=yy1-ones(n1,1)*mu1;

            mu2=mean(yy2);
            R2=yy2-ones(n2,1)*mu2;

            R=[R1;R2];

            if N > m
                pS=R'*R/(N-k);  %% pxp pooled covariance matrix
            else
                pS=R*R'/(N-k);
            end

            stat0=0;
            nni=0;
            for i=1:k
                ni=gsize(i);
                flag=(nni+1):(nni+ni);
                Ri=R(flag,:);

                % ith groups covariance
                if N > m
                    pSi=Ri'*Ri/(ni-1);  %% Vi: nixp
                    temp=trace((pSi-pS)^2);
                else
                    pSi=Ri*Ri'/(ni-1);
                    temp=trace(pSi^2)-2*trace(Ri*R'*R*Ri')/(N-k)/(ni-1)+trace(pS^2);
                end
                stat0=stat0+(ni-1)*temp;
                vstat_check(ii) = stat0;
                nni=nni+ni;
            end
            ts.progress()
        end
        ts.stop;ts.delete
        pvalue=mean(vstat>stat);

    case "Permutation-Test" %% Permutation  test
        ts = TimedProgressBar(self.N_permutations, 35, 'Running Permutation Test: ', '');
        vstat=nan(self.N_boot, 1);
        k=2;
        Y = [y1; y2];
        parfor ii=1:self.N_boot
            flag=randperm(N);
            Y_perm = Y(flag, :)
            yy1=Y(1:n1,:);
            yy2=Y(n1+1:end,:);

            %Computing pS : Boot-strap Pooled Covariance
            mu1=mean(yy1);
            R1=yy1-ones(n1,1)*mu1;

            mu2=mean(yy2);
            R2=yy2-ones(n2,1)*mu2;

            R=[R1;R2];

            if N > m
                pS=R'*R/(N-k);  %% pxp pooled covariance matrix
            else
                pS=R*R'/(N-k);
            end

            stat0=0;
            nni=0;
            for i=1:k
                ni=gsize(i);
                flag=(nni+1):(nni+ni);
                Ri=R(flag,:);

                % ith groups covariance
                if N > m
                    pSi=Ri'*Ri/(ni-1);  %% Vi: nixp
                    temp=trace((pSi-pS)^2);
                else
                    pSi=Ri*Ri'/(ni-1);
                    temp=trace(pSi^2)-2*trace(Ri*R'*R*Ri')/(N-k)/(ni-1)+trace(pS^2);
                end
                stat0=stat0+(ni-1)*temp;
                vstat_check(ii) = stat0;
                nni=nni+ni;
            end
            ts.progress()
        end
        ts.stop;ts.delete
        pvalue=mean(vstat>stat);

end

end




