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
%          "Bootstrap-Test", "Permutation  test"
%        - "L2-Simul" not supported yet.
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

%}

%% L2-norm based test

gsize = self.n_i;
N = self.N;
k = self.k_groups;
p = self.n_domain_points;

switch method
    case "L2-Simul"
        % TBD
        pvalue = nan;
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
            flag=fix(rand(N,1)*(N-1))+1;
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



