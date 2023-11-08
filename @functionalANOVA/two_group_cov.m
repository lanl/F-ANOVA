function pvalue = two_group_cov(self, method, y1, y2)
% TWO_GROUP_COV Test of equality of covariance functions from 2 groups
% PVALUE = TWO_GROUP_COV(SELF, METHOD, Y1, Y2) 2-Group Equality of Covariance
% Test of equality of covariance functions from 2 groups for functional data. All
% tests are versions of the L2-Norm based tests. Internally called through
% the funtionalANOVA class method: CovarianceTest_TwoSample
%
% <strong>Required Input</strong>
%
% method ([1x1] String)
%        - Options:  "L2-Simul", "L2-Naive", "L2-BiasReduced",
%          "Bootstrap-Test", "Permutation  test"
%        - "L2-Simul" not supported yet.
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
% See also FUNCTIONALANOVA, COVARIANCETEST_TWOSAMPLE

% History
% Initial   March 13, 2008 National University of Singapore, Singapore (Jin-Ting Zhang)
% Modified  May 26,   2023 Los Alamos National Laboratory, USA

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

n1 = self.n_i(1);
n2 = self.n_i(2);
Sigma1=cov(y1);
Sigma2=cov(y2);
N=n1+n2;
Sigma=((n1-1)*Sigma1+(n2-1)*Sigma2)/(N-2);
stat=n1*n2/N*trace((Sigma1-Sigma2)^2);

switch method
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
        y1=y1-ones(n1,1)*mean(y1);
        y2=y2-ones(n2,1)*mean(y2);
        vstat=nan(self.N_boot, 1);

        parfor ii=1:self.N_boot
            flag1=fix(rand(n1,1)*(n1-1))+1;
            flag2=fix(rand(n2,1)*(n2-1))+1;
            yy1=y1(flag1,:);
            yy2=y2(flag2,:);
            S1=cov(yy1);
            S2=cov(yy2);
            stat0=n1*n2/N*trace( ( (S1-S2)-(Sigma1-Sigma2) )^2 );
            vstat(ii)=stat0;
            ts.progress()
        end
        ts.stop;ts.delete
        pvalue=mean(vstat>stat);
%         pstat=[stat,pvalue];

    case "Permutation-Test" %% Permutation  test
        ts = TimedProgressBar(self.N_permutations, 35, 'Running Permutation Test: ', '');
        vstat=nan(self.N_permutations, 1);
        y1=y1-ones(n1,1)*mean(y1);
        y2=y2-ones(n2,1)*mean(y2);
        yy=[y1;y2];
        N=n1+n2;

        parfor ii=1:self.N_permutations
            flag=randperm(N);
            yy1=yy(flag(1:n1),:);
            yy2=yy(flag((n1+1):N),:);
            S1=cov(yy1);
            S2=cov(yy2);
            stat0=n1*n2/N*trace((S1-S2)^2);
            vstat(ii)=stat0;
            ts.progress()
        end
        ts.stop;ts.delete
        pvalue=mean(vstat>stat);
%         pstat=[stat,pvalue];

end

end

