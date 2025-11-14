function [pvalue, stat] = OneWay_BF(self, method, data, Contrast, c, varargin)
%ONEWAY_BF Heteroscedastic OneWay F-ANOVA
% [PVALUE, STAT]=ONEWAY_BF(SELF, METHOD, DATA, CONTRAST, C) 
% Heteroscedastic One-way ANOVA for functional data. Internally called
% through the funtionalANOVA class method: OneWayANOVA_BF.
%
% <strong>Required Input</strong>
%
%      method ([1x1] String)
%             - Options:  "L2-Simul", "L2-Naive", "L2-BiasReduced", "L2-Bootstrap",
%               "F-Simul", "F-Naive", "F-BiasReduced", "F-Bootstrap"
%        data ([Nxm] Numeric)
%             - N: Number of functional samples/realizations
%             - m: number of data points within a sample/realization
%    Contrast ([1xG] Numeric)
%             - Contrast vector for specific hypotheses.
%             - Let p and q be the number of levels within the primary
%               and secondary factors respectively.
%             - Contrast for Primary levels,   then G=p.
%             - Contrast for Secondary levels, then G=q.
%           c ([1x1] Numeric)
%             - Constant used for Hypothesis testing. Typically, c=0 to denote
%             that one expects a zero mean difference between factors and
%             levels within factors.
%
% <strong>Optional Inputs</strong>
%
% Indicator_A ([Nx1] Numeric)
%             - Custom Indicator array, generated in TwoWay_BF.m
%
% <strong>Output</strong>
%
%      pvalue ([1x1] Numeric)
%             - p value from running the test
%        stat ([1x1] Numeric)
%             - test statistic value from either the L2 or F-type
%               methods
%
% See also FUNCTIONALANOVA

%  History
% Initial  April 27, 2010 National University of Singapore (NUS), Singapore (Jin-Ting Zhang)
% Revised  Sept 29,  2010 
% Revised  May 04,   2012 Princeton University
% Revised  Jan 10,   2013 NUS, Singapore
% Modified May 26,   2023 Los Alamos National Laboratory, USA

ip = inputParser;

addRequired(ip, 'method', @(x) isstring(x)) 
addRequired(ip, 'data', @(x) isnumeric(x)) 
addRequired(ip, 'Contrast', @(x) isnumeric(x)) 
addRequired(ip, 'c', @(x) isnumeric(x)) 

addParameter(ip, 'Indicator_A', [], @(x)  isnumeric(x))  % 

parse(ip, method, data, Contrast, c, varargin{:});


N = self.N;
p = self.n_domain_points;

if isempty(ip.Results.Indicator_A)
    gsize = self.n_i';
    aflag = aflagMaker(gsize);
else
    aflag = ip.Results.Indicator_A;
end
yy=data;
aflag0=unique(aflag); %% Levels of Factor A
k=length(aflag0); %% Number of Factor A's levels

%% specifying the sample sizes of each cell
vmu=[];A=[];A2=[];B2=[]; gsize=nan(k,1); S_ii=cell(k, 1);
for i=1:k
    iflag=(aflag==aflag0(i));
    yi=yy(iflag,:);
    gsize(i,:)=size(yi,1);
    ni=gsize(i);
    mui=mean(yi);
    vmu=[vmu;mui];
    ri=yi-ones(ni,1)*mui;
    if N>p
        Si=ri'*ri/(ni-1);
    else
        Si=ri*ri'/(ni-1);
    end

    S_ii{i} = Si;

    Ai=trace(Si);
    A = [A;Ai];

    if contains(method, ["Naive", "BiasReduced"])  % Check for Naive or Bias-Reduced
        Bi=trace(Si^2);
        switch method
            case {"L2-Naive", "F-Naive"} %% the naive method
                A2i=Ai^2;
                B2i=Bi;
            case  {"L2-BiasReduced", "F-BiasReduced"} %% the bias-reduced method
                A2i=(ni-1)*ni/(ni-2)/(ni+1)*(Ai^2-2*Bi/ni);
                B2i=(ni-1)^2/ni/(ni+1)*(Bi-Ai^2/(ni-1));
        end
        A2=[A2;A2i];
        B2=[B2;B2i];
    end

end

%% Computing the statistic
D=diag(gsize.^(-1));
H=sqrtm(inv(Contrast*D*Contrast')); %% q x q matrix
stat0=trace(H*(Contrast*vmu-c)*(Contrast*vmu-c)'*H');  %%

if any(strcmp(method, ["L2-Naive", "L2-BiasReduced", "F-Naive", "F-BiasReduced"]))
    Dh=sqrt(D);
    W=Dh*Contrast'*H'*H*Contrast*Dh;  %% k x k matrix
    
    dd=diag(W);
    K1=sum(dd.*A);
    K2a=sum(dd.^2.*A2);
    K2b=sum(dd.^2.*B2);
    AB1=[];AB2=[];
    for i=1:(k-1)
        ni=gsize(i);
        iflag=(aflag==aflag0(i));
        yi=yy(iflag,:);
        ri=yi-ones(ni,1)*vmu(i,:);
        for j=(i+1):k
            nj=gsize(j);
            jflag=(aflag==aflag0(j));
            yj=yy(jflag,:);
            rj=yj-ones(nj,1)*vmu(j,:);
            if N>p
                temp=trace(ri'*ri*rj'*rj)/(ni-1)/(nj-1);
            else
                temp=trace(ri*rj'*rj*ri')/(ni-1)/(nj-1);
            end
            K2a=K2a+2*W(i,i)*W(j,j)*A(i)*A(j);
            AB1=[AB1,A(i)*A(j)];
            K2b=K2b+2*W(i,j)^2*temp;
            AB2=[AB2,temp];
        end
    end
end

if any(strcmp(method, ["F-Bootstrap", "F-Simul"]))
    Dh=sqrt(D);
    W=Dh*Contrast'*H'*H*Contrast*Dh;  %% k x k matrix

    dd=diag(W);
    K1=sum(dd.*A);

    f_stat = stat0/K1;

 % Adjust A_n_ii depending on contrast Matrix
    switch self.Hypothesis
        case {"FAMILY"}
            b_n = gsize.^(1/2);
            A_n = eye(k) - (b_n*b_n') / N;  % Equation 9.36
%           A_n_ii = ones(k, 1) - gsize/N;  % Equation 9.38
            A_n_ii = diag(A_n);
            mask = true(k, 1);
        case "PAIRWISE"
            mask = logical(abs(Contrast'));  % Find which Pairs
            g_n = gsize(mask);               % subset of cell size
            N_n = sum(g_n);                  % Total from pair cell-size
            k_n = numel(g_n);                % Dimension space

            b_n = g_n.^(1/2);
            A_n = eye(k_n) - (b_n*b_n') / N_n;  % Equation 9.36 Modified for PairWise test
%             A_n_ii = ones(k_n, 1) - g_n/N_n;
            A_n_ii = diag(A_n);
        case {'INTERACTION', 'PRIMARY', 'SECONDARY'}  % Could be combined with "FAMILY". 
            A_n = D^(1/2) * Contrast' * inv(Contrast*D*Contrast') * Contrast * D^(1/2);  % EQ 9.64 or 9.83
            A_n_ii = diag(A_n);
            mask = true(k, 1);
    end

% % Numerically Equivalent to above.
%     tr_gamma = [];
%     A_n_ii = ones(k, 1) - gsize/N;
%     for i=1:k
%         iflag=(aflag==aflag0(i));
%         yi=yy(iflag,:);
% 
%         tr_gamma_i = trace(cov(yi));
%         tr_gamma=[tr_gamma; tr_gamma_i];
%     end
% 
%     S_n = sum(A_n_ii .* tr_gamma);
%     stat0 = stat0 / S_n;
end

switch method
    case {"L2-Naive", "L2-BiasReduced"} %% L2-norm based test(s)
        beta=K2b/K1;
        df=K2a/K2b;
        stat=stat0/beta;
        pvalue=1-chi2cdf(stat,df);
        pstat=[stat0,pvalue];
        params=[beta,df,K1,K2a,2*K2b];
    case {"F-Naive", "F-BiasReduced"}  %% F-type test(s)
        f_stat=stat0/K1;
        K2c=sum((dd./gsize).^2.*B2./(gsize-1));
        df1=K2a/K2b;
        df2=K2a/K2c;
        pvalue=1-fcdf(f_stat,df1,df2);
        pstat=[f_stat,pvalue];
        params=[df1,df2,K2a,2*K2b,2*K2c];

    case "L2-Bootstrap" %% Bootstrap test (Pretty sure its the L2 type)
        Bstat=nan(self.N_boot, 1);

        ts = self.setUpTimeBar(method);

        parfor ii=1:self.N_boot

            Bmu=[];
            for i=1:k
                iflag=(aflag==aflag0(i));
                yi=yy(iflag,:);
                ni=gsize(i);
                % Bflag=fix(rand(ni,1)*(ni-1))+1;
                Bflag = randsample(ni, ni, true);
                Byi=yi(Bflag,:);
                Bmui=mean(Byi);
                Bmu=[Bmu;Bmui];
            end

            %% Computing the statistic
            temp=H*Contrast*(Bmu-vmu); %%
            temp=trace(temp*temp');  %%
            Bstat(ii) = temp;

            ts.progress();
        end
        ts.stop;ts.delete;

        pvalue=mean(Bstat>stat0);  % Confirmed that its the L-2 Norm Stat

        pstat=[stat0,pvalue];

    case "F-Bootstrap"
        Bstat=nan(self.N_boot, 1);
        ts = self.setUpTimeBar(method);

        parfor ii=1:self.N_boot

            Bmu=[];
            tr_gamma=[];
            for i=1:k  % Iterate over all k groups
                iflag=(aflag==aflag0(i));
                yi=yy(iflag,:);
                ni=gsize(i);
                % Bflag=fix(rand(ni,1)*(ni-1))+1;
                Bflag = randsample(ni, ni, true);
                Byi=yi(Bflag,:);
                Bmui=mean(Byi);
                Bmu=[Bmu;Bmui];


                if mask(i)  % Run only if Asked for said group
                    % Faster Than Calling Covariance
                    z_mean = Byi - Bmui;                    % demean of ith group of k
                    test_cov  =(z_mean * z_mean')/ (ni-1);  % covariance of ith group of k
                    tr_gamma_i = trace(test_cov);           % trace of covariance from ith group of k
    
    %                 tr_gamma_i = trace(cov(Byi));
                    tr_gamma=[tr_gamma; tr_gamma_i];
                end

            end

            %% Computing the statistic
            temp=H*Contrast*(Bmu-vmu); %%
            T_n=trace(temp*temp');  %%

            S_n = sum( A_n_ii .* tr_gamma );
            temp = T_n / S_n;

            Bstat(ii) = temp;

            ts.progress();
        end
        ts.stop;ts.delete;

        pvalue=mean(Bstat>f_stat);  

        pstat=[f_stat, pvalue];

    case "L2-Simul"  % Works for OneWay_BF and now TwoWay_BF
        build_Covar_star = zeros(self.n_domain_points, 0);
        mask = any(logical(Contrast'), 2);
        COV_Sum = 0;
        for i=1:k
            if mask(i)
                iflag=(aflag==aflag0(i));
                yi=yy(iflag,:);
                gsize(i,:)=size(yi,1);
                ni=gsize(i);
                mui=mean(yi);
                vmu=[vmu;mui];
                ri=yi-ones(ni,1)*mui;
                COV_Sum = COV_Sum +  ( cov(ri) * (ni-1) );
                build_Covar_star = [build_Covar_star; ri];
            end

        end
        g_n = gsize(mask);               % subset of cell size
        N_n = sum(g_n);                  % Total from pair cell-size
        k_n = numel(g_n);                % Dimension space

        COV_Sum = COV_Sum ./ ((N_n-k_n));    % Pooled Covariance

        eig_gamma_hat = eig(COV_Sum);
        eig_gamma_hat = eig_gamma_hat(eig_gamma_hat>0);


        switch self.Hypothesis
            case {'FAMILY', 'PAIRWISE'}
                q = k_n-1;
            case {'PRIMARY','SECONDARY','INTERACTION'}
                % key change: df equals the number of independent linear restrictions
                q = rank(Contrast);  % fix
        end
        
        T_null = functionalANOVA.Chi_sq_mixture(q, eig_gamma_hat, self.N_simul);

        T_NullFitted = fitdist(T_null, "Kernel");
        pvalue = 1 - T_NullFitted.cdf(stat0);

        pstat=[stat0, pvalue];


    case "F-Simul"
        Dh=sqrt(D);
        W=Dh*Contrast'*H'*H*Contrast*Dh;  %% k x k matrix

        dd=diag(W);
        K1=sum(dd.*A);
        f_stat=stat0/K1;


        build_Covar_star = zeros(self.n_domain_points, 0);
        COV_Sum = 0;
        for i=1:k
            if mask(i)
                iflag=(aflag==aflag0(i));
                yi=yy(iflag,:);
                gsize(i,:)=size(yi,1);
                ni=gsize(i);
                mui=mean(yi);
                vmu=[vmu;mui];
                ri=yi-ones(ni,1)*mui;
                COV_Sum = COV_Sum +  ( cov(ri) * (ni-1) );
                build_Covar_star = [build_Covar_star; ri];
            end
        end
        g_n = gsize(mask);               % subset of cell size
        N_n = sum(g_n);                  % Total from pair cell-size
        k_n = numel(g_n);                % Dimension space

        COV_Sum = COV_Sum ./ ((N_n-k_n));    % Pooled Covariance

        eig_gamma_hat = eig(COV_Sum);
        eig_gamma_hat = eig_gamma_hat(eig_gamma_hat>0);

        switch self.Hypothesis
            case {'FAMILY', 'PAIRWISE'}
                q = k_n-1;
            case {'PRIMARY','SECONDARY','INTERACTION'}
                % key change: df equals the number of independent linear restrictions
                q = rank(Contrast);  % fix
        end

        T_null = functionalANOVA.Chi_sq_mixture(q, eig_gamma_hat, self.N_simul);

        % calculate Denominator
        S_null = zeros(self.N_simul, 1);
        S_ii = S_ii(mask);
        for i=1:k_n
            eig_gamma_hat = eig(S_ii{i});
            eig_gamma_hat = eig_gamma_hat(eig_gamma_hat > 0);
            S_temp = functionalANOVA.Chi_sq_mixture( g_n(i) - 1, eig_gamma_hat, self.N_simul);
            S_temp = (S_temp * A_n_ii(i)) ./ (g_n(i) - 1); % EQ (9.53) and EQ (9.98)
            S_null = S_null + S_temp;
        end

        F_null = T_null ./ S_null;
        F_NullFitted = fitdist(F_null, "Kernel");
        pvalue = 1 - F_NullFitted.cdf(f_stat);

        pstat=[f_stat, pvalue];

        % Equation 9.53 how the numerator,S_n, is distributed: One-Way
        % Equation 9.98 how the numerator,S_n, is distributed: Two-Way
end

stat = pstat(1);
pvalue = pstat(2);

end


function aflag  = aflagMaker(n_i)

aflag=[];

for K = 1 :numel(n_i)
    idicator = repmat(K, n_i(K), 1);
    aflag = [aflag; idicator];
end

end


