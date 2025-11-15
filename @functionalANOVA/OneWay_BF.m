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

yy = data;
aflag0 = unique(aflag);                % Levels of Factor A
k = numel(aflag0);                     % Number of levels
m = size(yy,2);                        % = self.n_domain_points (alias)

% --- Method flags (avoid contains() in the loop)
isNaive = strcmp(method, "L2-Naive")  || strcmp(method, "F-Naive");
isBias  = strcmp(method, "L2-BiasReduced") || strcmp(method, "F-BiasReduced");
needAB  = isNaive || isBias;

% --- Preallocate (no growth inside loop)
vmu   = zeros(k, m, 'like', yy);
A     = zeros(k, 1);
if needAB
    A2 = zeros(k, 1);
    B2 = zeros(k, 1);
else
    A2 = []; 
    B2 = [];
end
gsize = zeros(k, 1);
S_ii  = cell(k, 1);

% --- Preserve original S_ii sizing rule (global)
storeAsH = (N > p);   % true => S_i is m x m (R'R); false => n_i x n_i (RR')

for i = 1:k
    mask = (aflag == aflag0(i));
    Yi   = yy(mask, :);
    ni   = size(Yi, 1);
    gsize(i) = ni;

    % Mean & residuals (implicit expansion; no ones(ni,1)*mui)
    mui     = mean(Yi, 1);
    vmu(i,:)= mui;
    Ri_raw  = Yi - mui;           % ni x m

    % Scale residuals once so S_i = (Ri_scaled' * Ri_scaled) or (Ri_scaled * Ri_scaled')
    denom       = ni - 1;
    inv_sqrt_dn = 1 / sqrt(denom);
    Ri          = Ri_raw .* inv_sqrt_dn;  % ni x m, scaled

    % A_i = trace(S_i) = ||Ri||_F^2 (works regardless of S_i shape)
    Ai = sum(Ri(:).^2);
    A(i) = Ai;

    % Build S_i once in the desired shape and reuse it for B_i and storage
    if storeAsH
        Si = Ri' * Ri;     % m x m
    else
        Si = Ri  * Ri';    % ni x ni
    end
    S_ii{i} = Si;

    if needAB
        % B_i = trace(S_i^2) = ||S_i||_F^2
        Bi = sum(Si(:).^2);

        if isNaive
            A2(i) = Ai^2;
            B2(i) = Bi;
        else
            % Bias-reduced (guard small ni)
            if ni > 2
                % Note: Ai and Bi already reflect scaling by (ni-1)
                A2(i) = (denom) * ni / (ni - 2) / (ni + 1) * (Ai^2 - 2*Bi/ni);
                B2(i) = (denom^2) / ni / (ni + 1) * (Bi - Ai^2/denom);
            else
                % Degenerate small-sample case (avoid division by zero)
                A2(i) = NaN;
                B2(i) = NaN;
            end
        end
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

    case "L2-Bootstrap"
        % ---- precompute things once (broadcast into parfor) ----
        HC = H * Contrast;                    % q x k
        m  = size(yy, 2);
    
        % group rows once; avoids recomputing logical masks inside parfor
        Ycell = cell(k,1);
        for i = 1:k
            Ycell{i} = yy(aflag == aflag0(i), :);   % ni x m for group i
        end
    
        Bstat = nan(self.N_boot, 1, 'like', yy);
        ts    = self.setUpTimeBar(method);
    
        % Tip: for reproducibility, set rng(seed,'combRecursive') before parfor
        parfor b = 1:self.N_boot
            % bootstrap means for all k groups
            Bmu = zeros(k, m, 'like', yy);
            for i = 1:k
                Yi = Ycell{i};
                ni = size(Yi,1);
    
                % uniform sampling with replacement, faster than randsample
                idx  = randi(ni, ni, 1);
                % mean of the resampled rows
                Bmu(i,:) = mean(Yi(idx,:), 1);
            end
    
            % T_n = || H*C*(Bmu - vmu) ||_F^2
            Bmu = Bmu - vmu;                  % k x m (in-place diff)
            T    = HC * Bmu;                  % q x m
            Bstat(b) = sum(T(:).^2);          % Frobenius norm squared
    
            ts.progress();
        end
        ts.stop; ts.delete;
    
        pvalue = mean(Bstat > stat0);         % L2 statistic
        pstat  = [stat0, pvalue];


    case "F-Bootstrap"
        % ---- precompute once (broadcast into parfor) ----
        HC = H * Contrast;                      % q x k
        m  = size(yy, 2);
    
        % group slices & sizes (avoid recomputing masks)
        Ycell   = cell(k,1);
        nvec    = zeros(k,1);
        invDen  = zeros(k,1);
        for i = 1:k
            Yi      = yy(aflag == aflag0(i), :);
            Ycell{i}= Yi;
            nvec(i) = size(Yi,1);
            invDen(i)= 1/(nvec(i)-1);
        end
    
        % mask → positions for tr_gamma
        mask_idx  = find(mask);
        k_mask    = numel(mask_idx);
        mask_pos  = zeros(k,1);
        mask_pos(mask_idx) = 1:k_mask;          % position of each masked group in tr_gamma
    
        Bstat = nan(self.N_boot, 1, 'like', yy);
        ts    = self.setUpTimeBar(method);
    
        % Tip: call rng(self.Seed,'combRecursive') before parfor for reproducibility
        parfor b = 1:self.N_boot
            % bootstrap means for all k groups
            Bmu = zeros(k, m, 'like', yy);
            % only for masked groups: traces of covariances
            tr_gamma = zeros(k_mask, 1);
    
            for i = 1:k
                Yi = Ycell{i};
                ni = nvec(i);
    
                % uniform bootstrap sample (faster than randsample)
                idx  = randi(ni, ni, 1);
                Byi  = Yi(idx, :);
                Bmui = mean(Byi, 1);
                Bmu(i, :) = Bmui;
    
                % trace(cov(Byi)) = ||Byi - mean||_F^2 / (ni-1)
                if mask(i)
                    sse = sum((Byi - Bmui).^2, 'all');
                    tr_gamma(mask_pos(i)) = sse * invDen(i);
                end
            end
    
            % T_n = || H*C*(Bmu - vmu) ||_F^2
            Tn = sum((HC * (Bmu - vmu)).^2, 'all');
    
            % S_n = sum_i A_n_ii(i) * tr_gamma(i) over masked groups
            S_n = sum(A_n_ii .* tr_gamma);
    
            Bstat(b) = Tn / S_n;
    
            ts.progress();
        end
        ts.stop; ts.delete;
    
        pvalue = mean(Bstat > f_stat);
        pstat  = [f_stat, pvalue];

    case "L2-Simul"  % Works for OneWay_BF and TwoWay_BF
        % groups used by the contrast
        mask = any(Contrast.' ~= 0, 2);        % k x 1 logical
        g_n  = gsize(mask);
        N_n  = sum(g_n);
        k_n  = nnz(mask);
    
        m = size(yy,2);                         % # domain points
    
        % pooled covariance numerator (sum_i R_i' R_i) without cov()
        COV_Sum = zeros(m, m, 'like', yy);
    
        if N > p
            % S_ii{i} = R_i'R_i/(n_i-1)  => multiply back by (n_i-1)
            for i = find(mask).'
                COV_Sum = COV_Sum + S_ii{i} * (gsize(i) - 1);
            end
        else
            % need R_i'R_i explicitly
            for i = find(mask).'
                Yi = yy(aflag == aflag0(i), :);
                Ri = Yi - vmu(i,:);            % demean with precomputed mean
                COV_Sum = COV_Sum + (Ri' * Ri);
            end
        end
    
        % pooled covariance
        COV_Sum = COV_Sum ./ (N_n - k_n);
    
        % spectrum (positive part only)
        eig_gamma_hat = eig(COV_Sum, 'vector');
        eig_gamma_hat = eig_gamma_hat(eig_gamma_hat > 0);
    
        % effective df for the mixture
        switch self.Hypothesis
            case {'FAMILY','PAIRWISE'}
                q = k_n - 1;
            case {'PRIMARY','SECONDARY','INTERACTION'}
                q = rank(Contrast);
        end
    
        % null draws and p-value
        T_null       = functionalANOVA.Chi_sq_mixture(q, eig_gamma_hat, self.N_simul);
    
        % (Fast option) empirical p-value:
        pvalue = mean(T_null > stat0);
        % (Keep original smoothing:)
        % T_NullFitted = fitdist(T_null, "Kernel");
        % pvalue       = 1 - T_NullFitted.cdf(stat0);
    
        pstat = [stat0, pvalue];

    case "F-Simul"
        % ---- W, dd, K1, f_stat (no full diag, reuse H) ----------------------
        sqrt_d = 1 ./ sqrt(gsize(:));          % k x 1
        B      = Contrast .* sqrt_d.';         % q x k  (so C*D*C' = B*B')
        HB     = H * B;                         % q x k
        W      = HB' * HB;                      % k x k
        dd     = diag(W);
        K1     = sum(dd .* A);
        f_stat = stat0 / K1;
    
        % ---- pooled covariance for masked groups (no cov(), no growth) ------
        % mask and A_n_ii must already be set in the shared block above
        idx_mask = find(mask).';                % 1 x k_n
        g_n      = gsize(mask);                 % k_n x 1
        N_n      = sum(g_n);
        k_n      = numel(idx_mask);
    
        m = size(yy, 2);
        COV_Sum = zeros(m, m, 'like', yy);
    
        if N > p
            % S_ii{i} = R_i'R_i / (n_i-1): multiply back by (n_i-1)
            for ii = idx_mask
                COV_Sum = COV_Sum + S_ii{ii} * (gsize(ii) - 1);
            end
        else
            % Rebuild R_i'R_i only for needed groups
            for ii = idx_mask
                rows = (aflag == aflag0(ii));
                Ri   = yy(rows, :) - vmu(ii, :);   % demean with existing mean
                COV_Sum = COV_Sum + (Ri' * Ri);
            end
        end
    
        COV_Sum = COV_Sum ./ (N_n - k_n);      % pooled covariance
    
        % ---- spectrum & numerator null --------------------------------------
        lam_num = eig(COV_Sum, 'vector');
        lam_num = lam_num(lam_num > 0);
    
        switch self.Hypothesis
            case {'FAMILY','PAIRWISE'}
                qeff = k_n - 1;
            otherwise  % {'PRIMARY','SECONDARY','INTERACTION'}
                qeff = rank(Contrast);
        end
    
        T_null = functionalANOVA.Chi_sq_mixture(qeff, lam_num, self.N_simul);
    
        % ---- denominator null S_null (only masked groups) --------------------
        S_null    = zeros(self.N_simul, 1);
        Sii_mask  = S_ii(mask);
        for i = 1:k_n
            lam_i = eig(Sii_mask{i}, 'vector');
            lam_i = lam_i(lam_i > 0);
            S_t   = functionalANOVA.Chi_sq_mixture(g_n(i) - 1, lam_i, self.N_simul);
            S_t   = (S_t .* A_n_ii(i)) ./ (g_n(i) - 1);   % EQ (9.53)/(9.98)
            S_null = S_null + S_t;
        end
    
        % ---- final p-value ---------------------------------------------------
        F_null       = T_null ./ S_null;
    
        % (Fast empirical alternative)  % 
        pvalue = mean(F_null > f_stat);
    
        pstat = [f_stat, pvalue];


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


