function [pvalue, stat] = TwoWay(self, method, data, Contrast)
% TWOWAY Perform Homogeneous Two-way F-ANOVA
% [PVALUE, STAT]=TWOWAY(SELF, METHOD, DATA,  CONTRAST)  Homogeneous Two-way
% ANOVA for functional data Internally called through the 
% funtionalANOVA class method: TwoWayANOVA
%
% <strong>Required Input</strong>
%
%   method ([1x1] String)
%          - Options:  "L2-Simul", "L2-Naive", "L2-BiasReduced", "L2-Bootstrap",
%            "F-Simul", "F-Naive", "F-BiasReduced", "F-Bootstrap"
%     data ([Nxm] Numeric)
%          - N: Number of functional samples/realizations
%          - m: number of data points within a sample/realization
% Contrast ([1xG] Numeric)
%          - Contrast vector for specific hypotheses.
%          - Let p and q be the number of levels within the primary
%            and secondary factors respectively.
%          - Contrast for Primary levels,   then G=p.
%          - Contrast for Secondary levels, then G=q.
% <strong>Output</strong>
%
%   pvalue ([1x1] Numeric)
%          - p value from running the test
%     stat ([1x1] Numeric)
%          - test statistic value from either the L2 or F-type
%            methods
%
% See also FUNCTIONALANOVA, TWOWAYANOVA

% History
% Inital   April 27, 2010, (Jin-Ting Zhang)
% Revised  April 18, 2011, National University of Singapore (NUS), Singapore
% Revised  Oct   14, 2011, Princeton University, USA
% Revised  Dec   15, 2011, Princeton University, USA
% Revised  Feb   14, 2012, Princeton University, USA
% Revised  Feb   23, 2012, Princeton University, USA
% Revised  Jan   09, 2013, NUS
% Modified June  06, 2023, Los Alamos National Laboratory, USA


N = self.N;
ddim = self.n_domain_points;
aflag  = functionalANOVA.aflagMaker(self.n_i');
bflag  = self.SubgroupIndicator;
aflag0 = unique(aflag); %% Levels of Factor A (Primary Factor)
p = length(aflag0); %% Number of Factor A's levels

yy=data;

bflag0=unique(self.SubgroupIndicator); %% levels of Factor B (Secondary Factor)
q=length(bflag0); %% Number of Factor B's level

k=p*q;  % also known as n in the textbook


%% specifying the sample sizes of each cell
gsize=[];vmu=[];V=[];
p_cell = cell(1, p);
error_string = {};
for i=1:p
    q_cell = nan(1, q);
    for j=1:q
        ijflag=(aflag==aflag0(i))&(bflag==bflag0(j));
        ni=sum(ijflag);

        if ni == 0
            error_string{end + 1} = sprintf('A missing combination of data occurs for Primary Label: %s and Secondary Label: %s', self.PrimaryLabels(i), self.SecondaryLabels(j));
        end

        gsize=[gsize;ni];
        yyi=yy(ijflag,:);
        cell_mui=mean(yyi);
        vmu=[vmu;cell_mui];
        V=[V;yyi-ones(ni,1)*cell_mui]; %% N x ddim
        q_cell(j) = ni;
    end
    p_cell{i} = q_cell;
end

if ~isempty(error_string)
    for K = 1:length(error_string)
        warning(error_string{K})
    end
    error('Missing Combinations Listed Above')
end

self.n_ii = p_cell;
n_ii = cell2mat(p_cell);

switch self.Hypothesis
    case {"PAIRWISE", "FAMILY"}
        % continue
    otherwise
        switch self.Weights
            case "UNIFORM"
                u=ones(p,1)/p;
                v=ones(q,1)/q;
            case "PROPORTIONAL"
                matsize=reshape(gsize, p, q); %% pxq matrix
                N=sum(gsize);
                u=sum(matsize, 2) / N;   %% row sums/N: px1
                v=sum(matsize, 1)' / N;  %% column sums/N: qx1
        end
        % Specify the constant coefficient matrix
        Hp=[eye(p-1),-ones(p-1,1)];
        Hq=[eye(q-1),-ones(q-1,1)];
        Ap=eye(p)-ones(p,1)*u';
        Aq=eye(q)-ones(q,1)*v';
end

switch self.Hypothesis
    case {"PAIRWISE", "FAMILY"}
        % continue
        r=size(Contrast,1);  % also known as q in the textbook
    case "INTERACTION"   %% Interaction test
        Contrast=kron(Hp,Hq)*kron(Ap,Aq);
        r=(p-1)*(q-1);
    case "PRIMARY"  %% Main effect A (Primnary) test
        Contrast=Hp*kron(Ap, v');
        r=p-1;
    case "SECONDARY" %% Main effect B (Secondary) test
        Contrast=Hq*kron(u',Aq);
        r=q-1;
    otherwise % Custom Hypothesis

        switch self.Contrast_Factor
            case 1  % Primary Factor
                Contrast = Contrast*kron(Ap, v');
            case 2  % Secondary Factor
                Contrast = Contrast*kron(u', Aq);
        end
        r=size(Contrast,1);  % also known as q in the textbook
end


%% computing pointwise SSH, SSE 


W=inv(Contrast*diag(gsize.^(-1))*Contrast');
SSH=diag((Contrast*vmu)'*W*(Contrast*vmu));
SSH0=sum(SSH);

SSE=diag(V'*V);
SSE0=sum(SSE);

pool_coeff = N-k;

if N>ddim
    Pooled_COVAR=(V'*V)./pool_coeff;  %% pxp pooled covariance matrix
else
    Pooled_COVAR=(V*V')./pool_coeff;  %% Nx N matrix whch has the same non-zero eigen values with Sigma
end

A=trace(Pooled_COVAR);
B=trace(Pooled_COVAR^2);

if contains(method,  "Naive") %% the naive method
    A2=A^2;
    B2=B;
end

if  contains(method,  "BiasReduced") %% the bias-reduced method
    A2=(N-k)*(N-k+1)/(N-k-1)/(N-k+2)*(A^2-2*B/(N-k+1));
    B2=(N-k)^2/(N-k-1)/(N-k+2)*(B-A^2/(N-k));
end

switch method

    case "L2-Simul"
        stat=SSH0;

        eig_gamma_hat = eig(Pooled_COVAR);
        eig_gamma_hat = eig_gamma_hat(eig_gamma_hat>0);   % Only Positive eigen values

        SSH_null = functionalANOVA.Chi_sq_mixture(r, eig_gamma_hat, self.N_simul);
        SSH_NullFitted = fitdist(SSH_null, "Kernel");
        pvalue =  1 - SSH_NullFitted.cdf(stat);
        pstat=[stat,pvalue];

    case {"L2-Naive", "L2-BiasReduced"} %% L2-norm based test
        stat=SSH0;
        beta=B2/A;
        kappa=A2/B2;
        pvalue=1-chi2cdf(stat/beta,r*kappa);
        pstat=[stat/beta,pvalue];

    case {"F-Naive", "F-BiasReduced"} %% L2-norm based test
        stat=SSH0/SSE0*(N-k)/r;
        kappa=A2/B2;
        pvalue=1-fcdf(stat,r*kappa, (N-k)*kappa);
        pstat=[stat,pvalue];

    case "F-Simul"
        stat=SSH0/SSE0*(N-k)/r;

        eig_gamma_hat = eig(Pooled_COVAR);
        eig_gamma_hat = eig_gamma_hat(eig_gamma_hat>0);   % Only Positive eigen values

        SSH_null = functionalANOVA.Chi_sq_mixture(r, eig_gamma_hat, self.N_simul);
        SSE_null = functionalANOVA.Chi_sq_mixture(N-k, eig_gamma_hat, self.N_simul); 

        ratio = (N-k) / r;

        F_Null = SSH_null./SSE_null .* ratio;
        F_NullFitted =  fitdist(F_Null, "Kernel");
        pvalue =  1 - F_NullFitted.cdf(stat);
        pstat=[stat,pvalue];

    case "L2-Bootstrap"
        stat=SSH0;
        Bstat=nan(self.N_boot, 1);

        ts = self.setUpTimeBar(method);
        
        parfor b = 1 : self.N_boot

            Bmu=[];V=[];
            counter = 1;

            for i=1:p
                for j=1:q
                    ijflag=(aflag==aflag0(i))&(bflag==bflag0(j));
                    ni = n_ii(counter);

                    yi=yy(ijflag,:);

                    Bootflag=randi([1, ni], ni, 1);

                    Byi=yi(Bootflag,:);
                    Bmui=mean(Byi);
                    Bmu=[Bmu;Bmui];

                    counter = counter + 1;
                end
            end
            W=inv(Contrast*diag(gsize.^(-1))*Contrast');
            SSH=diag( (Bmu-vmu)'*(Contrast)'*W*(Contrast*(Bmu-vmu)) );
            Bstat(b) = sum(SSH);
            ts.progress();
        end
        ts.stop;ts.delete;
        pvalue=mean(Bstat>stat);  % Confirmed that its the L-2 Norm Stat
        
    case  "F-Bootstrap"
        stat=SSH0/SSE0*(N-k)/r;

        ts = self.setUpTimeBar(method);

        ratio = (N-k) / r;

        parfor b = 1 : self.N_boot
            Bmu=[];V=[];
            counter = 1;

            for i=1:p
                for j=1:q
                    ijflag=(aflag==aflag0(i))&(bflag==bflag0(j));
                    ni = n_ii(counter);

                    yi=yy(ijflag,:);

                    Bootflag=randi([1, ni], ni, 1);

                    Byi=yi(Bootflag,:);
                    Bmui=mean(Byi);
                    Bmu=[Bmu;Bmui];

                    V=[V; Byi-ones(ni,1)*Bmui]; %% Nx ddim

                    counter = counter + 1;
                end
            end

            W=inv(Contrast*diag(gsize.^(-1))*Contrast');
            SSH=diag( (Bmu-vmu)'*(Contrast)'*W*(Contrast*(Bmu-vmu)) );
            SSE=diag(V'*V);


            Bstat(b) = sum(SSH)/sum(SSE) * ratio;
            ts.progress();
        end

        ts.stop;ts.delete;
        pvalue=mean(Bstat>stat); 

end

end