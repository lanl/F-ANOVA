function [pvalue, stat] = TwoWay_BF(self, method, data, Contrast, c)
% TWOWAY_BF Perform Heteroscedastic Two-way F-ANOVA
% [PVALUE, STAT] = TWOWAY_BF(SELF, METHOD, DATA, CONTRAST, C) 
% Heteroscedastic Two-way ANOVA for functional data. Internally called
% through the funtionalANOVA class method: TwoWayANOVA_BF
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
%       c ([1x1] Numeric)
%          - Constant used for Hypothesis testing. Typically, c=0 to denote
%          that one expects a zero mean difference between factors and
%          levels within factors.
% <strong>Output</strong>
%
%   pvalue ([1x1] Numeric)
%          - p value from running the test
%     stat ([1x1] Numeric)
%          - test statistic value from either the L2 or F-type
%            methods
%
% See also FUNCTIONALANOVA, TWOWAYANOVA_BF

% History
% Initial  April 27, 2010 National University of Singapore (NUS), Singapore (Jin-Ting Zhang)
% Revised  Sept 29,  2010
% Revised  Jan 10,   2013  NUS
% Modified May 26,   2023 Los Alamos National Laboratory, USA

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

   
bflag = self.SubgroupIndicator;
N = self.N;
aflag  = functionalANOVA.aflagMaker(self.n_i');
dim=self.n_domain_points; %% data dimension

aflag0=unique(aflag); %% Levels of Factor A
p=length(aflag0); %% Number of Factor A's levels
bflag0=unique(bflag); %% levels of Factor B
q=length(bflag0); %% Number of Factor B's level

%% specifying the sample sizes of each cell
gsize=zeros(p,q);
yy=[];
p_cell = cell(1, p);
ij=0;
for i=1:p
    q_cell = nan(1, q);
    for j=1:q
        ij=ij+1;
        ijflag=(aflag==aflag0(i))&(bflag==bflag0(j));
        nij=sum(ijflag);
        gsize(i,j)=nij;
        yij=data(ijflag,:);
        yy=[yy;[ij*ones(nij,1),yij]];
        q_cell(j) = nij;
    end
    p_cell{i} = q_cell;
end

self.n_ii = p_cell;

%% Specify the weights
switch self.Weights
    case "UNIFORM"
        u=ones(p,1)/p;
        v=ones(q,1)/q;
    case "PROPORTIONAL"
        u=sum(gsize, 2) / N;   %% row sums/N: px1
        v=sum(gsize, 1)' / N;  %% column sums/N: qx1

        %     u=sum(gsize')'/N; %% px1
        %     v=sum(gsize)'/N;  %% qx1
end

%% Specify the constant coefficient matrix
Ap=eye(p)-ones(p,1)*u';
Aq=eye(q)-ones(q,1)*v';

if isempty(Contrast)
    Hp=[eye(p-1),-ones(p-1,1)];
    Hq=[eye(q-1),-ones(q-1,1)];
    switch self.Hypothesis
        case "INTERACTION"   %% Interaction test
            H=kron(Hp,Hq);
        case "PRIMARY"  %% Main effect A (Primnary) test
            H=Hp;
        case "SECONDARY" %% Main effect B (Secondary) test
            H=Hq;
    end
else
    H=Contrast;
end


switch self.Hypothesis
    case "INTERACTION"   %% Interaction test
        Contrast=H*kron(Ap,Aq);
    case "PRIMARY"  %% Main effect A (Primnary) test
        Contrast=H*kron(Ap,v');
    case "SECONDARY" %% Main effect B (Secondary) test
        Contrast=H*kron(u',Aq);
    otherwise

end

pure_data = yy(:, 2:end);
A = yy(:, 1);
[pvalue, stat] = self.OneWay_BF(method,pure_data , Contrast, c, 'Indicator_A',  A) ;
end

