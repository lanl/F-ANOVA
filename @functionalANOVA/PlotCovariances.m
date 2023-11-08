function PlotCovariances(self, varargin)
%PLOTCOVARIANCES Plots Covariance Matrix Images for OneWay & TwoWay F-ANOVA 
% PLOTCOVARIANCES(self,...)
%
% Plotting tool for funtionalANOVA class. Default behavior is to plot the
% covariances for each group/level within a factor. OneWAY F-ANOVA
% utilizes the GroupLabels parameters while TwoWay F-ANOVA supports Primary
% and Secondary Labels for the respective factor groups/levels. TwoWay F-ANOVA
% plotting is specified by the "plotType" and requires the class to have a
% non-empty SubgroupIndicator property.
%
% <strong>Optional Inputs</strong>
%           plotType (String, default='default')
%                    - Type of plots to generate
%                    - Options are: 'default', 'primary', 'secondary', or 'interaction'
%                    - OneWay F-ANOVA only suppports: 'default'
%                    - For TwoWay F-ANOVA, 'primary' and 'default' are identical
%                      and shows the covariance from the primary factor levels
%                    - For TwoWay F-ANOVA, 'secondary' shows the covariance 
%                      from the secondary factor levels 
%                    - For TwoWay F-ANOVA, 'interaction' shows the covariances
%                      from of all combinations of the primary and secondary
%                      factor levels
%  SubgroupIndicator ([Nx1] numeric or [1xK] Cell Array, default=[])
%                    - Indicator array denoting the secondary factor levels for
%                      TwoWay F-ANOVA.
%                    - N represents the number of observations and must
%                      match the total number of concatenated observations
%                    - A represents the number of primary factor levels. Cell
%                      array of indicator arrays, one for for each primary
%                      factor level.
%        GroupLabels (String or Cell String Array [1xK], default=[])
%                     - For OneWay F-ANOVA, array labeling the
%                       different levels within the Main(primary) factor
%      PrimaryLabels (String or Cell String Array [1xK], default=[])
%                     - for TwoWay F-ANOVA, Array labeling the different
%                       Primary factor levels
%    SecondaryLabels (String or Cell String Array [1xB], default=[])
%                    - for TwoWay F-ANOVA, Array labeling the different
%                      Secondary factor levels
%   domainUnitsLabel ([1x1] string, default='')
%                    - Independent variable units
% responseUnitsLabel ([1x1] string, default='')
%                    - Dependent variable units
%           <>Scale (String, default='')
%                    - Replace '<>' with either 'x', 'y', or 'color'
%                    - Plot scaling for x, y, or color datum
%                    - Options are: '', 'linear', or 'log'.
%                    - Empty defaults to 'linear'
%        titleLabels ([1x1] cell, default={})
%                    - Labels to use for plot title.
%                    - If class was instantiated with Echo Records, then the
%                      label values from the records with those label names. 
%           savePath ([1x1] string, default='')
%                    - If empty, then figure won't be saved. Else, figure is
%                      saved in that file path 
%           position ([4x1] numeric, default= [90,257,2000,800])
%                    - Position and size of graphic.
%
% See also FUNCTIONALANOVA

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

p = inputParser;

addParameter(p, 'plotType', 'default', @(x) any(strcmpi(x, {'default', 'primary', 'secondary', 'interaction'})))
addParameter(p, 'SubgroupIndicator', self.SubgroupIndicator, @(x) isnumeric(x) || iscell(x)) % Indicator Array Denoting Subgroups
addParameter(p, 'GroupLabels', self.GroupLabels, @(x) ischar(x) || isstring(x));
addParameter(p, 'PrimaryLabels', self.PrimaryLabels, @(x) iscellstr(x) || isstring(x))
addParameter(p, 'SecondaryLabels', self.SecondaryLabels, @(x) iscellstr(x) || isstring(x))
addParameter(p, 'xScale', '', @(x)  any(strcmp(x, { '','log', 'linear'})))
addParameter(p, 'yScale', '', @(y) any(strcmp(y, {'','log', 'linear'})))
addParameter(p, 'colorScale', '', @(c) any(strcmp(c, {'','log', 'linear'})))
addParameter(p, 'domainUnitsLabel', self.domainUnitsLabel, @(x) ischar(x) || isstring(x));
addParameter(p, 'responseUnitsLabel', self.responseUnitsLabel, @(x) ischar(x) || isstring(x));
addParameter(p, 'titleLabels', {}, @(x) iscell(x));
addParameter(p, 'savePath', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'position', [90,257,2000,800], @(x) isnumeric(x) || numel(x)==4)


parse(p, varargin{:})
self.checkForNewLabels(p)


plotType = upper(p.Results.plotType);
self.SubgroupIndicator = p.Results.SubgroupIndicator;
self.GroupLabels = p.Results.GroupLabels;
self.PrimaryLabels = p.Results.PrimaryLabels;
self.SecondaryLabels = p.Results.SecondaryLabels;
self.domainUnitsLabel = p.Results.domainUnitsLabel;
self.responseUnitsLabel = p.Results.responseUnitsLabel;
xScale = p.Results.xScale;
yScale = p.Results.yScale;
colorScale = p.Results.colorScale;
titleLabels = p.Results.titleLabels;
savePath = p.Results.savePath;
position = p.Results.position;

fig_label = 'Group';
tempLabel = '';

if isempty(self.SubgroupIndicator)
    assert(strcmp(plotType, 'DEFAULT'), 'TwoWay plotting options require a SubgroupIndicator argument')
    
    if self.genericGroupLabels
        display_label = "Group " + string(1:self.k_groups);
    else
        display_label = self.GroupLabels;
    end

else
    self.setUpTwoWay()  % Creates Indicator Matrices and Labels
    switch plotType
        case {'DEFAULT', 'PRIMARY'}
            plotType = 'PRIMARY';
            fig_label = 'Primary Factor';
            display_label = self.PrimaryLabels;
            if self.genericGroupLabels
                display_label = "Group " + display_label;
            end
        case "SECONDARY"
            fig_label = 'Secondary Factor';
            display_label = self.SecondaryLabels;
            if self.genericGroupLabels
                display_label = "Group " + display_label;
            end
        case "INTERACTION"
            fig_label = 'Primary & Secondary Factor';
            % Prepare Labels
            combinations = generateTwoWayComb(self);
            display_label = combinations;
    end
end

if ~isempty(self.EchoEnsembleRecs)  % If Echo Ensemble Records (get MetaData)

    if ~isempty(titleLabels)
        Title_Labels=self.EchoEnsembleRecs.makeSummaryString(titleLabels, true);
        saveLabels =self.EchoEnsembleRecs.makeSummaryString(titleLabels, true, 'sanitizeString', true);
    else
        Title_Labels = '';
        saveLabels = '';
    end
else
    if ~isempty(titleLabels)
        Title_Labels = titleLabels;
        saveLabels = Title_Labels;
    else
        Title_Labels = '';
        saveLabels = '';
    end
end

if isempty(self.responseUnitsLabel)
    color_barLabel = {'(Response)^2'};
else
    color_barLabel = {sprintf("(%s)^2", self.responseUnitsLabel)};
end

switch plotType
    case {'DEFAULT', 'PRIMARY'}
        [gamma_hat_i, pooled_covar, n_ii] = MakeCovariances(self.data, self.k_groups, self.n_domain_points);
        K = self.k_groups;

    case 'SECONDARY'
        yy = [];
        for K = 1:self.k_groups
            yy = [yy; self.data{K}'];
        end

        bflag  = self.SubgroupIndicator;
        bflag0=unique(self.SubgroupIndicator); %% levels of Factor B
        N_Secondary = numel(bflag0);

        subData = cell(1, N_Secondary);

        for K = 1:N_Secondary
            flag = bflag==bflag0(K);
            subData{K} = yy(flag, :)';
        end

        [gamma_hat_i, pooled_covar, n_ii] = MakeCovariances(subData, N_Secondary, self.n_domain_points);
        K = N_Secondary;

    case 'INTERACTION'

        yy = [];
        for K = 1:self.k_groups
            yy = [yy; self.data{K}'];
        end

        aflag  = functionalANOVA.aflagMaker(self.n_i');
        bflag  = self.SubgroupIndicator;
        aflag0 = unique(aflag); %% Levels of Factor A
        p = length(aflag0); %% Number of Factor A's levels

        bflag0=unique(self.SubgroupIndicator); %% levels of Factor B
        q=length(bflag0); %% Number of Factor B's level

        ab=p*q;  % also known as n in the textbook
        self.AB_groups = ab;

        subData = cell(1, ab);
        counter = 1;
        for i=1:p
            for j=1:q
                ijflag=(aflag==aflag0(i))&(bflag==bflag0(j));
                yyi=yy(ijflag,:);
                subData{counter} = yyi';
                counter = counter + 1;
            end
        end

        [gamma_hat_i, pooled_covar, n_ii] = MakeCovariances(subData, ab, self.n_domain_points);
        K = ab;

end





% for k = 1 : size(gamma_hat_i, 3)
%     cmin = min(gamma_hat_i(:, :, k), [], 'all');
%     cmax = max(gamma_hat_i(:, :, k), [], 'all');
%     gamma_hat_i(:, :, k) = normalize(gamma_hat_i(:, :, 1), "range", [cmin, cmax]);
% end



% Visualize Kernel over entire Domain
FIG1 = figure('Name', sprintf('%s Covariances Visualized', fig_label), 'units', 'points', 'position', position);


[M, N] = GridGenerator(K);

M_exta = M + 1;
plot_counter = 1;
for ii = 1:K
    cmin = min(gamma_hat_i(:, :, ii), [], 'all');
    cmax = max(gamma_hat_i(:, :, ii), [], 'all');

    if mod(plot_counter, M_exta) == 0
        plot_counter = plot_counter + 1;
        subh = subplot(N, M_exta, plot_counter);
    else
        subh = subplot(N, M_exta, plot_counter);
    end

    image(self.d_grid, self.d_grid, gamma_hat_i(:,:,ii), "CDataMapping","scaled")
    clim(subh,[cmin,cmax])
    title(display_label(ii) + " | n = " + n_ii(ii))%, "interpreter", 'latex', 'fontSize', 28)
    C = colorbar;
    if ~isempty(self.EchoEnsembleRecs)  % If Echo Ensemble Records (use MetaData)
        if isempty(self.domainUnitsLabel)
            tempLabel = self.domainLabel;
        else
            tempLabel = self.domainLabel + sprintf(' (%s)', self.domainUnitsLabel);
        end

    else
        if ~isempty(self.domainUnitsLabel)
            tempLabel = sprintf('(%s)', self.domainUnitsLabel);
        end
    end

    xlabel(tempLabel)
    ylabel(tempLabel)

    if isempty(xScale)
        set(gca,'Xscale','linear')
    else
        set(gca,'Xscale', xScale)
    end

    if isempty(yScale)
        set(gca,'Yscale', 'linear')
    else
        set(gca,'Yscale', yScale)
    end

    plot_counter = plot_counter + 1;
end

% Pooled Options

cmin = min(pooled_covar, [], 'all');
cmax = max(pooled_covar, [], 'all');

subh = subplot(1, K + 1, K + 1);
image(self.d_grid, self.d_grid, pooled_covar, "CDataMapping","scaled")
clim(subh,[cmin,cmax])
title("Pooled" + " | n = " + self.N)
xlabel(tempLabel)
ylabel(tempLabel)
sgtitle({sprintf('%s Covariances', fig_label),Title_Labels})

C = colorbar;

title(C, color_barLabel)
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 18)

saveLabels = "Covariance_" + saveLabels;


if isfolder(savePath)
    saveas(FIG1, fullfile(savePath, saveLabels), 'png')
    close(FIG1)
end


if isempty(xScale)
    set(gca,'Xscale','linear')
else
    set(gca,'Xscale', xScale)
end

if isempty(yScale)
    set(gca,'Yscale', 'linear')
else
    set(gca,'Yscale', yScale)
end

if isempty(colorScale)
    set(gca,'ColorScale', 'linear')
else
    set(gca,'ColorScale', colorScale)
end

end

function [gamma_hat_i, pooled_covar, n_ii] = MakeCovariances(DATA, K, domainPoints)
% Compute basic group means and covariances
eta_i = zeros(domainPoints, K); % estimated means
v_hat_i = cell(K); %residuals matrix
gamma_hat_i = zeros(domainPoints, domainPoints,K); % groupwise covariance
n_ii = zeros(1, K);
for kk = 1:K
    n_i = size(DATA{kk}, 2);
    n_ii(kk) = n_i;
    eta_i(:, kk) = mean(DATA{kk}, 2); %estimated means
    zeroMeanData_k_subset = DATA{kk} - eta_i(:, kk); %residuals
    v_hat_i{kk} = zeroMeanData_k_subset; %residuals
    gamma_hat_i(:,:,kk) = 1 / (n_i - 1) .* (zeroMeanData_k_subset * zeroMeanData_k_subset'); % groupwise covariance
end


N = sum(n_ii);
% Compute pooled covariance
pooled_covar_terms = zeros(domainPoints, domainPoints, K);
for kk = 1:K
    pooled_covar_terms(:,:,kk) = (n_ii(kk) - 1) .* gamma_hat_i(:,:,kk);
end

pooled_covar = sum(pooled_covar_terms, 3) ./ (N - K);

end

function [M, N] = GridGenerator(K)

%Check if Odd or Even
if mod(K, 2) == 0
    EvenGroups = true;
else
    EvenGroups = false;
end

if K <= 3
    N = 1;
    M = K;
elseif K == 4
    N = 2;
    M = 2;
elseif K == 5
    N = 3;
    M = 2;

elseif K == 6
    N = 3;
    M = 2;
elseif K == 7
    N = 4;
    M = 3;
elseif K >= 8 && EvenGroups
    N = 4;
    M = ceil(K/N);
elseif K >= 9 && ~EvenGroups
    if floor(sqrt(K))^2 == K
        N = sqrt(K);
        M = N;
    else
        N = 4;
        M = ceil(K/N);
    end
end

end