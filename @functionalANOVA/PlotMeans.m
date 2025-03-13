function PlotMeans(self, varargin)
%PLOTMEANS Plot Means of data groups for either OneWay and TwoWay F-ANOVA
% PLOTMEANS(self,...)
%
% Plotting tool for funtionalANOVA class. Default behavior is to overlay
% the group means over the realizations of the data. OneWAY F-ANOVA
% utilizes the GroupLabels parameters while TwoWay F-ANOVA supports Primary
% and Secondary Labels for the respective factor groups/levels. TwoWay F-ANOVA
% plotting is specified by the "plotType" and requires the class to have a
% non-empty SubgroupIndicator property.
%
% <strong>Optional Inputs</strong>
%             plotType (String, default='default')
%                      - Type of plots to generate
%                      - Options are: 'default', 'primary', 'secondary', or 'interaction'
%                      - OneWay F-ANOVA only suppports: 'default'
%                      - For TwoWay F-ANOVA, 'primary' and 'default' are identical
%                        and show the mean response and data from the primary
%                        factor levels
%                      - For TwoWay F-ANOVA, 'secondary' shows the mean response 
%                        and data from the secondary factor levels 
%                      - For TwoWay F-ANOVA, 'interaction' shows the mean response and
%                        data of all combinations of the primary and secondary
%                        factor levels
%    subgroupIndicator ([Nx1] numeric or [Mx1] Cell Array, default=[])
%                      - Indicator array denoting the secondary factor levels for
%                        TwoWay F-ANOVA.
%                      - N represents the number of observations and must
%                        match the total number of concatenated observations
%                      - A represents the number of primary factor levels. Cell
%                        array of indicator arrays, one for for each primary
%                        factor level.  
% observationSizeLabel ([1x1] logical, default=true)
%                      - Provides observation count for each level in
%                        plot legend
%          GroupLabels (String or Cell String Array [1xK], default=[])
%                       - For OneWay F-ANOVA, array labeling the
%                         different levels within the Main(primary) factor
%        PrimaryLabels (String or Cell String Array [1xK], default=[])
%                       - for TwoWay F-ANOVA, Array labeling the different
%                         Primary factor levels
%      SecondaryLabels (String or Cell String Array [1xB], default=[])
%                      - for TwoWay F-ANOVA, Array labeling the different
%                        Secondary factor levels
%     domainUnitsLabel ([1x1] string, default='')
%                      - Independent variable units
%   responseUnitsLabel ([1x1] string, default='')
%                      - Dependent variable units
%       xScale/yScale (String, default='')
%                      - Plot scaling for x or y datum
%                      - Options are: '', 'linear', or 'log'.
%                      - Empty defaults to 'linear'
%     dataTransparency ([1x1] numeric, default=0.1)
%                      - Line transparency for data
%   legendTransparency ([1x1] numeric, default=1/3)
%                      - Legend line transparency for data
%        dataLineWidth ([1x1] numeric, default=1.75)
%                      - Data line width
%        meanLineWidth ([1x1] numeric, default=5)
%                      - Mean line width
%             fontSize ([1x1] numeric, default=18)
%                       - Font Size for entire graphic
%          titleLabels ([1x1] cell, default={})
%                      - Labels to use for plot title.
%                      - If class was instantiated with Echo Records, then the
%                        label values from the records with those label names. 
%             savePath ([1x1] string, default='')
%                      - If empty, then figure won't be saved. Else, figure is
%                        saved in that file path 
%       legendLocation ([1x1] string, default='best')
%                      - Legend location using MATLABs syntax.
%           numColumns ([1x1] numeric, default=1)
%                      - Number of columns to for Legend
%          legendTitle ([1x1] string, default='')
%                      - Custom labels to use for legend title. Will overwrite
%                        default generated one depending on the type of F-ANOVA used
%                        and the plotType specified.
%            newcolors ([Kx3] numeric, default=[])
%                      - Custom RGB color array to plot levels or groups. 
%                      - Depending on the plotType, K must match the number of
%                       primary and/or secondary factor levels.
%             position ([4x1] numeric, default=[90, 90, 1400, 800])
%                      - Position and size of graphic.
% 
% See also FUNCTIONALANOVA

p=inputParser;
addParameter(p, 'plotType', 'default', @(x) any(strcmpi(x, {'default', 'primary', 'secondary', 'interaction'})))
addParameter(p, 'subgroupIndicator', self.SubgroupIndicator, @(x) isnumeric(x) || iscell(x)) % Indicator Array Denoting Subgroups
addParameter(p, 'observationSizeLabel', true, @(x) islogical(x))
addParameter(p, 'GroupLabels', self.GroupLabels, @(x) iscellstr(x) || isstring(x))
addParameter(p, 'PrimaryLabels', self.PrimaryLabels, @(x) iscellstr(x) || isstring(x))
addParameter(p, 'SecondaryLabels', self.SecondaryLabels, @(x) iscellstr(x) || isstring(x))
addParameter(p, 'xScale', '', @(x)  any(strcmp(x, { '','log', 'linear'})))
addParameter(p, 'yScale', '', @(y) any(strcmp(y, {'','log', 'linear'})))
addParameter(p, 'domainUnitsLabel', self.domainUnitsLabel, @(x) ischar(x) || isstring(x));
addParameter(p, 'responseUnitsLabel', self.responseUnitsLabel, @(x) ischar(x) || isstring(x));
addParameter(p, 'dataTransparency', 0.1, @(x) x<=1 || x >= 0)
addParameter(p, 'legendTransparency', 0.3333, @(x) x<=1 || x >= 0)
addParameter(p, 'dataLineWidth', 1.75, @(x)  x >= 0)
addParameter(p, 'meanLineWidth', 5, @(x)  x >= 0)
addParameter(p, 'fontSize', 18, @(x)  x >= 0)
addParameter(p, 'titleLabels', {}, @(x) iscell(x));
addParameter(p, 'savePath', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'legendLocation', 'best', @(x) any(strcmpi(x, LegendLocations())))
addParameter(p, 'numColumns', 1, @(x) isnumeric(x) && x>=1)
addParameter(p, 'legendTitle', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'newcolors', [], @(x) isnumeric(x))
addParameter(p, 'position', [90, 90, 1400, 800], @(x) isnumeric(x) || numel(x)==4)


parse(p, varargin{:})
self.checkForNewLabels(p)

self.SubgroupIndicator = p.Results.subgroupIndicator;
self.GroupLabels = p.Results.GroupLabels;
self.PrimaryLabels = p.Results.PrimaryLabels;
self.SecondaryLabels = p.Results.SecondaryLabels;
plotType = upper(p.Results.plotType);
xScale = p.Results.xScale;
yScale = p.Results.yScale;
dataTransparency = p.Results.dataTransparency;
legendTransparency = p.Results.legendTransparency;
dataLineWidth =  p.Results.dataLineWidth;
meanLineWidth = p.Results.meanLineWidth;
fontSize = p.Results.fontSize;
titleLabels = p.Results.titleLabels;
savePath = p.Results.savePath;
legendLoc = p.Results.legendLocation;
legendNumColumns = p.Results.numColumns;
legendTitle = p.Results.legendTitle;
newcolors = p.Results.newcolors;
position = p.Results.position;
self.domainUnitsLabel = p.Results.domainUnitsLabel;
self.responseUnitsLabel = p.Results.responseUnitsLabel;
observationSizeLabel = p.Results.observationSizeLabel;


domainLabel = self.domainLabel;
responseLabel = self.responseLabel;



if isempty(self.SubgroupIndicator)
    assert(strcmp(plotType, 'DEFAULT'), 'TwoWay plotting options require a SubgroupIndicator argument')
    TheLabels = self.GroupLabels;
    nLabels = self.n_i';

else
    self.setUpTwoWay()  % Creates Indicator Matrices and Labels
    switch plotType
        case {'DEFAULT', 'PRIMARY'}
            plotType = 'PRIMARY';
            TheLabels = self.PrimaryLabels;

            nLabels = self.n_i';

        case "SECONDARY"
            TheLabels = self.SecondaryLabels;
            nLabels = zeros(self.B_groups, 1);
            for K = 1 :self.A_groups
                for KK = 1 : self.B_groups
                    nLabels(KK) = nLabels(KK) + self.n_ii{K}(KK) ;
                end
            end
        case "INTERACTION"
            % Prepare Labels
            TheLabels = generateTwoWayComb(self);

            DataLabels_interact = cellfun(@(r) string(r), self.n_ii, 'UniformOutput', false);
            nLabels = [DataLabels_interact{:}]';
    end
end


if observationSizeLabel
    if self.genericGroupLabels
        TheDataLabels = ": Group Data " + "(" + string(nLabels) + ")";
    else
        TheDataLabels = ": Data " + "(" + string(nLabels) + ")";
    end
    TheDataLabels = TheDataLabels(:)';
end


if self.genericGroupLabels
    TheLabels = TheLabels(:)';
    MeangroupLabels = TheLabels + ": Group Mean";
    Data_groupLabels = TheLabels + TheDataLabels;

else
    TheLabels = TheLabels(:)';
    MeangroupLabels = TheLabels + ": Mean";
    Data_groupLabels = TheLabels  + TheDataLabels;

end

if ~isempty(self.EchoEnsembleRecs)  % If Echo Ensemble Records (get MetaData)
    if ~isempty(titleLabels)
        Title_Labels=self.EchoEnsembleRecs.makeSummaryString(titleLabels, true);
        saveLabels =self.EchoEnsembleRecs.makeSummaryString(titleLabels, true, 'sanitizeString', true);
    else
        Title_Labels = '';
        saveLabels = '';
    end

    classType = class(self.EchoEnsembleRecs.pullSample.sources);

else
    if ~isempty(titleLabels)
        Title_Labels = titleLabels;
        saveLabels = Title_Labels;
    else
        Title_Labels = '';
        saveLabels = '';
    end

    classType = '';
end

switch plotType
    case 'DEFAULT'

        mean_groups = cellfun(@(x) mean(x, 2), self.data, 'UniformOutput',false);

        if isempty(newcolors)
            if self.k_groups < 7
                colorList = lines(self.k_groups);
            else
                colorList = turbo(self.k_groups);
            end

        else
            assert(size(colororder, 2) == 3, 'Color order Matrix must have 3 columns reperesenting RGB values from 0 to 1')
            mesg = sprintf('Color order Matrix must have at least %i rows reperesenting all the groups', self.k_groups);
            assert(size(colororder, 1)  >= self.k_groups, mesg)
            colorList = newcolors;
        end

        FIG1 = figure('Name','GroupMeans and Realizations', 'units', 'points', 'position', position);
        hold on

        PlotReals = cell(1, self.k_groups);
        legendExampleLines = cell(1, self.k_groups);
        for K = 1:self.k_groups
            PlotReals{K} = plot(self.d_grid, self.data{K}, "Color", colorList(K, :), 'lineWidth', dataLineWidth);

            for KK = 1 :numel(PlotReals{K})
                PlotReals{K}(KK).Color(4)=dataTransparency;
            end
            legendExampleLines{K} = copy(PlotReals{K}(1));
            legendExampleLines{K}.Color(4)=legendTransparency;
        end

        PlotMeans = cell(1, self.k_groups);
        for K = 1:self.k_groups
            PlotMeans{K} = plot(self.d_grid, mean_groups{K}, "Color", colorList(K, :), 'lineWidth',meanLineWidth,'lineStyle','--');
        end

        title({'F-ANOVA Group Means and Realizations', Title_Labels})
        LG = legend([PlotMeans{:}, legendExampleLines{:}], [MeangroupLabels, Data_groupLabels], 'location', legendLoc, 'NumColumns', legendNumColumns);

        if ~isempty(legendTitle)
            title(LG, legendTitle);
        end

    case 'PRIMARY'

        if isempty(newcolors)
            if self.k_groups < 7
                colorList = lines(self.k_groups);
            else
                colorList = turbo(self.k_groups);
            end

        else
            assert(size(colororder, 2) == 3, 'Color order Matrix must have 3 columns reperesenting RGB values from 0 to 1')
            mesg = sprintf('Color order matrix must have at least %i rows reperesenting all the Primary factor levels', self.k_groups);
            assert(size(colororder, 1)  >= self.k_groups, mesg)
            colorList = newcolors;
        end


%         MeangroupLabels = string(self.PrimaryLabels) + ": Mean";
%         Data_groupLabels = self.PrimaryLabels  + ": Data";

        mean_groups = cellfun(@(x) mean(x, 2), self.data, 'UniformOutput',false);


        FIG1 = figure('Name','Primary Factor Means and Realizations', 'units', 'points', 'position', position);
        hold on

        PlotReals = cell(1, self.k_groups);
        legendExampleLines = cell(1, self.k_groups);
        for K = 1:self.k_groups
            PlotReals{K} = plot(self.d_grid, self.data{K}, "Color", colorList(K, :), 'lineWidth', dataLineWidth);

            for KK = 1 :numel(PlotReals{K})
                PlotReals{K}(KK).Color(4)=dataTransparency;
            end
            legendExampleLines{K} = copy(PlotReals{K}(1));
            legendExampleLines{K}.Color(4)=legendTransparency;
        end

        PlotMeans = cell(1, self.k_groups);
        for K = 1:self.k_groups
            PlotMeans{K} = plot(self.d_grid, mean_groups{K}, "Color", colorList(K, :), 'lineWidth',meanLineWidth,'lineStyle','--');
        end

        title({'F-ANOVA Primary Factor Means and Realizations', Title_Labels})
        LG = legend([PlotMeans{:}, legendExampleLines{:}], [MeangroupLabels, Data_groupLabels], 'location', legendLoc, 'NumColumns', legendNumColumns);


        if isempty(legendTitle)
            title(LG, "Primary Factor Levels")
        else
            title(LG, legendTitle);
        end

    case 'SECONDARY'
%         MeangroupLabels = string(self.SecondaryLabels) + ": Mean";
%         Data_groupLabels = self.SecondaryLabels  + ": Data";

        FIG1 = figure('Name','Secondary Factor Means and Realizations', 'units', 'points', 'position', position);
        hold on

        yy = [];
        for K = 1:self.k_groups
            yy = [yy; self.data{K}'];
        end

        bflag  = self.SubgroupIndicator;
        bflag0=unique(self.SubgroupIndicator); %% levels of Factor B
        N_Secondary = numel(bflag0);

        if isempty(newcolors)
            if N_Secondary < 7
                colorList = lines(N_Secondary);
            else
                colorList = turbo(N_Secondary);
            end
        else
            assert(size(colororder, 2) == 3, 'Color order Matrix must have 3 columns reperesenting RGB values from 0 to 1')
            mesg = sprintf('Color order matrix must have at least %i rows reperesenting all the Secondary factor levels', N_Secondary);
            assert(size(colororder, 1)  >= N_Secondary, mesg)
            colorList = newcolors;
        end

        PlotMeans = cell(1, N_Secondary);
        mean_groups = cell(1, N_Secondary);
        subData = cell(1, N_Secondary);

        for K = 1:N_Secondary
            flag = bflag==bflag0(K);
            subData{K} = yy(flag, :)';
            mean_groups{K} = mean(subData{K}, 2);
        end

        PlotReals = cell(1, N_Secondary);
        legendExampleLines = cell(1, N_Secondary);

        for K = 1:N_Secondary
            PlotReals{K} = plot(self.d_grid, subData{K}, "Color", colorList(K, :), 'lineWidth', dataLineWidth);

            for KK = 1 :numel(PlotReals{K})
                PlotReals{K}(KK).Color(4)=dataTransparency;
            end
            legendExampleLines{K} = copy(PlotReals{K}(1));
            legendExampleLines{K}.Color(4)=legendTransparency;
        end


        for K = 1:N_Secondary
            PlotMeans{K} = plot(self.d_grid, mean_groups{K}, "Color", colorList(K, :), 'lineWidth',meanLineWidth,'lineStyle','--');
        end


        title({'F-ANOVA Secondary Factor Means and Realizations', Title_Labels})
        LG = legend([PlotMeans{:}, legendExampleLines{:}], [MeangroupLabels, Data_groupLabels], 'location', legendLoc, 'NumColumns', legendNumColumns);

        if isempty(legendTitle)
            title(LG, "Secondary Factor Levels")
        else
            title(LG, legendTitle);
        end
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

        % Prepare Labels
        combinations = generateTwoWayComb(self);
        MeanLabels = strings(ab, 1);
        DataLabels = strings(ab, 1);
        for K = 1:numel(combinations)
            MeanLabels(K) = combinations(K) + ": Mean";

            if observationSizeLabel
                DataLabels(K) = combinations(K)  + ": Data " + "(" + string(nLabels(K)) + ")";
            else
                DataLabels(K) = combinations(K)  + ": Data";
            end
        end


        %% specifying the sample sizes of each cell
        if isempty(newcolors)
            if ab < 7
                colorList = lines(ab);
            else
                colorList = turbo(ab);
            end
        else
            assert(size(colororder, 2) == 3, 'Color order Matrix must have 3 columns reperesenting RGB values from 0 to 1')
            mesg = sprintf('Color order matrix must have at least %i rows reperesenting all the TwoWay factor levels', ab);
            assert(size(colororder, 1)  >= ab, mesg)
            colorList = newcolors;
        end


        PlotMeans = cell(1, ab);
        mean_groups=cell(ab, 1);
        subData = cell(1, ab);
        counter = 1;
        for i=1:p
            for j=1:q
                ijflag=(aflag==aflag0(i))&(bflag==bflag0(j));
                yyi=yy(ijflag,:);
                subData{counter} = yyi';
                cell_mui=mean(yyi);
                mean_groups{counter}=cell_mui';
                counter = counter + 1;
            end
        end

        PlotReals = cell(1, ab);
        legendExampleLines = cell(1, ab);

        FIG1 = figure('Name','Primary & Secondary Combinatorial Means and Realizations', 'units', 'points', 'position', position);
        hold on

        for K = 1:ab
            PlotReals{K} = plot(self.d_grid, subData{K}, "Color", colorList(K, :), 'lineWidth', dataLineWidth);

            for KK = 1 :numel(PlotReals{K})
                PlotReals{K}(KK).Color(4)=dataTransparency;
            end
            legendExampleLines{K} = copy(PlotReals{K}(1));
            legendExampleLines{K}.Color(4)=legendTransparency;
        end
       
        for K = 1:ab
            PlotMeans{K} = plot(self.d_grid, mean_groups{K}, "Color", colorList(K, :), 'lineWidth',meanLineWidth,'lineStyle','--');
        end

        title({'F-ANOVA Primary and Secondary Factor Combinatorial Means and Realizations', Title_Labels})
        LG = legend([PlotMeans{:}, legendExampleLines{:}], [MeanLabels, DataLabels], 'location', legendLoc, 'NumColumns', legendNumColumns);

        if isempty(legendTitle)
            title(LG, "TwoWay Factor Levels")
        else
            title(LG, legendTitle);
        end


end

if ~isempty(self.EchoEnsembleRecs)  % If Echo Ensemble Records (use MetaData)
    if isempty(self.domainUnitsLabel)
        xlabel(domainLabel)
    else
        xlabel(domainLabel + sprintf(' (%s)', self.domainUnitsLabel))
    end

    if isempty(self.responseUnitsLabel)
        ylabel(responseLabel)
    else
        ylabel(responseLabel + sprintf(' (%s)', self.responseUnitsLabel))
    end
else
    if ~isempty(self.domainUnitsLabel)
        xlabel(sprintf('(%s)', self.domainUnitsLabel))
    end

    if ~isempty(self.responseUnitsLabel)
        ylabel(sprintf('(%s)', self.responseUnitsLabel))
    end
end


if any(strcmpi(classType, {'PSDRecord', 'SRSRecord', 'SpectralRecord'})) || contains(responseLabel, 'SRS', 'IgnoreCase',true) && (isempty(xScale) ||  isempty(yScale))
    if isempty(xScale)
        set(gca,'Xscale','log')
    else
        set(gca,'Xscale', xScale)
    end

    if isempty(yScale)
        set(gca,'Yscale','log')
    else
        set(gca,'Yscale', yScale)
    end
else
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
end

set(findall(gcf, '-property', 'FontSize'), 'FontSize', fontSize)

saveLabels = "GroupMeans_" + saveLabels;
if isfolder(savePath)
    saveas(FIG1, fullfile(savePath, saveLabels), 'png')
    close(FIG1)
end





end


function LG_loc = LegendLocations()
LG_loc  = {...
    'north',...
    'south',...
    'east',...
    'west',...
    'northeast',...
    'northwest',...
    'southeast',...
    'southwest',...
    'northoutside',...
    'southoutside',...
    'eastoutside',...
    'westoutside',...
    'northeastoutside',...
    'northwestoutside',...
    'southeastoutside',...
    'southwestoutside',...
    'best',...
    'bestoutside',...
    'none'};
end