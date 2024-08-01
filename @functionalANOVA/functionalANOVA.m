classdef functionalANOVA < handle & matlab.mixin.Copyable
    %functionalANOVA a class for performing functional analysis of variance (F-ANOVA)
    %
    %   Functional analyis of variance investigates whether categorical
    %   variables affect the mean response. OneWay ANOVA investigates how
    %   a single categorical variable with different levels affects the
    %   mean response. TwoWay ANOVA investigates how a two categorical
    %   variables, denoted Primary and Secondary Factors, each with different
    %   levels affects the mean response.
    %
    %   The Null hypothesis for all the different hypotheses is that the
    %   affect, or hypothesis, is <strong>not</strong> statistically significant.
    %   The alternative hypothesis is that the affect, or hypothesis, is
    %   statistically significant.
    %
    %   F-ANOVA is broken down into OneWay or TwoWay depending on the
    %   heirariachal nature of the experiment/date set. Run the
    %   "CovarianceTest" class method to verify the assumption that the
    %   covariance function between the factors and levels are all
    %   equaivalent. If the the test returns statistically signinficant,
    %   then the Behrens-Fischer (BF) ANOVA methods need to be used.
    %
    %   OneWay ANOVA allows for two different hypotheses.
    %   (1) "Family"
    %       Investigates whether at least one of the levels statistically affects
    %       the mean response. However, it doesnt provide which level is
    %       statistically significant.
    %   (2) "PairWise"
    %       Investigates all combinatorial pairs of the levels which allows for
    %       further analysis as to which pairs lead to a statistically
    %       signficant mean response.
    %
    %   TwoWay ANOVA allows for various different hypotheses.
    %   (1) "Family"
    %       Investigates whether at least one of the levels in either the
    %       primary or secondary factor statistically affects the mean
    %       response. Doesnt provide which level is statistically significant.
    %   (2) "PairWise"
    %       Investigates all combinatorial pairs of both primary and
    %       secondary factors which allows for further analysis as to
    %       which pairs lead to a statistically signficant mean response.
    %   (3) "Primary" or "Secondary"
    %       Investigates whether the primary or secondary factor,
    %       respectively, statistically affects the mean response.
    %   (4) "Interaction"
    %       Investigates whether there exists a statistically significant
    %       interaction affect that alters the mean response between
    %       at least one of the combinatorial pairs from primary and
    %       secondary factor levels.
    %   (4) "Custom"
    %       When provided with a user specified, contrast vector, a custom
    %       hypothesis. Contrast vector is a linear combination of unit
    %       vectors that sum to 0, in this case.
    %
    % <strong>Required Input</strong>
    %
    %          dataArray ([Mx1] Cell Array or Echo EnsembleRecord)
    %                    - Each cell or EnsembleRecord contains the data matrix for
    %                      each level within the Primary/Main factor
    %                    - Data matrix is assumed to be in long format. Where each
    %                      column represents a functional sample and rows are the
    %                      functional response over the common domain of interest
    %        boundsArray ([2x1] Numeric)
    %                    - The lower and upper bounds over the common domain of interest.
    %                    - If the analysis if over the entire domain: [-inf, inf]
    %
    % <strong>Optional Inputs</strong>
    %
    %  SubgroupIndicator ([Nx1] numeric or [Mx1] Cell Array, default=[])
    %                    - Indicator array denoting the secondary factor levels for
    %                      TwoWay F-ANOVA.
    %                    - N represents the number of observations and must
    %                      match the total number of concatenated observations
    %                    - A represents the number of primary factor levels. Cell
    %                      array of indicator arrays, one for for each primary
    %                      factor level.
    %        GroupLabels (String or Cell String Array [Mx1], default=[])
    %                    - For OneWay F-ANOVA, array labeling the different groups
    %      PrimaryLabels (String or Cell String Array [Mx1], default=[])
    %                    - for TwoWay F-ANOVA, Array labeling the different
    %                      Primary factor levels
    %    SecondaryLabels (String or Cell String Array [Bx1], default=[])
    %                    - for TwoWay F-ANOVA, Array labeling the different
    %                      Secondary factor levels

    % Read Only but Public
    properties(SetAccess=private)
        ANOVA_Methods = [...
            "L2-Simul", "L2-Naive", "L2-BiasReduced", "L2-Bootstrap",...
            "F-Simul", "F-Naive", "F-BiasReduced", "F-Bootstrap",...
            ];

        ANOVA_Methods_Used = [];

        COVAR_Methods = ["L2-Simul", "L2-Naive", "L2-BiasReduced", "Permutation-Test", "Bootstrap-Test"];
        COVAR_Methods_Used = [];


        OneWay_P_Table    % OneWay Homoscedastic   F-ANOVA
        OneWay_BF_P_Table % OneWay Heteroscedastic F-ANOVA

        TwoWay_P_Table    % TwoWay Homoscedastic   F-ANOVA
        TwoWay_BF_P_Table % TwoWay Heteroscedastic F-ANOVA

        COVAR_P_Table     % Eequality of Covariance Tests

        SW_P_Table % shapiro-wilkes table from PointwiseDistributionCheck Method

        Weights = "proportional"
        Hypothesis = "FAMILY";
    end

    % Public and Modifiable
    properties
        N_boot = 10000; % Number of Bootstrap Samples
        N_permutations = 10000 % % Number of Permutations for Covariance Equality Tests
        N_simul = 100000;  % Number of Random Variables to simulate to Estimate KDE From
        alpha = 0.05;  % Significance Level

        verbose = true;  % Show Summary and Results.

        genericGroupLabels    % true if generic labels are produce, false if user input them

        domainUnitsLabel = '';   % independent variable unit
        responseUnitsLabel = '';  % Response variable (data) unit

        domainLabel = '';   % independent variable unit
        responseLabel = '';  % Response variable (data) unit


    end

    % Hidden but Modifiable
    properties(Hidden, SetAccess=public)
        showSimulPlot = false;  % Shows Null distribution plots for "Simul" Methods
        consistantTable = true; % Modifies OneWay P-table to match the others
    end

    % not Public and not Modifiable
    properties(Hidden, SetAccess=private)

        H0_OneWay = ["FAMILY", "PAIRWISE"]; % Hypothesis Choices
        H0_TwoWay = ["FAMILY", "PAIRWISE", "INTERACTION", "PRIMARY", "SECONDARY", "CUSTOM"]; % Hypothesis Choices

        d_grid = [];    % Independent array
        k_groups        % Used for OneWay ANOVA. # of Groups


        SubgroupIndicator = []; % indicator Array for B

        A_groups        % Used for twoWay ANOVA (primary level) Unique # levels
        B_groups        % Used for twoWay ANOVA (secondary level) Unique # levels
        AB_groups       % total groups for twoWay ANOVA (primary & Secondary)

        Contrast = [];         % User Specified Constrast vector
        Contrast_Factor = [];  % Either 1 for Primary, 2 for Secondary

        GroupLabels = strings(1,0); % Group Labels for One-Way
        PrimaryLabels = strings(1,0); % Primary Factor labels for TwoWay
        SecondaryLabels = strings(1,0); % Secondary Factor Labels  TwoWay


        hypothesis_LABEL = '';

        data
        n_i            % Used for OneWay ANOVA for group count
        n_ii           % Nested cell for TwoWay ANOVA
        N              % Total number of observations

        n_domain_points
        EchoEnsembleRecs
        boundsArray

        lb_index
        ub_index

        table_sigFigs = 4;

        CriticalValues
        omega_hat = []; % used for covariance test

    end

    methods
        function self = functionalANOVA(dataArray, boundsArray, varargin)
            %FUNCTIONALANOVA Constructor for the FUNCTIONALANOVA Class
            % SELF = FUNCTIONALANOVA(DATAARRAY, BOUNDSARRAY, VARARGIN)
            %
            % <strong>Required Input</strong>
            %
            %          dataArray ([1xK] Cell Array or [1xK] Echo EnsembleRecords)
            %                    - Cell Array of Data consists of Matrix response data.
            %                    - One cell for each group/level within the main(primary) factor
            %                    - Data must be in long format. Columns are samples,
            %                      and rows represent values of the functional data
            %                    - All functional data must be sampled at the same
            %                      independent datum and are the same length
            %        boundsArray ([1x2] Numeric)
            %                    - Bounds to subset the funtional data response by
            %                      the values from d_grid.
            %                    - To utilize all the functional data, it is
            %                      recommended to use boundsArray=[-inf, inf]
            %
            % <strong>Optional Inputs</strong>
            %
            %             d_grid ([mx1] Numeric, default=[])
            %                    - Domain grid consisting of an array of the independent
            %                      variable datums.
            %                    - Required when passing in a cell array of functional
            %                      response data. Not Required for Echo Records.
            %             N_boot ([1x1] Numeric, default=10,000)
            %                    - Number of bootstrap replicate.
            %                    - Used for any methods that have "Bootstrap" in them.
            %            N_simul ([1x1] Numeric, default=10,000)
            %                    - Number of data points to simulate null distribution
            %                    - Used for any methods that have "Simul" in them.
            %      ANOVA_Methods ([1x1] string Array or Cell String, default=self.ANOVA_Methods)
            %                    - ANOVA Methods used to calculate P-values
            %              alpha ([1x1] Numeric, default=0.05)
            %                    - Critical value for assigning statistical significance
            %  SubgroupIndicator ([Nx1] numeric or [1xK] Cell Array, default=[])
            %                    - Indicator array denoting the secondary factor levels for
            %                      TwoWay F-ANOVA.
            %                    - N represents the number of observations and must
            %                      match the total number of concatenated observations
            %                    - K represents the number of primary factor levels. Cell
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


            p = inputParser;

            addRequired(p, 'dataArray', @(x)  iscell(x) || isa(x, 'EnsembleRecord'))  % function values are rows and realizations are columns
            addRequired(p, 'boundsArray', @(x) isnumeric(x)) % Domain for which to compute the L2 Statistic Over

            addParameter(p, 'd_grid', self.d_grid, @(x)  isnumeric(x) & min(size((x))))  % Entire Time Domain of Function. Must be mx1 or 1xm
            addParameter(p, 'N_boot', self.N_boot, @(x) isinteger(int64(x)));
            addParameter(p, 'N_simul', self.N_simul, @(x) isinteger(int64(x)))
            addParameter(p, 'ANOVA_Methods', self.ANOVA_Methods, @(x) iscellstr(x) || isstring(x))
            addParameter(p, 'alpha', self.alpha, @(x) x<=1 && x>= 0)
            addParameter(p, 'SubgroupIndicator', self.SubgroupIndicator, @(x) isnumeric(x) || iscell(x)) % Indicator Array Denoting Subgroups
            addParameter(p, 'GroupLabels', self.GroupLabels, @(x) iscellstr(x) || isstring(x))
            addParameter(p, 'PrimaryLabels', self.PrimaryLabels, @(x) iscellstr(x) || isstring(x))     % If twoWayANOVA
            addParameter(p, 'SecondaryLabels', self.SecondaryLabels, @(x) iscellstr(x) || isstring(x)) % If twoWayANOVA
            addParameter(p, 'domainUnitsLabel', self.domainUnitsLabel, @(x) ischar(x) || isstring(x));
            addParameter(p, 'responseUnitsLabel', self.responseUnitsLabel, @(x) ischar(x) || isstring(x));

            parse(p, dataArray, boundsArray, varargin{:});
            self.d_grid = p.Results.d_grid;
            self.N_boot = int64(p.Results.N_boot);
            self.N_simul = int64(p.Results.N_simul);
            self.ANOVA_Methods_Used = p.Results.ANOVA_Methods;
            self.alpha = p.Results.alpha;
            self.SubgroupIndicator = p.Results.SubgroupIndicator;
            self.GroupLabels = p.Results.GroupLabels;
            self.PrimaryLabels = p.Results.PrimaryLabels;
            self.SecondaryLabels = p.Results.SecondaryLabels;
            self.domainUnitsLabel = p.Results.domainUnitsLabel;
            self.responseUnitsLabel = p.Results.responseUnitsLabel;

            self.k_groups = length(dataArray);

            
            self.boundsArray = boundsArray;
            assert(numel(self.boundsArray) == 2, 'Requires 2 Bounds')
            assert(numel(unique(self.boundsArray)) == 2, 'Lower and Upper Bounds Must be Different from each Other')

            assert(self.k_groups >=2 , 'Must have at least 2 categories/Groups for F-ANOVA to work')

            % Check if EnsembleRecord, if Ensemble Generate dataMatrixCell
            if isa(dataArray, 'EnsembleRecord')
                self.EchoEnsembleRecs = dataArray;
                dataArray = cell(size(dataArray));
                d_grids = cell(size(dataArray));
                for K = 1 : self.k_groups
                    dataArray{K} = self.EchoEnsembleRecs(K).data;
                    d_grids{K} = self.EchoEnsembleRecs(K).y;
                end

                % Check if all the same size
                first_dim = cellfun(@(x) size(x, 1), d_grids);
                assert(length(unique(first_dim))==1, 'Ensembles have mismatching domain vectors')

                % Check all domains are equal
                domainData = horzcat(d_grids{:});
                assert(~any((diff(domainData,1, 2)), 'all'), 'The domain vectors from the Ensemble records don''t  have identical datumns')
                self.d_grid = domainData(:, 1);  % Just Pick First domain vector since all are equal

                self.n_i = cellfun(@(x) size(x, 2), dataArray)';  % Transpose,  number within each Group

                % Get Info
                self.domainUnitsLabel = unique(string(self.EchoEnsembleRecs.arrayEval(@(r) r.getUnits('y').symbol)));
                self.responseUnitsLabel = unique(string(self.EchoEnsembleRecs.arrayEval(@(r) r.getUnits('data').symbol)));

                assert(length(self.domainUnitsLabel) == 1, "Ensembles do not share a consistant unit across their domains. Convert to a consistant unit")
                assert(length(self.responseUnitsLabel) == 1, "Ensembles do not share a consistant unit across their response. Convert to a consistant unit")

                if (self.domainUnitsLabel == "")
                   self.domainUnitsLabel = strings(0,1);
                end

                if (self.responseUnitsLabel == "")
                    self.responseUnitsLabel = strings(0,1);
                end

                self.responseLabel = string(self.EchoEnsembleRecs.pullSample.getQuantity('data'));
                self.domainLabel = string(self.EchoEnsembleRecs.pullSample.getQuantity('y'));

                if isempty(self.domainLabel) || strcmpi(self.domainLabel, 'y')
                    self.domainLabel = string(self.EchoEnsembleRecs(1).sources.pullSample.getQuantity('x'));
                end

                if isempty(self.responseLabel) || strcmpi(self.responseLabel, 'data')
                    self.responseLabel = string(self.EchoEnsembleRecs(1).sources.pullSample.getQuantity('y'));
                end

            else
                self.EchoEnsembleRecs = [];
                assert(~isempty(self.d_grid), 'Must include parameter called ''d_grid'' when using a cell array of matrices')
                self.n_i = cellfun(@(x) size(x, 2), dataArray);  % number within each Group
            end
            self.N = sum(self.n_i); % Total Samples combined
            % Check All groups and Reps have same vector lengths as the domain
            n_rows_per_group = cellfun(@(x) size(x, 1), dataArray);
            assert(all(n_rows_per_group == length(self.d_grid)), 'All Groups and Replicates must have the same vector length as the domain')

            self.function_Subsetter();
            self.n_domain_points = numel(self.lb_index : self.ub_index);
            self.d_grid = self.d_grid(self.lb_index : self.ub_index); % subset Domain
            %% Trim and plot data
            for k = 1:self.k_groups
                self.data{k} = dataArray{k}(self.lb_index : self.ub_index, :); %trimming
            end

            if isempty(self.SubgroupIndicator) % One Way ANOVA set up
                if ~isempty(self.GroupLabels)
                    assert(numel(self.GroupLabels) == self.k_groups, "Each Group Must Have Exactly One Label Associated to it")
                    self.genericGroupLabels = false;
                else
                    self.GroupLabels = string(1:self.k_groups); % Automatically Assign Group Labels
                    self.genericGroupLabels = true;
                end
            else
                self.setUpTwoWay()  % Creates Indicator Matrices and Labels
            end


        end

        % OneWay Homoscedastic (Common Covariance)
        function OneWayANOVA(self, varargin)
            p = inputParser;
            % Keep ability to change inputs
            addParameter(p, 'N_boot', self.N_boot, @(x) isinteger(int64(x)));
            addParameter(p, 'N_simul', self.N_simul, @(x) isinteger(int64(x)))
            addParameter(p, 'ANOVA_Methods', self.ANOVA_Methods, @(x) iscellstr(x) || isstring(x))
            addParameter(p, 'alpha', self.alpha, @(x) x<=1 && x>= 0)
            addParameter(p, 'GroupLabels', self.GroupLabels, @(x) iscellstr(x) || isstring(x))
            addParameter(p, 'showSimulPlot', self.showSimulPlot, @(x) islogical(x));
            addParameter(p, 'Hypothesis', self.Hypothesis, @(x)  (iscellstr(x) || isstring(x) || ischar(x)) && any(strcmpi(x, self.H0_OneWay)))  % Hypothesis of interest

            parse(p, varargin{:});
            self.N_boot = int64(p.Results.N_boot);
            self.N_simul = int64(p.Results.N_simul);
            self.ANOVA_Methods_Used = p.Results.ANOVA_Methods;
            self.alpha = p.Results.alpha;
            self.GroupLabels = p.Results.GroupLabels;
            self.showSimulPlot = p.Results.showSimulPlot;
            self.Hypothesis = upper(p.Results.Hypothesis);  % Cast to Uppe
            self.ANOVA_Methods_Used = sort(self.ANOVA_Methods_Used, "descend");  % Maintain L2 first than F-test

            self.castANOVAMethods()
            %% Estimate K Group Mean, Grand Mean, Pooled Covariance Function, between-subject and  within-subject variations

            eta_i = zeros(self.n_domain_points, self.k_groups);
            % zeroMeanData_k_subset = cell(k_groups);

            build_Covar_star = zeros(self.n_domain_points, 0);

            for k = 1:self.k_groups
                eta_i(:, k) = mean(self.data{k}, 2); %estimated means
                zeroMeanData_k_subset = self.data{k} - eta_i(:, k);
                build_Covar_star = [build_Covar_star, zeroMeanData_k_subset];
            end

            eta_grand = sum(eta_i .* self.n_i, 2) / self.N;

            % Compute SSH(t)
            pair_vec = string();

            switch self.Hypothesis

                case "FAMILY"
                    q = self.k_groups - 1; %rank of contrast matrix
                    n_tests = 1;
                    SSH_k = zeros(self.n_domain_points, self.k_groups);
                    for k = 1:self.k_groups
                        SSH_k(:,k) = self.n_i(k) .* (eta_i(:,k) - eta_grand).^2;
                    end
                    SSH_t = sum(SSH_k,2);
                    pair_vec(1) = "FAMILY";
                    C = [];
                    D = [];
                case "PAIRWISE"
                    C = self.construct_pairwise_contrast_matrix(self.k_groups);
                    n_tests = size(C,1); %k choose 2

                    assert(numel(self.GroupLabels) == self.k_groups, "Each Group Must Have Exactly One Label Associated to it")

                    c = zeros(1, self.n_domain_points);
                    D = diag(1./self.n_i);
                    q = 1; %rank of PAIRWISE contrasts
                    SSH_t = zeros(self.n_domain_points, n_tests);
                    for cc = 1:n_tests
                        Ct = C(cc,:);
                        part12 = (Ct*eta_i' - c) ;
                        SSH_t(:,cc) = part12.^2 .* inv(Ct*D*Ct'); %done in this method because functions don't lend themselves to simple matrix algebra

                        %The following acquires label pairs
                        [t1, t2] = self.GroupLabels{logical(C(cc,:))};
                        pair_vec(cc) = convertCharsToStrings([t1 ' & ' t2]);
                    end

            end

            if (SSH_t < eps)
                warning('Pointwise between-subject variation is 0. Check for Duplicated Data Matrices')
            end

            %Compute estimated covariance
            gamma_hat = (1 / (self.N-self.k_groups)) * (build_Covar_star * build_Covar_star'); %changed transpose order

            % Compute SSE(t)
            SSE_t = (self.N - self.k_groups) .* diag(gamma_hat); % Pointwise within-subject  (Group/Categorical) variations


            %% Methods
            n_methods = numel(self.ANOVA_Methods_Used);

            switch self.Hypothesis
                case 'FAMILY'
                    self.OneWay_P_Table =  table(self.ANOVA_Methods_Used',  nan(n_methods, 1), nan(n_methods, 1), strings(n_methods, 1), strings(n_methods, 1), cell(n_methods, 1), ...
                        strings(n_methods, 1), cell(n_methods, 1), ...
                        'VariableNames', ...
                        {'Family-Wise Method', 'Test-Statistic' ,'P-Value', 'Verdict', 'Parameter 1 Name', 'Parameter 1 Value', 'Parameter 2 Name', 'Parameter 2 Value'});

                    % Generate ResultTable, Test Statistics, estimates used for Null Dist Approximations
                    [T_n, F_n, beta_hat, kappa_hat, beta_hat_unbias, kappa_hat_unbias] = self.SetUpANOVA(SSH_t, SSE_t, gamma_hat, q);
                case 'PAIRWISE'
                    [T_n, F_n, beta_hat, kappa_hat, beta_hat_unbias, kappa_hat_unbias] = self.SetUpANOVA(SSH_t, SSE_t, gamma_hat, q);

                    T_hypothesis = table(pair_vec', 'VariableNames', {'Hypothesis'});
            end

            % Set up param Input
            params = {T_n, F_n, beta_hat, kappa_hat, beta_hat_unbias, kappa_hat_unbias, C, D};

            %             n_methods = numel(self.ANOVA_Methods_Used);
            %             p_value_matrix = zeros(n_tests, n_methods);
            self.CriticalValues = cell(n_methods, 3);

            eig_gamma_hat = eig(gamma_hat);
            eig_gamma_hat = eig_gamma_hat(eig_gamma_hat>0);   % Only Positive eigen values

            % Perform elsewhere to clean up main class code
            p_value_matrix = self.OneWay(n_tests, q, eig_gamma_hat, eta_i, eta_grand, params, pair_vec);

            switch self.Hypothesis
                case "PAIRWISE"
                    T_p_value = array2table(p_value_matrix, 'VariableNames', cellstr(self.ANOVA_Methods_Used));
                    self.OneWay_P_Table = [T_hypothesis, T_p_value];
                case "FAMILY"
                    self.OneWay_P_Table.("P-Value") = p_value_matrix';
                    Signif_results = self.OneWay_P_Table.("P-Value") < self.alpha;
                    self.OneWay_P_Table{Signif_results, "Verdict"} = "Reject Null Hypothesis for Alternative Hypothesis";
                    self.OneWay_P_Table{~Signif_results, "Verdict"} = "Fail to Reject Null Hypothesis";

                    if self.consistantTable
                        self.OneWay_P_Table = self.OneWay_P_Table(:, 1:4);
                    end
            end

            if self.verbose
                self.DataSummaryReportOneWay("homoskedastic")
                self.showTable(self.OneWay_P_Table)
            end
        end

        % OneWay Heteroscedastic (Different Covariances)
        function OneWayANOVA_BF(self, varargin)

            p = inputParser;
            % Keep ability to change inputs
            addParameter(p, 'N_boot', self.N_boot, @(x) isinteger(int64(x)));
            addParameter(p, 'N_simul', self.N_simul, @(x) isinteger(int64(x)))
            addParameter(p, 'ANOVA_Methods', self.ANOVA_Methods, @(x) iscellstr(x) || isstring(x))
            addParameter(p, 'alpha', self.alpha, @(x) x<=1 && x>= 0)
            addParameter(p, 'GroupLabels', self.GroupLabels, @(x) iscellstr(x) || isstring(x))
            addParameter(p, 'Hypothesis', self.Hypothesis, @(x)  (iscellstr(x) || isstring(x) || ischar(x)) && any(strcmpi(x, self.H0_OneWay)))  % Hypothesis of interest

            parse(p, varargin{:});
            self.N_boot = int64(p.Results.N_boot);
            self.N_simul = int64(p.Results.N_simul);
            self.ANOVA_Methods_Used = p.Results.ANOVA_Methods;
            self.alpha = p.Results.alpha;
            self.GroupLabels = p.Results.GroupLabels;
            self.Hypothesis = upper(p.Results.Hypothesis);  % Cast to Uppe
            self.ANOVA_Methods_Used = sort(self.ANOVA_Methods_Used, "descend");  % Maintain L2 first than F-test
            self.castANOVAMethods()

            pair_vec = string();

            yy = [];
            for K = 1:self.k_groups
                yy = [yy; self.data{K}'];
            end

            switch self.Hypothesis
                case "PAIRWISE"
                    C = self.construct_pairwise_contrast_matrix(self.k_groups);
                    n_tests = size(C,1); %k choose 2
                    assert(numel(self.GroupLabels) == self.k_groups, "Each Group Must Have Exactly One Label Associated to it")

                case "FAMILY"
                    n_tests = 1;
                    C = [eye(self.k_groups-1),-ones(self.k_groups-1,1)];
                    pair_vec(1) = "FAMILY";
            end


            %% Methods
            n_methods = numel(self.ANOVA_Methods_Used);

            switch self.Hypothesis
                case 'FAMILY'
                    self.OneWay_BF_P_Table =  table(self.ANOVA_Methods_Used',  nan(n_methods, 1), nan(n_methods, 1), strings(n_methods, 1),  ...
                        'VariableNames', ...
                        {'Family-Wise Method', 'Test-Statistic' ,'P-Value', 'Verdict'});

                case 'PAIRWISE'
                    for cc = 1:n_tests
                        %The following acquires label pairs
                        [t1, t2] = self.GroupLabels{logical(C(cc,:))};
                        pair_vec(cc) = convertCharsToStrings([t1 ' & ' t2]);
                    end
            end

            T_hypothesis = table(pair_vec', 'VariableNames', {'Hypothesis'});
            p_value_matrix = nan(n_tests, n_methods);
            test_stat = nan(1, n_methods);

            counter = 0;
            c = 0;  % Assuming Equality to each other and not a generic constant
            for method = self.ANOVA_Methods_Used   % Generate Test Statistic
                counter = counter + 1;

                p_value = zeros(n_tests, 1);
                for ii = 1:n_tests
                    switch self.Hypothesis
                        case 'FAMILY'
                            C_input = C;
                        case 'PAIRWISE'
                            C_input = C(ii, :);
                            self.hypothesis_LABEL = pair_vec(ii);
                    end
                    [p_value(ii), statistic] = self.OneWay_BF(method, yy, C_input, c);
                end
                p_value_matrix(:, counter) = p_value;
                test_stat(counter) = statistic;

            end

            NotNan_mask = ~isnan(p_value_matrix);
            NotNan_mask = all(NotNan_mask, 1);

            if sum(NotNan_mask) ~= length(self.ANOVA_Methods_Used)
                if self.verbose
                    warning('Removing ANOVA results from methods that arent implemented yet.')
                end
            end

            switch self.Hypothesis
                case "PAIRWISE"
                    T_p_value = array2table(p_value_matrix(:, NotNan_mask), 'VariableNames', cellstr(self.ANOVA_Methods_Used(:, NotNan_mask)));
                    self.OneWay_BF_P_Table = [T_hypothesis, T_p_value];
                case "FAMILY"
                    self.OneWay_BF_P_Table = self.OneWay_BF_P_Table(NotNan_mask', :);
                    self.OneWay_BF_P_Table.("P-Value") = p_value_matrix(:, NotNan_mask)';
                    self.OneWay_BF_P_Table.('Test-Statistic') = test_stat(:, NotNan_mask)';
                    Signif_results = self.OneWay_BF_P_Table.("P-Value") < self.alpha;
                    self.OneWay_BF_P_Table{Signif_results, "Verdict"} = "Reject Null Hypothesis for Alternative Hypothesis";
                    self.OneWay_BF_P_Table{~Signif_results, "Verdict"} = "Fail to Reject Null Hypothesis";
            end

            if self.verbose
                self.DataSummaryReportOneWay("heteroskedastic")
                self.showTable(self.OneWay_BF_P_Table)
            end

        end

        % TwoWay Homoscedastic (Common Covariance)
        function TwoWayANOVA(self, varargin)
            self.Contrast = []; % Emptyat start

            p = inputParser;
            % Keep ability to change inputs
            addParameter(p, 'SubgroupIndicator', self.SubgroupIndicator, @(x) isnumeric(x)) % Indicator Array Denoting Subgroups
            addParameter(p, 'N_boot', self.N_boot, @(x) isinteger(int64(x)));
            addParameter(p, 'N_simul', self.N_simul, @(x) isinteger(int64(x)))
            addParameter(p, 'ANOVA_Methods', self.ANOVA_Methods, @(x) iscellstr(x) || isstring(x))
            addParameter(p, 'alpha', self.alpha, @(x) x<=1 && x>= 0)
            addParameter(p, 'Contrast', self.Contrast, @(x) isvector(x) || isempty(x))
            addParameter(p, 'PrimaryLabels', self.PrimaryLabels, @(x) iscellstr(x) || isstring(x))
            addParameter(p, 'SecondaryLabels', self.SecondaryLabels, @(x) iscellstr(x) || isstring(x))
            addParameter(p, 'Hypothesis', self.Hypothesis, @(x)  (iscellstr(x) || isstring(x) || ischar)  && any(strcmpi(x, self.H0_TwoWay)))  % Hypothesis of interest
            addParameter(p, 'Weights', self.Weights, @(x)  (iscellstr(x) || isstring(x) || ischar(x)) && any(strcmpi(x, ["uniform", "proportional"])))  % weighting style

            parse(p, varargin{:});
            self.SubgroupIndicator = p.Results.SubgroupIndicator;
            self.N_boot = int64(p.Results.N_boot);
            self.N_simul = int64(p.Results.N_simul);
            self.ANOVA_Methods_Used = p.Results.ANOVA_Methods;
            self.alpha = p.Results.alpha;
            self.Contrast = p.Results.Contrast;
            self.PrimaryLabels = p.Results.PrimaryLabels;
            self.SecondaryLabels = p.Results.SecondaryLabels;
            self.Weights = upper(p.Results.Weights);  % Cast to Upper
            self.Hypothesis = upper(p.Results.Hypothesis);  % Cast to Upper
            self.ANOVA_Methods_Used = sort(self.ANOVA_Methods_Used, "descend");  % Maintain L2 first than F-test

            self.castANOVAMethods() % Casts Methods to write case


            yy = [];
            for K = 1:self.k_groups
                yy = [yy; self.data{K}'];
            end

            if isempty(self.SubgroupIndicator)
                error('SubgroupIndicator must be provided for TwoWay ANOVAs')
            else
                self.setUpTwoWay()  % Creates Indicator Matrices and Labels
                self.setUpTwoWayHypothesis()  % For Custom Hypothesis
            end

            % Define C and number of tests
            switch self.Hypothesis
                case "PAIRWISE"
                    C = self.construct_pairwise_contrast_matrix(self.AB_groups);
                    n_tests = size(C,1); % AB_groups choose 2
                case "FAMILY"
                    n_tests = 1;
                    C = [eye(self.AB_groups-1),-ones(self.AB_groups-1,1)];
                case {"INTERACTION", "PRIMARY", "SECONDARY"}
                    n_tests = 1;
                    C = [];
                otherwise % Custom Hypothesis
                    n_tests = 1;
                    C = self.Contrast;
            end


            %% Methods
            n_methods = numel(self.ANOVA_Methods_Used);

            switch self.Hypothesis
                case 'FAMILY'
                    self.TwoWay_P_Table =  table(self.ANOVA_Methods_Used',  nan(n_methods, 1), nan(n_methods, 1), strings(n_methods, 1),  ...
                        'VariableNames', ...
                        {'Family-Wise Method', 'Test-Statistic' ,'P-Value', 'Verdict'});

                    self.hypothesis_LABEL = '';
                case 'PAIRWISE'
                    pair_vec = string();

                    all_labels = generateTwoWayComb(self);
                    for cc = 1:n_tests
                        %  The following acquires label pairs
                        [t1, t2] = all_labels{logical(C(cc,:))};
                        pair_vec(cc) = convertCharsToStrings([t1 ' & ' t2]);
                    end

                    self.hypothesis_LABEL = pair_vec';

                case {"INTERACTION", "PRIMARY", "SECONDARY"}
                    my_str = char(lower(self.Hypothesis));
                    my_str(1) = string(upper(my_str(1)));
                    self.hypothesis_LABEL = my_str;

                    self.TwoWay_P_Table =  table(self.ANOVA_Methods_Used',  nan(n_methods, 1), nan(n_methods, 1), strings(n_methods, 1),  ...
                        'VariableNames', ...
                        {sprintf('%s Effect', my_str), 'Test-Statistic' ,'P-Value', 'Verdict'});

                otherwise % Custom Hypothesis

                    switch self.Contrast_Factor
                        case 1
                            mask_combine_Subset1 = strjoin(self.PrimaryLabels(self.Contrast >= 1), '+');
                            mask_combine_Subset2 = strjoin(self.PrimaryLabels(self.Contrast <= -1), '+');

                        case 2
                            mask_combine_Subset1 = strjoin(self.SecondaryLabels(self.Contrast >= 1), '+');
                            mask_combine_Subset2 = strjoin(self.SecondaryLabels(self.Contrast <= -1), '+');
                    end

                    self.hypothesis_LABEL = mask_combine_Subset1 + " vs. " + mask_combine_Subset2;
                    self.TwoWay_P_Table =  table(self.ANOVA_Methods_Used',  nan(n_methods, 1), nan(n_methods, 1), strings(n_methods, 1),  ...
                        'VariableNames', ...
                        {sprintf('%s Effect', self.hypothesis_LABEL), 'Test-Statistic' ,'P-Value', 'Verdict'});
            end

            p_value_matrix = nan(n_tests, n_methods);
            test_stat = nan(1, n_methods);

            counter = 0;
            c = 0;  % Assuming Equality to each other and not a generic constant

            for method = self.ANOVA_Methods_Used   % Generate Test Statistic
                counter = counter + 1;

                p_value = zeros(n_tests, 1);
                for ii = 1:n_tests
                    switch self.Hypothesis
                        case 'PAIRWISE'
                            C_input = C(ii, :);
                            self.hypothesis_LABEL = pair_vec(ii);
                        case {"INTERACTION", "PRIMARY", "SECONDARY"}
                            C_input = C;
                        otherwise
                            C_input = C;
                    end
                    [p_value(ii), statistic] = self.TwoWay(method, yy, C_input);
                end
                p_value_matrix(:, counter) = p_value;
                test_stat(counter) = statistic;

            end

            NotNan_mask = ~isnan(p_value_matrix);
            NotNan_mask = all(NotNan_mask, 1);

            %             if sum(NotNan_mask) ~= length(self.ANOVA_Methods_Used)
            %                 warning('Removing ANOVA results from methods that arent implemented yet.')
            %             end

            switch self.Hypothesis
                case "PAIRWISE"
                    T_hypothesis = table(pair_vec', 'VariableNames', {'Hypothesis'});
                    T_p_value = array2table(p_value_matrix(:, NotNan_mask), 'VariableNames', cellstr(self.ANOVA_Methods_Used(:, NotNan_mask)));
                    self.TwoWay_P_Table = [T_hypothesis, T_p_value];
                otherwise
                    self.TwoWay_P_Table = self.TwoWay_P_Table(NotNan_mask', :);
                    self.TwoWay_P_Table.("P-Value") = p_value_matrix(:, NotNan_mask)';
                    self.TwoWay_P_Table.('Test-Statistic') = test_stat(:, NotNan_mask)';
                    Signif_results = self.TwoWay_P_Table.("P-Value") < self.alpha;
                    self.TwoWay_P_Table{Signif_results, "Verdict"} = "Reject Null Hypothesis for Alternative Hypothesis";
                    self.TwoWay_P_Table{~Signif_results, "Verdict"} = "Fail to Reject Null Hypothesis";
            end

            if self.verbose
                self.DataSummaryReportTwoWay("homoskedastic")
                self.showTable(self.TwoWay_P_Table)
            end
        end

        % TwoWay Heteroscedastic (Different Covariances)
        function TwoWayANOVA_BF(self, varargin)
            self.Contrast = []; % Reset to Empty at start

            p = inputParser;
            % Keep ability to change inputs
            addParameter(p, 'SubgroupIndicator', self.SubgroupIndicator, @(x) isnumeric(x)) % Indicator Array Denoting Subgroups
            addParameter(p, 'N_boot', self.N_boot, @(x) isinteger(int64(x)));
            addParameter(p, 'N_simul', self.N_simul, @(x) isinteger(int64(x)))
            addParameter(p, 'ANOVA_Methods', self.ANOVA_Methods, @(x) iscellstr(x) || isstring(x))
            addParameter(p, 'alpha', self.alpha, @(x) x<=1 && x>= 0)
            addParameter(p, 'Contrast', self.Contrast, @(x) isvector(x) || isempty(x))
            addParameter(p, 'PrimaryLabels', self.PrimaryLabels, @(x) iscellstr(x) || isstring(x))
            addParameter(p, 'SecondaryLabels', self.SecondaryLabels, @(x) iscellstr(x) || isstring(x))
            addParameter(p, 'Weights', self.Weights, @(x)  (iscellstr(x) || isstring(x) || ischar) && any(strcmpi(x, ["uniform", "proportional"])))  % weighting style
            addParameter(p, 'Hypothesis', self.Hypothesis, @(x)  (iscellstr(x) || isstring(x) || ischar(x)) && any(strcmpi(x, self.H0_TwoWay)))  % Hypothesis of interest

            parse(p, varargin{:});
            self.SubgroupIndicator = p.Results.SubgroupIndicator;
            self.N_boot = int64(p.Results.N_boot);
            self.N_simul = int64(p.Results.N_simul);
            self.ANOVA_Methods_Used = p.Results.ANOVA_Methods;
            self.alpha = p.Results.alpha;
            self.PrimaryLabels = p.Results.PrimaryLabels;
            self.SecondaryLabels = p.Results.SecondaryLabels;
            self.Weights = upper(p.Results.Weights);  % Cast to Upper
            self.Hypothesis = upper(p.Results.Hypothesis);  % Cast to Upper
            self.ANOVA_Methods_Used = sort(self.ANOVA_Methods_Used, "descend");  % Maintain L2 first than F-test
            self.castANOVAMethods()

            yy = [];
            for K = 1:self.k_groups
                yy = [yy; self.data{K}'];
            end

            if isempty(self.SubgroupIndicator)
                error('SubgroupIndicator must be provided for TwoWay ANOVAs')
            else
                self.setUpTwoWay()  % Creates Indicator Matrices and Labels
                self.setUpTwoWayHypothesis()  % For Custom Hypothesis
            end

            % Define C and number of tests
            switch self.Hypothesis
                case "PAIRWISE"
                    C = self.construct_pairwise_contrast_matrix(self.AB_groups);
                    n_tests = size(C,1); % AB_groups choose 2
                case "FAMILY"
                    n_tests = 1;
                    C = [eye(self.AB_groups-1),-ones(self.AB_groups-1,1)];
                case {"INTERACTION", "PRIMARY", "SECONDARY"}
                    n_tests = 1;
                    C = [];
                otherwise % Custom Hypothesis
                    n_tests = 1;
                    C = self.Contrast;
            end

            %% Methods
            n_methods = numel(self.ANOVA_Methods_Used);

            switch self.Hypothesis
                case 'FAMILY'
                    self.TwoWay_BF_P_Table =  table(self.ANOVA_Methods_Used',  nan(n_methods, 1), nan(n_methods, 1), strings(n_methods, 1),  ...
                        'VariableNames', ...
                        {'Family-Wise Method', 'Test-Statistic' ,'P-Value', 'Verdict'});

                    self.hypothesis_LABEL = '';
                case 'PAIRWISE'
                    pair_vec = string();

                    all_labels = generateTwoWayComb(self);
                    for cc = 1:n_tests
                        %  The following acquires label pairs
                        [t1, t2] = all_labels{logical(C(cc,:))};
                        pair_vec(cc) = convertCharsToStrings([t1 ' & ' t2]);
                    end

                    self.hypothesis_LABEL = pair_vec';

                case {"INTERACTION", "PRIMARY", "SECONDARY"}
                    my_str = char(lower(self.Hypothesis));
                    my_str(1) = string(upper(my_str(1)));
                    self.hypothesis_LABEL = my_str;

                    self.TwoWay_BF_P_Table =  table(self.ANOVA_Methods_Used',  nan(n_methods, 1), nan(n_methods, 1), strings(n_methods, 1),  ...
                        'VariableNames', ...
                        {sprintf('%s Effect', my_str), 'Test-Statistic' ,'P-Value', 'Verdict'});

                otherwise % Custom Hypothesis

                    switch self.Contrast_Factor
                        case 1
                            mask_combine_Subset1 = strjoin(self.PrimaryLabels(self.Contrast == 1), '+');
                            mask_combine_Subset2 = strjoin(self.PrimaryLabels(self.Contrast == -1), '+');

                        case 2
                            mask_combine_Subset1 = strjoin(self.SecondaryLabels(self.Contrast == 1), '+');
                            mask_combine_Subset2 = strjoin(self.SecondaryLabels(self.Contrast == -1), '+');
                    end

                    self.hypothesis_LABEL = mask_combine_Subset1 + " vs. " + mask_combine_Subset2;
                    self.TwoWay_P_Table =  table(self.ANOVA_Methods_Used',  nan(n_methods, 1), nan(n_methods, 1), strings(n_methods, 1),  ...
                        'VariableNames', ...
                        {sprintf('%s Effect', self.hypothesis_LABEL), 'Test-Statistic' ,'P-Value', 'Verdict'});
            end

            p_value_matrix = nan(n_tests, n_methods);
            test_stat = nan(1, n_methods);

            counter = 0;
            c = 0;  % Assuming Equality to each other and not a generic constant

            for method = self.ANOVA_Methods_Used   % Generate Test Statistic
                counter = counter + 1;

                p_value = zeros(n_tests, 1);
                for ii = 1:n_tests
                    switch self.Hypothesis
                        case 'PAIRWISE'
                            C_input = C(ii, :);
                            self.hypothesis_LABEL = pair_vec(ii);
                        case {"INTERACTION", "PRIMARY", "SECONDARY"}
                            C_input = C;
                        otherwise
                            C_input = C;
                    end
                    [p_value(ii), statistic] = self.TwoWay_BF(method, yy, C_input, c);
                end
                p_value_matrix(:, counter) = p_value;
                test_stat(counter) = statistic;

            end

            NotNan_mask = ~isnan(p_value_matrix);
            NotNan_mask = all(NotNan_mask, 1);

            if sum(NotNan_mask) ~= length(self.ANOVA_Methods_Used)
                if self.verbose
                    warning('Removing ANOVA results from methods that arent implemented yet.')
                end
            end

            switch self.Hypothesis
                case "PAIRWISE"
                    T_hypothesis = table(pair_vec', 'VariableNames', {'Hypothesis'});
                    T_p_value = array2table(p_value_matrix(:, NotNan_mask), 'VariableNames', cellstr(self.ANOVA_Methods_Used(:, NotNan_mask)));
                    self.TwoWay_BF_P_Table = [T_hypothesis, T_p_value];
                otherwise
                    self.TwoWay_BF_P_Table = self.TwoWay_BF_P_Table(NotNan_mask', :);
                    self.TwoWay_BF_P_Table.("P-Value") = p_value_matrix(:, NotNan_mask)';
                    self.TwoWay_BF_P_Table.('Test-Statistic') = test_stat(:, NotNan_mask)';
                    Signif_results = self.TwoWay_BF_P_Table.("P-Value") < self.alpha;
                    self.TwoWay_BF_P_Table{Signif_results, "Verdict"} = "Reject Null Hypothesis for Alternative Hypothesis";
                    self.TwoWay_BF_P_Table{~Signif_results, "Verdict"} = "Fail to Reject Null Hypothesis";
            end

            if self.verbose
                self.DataSummaryReportTwoWay("heteroskedastic")
                self.showTable(self.TwoWay_BF_P_Table)
            end
        end

        % K-Sample/Group Check for Equality of Covariances
        function CovarianceTest(self, varargin)
            p = inputParser;
            addParameter(p, 'N_permutations', self.N_permutations,  @(x) isinteger(int64(x)))
            addParameter(p, 'N_boot', self.N_boot, @(x) isinteger(int64(x)));
            addParameter(p, 'COVAR_Methods', self.COVAR_Methods)

            parse(p, varargin{:});
            self.N_permutations = p.Results.N_permutations;
            self.N_boot = p.Results.N_boot;
            self.COVAR_Methods_Used = p.Results.COVAR_Methods;

            % The goal of the equal covariance test is to determine if all
            % groups have equal covariance. Always a family-wise test.
            % Subsetting should be done outside of this function, that is not
            % the goal of this function
            % Uses domain information for subsetting only

            % Compute basic group means and covariances


            %% Set up Data Matrix
            vmu=[]; V=[];
            for ii=1:self.k_groups
                yyi=self.data{ii}';
                mui=mean(yyi);
                Vi=yyi-ones(self.n_i(ii),1)*mui;
                vmu=[vmu;mui];
                V=[V;Vi];
            end

            m = self.n_domain_points;

            if self.N > m
                Sigma=V'*V/(self.N-self.k_groups);  %% [m x m] pooled covariance matrix
            else
                Sigma=V*V'/(self.N-self.k_groups);  %% [N x N] pooled covariance matrix
            end

            %% Computing the test statistic
            stat=0;  % T_n test statistic
            nni=0;
            for ii=1:self.k_groups
                ni=self.n_i(ii);
                flag=(nni+1):(nni+ni);
                Vi=V(flag,:);
                if self.N>m
                    Si=Vi'*Vi/(ni-1);  %% Vi: nixp
                    temp=trace((Si-Sigma)^2);
                else
                    Si=Vi*Vi'/(ni-1);
                    temp=trace(Si^2)-2*trace(Vi*V'*V*Vi')/(self.N-self.k_groups)/(ni-1)+trace(Sigma^2);
                end

                stat=stat+(ni-1)*temp;
                nni=nni+ni;
            end

            n_methods = numel(self.COVAR_Methods_Used);
            p_values = nan(n_methods, 1);


            %% Start Methods
            for K = 1 : numel(self.COVAR_Methods_Used)
                method = self.COVAR_Methods_Used(K);
                p_values(K) = self.k_group_cov(method, stat, Sigma, V);
            end

            NotNan_mask = ~isnan(p_values);
            NotNan_mask = all(NotNan_mask, 2);

            if sum(NotNan_mask) ~= length(self.COVAR_Methods_Used)
                if self.verbose
                    warning('Removing equality of covariance results for methods that arent implemented yet.')
                end
            end

            self.COVAR_P_Table = table(self.COVAR_Methods_Used', p_values, 'VariableNames', {'Method', 'P-Value'});
            self.COVAR_P_Table = self.COVAR_P_Table(NotNan_mask, :);
            Signif_results = self.COVAR_P_Table.("P-Value") < self.alpha;
            self.COVAR_P_Table{:, "Verdict"} = strings(height(self.COVAR_P_Table), 1); % Place Holders
            self.COVAR_P_Table{Signif_results, "Verdict"} = "Reject Null Hypothesis for Alternative Hypothesis";
            self.COVAR_P_Table{~Signif_results, "Verdict"} = "Fail to Reject Null Hypothesis";
            %             self.COVAR_P_Table{isnan(p_values), "Verdict"} = "";

            if self.verbose
                disp(self.COVAR_P_Table)
            end

        end

        % 2-Sample/Group Check for Equality of Covariances (slightly different statistics)
        function CovarianceTest_TwoSample(self, varargin)
            assert(self.k_groups == 2, 'Only supports 2 Groups')

            p = inputParser;
            addParameter(p, 'N_permutations', self.N_permutations,  @(x) isinteger(int64(x)))
            addParameter(p, 'N_boot', self.N_boot, @(x) isinteger(int64(x)));
            addParameter(p, 'COVAR_Methods', self.COVAR_Methods)

            parse(p, varargin{:});
            self.N_permutations = p.Results.N_permutations;
            self.N_boot = p.Results.N_boot;
            self.COVAR_Methods_Used = p.Results.COVAR_Methods;


            n_methods = numel(self.COVAR_Methods_Used);
            p_values = nan(n_methods, 1);

            %% Start Methods
            for K = 1 : numel(self.COVAR_Methods_Used)
                method = self.COVAR_Methods_Used(K);

                switch method
                    case "L2-Simul"
                        % Simulating the null distribution
                        %                         p_values(K) = 1 - T_NullFitted.cdf(T_n);
                         p_values(K) = self.two_group_cov(method, self.data{1}', self.data{2}');
                    case "L2-Naive"
                        p_values(K) = self.two_group_cov(method, self.data{1}', self.data{2}');
                    case "L2-BiasReduced"
                        p_values(K) = self.two_group_cov(method, self.data{1}', self.data{2}');
                    case "Permutation-Test"
                        p_values(K) = self.two_group_cov(method, self.data{1}', self.data{2}');
                    case "Bootstrap-Test"
                        p_values(K) = self.two_group_cov(method, self.data{1}', self.data{2}');
                end
            end

            self.COVAR_P_Table = table(self.COVAR_Methods_Used', p_values, 'VariableNames', {'Method', 'P-Value'});
            Signif_results = self.COVAR_P_Table.("P-Value") < self.alpha;
            self.COVAR_P_Table{:, "Verdict"} = strings(height(self.COVAR_P_Table), 1); % Place Holders
            self.COVAR_P_Table{Signif_results, "Verdict"} = "Reject Null Hypothesis for Alternative Hypothesis";
            self.COVAR_P_Table{~Signif_results, "Verdict"} = "Fail to Reject Null Hypothesis";
            self.COVAR_P_Table{isnan(p_values), "Verdict"} = "";

            if self.verbose
                disp(self.COVAR_P_Table)
            end


        end

        function PointwiseDistributionCheck(self, checkAt, checkType, varargin)
            %POINTWISEDISTRIBUTIONCHECK Examine the distribution of the residuals.
            % POINTWISEDISTRIBUTIONCHECK(SELF, CHECKAT, CHECKTYPE, …)
            % Currently only supports checking residuals for OneWay ANOVA.
            %
            % <strong>Required Input</strong>
            %
            %       checkAt ([Nx1] Numeric)
            %               - Array of points in the domain to check the
            %                 residuals at.
            %     checkType ([1x1] String)
            %               - Options: "boxplot", "qq", "shairo-wilkes"
            %               - "boxplot" constructs simple boxplots of residuals by group
            %               - "qq" constructs normal QQ plots by group
            %               - "shairo-wilkes" calculates SW test by group
            %   domainUnitsLabel ([1x1] string, default='')
            %                    - Independent variable units
            % responseUnitsLabel ([1x1] string, default='')
            %                    - Dependent variable units

            p = inputParser;
            addRequired(p, 'checkAt', @(x) isvector(x))
            addRequired(p, 'checkType', @(x) any(strcmpi(x, {'boxplot', 'qq', 'shapiro-wilkes'})))
            addParameter(p, 'domainUnitsLabel', self.domainUnitsLabel, @(x) ischar(x) || isstring(x));
            addParameter(p, 'responseUnitsLabel', self.responseUnitsLabel, @(x) ischar(x) || isstring(x));

            parse(p, checkAt, checkType, varargin{:});
            self.domainUnitsLabel = p.Results.domainUnitsLabel;
            self.responseUnitsLabel = p.Results.responseUnitsLabel;

            assert(isempty(self.SubgroupIndicator), 'PointwiseDistributionCheck Currently only supports OneWAY ANOVA Analysis')

            % Compute basic group means and covariances
            eta_i = zeros(self.n_domain_points, self.k_groups); % estimated means
            v_hat_i = cell(1, self.k_groups); %residuals matrix
            for kk = 1:self.k_groups
                eta_i(:, kk) = mean(self.data{kk}, 2); %estimated means
                v_hat_i{kk} = self.data{kk} - eta_i(:, kk); %residuals
            end

            [S, L] = bounds(self.d_grid);

            n_check = numel(checkAt);
            check_ind = nan(n_check, 1);

            for K = 1:n_check
                [~, check_idx] = min(abs(self.d_grid - checkAt(K))); %find index to check
                check_ind(K) = check_idx;
            end

            assert(all(checkAt <= L), 'Requested Value Exceeds the Domain %sDomain Upper Limit = %0.2f',newline, L)
            assert(all(checkAt >= S), 'Requested Value is Less than the Domain %sDomain Lower Limit = %0.2f',newline, S)


            n_checks = numel(check_ind); % number of checks

            res_kk_ii = cell(n_checks, self.k_groups);
            res_ii = cell(n_checks,1);
            sw_p_vals = zeros(n_checks, self.k_groups);
            for ii = 1:n_checks
                for kk = 1:self.k_groups
                    res_kk_ii{ii,kk} = v_hat_i{kk}(check_ind(ii), :);
                    [~, sw_p_vals(ii,kk), ~] = swtest(res_kk_ii{ii, kk}, 0.05);
                end
                res_ii{ii} = [res_kk_ii{ii,:}];
            end


            if self.genericGroupLabels
                display_label = "Group " + self.GroupLabels;
            end

            grps = repelem(display_label, self.n_i);


            switch checkType
                case "boxplot"
                    % Construct Boxplots of residuals
                    figure
                    if isempty(self.responseUnitsLabel)
                        y_label = "Residuals";
                    else
                        y_label = sprintf("Residuals (%s)", self.responseUnitsLabel);
                    end

                    for ii = 1:n_checks
                        subplot(1, n_checks, ii)
                        boxplot(res_ii{ii}, grps)
                        ylim([-max(abs(res_ii{ii})) max(abs(res_ii{ii}))])
                        title(checkAt(ii) + " " + self.domainUnitsLabel)
                        xlabel("Group Label")
                        ylabel(y_label);
                        fontsize(gca, 20, "points")

                    end

                case "qq"
                    % Construct QQ plots figure
                    fig = figure('Name', 'QQ-GRID','units', 'points', 'position', [90, 90, 1000, 800]);
                    counter = 0;
                    for ii = 1:n_checks
                        for kk = 1:self.k_groups
                            counter = counter + 1;
                            subplot(n_checks, self.k_groups, counter)
                            qqplot(res_kk_ii{ii,kk})


                            title(checkAt(ii) + " " + self.domainUnitsLabel + " | " + display_label(kk))


                            fontsize(gca,20, "points")
                            xlabel('')
                            ylabel('')
                        end
                    end
                    % Give common xlabel, ylabel and title to your figure
                    han=axes(fig,'visible','off');
                    han.Title.Visible='on';
                    han.XLabel.Visible='on';
                    han.YLabel.Visible='on';
                    xlabel(han,{'', 'Standard Normal Quantiles'}, 'FontSize', 24);
                    y_label = 'Quantiles of Input Data';
                    ylabel(han,{y_label, ''}, 'FontSize', 24);
                    %                     title(han,'yourTitle');


                case "shapiro-wilkes"
                    self.SW_P_Table = array2table(sw_p_vals, 'VariableNames', display_label, 'RowNames', string(checkAt) + " " + self.domainUnitsLabel);

                    if self.verbose
                        disp(self.SW_P_Table)
                    end
            end
        end

        % Plot Means for OneWay or TwoWay Data set
        PlotMeans(self, varargin)

        % Plot Means for OneWay... only
        PlotCovariances(self, varargin)

    end

    methods(Hidden=true)  % In-Progress
        function VerifyOmegaCalculation(self)
            assert(self.k_groups == 2, 'Verification on 2 groups')

            %% Set up Data Matrix
            vmu=[]; V=[];
            for ii=1:self.k_groups
                yyi=self.data{ii}';
                mui=mean(yyi);
                Vi=yyi-ones(self.n_i(ii),1)*mui;
                vmu=[vmu;mui];
                V=[V;Vi];
            end

            p = self.n_domain_points;

            if self.N>p
                Sigma=V'*V/(self.N-self.k_groups);  %% [p x p] pooled covariance matrix
            else
                Sigma=V*V'/(self.N-self.k_groups);
            end

            %% Calculate Omega
            % STILL NOT CALCULATING CORRECTLY
            % Code for constructing omega, the covariance of the covariance
            % ii is the ith group
            % jj is jth observation within the ith group
            %
            n = self.N;
            n_ii = self.n_i;
            eta_i = zeros(self.n_domain_points, self.k_groups); % estimated means
            v_hat_i = cell(self.k_groups, 1); %residuals matrix
            gamma_hat_i = zeros(self.n_domain_points, self.n_domain_points, self.k_groups); % groupwise covariance
            for kk = 1:self.k_groups
                eta_i(:, kk) = mean(self.data{kk}, 2); %estimated means
                zeroMeanData_k_subset = self.data{kk} - eta_i(:, kk); %residuals
                v_hat_i{kk} = zeroMeanData_k_subset; %residuals
                gamma_hat_i(:,:,kk) = 1 / (self.n_i(kk) - 1) .* (zeroMeanData_k_subset * zeroMeanData_k_subset'); % groupwise covariance
            end

            % Compute pooled covariance
            pooled_covar_terms = zeros(self.n_domain_points, self.n_domain_points, self.k_groups);
            for kk = 1:self.k_groups
                pooled_covar_terms(:,:,kk) = (self.n_i(kk) - 1) .* gamma_hat_i(:,:,kk);
            end
            pooled_covar = sum(pooled_covar_terms, 3) ./ (n - self.k_groups);

            RHS  = (pooled_covar * pooled_covar); % RHS of minus sign in Eq 10.29  %trace(RHS) == trace(Sigma * Sigma)
            LHS  = zeros(size(RHS)); % LHS of minus sign in Eq 10.29  % 3-Dimensions, parallize is fastest on the inner loop. Slower than serial on outerloop

            T = TimedProgressBar(n, 35, ...
                'Calculating Omega for Equality of Covariance in: ', ...
                'Finished Calculating Omega in: ');

            % Summation portion of Eq 10.29
            for ii = 1:self.k_groups
                temp_hat = v_hat_i{ii};
                temp_LHS = zeros(size(RHS));

                parfor jj = 1:n_ii(ii)
                    % v_ij = v_hat_i{ii}(:, jj);
                    v_ij =  temp_hat(:, jj);
                    temp = (v_ij * v_ij') * (v_ij * v_ij');
                    temp_LHS = temp_LHS + temp;
                    T.progress;
                end
                LHS = LHS + temp_LHS;
            end

            %                         LHS = sum(LHS, 3); % Sum over all Groups
            LHS = LHS ./ (n);  % Coefficient outside summation in Eq 10.29

            self.omega_hat = LHS - RHS;
            T.stop;T.delete;
            trace_omega_sum = trace(Sigma*Sigma + Sigma*Sigma);
            fprintf('%sTrace(Omega) = %0.2f : Ground Truth using Equality from Equation 10.15', newline, trace_omega_sum)
            fprintf('%sTrace(Omega) = %0.2f : Estimated Using Equation 10.7', newline, trace(self.omega_hat))
        end
    end

    methods (Access = private, Hidden=true)

        % Runs K-group Covariance Tests
        pvalue = k_group_cov(self, method, stat, Sigma, V)

        % Runs 2-group Covariance Tests (slightly different test statistics
        % than the generalized K-group)
        pvalue = two_group_cov(self, method, y1, y2)

        % Runs OneWay homogeneous tests
        p_value_matrix =       OneWay(self, n_tests, q, eig_gamma_hat, eta_i, eta_grand, params, pair_vec)

        % Runs OneWay Heteroscedastic tests
        [p_value, statistic] = OneWay_BF(self, method, data, Contrast, c, varargin)

        % Runs TwoWay homogeneous tests
        [pvalue, stat]  =      TwoWay(self, method, data, Contrast)

        % Runs TwoWay Heteroscedastic tests
        [pvalue, stat]  =      TwoWay_BF(self, method, data, Contrast,c)

        % Subsets Data into smaller Domain
        self = function_Subsetter(self)

        % Helper Function for OneWay Homoscedastic F-ANOVA F-type bootstrap
        function [p_value, stat, btstat]= FBootstrap(self, yy)
            % F-type test for k-sample problem with a common covariance function.
            % One-way ANOVA for functional data: Main-effects test
            % yy=Matrix [y1,y2,...,yp]: nxp data matrix,  each row: discretization of a func.
            % stat=test statistic, btstat= a bootstrap sample of statistic
            % Jin-Ting Zhang
            % March 13, 2008; 
            % Revised July 11, 2008  NUS, Singapore
            % Revised July 31, 2024 Los Alamos National Laboratory

            [n,p]=size(yy);
            % p=p-1;

            gsize = self.n_i';
            aflag = functionalANOVA.aflagMaker(gsize);
            % aflag=yy(:,1);
            aflag0=unique(aflag);
            k=length(aflag0); %% Level of Factor A
            % yy=yy(:,2:(p+1)); %% nxp data matrix,  each row: discretization of a func.

            mu0=mean(yy); %% pooled sample mean function
            gsize=[];
            vmu=[];
            z=[];
            SSR=0;

            for i=1:k
                iflag=(aflag==aflag0(i));
                yi=yy(iflag,:); %%Samle i
                ni=size(yi,1);
                mui=mean(yi);
                zi=yi-ones(ni,1)*mui;
                gsize=[gsize;ni];
                vmu=[vmu;mui]; %% each row is a group mean vector
                z=[z;zi];
                SSR=SSR+ni*(mui-mu0).^2; %% 1xp vector
            end

            if n>p
                Sigma=z'*z/(n-k);  %% pxp pooled covar. matrix
            else
                Sigma=z*z'/(n-k); %% nxn matrix, having the same eigenvalues with pooled cov. matrix
            end

            A=trace(Sigma);
            % B=trace(Sigma^2);
            stat=sum(SSR)/A/(k-1);

            btstat = zeros(1, self.N_boot);
            parfor ii=1:self.N_boot
                btvmu=[];btz=[];
                for i=1:k
                    iflag=(aflag==aflag0(i));
                    yi=yy(iflag,:);
                    ni=gsize(i);
                    % btflag=fix(rand(ni,1)*(ni-1))+1;
                    btflag = randsample(ni, ni, true); 
                    
                    btyi=yi(btflag,:);
                    btmui=mean(btyi);
                    btzi=btyi-ones(ni,1)*btmui;btz=[btz;btzi];
                    btmui=btmui-vmu(i,:);
                    btvmu=[btvmu;btmui];
                end
                btmu0=gsize'*btvmu/n;
                btSSR=0;
                for i=1:k
                    btSSR=btSSR+gsize(i)*(btvmu(i,:)-btmu0).^2;
                end

                if n>p
                    btSigma=btz'*btz/(n-k);  %% pxp pooled covar. matrix
                else
                    btSigma=btz*btz'/(n-k); %% nxn matrix, having the same eigenvalues with pooled cov. matrix
                end
                btA=trace(btSigma);
                btstat(ii)=sum(btSSR)/btA/(k-1);

            end
            p_value=mean(btstat>=stat);         
        end

        % Helper Function for OneWay Homoscedastic F-ANOVA L2-norm bootstrap
        function [p_value, stat, btstat]= L2Bootstrap(self, yy)
            % L2-norm based  test for k-sample problem with a common covariance function.
            % One-way ANOVA for functional data: Main-effects test
            % yy=Matrix [y1,y2,...,yp]: nxp data matrix,  each row: discretization of a func.
            % stat=test statistic, btstat= a bootstrap sample of statistic
            % Jin-Ting Zhang
            % March 13, 2008;
            % Revised July 11, 2008  NUS, Singapore
            % Revised Oct 19, 2011 Princeton University
            % Revised July 31, 2024 Los Alamos National Laboratory


            [n,p]=size(yy);
            % p=p-1;

            gsize = self.n_i';
            aflag = functionalANOVA.aflagMaker(gsize);
            % aflag=yy(:,1);
            aflag0=unique(aflag);
            k=length(aflag0); %% Level of Factor A
            % yy=yy(:,2:(p+1)); %% nxp data matrix,  each row: discretization of a func.

            mu0=mean(yy); %% pooled sample mean function
            gsize=[];
            vmu=[];
            z=[];
            SSR=0;

            for i=1:k
                iflag=(aflag==aflag0(i));
                yi=yy(iflag,:); %%Samle i
                ni=size(yi,1);
                mui=mean(yi);
                zi=yi-ones(ni,1)*mui;
                gsize=[gsize;ni];
                vmu=[vmu;mui]; %% each row is a group mean vector
                z=[z;zi];
                SSR=SSR+ni*(mui-mu0).^2; %% 1xp vector
            end

            stat=sum(SSR);

            btstat = zeros(1, self.N_boot);
            parfor ii=1:self.N_boot
                btvmu=[];
                for i=1:k
                    iflag=(aflag==aflag0(i));
                    yi=yy(iflag,:);
                    ni=gsize(i);
                    % btflag=fix(rand(ni,1)*(ni-1))+1;
                    btflag = randsample(ni, ni, true);
                    btyi=yi(btflag,:);
                    btmui=mean(btyi)-vmu(i,:);
                    btvmu=[btvmu;btmui];
                end
                btmu0=gsize'*btvmu/n;
                btSSR=0;
                for i=1:k
                    btSSR=btSSR+gsize(i)*(btvmu(i,:)-btmu0).^2;
                end
                btstat(ii)=sum(btSSR);
            end
            p_value=mean(btstat>=stat);
  
            
        end

        % Helper Function to set up OneWay Homoscedastic F-ANOVA
        function [T_n, F_n, beta_hat, kappa_hat, beta_hat_unbias, kappa_hat_unbias] = SetUpANOVA(self, SSH_t, SSE_t, gamma_hat, q)

            FamilyFlag = strcmpi(self.Hypothesis, 'FAMILY');

            F_n = 0;
            beta_hat = 0;
            kappa_hat = 0;
            beta_hat_unbias = 0;
            kappa_hat_unbias = 0;

            T_n = sum(SSH_t, 1); %dimension may be wrong

            if FamilyFlag
                if any(contains(self.ANOVA_Methods_Used, 'L2'))
                    mask  = contains(self.ANOVA_Methods_Used, 'L2');
                    self.OneWay_P_Table.("Test-Statistic")(mask) = T_n;
                end
            end


            if any(contains(self.ANOVA_Methods_Used, 'F'))   % Check for F-test Methods
                numerator = T_n ./ q;
                denomenator = sum(SSE_t) / (self.N - self.k_groups);
                F_n = numerator ./ denomenator;

                if FamilyFlag
                    mask  = contains(self.ANOVA_Methods_Used, 'F');
                    self.OneWay_P_Table.("Test-Statistic")(mask) = F_n;
                end
            end

            if any(contains(self.ANOVA_Methods_Used, 'Naive'))  % Check for Naive Approx
                beta_hat = functionalANOVA.BetaHat(gamma_hat);
                kappa_hat = functionalANOVA.KappaHat(gamma_hat);

                if FamilyFlag
                    mask_1  = contains(self.ANOVA_Methods_Used, 'Naive');
                    mask_2  =  contains(self.ANOVA_Methods_Used, 'F');
                    mask = and(mask_1, ~mask_2);
                    if any(mask)
                        self.OneWay_P_Table.("Parameter 1 Name")(mask) = 'beta';
                        self.OneWay_P_Table.("Parameter 2 Name")(mask) = 'd';
                        self.OneWay_P_Table.("Parameter 1 Value"){mask} = beta_hat;
                        self.OneWay_P_Table.("Parameter 2 Value"){mask} = q * kappa_hat;
                    end

                    mask = and(mask_1, mask_2);
                    if any(mask)  % Unbias and F-test
                        d_1 = q * kappa_hat;
                        d_2 = (self.N - self.k_groups) * kappa_hat;

                        self.OneWay_P_Table.("Parameter 1 Name")(mask) = 'd1';
                        self.OneWay_P_Table.("Parameter 2 Name")(mask) = 'd2';
                        self.OneWay_P_Table.("Parameter 1 Value"){mask} = d_1;
                        self.OneWay_P_Table.("Parameter 2 Value"){mask} = d_2;
                    end
                end

            end

            if any(contains(self.ANOVA_Methods_Used, 'BiasReduced')) % Check for BiasReduced Approx
                beta_hat_unbias = functionalANOVA.BetaHatUnbias(self.N, self.k_groups, gamma_hat);
                kappa_hat_unbias = functionalANOVA.KappaHatUnbias(self.N, self.k_groups, gamma_hat);

                if FamilyFlag
                    mask_1  = contains(self.ANOVA_Methods_Used, 'BiasReduced');
                    mask_2  =  contains(self.ANOVA_Methods_Used, 'F');
                    mask = and(mask_1, ~mask_2);
                    if any(mask)
                        self.OneWay_P_Table.("Parameter 1 Name")(mask) = 'beta';
                        self.OneWay_P_Table.("Parameter 2 Name")(mask) = 'd';
                        self.OneWay_P_Table.("Parameter 1 Value"){mask} = beta_hat_unbias;
                        self.OneWay_P_Table.("Parameter 2 Value"){mask} = q * kappa_hat_unbias;
                    end

                    mask = and(mask_1, mask_2);
                    if any(mask)  % Unbias and F-test
                        d_1 = q * kappa_hat_unbias;
                        d_2 = (self.N - self.k_groups) * kappa_hat_unbias;

                        self.OneWay_P_Table.("Parameter 1 Name")(mask) = 'd1';
                        self.OneWay_P_Table.("Parameter 2 Name")(mask) = 'd2';
                        self.OneWay_P_Table.("Parameter 1 Value"){mask} = d_1;
                        self.OneWay_P_Table.("Parameter 2 Value"){mask} = d_2;
                    end
                end
            end

        end

        % Updates OneWay homogeneous P-value table
        function updateFamilyTable(self, A_method, params)

            switch A_method

                case {"L2-Simul", "F-Simul"}

                    mask = self.OneWay_P_Table{:, 1}== A_method;
                    self.OneWay_P_Table.("Parameter 1 Name")(mask) = 'KDE: Kernel';
                    self.OneWay_P_Table.("Parameter 2 Name")(mask) = 'KDE: BandWidth';
                    self.OneWay_P_Table.("Parameter 1 Value"){mask} = params{1}.Kernel;
                    self.OneWay_P_Table.("Parameter 2 Value"){mask} = params{1}.Bandwidth;

                case {"L2-Bootstrap", "F-Bootstrap"}

                    mask = self.OneWay_P_Table{:, 1}== A_method;
                    self.OneWay_P_Table.("Parameter 1 Name")(mask) = 'Bootstrap: Resamples';
                    self.OneWay_P_Table.("Parameter 2 Name")(mask) = 'Bootstrap: Type';
                    self.OneWay_P_Table.("Parameter 1 Value"){mask} = params{1};
                    self.OneWay_P_Table.("Parameter 2 Value"){mask} = "nonparametric";
            end


        end

        % Creates String for generic data summary for OneWay F-ANOVA, similar to the R library version
        function DataSummaryReportOneWay(self, ANOVA_TYPE)

            n_groups = numel(self.n_i);

            Mystring = sprintf("\nOne-Way %s F-ANOVA Data Summary:\n\n", ANOVA_TYPE);
            Mystring = Mystring + sprintf("Confidence Level = %0.3f %s\n", (1-self.alpha) * 100, '%%');
            Mystring = Mystring  + sprintf("Number of Observations Total = %i\n", sum(self.n_i));
            Mystring = Mystring + sprintf("Number of Points in Domain = %i\n", sum(self.n_domain_points));
            Mystring = Mystring + sprintf("Number of Groups = %i\n", n_groups);
            Mystring = Mystring + sprintf("Domain Range = [%0.3f, %0.3f]\n", self.d_grid(1), self.d_grid(end));
            Mystring = Mystring + sprintf("Domain Subset = [%0.3f, %0.3f]\n", self.boundsArray(1), self.boundsArray(end));
            Mystring = Mystring + sprintf("Group Observation Size: [%s]\n", strjoin(string(self.n_i), ', '));


            if ~isempty(self.GroupLabels)
                Mystring = Mystring + sprintf("Group Labels: [%s]\n", strjoin(self.GroupLabels, ', '));
            end

            Mystring = Mystring + sprintf('%s', newline);
            fprintf(Mystring)

        end

        % A function for constructing all PAIRWISE contrasts coefficients
        function C = construct_pairwise_contrast_matrix(~, total_groups)
            k = total_groups;

            if k == 2
                C = [eye(total_groups-1),-ones(total_groups-1,1)];
            else
                C = {1,k-1};
                for mm = 1:(k-1)
                    C{mm} = zeros(k - mm, k - 1);
                    for cc = 1:(k - mm)
                        C{mm}(cc,mm) = 1;
                        C{mm}(cc,cc+mm) = -1;
                    end
                end
                C = cell2mat(C');

            end
        end

        % Generates All the labels for the interactions between Factors
        function combinations = generateTwoWayComb(self)

            combinations = strings(self.AB_groups, 1);
            counter = 1;
            for K = 1:self.A_groups
                for KK = 1:self.B_groups
                    combinations(counter) = strjoin([string(self.PrimaryLabels(K)),  string(self.SecondaryLabels(KK))], '-' );
                    counter = counter + 1;
                end
            end

        end

        % Creates String for generic data summary for OneWay F-ANOVA, similar to the R library version
        function DataSummaryReportTwoWay(self, ANOVA_TYPE)


            % Set up Nested String for nested observations
            B_obs_string = [];
            for k = 1 :self.A_groups
                B_obs_string = [B_obs_string,  string(sprintf('[%s]', strjoin(string(self.n_ii{k}), ' ')))];
            end
            B_obs_string =  strjoin(B_obs_string, ', ');


            Mystring = sprintf("\nTwo-Way %s F-ANOVA Data Summary:\n\n", ANOVA_TYPE);
            Mystring = Mystring + sprintf("Confidence Level = %0.3f %s\n", (1-self.alpha) * 100, '%%');
            Mystring = Mystring  + sprintf("Number of Observations Total = %i\n", sum(self.n_i));
            Mystring = Mystring + sprintf("Number of Points in Domain = %i\n", sum(self.n_domain_points));
            Mystring = Mystring + sprintf("Number of Groups within Primary Factor = %i\n", self.A_groups);
            Mystring = Mystring + sprintf("Number of Groups within Secondary Factor= %i\n", self.B_groups);
            Mystring = Mystring + sprintf("Number of Total Groups = %i\n", self.AB_groups);
            Mystring = Mystring + sprintf("Domain Range = [%0.3f, %0.3f]\n", self.d_grid(1), self.d_grid(end));
            Mystring = Mystring + sprintf("Domain Subset = [%0.3f, %0.3f]\n", self.boundsArray(1), self.boundsArray(end));
            Mystring = Mystring + sprintf("Primary Factor Observation Size: [%s]\n", strjoin(string(self.n_i), ', '));
            Mystring = Mystring + sprintf("Secondary Factor Observation Size: [%s]\n", B_obs_string);

            if ~isempty(self.PrimaryLabels)
                Mystring = Mystring + sprintf("Primary Factor Labels: [%s]\n", strjoin(self.PrimaryLabels, ', '));
            end

            if ~isempty(self.SecondaryLabels)
                Mystring = Mystring + sprintf("Secondary Factor Labels: [%s]\n", strjoin(self.SecondaryLabels, ', '));
            end

            Mystring = Mystring + sprintf('%s', newline);
            fprintf(Mystring)


        end

        % Turns Numbers into Characters as default secondary labels factors
        function letters = numbersToLetters(~, numbers)
            Alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
            N_nums = numel(numbers);

            letters = strings(0, N_nums);
            for K = 1: N_nums
                letters(K) = Alphabet(numbers(K));
            end
        end

        function delta_y = TwoSampleDeltaCalcultor(T_sig, n, n_1, n_2, M)

            delta_y = sqrt( (T_sig * n) / (n_1 * n_2 * M));
        end

        % Cleans up table results with fewer digits for displaying to Window
        function showTable(self, table_to_show)
            temp_table = table_to_show;

            % Check whether Family or PairWise
            switch self.Hypothesis
                case "PAIRWISE"
                    all_varNames = temp_table.Properties.VariableNames;
                    n_items = height(temp_table);

                    for K = 2 : numel(all_varNames)
                        data_str = strings(n_items, 1);
                        for KK = 1 : n_items
                            data_str(KK) = sprintf('%*.*f',6, self.table_sigFigs, temp_table.(all_varNames{K})(KK,1));
                        end
                        temp_table.(all_varNames{K}) = data_str;
                    end

                    [~,Y] = ismember(all_varNames(2:end), self.ANOVA_Methods);
                    [~, Z] = sort(Y);

                    subTable = temp_table(:, 2:end);
                    temp_table = [temp_table(:,1), subTable(:, Z)];
                otherwise
                    n_items = height(temp_table);

                    data_str = strings(n_items, 1);
                    %round  decimal places
                    for i = 1:n_items
                        data_str(i) = sprintf('%*.*f',6, self.table_sigFigs, temp_table.("Test-Statistic")(i,1));
                    end

                    temp_table.("Test-Statistic") = data_str;

                    data_str = strings(n_items, 1);
                    %round  decimal places
                    for i = 1:n_items
                        data_str(i) = sprintf('%*.*f',6, self.table_sigFigs, temp_table.("P-Value")(i,1));
                    end

                    temp_table.("P-Value") = data_str;

                    % Sort in a Specific Manner

                    [~,Y] = ismember(temp_table{:, 1}, self.ANOVA_Methods);
                    [~, Z] = sort(Y);

                    temp_table = temp_table(Z, :);

            end
            disp(temp_table)
        end

        % Casts ANOVA methods into specific string format
        function castANOVAMethods(self)
            mask =  contains(self.ANOVA_Methods, self.ANOVA_Methods_Used, "IgnoreCase", true);
            self.ANOVA_Methods_Used = self.ANOVA_Methods(mask);

            assert(sum(mask) >=1, 'No ANOVA Methods were selected!%sMust be at least one of the following: %s.', newline, strjoin(self.ANOVA_Methods, ', '))
        end

        % Verifies Indicator Array to be used in TwoWay Analysis
        function verifyIndicator(self)
            % Verify
            if isa(self.SubgroupIndicator, 'cell')
                assert(numel(self.SubgroupIndicator) == self.k_groups, 'Cell array must match the number of Primary Factors')
                self.SubgroupIndicator = cell2mat(self.SubgroupIndicator); % Cast to one numeric Array
                self.SubgroupIndicator = reshape(self.SubgroupIndicator, self.N, 1); % Reshape
            end

            N_subgroup = numel(self.SubgroupIndicator);
            assert(N_subgroup== self.N, 'Number of elements in SubgroupIndicator=%i doesnt match the total number of samples N=%i',N_subgroup, self.N )
        end

        % Sets up TwoWay ANOVA with labels for each level within the factors
        function setUpTwoWay(self)

            self.verifyIndicator()

            % Generate More info about 2-Way
            self.B_groups = numel(unique(self.SubgroupIndicator));
            self.A_groups = self.k_groups;
            self.AB_groups = self.A_groups * self.B_groups;

            if isempty(self.PrimaryLabels) && isempty(self.SecondaryLabels)
                self.genericGroupLabels = true;
            end

            % Provide Labels if empty
            if isempty(self.PrimaryLabels)
                self.PrimaryLabels = string(1:numel(self.data));
            end

            if isempty(self.SecondaryLabels)
                self.SecondaryLabels = self.numbersToLetters(unique(self.SubgroupIndicator));
            end

            % Check that each level label is one-to-one
            assert(numel(self.PrimaryLabels)   == self.A_groups, "Labels for each Primary factor level must have a one-to-one correspondance to each level")
            assert(numel(self.SecondaryLabels) == self.B_groups, "Labels for each Secondary factor level must have a one-to-one correspondance to each level")


            assert(isempty(self.GroupLabels), ...
                'TwoWay ANOVA requires using "PrimaryLabels" and "SecondaryLabels" as input arguments.%sIt doesnt support the "GroupLabels" argument due to ambiguity.', newline)

        end

        % Sets up TwoWay ANOVA for custom Hypothesis Testing
        function setUpTwoWayHypothesis(self)
            % If custom Hypothesis specified and verify size.
            if strcmp(self.Hypothesis, "CUSTOM")
                assert(~isempty(self.Contrast), 'A user-specified Hypothesis requires a user specified, non-empty, Contrast vector')
                assert(sum(self.Contrast)==0,   'Contrast Vector must sum to zero')

                if self.B_groups == self.A_groups

                    assert(numel(self.Contrast) == self.A_groups, ...
                        'Size of Contrast Vector=%i doesnt match the number of levels in either the Primary or Secondary Factors, n=%i.', numel(self.Contrast), self.A_groups)

                    % Prompt contrast answer for Primary or Secondary Terms
                    list = {'Primary', 'Secondary'};
                    [indx,tf] = listdlg('PromptString',{'Select a Factor.'},...
                        'SelectionMode','single','ListString',list, 'Name','Factor Selector', 'ListSize',[250,250]);

                    if tf
                        switch indx
                            case 1
                                assert(numel(self.Contrast) == self.A_groups, ...
                                    'Size of Contrast Vector=%i doesnt match the number of levels for the Primary=%i Factor', numel(self.Contrast), self.A_groups)
                                disp('Contrast Vector Applied to Primary Factor')
                                self.Contrast_Factor = 1;
                            case 2
                                assert(numel(self.Contrast) == self.A_groups, ...
                                    'Size of Contrast Vector=%i doesnt match the number of levels for the Secondary=%i Factor', numel(self.Contrast), self.B_groups)
                                disp('Contrast Vector Applied to Secondary Factor')
                                self.Contrast_Factor = 2;
                        end

                    else
                        exit('No Selection Made')
                    end

                else % Different sizes Check size of contrast
                    if numel(self.Contrast) == self.A_groups
                        disp('Contrast Vector Applied to Primary Factor')
                        self.Contrast_Factor = 1;

                    elseif numel(self.Contrast) == self.B_groups
                        disp('Contrast Vector Applied to Secondary Factor')
                        self.Contrast_Factor = 2;

                    else
                        error('Size of Contrast Vector=%i doesnt match the number of levels in either the Primary=%i or Secondary=%i Factors.', numel(self.Contrast), self.A_groups, self.B_groups)
                    end
                end
            elseif ~isempty(self.Contrast)
                warning('Ignoring Contrast Vector. Contrast Vector only used when Hypothesis="CUSTOM"')
            end
        end

        % Sets up TimeBar display outputs
        function ts = setUpTimeBar(self, method)

            if self.verbose
                switch method

                    case "L2-Bootstrap"
                        switch self.Hypothesis
                            case 'FAMILY'
                                ts = TimedProgressBar(self.N_boot, 35, ...
                                    'Calculating, Family-wise, Bootstrap L2 test in: ', ...
                                    'Finished, Family-wise, Bootstrap L2 test in: ');
                            case "PAIRWISE"
                                ts = TimedProgressBar(self.N_boot, 35, ...
                                    sprintf('Calculating, Pair-wise (%s), Bootstrap L2 test in: ', self.hypothesis_LABEL), ...
                                    sprintf('Finished, Pair-wise (%s), Bootstrap L2 test in: ', self.hypothesis_LABEL));
                            otherwise
                                ts = TimedProgressBar(self.N_boot, 35, ...
                                    sprintf('Calculating, %s Effect, Bootstrap L2 test in: ', self.hypothesis_LABEL), ...
                                    sprintf('Finished, %s Effect, Bootstrap L2 test in: ', self.hypothesis_LABEL));
                        end

                    case "F-Bootstrap"

                        switch self.Hypothesis
                            case 'FAMILY'
                                ts = TimedProgressBar(self.N_boot, 35, ...
                                    'Calculating, Family-wise, Bootstrap F-type test in: ', ...
                                    'Finished, Family-wise, Bootstrap F-type test in: ');
                            case "PAIRWISE"
                                ts = TimedProgressBar(self.N_boot, 35, ...
                                    sprintf('Calculating, Pair-wise (%s), Bootstrap F-type test in: ', self.hypothesis_LABEL), ...
                                    sprintf('Finished, Pair-wise (%s), Bootstrap F-type test in: ', self.hypothesis_LABEL));
                            otherwise
                                ts = TimedProgressBar(self.N_boot, 35, ...
                                    sprintf('Calculating, %s Effect, Bootstrap F-type test in: ', self.hypothesis_LABEL), ...
                                    sprintf('Finished, %s Effect, Bootstrap F-type test in: ', self.hypothesis_LABEL));
                        end

                end
            else
                ts = struct('progress', [], 'stop', [], 'delete', []);
            end


        end

        % Checks to see if newLabels are passed in
        function checkForNewLabels(self, p)
            if isempty(self.SubgroupIndicator)  % OneWay
                if ~all(contains(self.GroupLabels, p.Results.GroupLabels))
                    self.genericGroupLabels = false;
                end
            else % TwoWay
               if ~all(contains( self.PrimaryLabels, p.Results.GroupLabels))
                   self.genericGroupLabels = false;
                   return
               end

               if ~all(contains(self.SecondaryLabels, p.Results.GroupLabels))
                   self.genericGroupLabels = false;
               end
            end
        end

    end

    methods (Static, Access = private,  Hidden=true)

        % Plots Null Distribution Results for One-Way Homogeneous Simulation 
        PlotTestStats(p_value, alpha, NullDist, testStat, testName, Hypothesis, HypothesisLabel)

        % Function used for BootStrapping
        function [eta_i_star, eta_grand_star, gamma_hat_star] = GroupBooter(DataMatrixCell, n_domain_points, k_groups, n_i, n)

            eta_i_star = zeros(n_domain_points, k_groups);
            build_Covar_star = zeros(n_domain_points, 0);

            for k = 1:k_groups
                data_subset_boot = datasample(DataMatrixCell{k}, n_i(k), 2, "Replace", true);
                eta_i_star(:, k) = mean(data_subset_boot, 2);

                zeroMeanData_k_subset = data_subset_boot -  eta_i_star(:, k);
                build_Covar_star = [build_Covar_star, zeroMeanData_k_subset];
            end

            gamma_hat_star = (1 / (n-k_groups)) * (build_Covar_star' * build_Covar_star);
            eta_grand_star = sum(eta_i_star .* n_i, 2) / n;


        end

        % SSH_n Bootstrap Test Statistic
        function T_n_star = SSH_n_boot(eta_i, eta_grand, eta_i_star, eta_grand_star, n_i)

            SSH_n = n_i .* ((eta_i_star - eta_grand_star) - (eta_i - eta_grand)).^2;
            SSH_n = sum(SSH_n, 2);
            T_n_star = sum(SSH_n);

        end


        %Behrens-Fischer Booter
        function [T_n_star, F_n_star] = BF_boot(method, eta_i, DataMatrixCell, n_i, N, n_domain_points)
            n1 = n_i(1);
            n2 = n_i(2);
            eta_i_star = zeros(n_domain_points, 2);
            gamma_i = cell(1,2);

            for k = 1:2
                data_subset_boot = datasample(DataMatrixCell{k}, n_i(k), 2, "Replace", true);
                eta_i_star(:, k) = mean(data_subset_boot, 2);
                zeroMeanData_k_subset = data_subset_boot -  eta_i_star(:, k);
                gamma_i{k} = (n_i(k)-1)^-1 .* (zeroMeanData_k_subset * zeroMeanData_k_subset');
            end

            % Fixed Synax Mistake as it wasnt using the entire column vector
            eta_delta = eta_i(:, 1) - eta_i(:, 2);
            eta_delta_star = eta_i_star(:, 1) - eta_i_star(:, 2);
            delta_N_star = sqrt(n1*n2 / N).*(eta_delta_star - eta_delta);
            T_n_star = sum(delta_N_star.^2);

            switch method
                case "F-Bootstrap"
                    gamma_deltaN = (n2 / N) .* gamma_i{1} + (n1 / N) .* gamma_i{2}; % gamma_deltaN
                    F_n_star = T_n_star / trace(gamma_deltaN);
                otherwise
                    F_n_star = [];
            end
        end

        % Generate L2 Null Distribution
        function T_null = Chi_sq_mixture(df, coefficients, N_samples)
            n_eigs = length(coefficients); %n eigvals
            chi2rvs = chi2rnd(df, n_eigs, N_samples); %get chi2 rand vars w/ q dof
            T_null = (coefficients' * chi2rvs)'; %compute T null dist
            % gamma_hat = (1 / (n-k_groups)) * (COV * COV');
            % T_null = zeros(N_size, 1);
            % parfor K = 1:N_size
            %     T_null(K) = sum(lambda .* chi2rnd(df, n_eigs, 1));
            % end

        end

        % Calculates Parameters for estimating Null Distribution
        function traceSquared = unbiasedEstimator_traceSquared(n, k, COV)
            n_adj = ((n - k) * (n - k + 1)) / ( (n-k-1) * (n-k+2) );

            factor = trace(COV)^2 - ( (2 * trace(COV^2) )  / (n-k+1) );

            % e = eig(COV);
            % factor = sum(e)^2 - ( (2 * sum(e.^2) ) / n);

            traceSquared = n_adj * factor;
        end

        function traceCovSquared = unbiasedEstimator_traceCovSquared(n, k, COV)
            n_adj = ((n-k)^2) / ((n-k-1)*(n-k+2));

            % e = eig(COV);
            % factor = sum(e.^2) - ((sum(e)^2 / (n-1)));

            factor = trace(COV^2) - ( trace(COV)^2 / (n-k));

            traceCovSquared =  n_adj * factor;
        end

        function beta_hat = BetaHat(COV)
            %%% Slower O(n^3)
            % e = eig(COV);
            % beta_hat = sum(e.^2) / sum(e);

            % Faster O(n^2)
            beta_hat = trace(COV^2) / trace(COV);
        end

        function beta_hat_unbias = BetaHatUnbias(n, k, COV)
            traceCovSquared = functionalANOVA.unbiasedEstimator_traceCovSquared(n, k, COV);
            beta_hat_unbias = traceCovSquared / trace(COV);
        end

        function kappa_hat = KappaHat(COV)
            % e = eig(COV);
            % d_hat = sum(e)^2 / sum(e.^2);
            kappa_hat = trace(COV)^2 / trace(COV^2);
        end

        function kappa_hat_unbias = KappaHatUnbias(n, k, COV)

            kappa_hat_unbias = functionalANOVA.unbiasedEstimator_traceSquared(n, k, COV) / functionalANOVA.unbiasedEstimator_traceCovSquared(n, k, COV);

        end

        % Fuction to generate Indicator array
        function aflag  = aflagMaker(n_i)

            aflag=[];

            for K = 1 :numel(n_i)
                idicator = repmat(K, n_i(K), 1);
                aflag = [aflag; idicator];
            end

        end



    end

    methods (Static)

        function F_ANOVA_Obj = EchoWrapper(R, ANOVA_Labels, boundsArray, varargin)
            % EchoWrapper Returns F-ANOVA Object using Echo Labels
            % EchoWrapper(R, ANOVA_Labels, boundsArray, varargin)
            %
            % Echo Wrapper using the label names in the record array to help
            % group the records for either one-way or two-way ANOVA.
            % Automatically creates ensembles and correct grouping labels.
            %
            % <strong>Required Inputs</strong>
            %                    R (OneD Record array)
            %                      - Record array of OneD Echo Records
            %         ANOVA_Labels ([1xm] String or Cell String Array)
            %                      - Echo label names for grouping records.
            %                      - if m=1 then it prepares a OneWay Analysis
            %                      - if m=2 then it prepares a Two-Way Analysis.
            %                        The first label name is used to generate the
            %                        primary factor lavels and the second label
            %                        name is used to create the secondary factor levels.
            %                      - Only supports 1 or 2 label names.
            %          boundsArray ([1x2] Numeric)
            %                      - Bounds to subset the funtional data response by
            %                        the values from d_grid.
            %                      - To utilize all the functional data, it is
            %                        recommended to use boundsArray=[-inf, inf]
            %
            % <strong>Optional Inputs</strong>
            %               N_boot ([1x1] Numeric, default=10,000)
            %                      - Number of bootstrap replicate.
            %                      - Used for any methods that have "Bootstrap" in them.
            %              N_simul ([1x1] Numeric, default=10,000)
            %                      - Number of data points to simulate null distribution
            %                      - Used for any methods that have "Simul" in them.
            %        ANOVA_Methods ([1x1] string Array or Cell String, default=self.ANOVA_Methods)
            %                      - ANOVA Methods used to calculate P-values
            %                alpha ([1x1] Numeric, default=0.05)
            %                      - Critical value for assigning statistical significance
            %  FirstLabelSortOrder (1x1 String, default={})
            %                      - Custom sort order of the label values
            %                        corresponding to the first label name
            %                        in ANOVA_Labels
            % SecondLabelSortOrder (1x1 String, default={})
            %                      - Custom sort order of the label values
            %                        corresponding to the second label name
            %                        in ANOVA_Labels
            %        LabelSortMode (1x1 String, default='ascend')
            %                      - SortMode in Echos sortBy() method
            %                      - Options are 'ascend' or 'descend'.
            %
            % See also FUNCTIONALANOVA
            
            p = inputParser;
            classInspect = ?functionalANOVA;

            addRequired(p, 'R', @(x)   ismember('OneD', superclasses(R)))                       % OneD Records to be used in F-Anova
            addRequired(p, 'ANOVA_Labels', @(x) iscellstr(x) || isstring(x) || ischar(x))
            addRequired(p, 'boundsArray', @(x) isnumeric(x)) % Domain for which to compute the L2 Statistic Over

            addParameter(p, 'N_boot', classInspect.PropertyList(13, 1).DefaultValue  , @(x) isinteger(int64(x)));
            addParameter(p, 'N_simul',classInspect.PropertyList(15, 1).DefaultValue  , @(x) isinteger(int64(x)))
            addParameter(p, 'ANOVA_Methods', classInspect.PropertyList(1).DefaultValue, @(x) iscellstr(x) || isstring(x))
            addParameter(p, 'alpha', classInspect.PropertyList(16, 1).DefaultValue  , @(x) x<=1 && x>= 0)
            addParameter(p, 'FirstLabelSortOrder', {})
            addParameter(p, 'SecondLabelSortOrder', {})
            addParameter(p, 'LabelSortMode', 'ascend')


            parse(p, R, ANOVA_Labels, boundsArray, varargin{:});

            if iscellstr(ANOVA_Labels) || ischar(ANOVA_Labels)
                ANOVA_Labels = string(ANOVA_Labels);
            end

            if iscellstr(p.Results.LabelSortMode) || ischar(p.Results.LabelSortMode)
                LabelSortMode = string(p.Results.LabelSortMode);
                assert(any(strcmpi(LabelSortMode, ["ascend", "descend"])), 'LabelSortMode must be either ''ascend'' or ''descend''')
            end

            %Verify label(1) exist for all Records
            assert(R.isLabelSet(ANOVA_Labels(1)), '"%s" Doesnt exist for all objects in array', ANOVA_Labels(1))
            % Create Cell for first Label
            label_one = ANOVA_Labels(1);
            unique_values_1 = R.uniqueLabelValues(label_one);
            n_values_1 = numel(unique_values_1);

            if numel(ANOVA_Labels) == 1
                F_WAY = 1;

                if ~isempty(p.Results.FirstLabelSortOrder)
                    assert(any(ismember(p.Results.FirstLabelSortOrder, string(unique_values_1))), 'None of the provide label values exist for label Name: %s', label_one)
                end

                R = R.sortBy(label_one, p.Results.FirstLabelSortOrder, sortMode=LabelSortMode);
            elseif numel(ANOVA_Labels) == 2
                F_WAY = 2;

                %Verify label(2) exist for all Records
                assert(R.isLabelSet(ANOVA_Labels(2)), '"%s" Doesnt exist for all objects in array', ANOVA_Labels(2))
                label_two = ANOVA_Labels(2);

                % Level 2
                unique_values_2 = R.uniqueLabelValues(label_two);
                n_values_2 = numel(unique_values_2);


                if ~isempty(p.Results.FirstLabelSortOrder)
                    assert(any(ismember(p.Results.FirstLabelSortOrder, string(unique_values_1))), 'None of the provide label values exist for label Name: %s', label_one)
                end


                if ~isempty(p.Results.SecondLabelSortOrder)
                    assert(any(ismember(p.Results.SecondLabelSortOrder, string(unique_values_2))), 'None of the provide label values exist for label Name: %s', label_two)
                end

                
                R = R.sortBy(label_one, p.Results.FirstLabelSortOrder, label_two, p.Results.SecondLabelSortOrder, sortMode=LabelSortMode);

                % Re-do after sorting
                
                % Level 1
                unique_values_1 = R.uniqueLabelValues(label_one);
                n_values_1 = numel(unique_values_1);


                % Level 2
                unique_values_2 = R.uniqueLabelValues(label_two);
                n_values_2 = numel(unique_values_2);

                N_indicators = cell(1, n_values_1);

            else
                error('F-ANOVA can only support 1 or 2 Label Names utilized for One-WAY or Two-Way Analysis')
            end

            if F_WAY == 1
                G1 = R.groupBy(label_one);
            else
                G1 = R.groupBy(label_one);

                G2_levels = cell(1, n_values_1);
                for K = 1:n_values_1
                    G2_levels{K} = G1{K}.groupBy(label_two);
                end

            end



            try
                GroupLevel1 = R.evalBy(label_one, @(r) r.toEnsemble);
            catch ME
                fprintf('%sError in converting OneD Records to Ensemble Records%s', newline, newline)
                rethrow(ME)
            end

            G1_Labels = GroupLevel1.labelTable("-unsorted", "labelNames", label_one).Variables;


            if F_WAY == 2
                for K = 1: n_values_1 % Primary Group
                    subIndicator = nan(numel(G1{K}), 1);
                    % GroupLevel1(K) = GroupLevel1(K).sortBy(level_two, {}); % Sort Before Masking
                    for KK = 1 : n_values_2 % Secondary Group
                        ReturnedValues = string(G1{K}.arrayEval(@(r) r.labelValue(label_two)));
                        mask = strcmp(ReturnedValues, string(unique_values_2{KK}));
                        subIndicator(mask) = KK;
                        N_indicators{K} = subIndicator;
                    end
                end
            end


            if F_WAY == 1
                F_ANOVA_Obj = functionalANOVA(GroupLevel1, boundsArray, ...
                    GroupLabels=string(G1_Labels));
            else
                F_ANOVA_Obj = functionalANOVA(GroupLevel1, boundsArray, ...
                    PrimaryLabels=string(G1_Labels), ...
                    SecondaryLabels=string(unique_values_2), ...
                    SubgroupIndicator=N_indicators);
            end

        end

        function CovTrace =  fastCovarianceTrace(Y, format)
            % Assumes Number of domain points is greater than number of samples

            switch format
                case 'long'
                    % Rows are elements of data and columns are realizations
                case 'wide'
                    % Rows are realizations and columns are elements of data
                    Y = Y';  % Covert to Long Format
            end

            [p, n]  = size(Y);
            dim = 2;

            if p > n
                z_mean = Y - mean(Y, dim);
                COV  =(z_mean' * z_mean)/ (n-1);
            else
                z_mean = Y - mean(Y, dim);
                COV  =(z_mean * z_mean')/ (n-1);
            end
            CovTrace = trace(COV);
        end

    end




end