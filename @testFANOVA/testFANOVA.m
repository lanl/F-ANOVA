classdef testFANOVA < matlab.unittest.TestCase
    properties
        eta_t
        COV_MATRIX = cell(3,1);
        t_grid

        genericANOVA
        genericTwoWayANOVA
    end


    methods(TestClassSetup)

        function setUpClass(tc)
            % Shared setup for the entire test class
    
            tc.t_grid = 0:0.01:25;

            % Logistic Params

            L = 7;
            x_0 = 10;
            k = 0.75;

            tc.eta_t = testFANOVA.LogisticFunction(L, k, x_0, tc.t_grid);

            % Kernel Params
            l = [0.25, 2.5, 25];
            for ii = 1:3
                Kernel_Handle = @(s,t) testFANOVA.radial_basis_kernel(s, t, l(ii));
                tc.COV_MATRIX{ii} = testFANOVA.CovarianceMatrix(Kernel_Handle, tc.t_grid);
            end

            % Generic Test 3 Groups with Same Mean Response
            rng(9000, "simdTwister")
            n_reals = 25;
            n_groups = 3;

            COV = tc.COV_MATRIX{2};
            bundle = cell(1, n_groups);
            for K =1:n_groups
                bundle{K} =  testFANOVA.Realizations(tc.eta_t, tc.t_grid, COV, n_reals);
            end

            tc.genericANOVA = functionalANOVA(bundle, [-inf, inf], d_grid=tc.t_grid);
            tc.genericANOVA.verbose = false;


            % Number of Samples within a group and Number of groups
            rng(9001, "simdTwister")
            n_groups_A = 2;
            n_groups_B = 4;

            n_sub_reals = 15;

            bundle = cell(1, n_groups_A);
            COV = tc.COV_MATRIX{2};
            Indicator_B = [];
            for K = 1:n_groups_A
                sub_bundle = [];
                for KK = 1 :n_groups_B
                    temp_bundle = testFANOVA.Realizations(tc.eta_t, tc.t_grid, COV, n_sub_reals);
                    sub_bundle = [sub_bundle, temp_bundle];

                    temp_indicator = repmat(KK, n_sub_reals, 1);
                    Indicator_B = [Indicator_B; temp_indicator];
                end
                bundle{K} =  sub_bundle;
            end

            tc.genericTwoWayANOVA = functionalANOVA(bundle, [-inf, inf], ...
                d_grid=tc.t_grid, ...
                SubgroupIndicator=Indicator_B);

            tc.genericTwoWayANOVA.verbose = false;

        end

    end

    methods(TestMethodSetup)
        % Setup for each test
    end

    methods(Test)
        % Test methods
        function R_verification_test(tc)
            group1 = [37	47	46	37	20	57	46	46	46	35	38	35	34
                36	46	44	36	18	48	38	46	42	34	37	31	31
                33	42	39	27	11	44	33	43	37	29	33	29	27
                29	34	34	20	8	35	25	40	34	28	29	26	23
                23	27	33	15	7	31	18	36	31	19	26	22	19
                18	21	27	15	5	27	15	30	25	15	20	19	15
                15	16	21	10	1	21	14	24	19	13	14	15	13
                12	12	15	6	0	18	10	17	17	6	9	11	10
                9	4	7	3	-2	14	8	12	13	0	3	6	6
                6	-1	1	0	-4	10	4	7	10	-1	-2	3	1
                4	-2	-3	3	-4	8	0	1	5	2	-1	0	0
                6	2	-1	4	1	8	3	3	7	11	5	4	3
                13	12	10	10	10	11	9	12	13	20	15	16	8
                22	24	22	14	18	21	19	24	21	30	24	27	16
                31	34	34	23	27	32	31	35	29	39	36	35	27
                38	42	43	31	36	40	37	44	39	45	44	41	34
                43	47	51	33	34	49	46	49	47	48	47	44	40
                44	48	53	35	33	55	45	49	51	48	45	45	41
                40	45	50	34	32	56	46	47	52	43	42	40	39
                35	43	49	32	28	55	41	43	47	39	40	39	36];

            group2 = [43	43	40	51	52	36	35	46	43	55	39	37	36
                41	37	41	49	46	33	37	38	41	51	38	34	33
                36	35	36	45	41	28	33	30	37	47	31	30	28
                31	28	32	39	35	22	27	23	30	41	27	27	22
                26	26	27	31	31	18	22	17	24	35	21	26	18
                20	21	20	23	24	13	14	13	16	30	14	19	13
                16	20	17	16	18	11	9	7	12	26	9	15	11
                13	15	10	9	12	7	6	3	8	23	9	10	7
                10	13	5	4	6	1	4	0	5	20	6	4	1
                8	7	-2	0	3	-2	-1	-3	2	11	4	0	-2
                4	3	-4	-5	2	-3	1	-3	1	8	2	0	-3
                12	8	-2	1	5	-3	2	1	8	10	3	4	-3
                19	16	7	10	12	2	11	11	16	19	13	12	2
                27	25	24	21	22	13	17	20	24	35	30	22	13
                37	35	38	33	32	22	34	36	33	41	37	31	22
                41	43	49	44	42	32	35	48	42	52	47	40	32
                44	49	57	51	46	39	40	55	48	57	53	43	39
                44	50	59	55	48	41	43	57	48	61	53	43	41
                41	45	54	56	49	38	43	56	48	63	49	38	38
                37	46	46	51	46	34	42	50	46	58	44	38	34];

            group3 = [36	42	38	46	54	52	32	46	46	48	44	55	48
                33	40	34	47	48	44	28	41	44	42	41	56	50
                30	40	30	44	44	44	26	38	40	42	38	51	47
                28	34	23	37	37	33	22	31	35	35	32	46	42
                21	23	17	29	30	28	19	25	31	30	24	41	37
                15	15	12	23	27	27	13	20	25	23	18	36	29
                8	11	8	19	21	23	8	13	19	19	10	30	22
                1	7	4	14	18	24	5	7	15	14	6	25	14
                -5	5	1	8	15	19	-1	1	10	9	3	21	8
                -11	6	0	3	11	15	-6	-4	5	4	0	15	5
                -12	8	-4	-2	12	15	-5	-3	3	-1	-2	8	8
                -7	12	-4	1	20	16	0	-2	6	3	-1	5	15
                4	22	1	12	30	24	12	8	14	9	6	9	24
                16	33	10	26	41	32	22	22	27	25	20	19	36
                26	43	22	39	50	43	30	34	40	37	36	31	51
                37	51	32	48	56	52	36	45	48	45	44	43	59
                44	57	38	52	61	58	39	53	53	52	49	52	63
                47	58	41	48	59	59	36	57	53	53	46	56	64
                44	54	41	43	57	57	30	55	50	52	38	59	61
                37	46	40	42	58	52	29	43	47	46	35	59	55];

            gaitData_T = [0.0250000000000000
                0.0750000000000000
                0.125000000000000
                0.175000000000000
                0.225000000000000
                0.275000000000000
                0.325000000000000
                0.375000000000000
                0.425000000000000
                0.475000000000000
                0.525000000000000
                0.575000000000000
                0.625000000000000
                0.675000000000000
                0.725000000000000
                0.775000000000000
                0.825000000000000
                0.875000000000000
                0.925000000000000
                0.975000000000000];

            groupData = {group1, group2, group3};

            f_verify = functionalANOVA(groupData,[-inf, inf], GroupLabels=["1", "2", "3"], d_grid=gaitData_T, N_boot=100000);
            f_verify.verbose = false;
            f_verify.OneWayANOVA

            L2_stat_R = 2637.1282051282050816;
            f_stat_R = 1.4682175466420953;

            % Verify Test statistics are within 0.1% of R-fdANOVA results
            tc.verifyEqual(f_verify.OneWay_P_Table.("Test-Statistic")(1), L2_stat_R, RelTol= 1e-3)  %  accurate within 0.1%.
            tc.verifyEqual(f_verify.OneWay_P_Table.("Test-Statistic")(end), f_stat_R, RelTol= 1e-3) %  accurate within 0.1%.

            % fdANOVA from documentation (non-boostrap methods)
            R_values = table(["L2-Naive";  "L2-BiasReduced"; "F-Naive";"F-BiasReduced"], ...
                [0.2106562;0.1957646;0.2226683; 0.2198691], VariableNames={'Method', 'p_value'});

            % Verify Test p-values are within 0.1% of R-fdANOVA results (non-boostrap methods)
            tc.verifyEqual(f_verify.OneWay_P_Table{[2,3, 6,7], 'P-Value'}, ...
                 R_values.p_value, RelTol= 1e-3) %  accurate within 0.1%.

            % % fdANOVA from documentation (non-boostrap methods)
            % R_boot_values = table(["L2-Bootstrap"; "F-Bootstrap"], ...
            %     [0.2151670000000000; 0.2728410000000000], VariableNames={'Method', 'p_value'});
            % 
            % % Verify Test p-values are within 1% of R-fdANOVA results (boostrap methods)
            % tc.verifyEqual(f_verify.OneWay_P_Table{[4, 8], 'P-Value'}, ...
            %      R_boot_values.p_value, RelTol= 1e-2) %  accurate within 1%.

        end

        function OneWay_Test(tc)

            % Run Default Test
            OneWay = tc.genericANOVA.copy();
            OneWay.OneWayANOVA();
            tc.verifyEqual(OneWay.OneWay_P_Table.("P-Value"), ones(height(OneWay.OneWay_P_Table),1)*0.78, 'AbsTol', 0.02)


            % Test 3 Groups with one having a different Response
            rng(9001, "simdTwister")
            n_reals = 15;
            n_groups = 3;
            COV = tc.COV_MATRIX{2};

            bundle = cell(1, n_groups);
            for K =1:n_groups
                bundle{K} =  testFANOVA.Realizations(tc.eta_t, tc.t_grid, COV, n_reals);
                if K == 1
                    bundle{K} = randi([1 5], 1) + bundle{K};  % Shift up 1 to 5 units
                end
            end

            Test1_3Groups = functionalANOVA(bundle, [-inf, inf], d_grid=tc.t_grid);
            Test1_3Groups.verbose=false;
            Test1_3Groups.OneWayANOVA

            tc.verifyLessThan(Test1_3Groups.OneWay_P_Table.("P-Value"), 0.05)

            Test1_3Groups.OneWayANOVA(Hypothesis='pairwise')

            tc.verifyLessThan(Test1_3Groups.OneWay_P_Table{1, 2:end}, 0.05)    % One and Two are different
            tc.verifyLessThan(Test1_3Groups.OneWay_P_Table{2, 2:end}, 0.05)    % One and Three are different
            tc.verifyGreaterThan(Test1_3Groups.OneWay_P_Table{3, 2:end}, 0.05) % Two and Three are equiv


            % Test 2 Groups with Different Mean Response
            rng(9000, "simdTwister")
            n_reals = 25;
            n_groups = 2;

            bundle = cell(1, n_groups);
            for K =1:n_groups
                bundle{K} =  testFANOVA.Realizations(tc.eta_t, tc.t_grid, COV, n_reals);
                if K > 1
                    bundle{K} = randi([1 5], 1) + bundle{K};  % Shift up 1 to 5 units
                end
            end

            Test1_2Groups = functionalANOVA(bundle, [-inf, inf], d_grid=tc.t_grid);
            Test1_2Groups.verbose=false;
            Test1_2Groups.OneWayANOVA

            tc.verifyLessThan(Test1_2Groups.OneWay_P_Table.("P-Value"), 0.05)
        end

        function OneWayBF_Test(tc)
            % Test 3 Groups with Same Mean Response
            rng(9001, "simdTwister")
            n_reals = 30;
            n_groups = 3;

            bundle = cell(1, n_groups);
            for K =1:n_groups
                COV = tc.COV_MATRIX{K};
                bundle{K} =  testFANOVA.Realizations(tc.eta_t, tc.t_grid, COV, n_reals);
            end

            Test2_3Groups = functionalANOVA(bundle, [-inf, inf], d_grid=tc.t_grid);
            Test2_3Groups.verbose=false;
            Test2_3Groups.OneWayANOVA_BF

            tc.verifyGreaterThan(Test2_3Groups.OneWay_BF_P_Table.("P-Value"), 0.05) 

            %  Shift last group
            bundle{end} = bundle{end} + 1;
            Test2_3Groups = functionalANOVA(bundle, [-inf, inf], d_grid=tc.t_grid);
            Test2_3Groups.verbose=false;
            Test2_3Groups.OneWayANOVA_BF

             tc.verifyLessThan(Test2_3Groups.OneWay_BF_P_Table.("P-Value"), 0.05)

             % Verify PairWise with 3 being different
             Test2_3Groups.OneWayANOVA_BF(Hypothesis='pairwise')
             tc.assertGreaterThan(Test2_3Groups.OneWay_BF_P_Table{1, 2:end}, 0.05)    % One and Two are equiv
             tc.verifyLessThan(Test2_3Groups.OneWay_BF_P_Table{2, 2:end}, 0.05)       % One and Three are different
             tc.verifyLessThan(Test2_3Groups.OneWay_BF_P_Table{3, 2:end}, 0.05)       % Two and Three are different
        end

        function TwoWay_Test(tc)

            test2Way = tc.genericTwoWayANOVA.copy();
            
            test2Way.TwoWayANOVA(hypothesis="FAMILY");
            tc.verifyGreaterThan(test2Way.TwoWay_P_Table.("P-Value"), 0.05)

            test2Way.TwoWayANOVA(hypothesis="interaction")
            tc.verifyGreaterThan(test2Way.TwoWay_P_Table.("P-Value"), 0.05)

            tc.genericTwoWayANOVA.TwoWayANOVA(hypothesis="primary")
            tc.verifyGreaterThan(test2Way.TwoWay_P_Table.("P-Value"), 0.05)

            tc.genericTwoWayANOVA.TwoWayANOVA(hypothesis="secondary")
            tc.verifyGreaterThan(test2Way.TwoWay_P_Table.("P-Value"), 0.05)

            test2Way.TwoWayANOVA(hypothesis="PAIRWISE", N_boot=2500);
            tc.verifyGreaterThan(test2Way.TwoWay_P_Table{:, 2:end}, 0.05)

            % Modify for significant diff
            n_groups_A = 2;
            n_groups_B = 4;

            n_sub_reals = 15;
            COV = tc.COV_MATRIX{2};
            bundle = cell(1, n_groups_A);
            Indicator_B = [];
            for K = 1:n_groups_A
                sub_bundle = [];
                for KK = 1 :n_groups_B
                    temp_bundle = testFANOVA.Realizations(tc.eta_t, tc.t_grid, COV, n_sub_reals);
                    
                    if KK == 1
                        temp_bundle = temp_bundle + ones(size(temp_bundle));
                    end

                    sub_bundle = [sub_bundle, temp_bundle];

                    temp_indicator = repmat(KK, n_sub_reals, 1);
                    Indicator_B = [Indicator_B; temp_indicator];
                end
                bundle{K} =  sub_bundle;
            end

            test2Way_pants = functionalANOVA(bundle, [-inf, inf], d_grid=tc.t_grid, ...
                SubgroupIndicator=Indicator_B,...
                PrimaryLabels=["Bladder", "Bottle"],...
                SecondaryLabels=["Pants", "Capri", "Skirt", "Shorts"]);

            test2Way_pants.verbose = false;

            test2Way_pants.TwoWayANOVA(hypothesis="FAMILY");
            tc.verifyLessThan(test2Way_pants.TwoWay_P_Table.("P-Value"), 0.05)  % Should return that something is causing a significant affect

            test2Way_pants.TwoWayANOVA(hypothesis="primary");     % should return as not-significant. Constant was added to first level of secondary factor
            tc.verifyGreaterThan(test2Way_pants.TwoWay_P_Table.("P-Value"),  0.05)  % Not signficant

            test2Way_pants.TwoWayANOVA(hypothesis="secondary");    % should return as significant. Constant was added to first level of secondary factor
            tc.verifyLessThan(test2Way_pants.TwoWay_P_Table.("P-Value"), 0.05)  % Should return that something is causing a significant affect

            test2Way_pants.TwoWayANOVA(hypothesis="interaction");
            tc.verifyGreaterThan(test2Way_pants.TwoWay_P_Table.("P-Value"),  0.05)  % Not signficant
% 
            test2Way_pants.TwoWayANOVA(hypothesis="pairwise", N_boot=2500);
            pantsCountMask = count(test2Way_pants.TwoWay_P_Table.Hypothesis, 'Pants') == 1;  % Only check where pants occurs once in pairwise
            subTable = test2Way_pants.TwoWay_P_Table(pantsCountMask, :);

            tc.verifyLessThan(subTable{:, 2:end}, 0.05)  % Any pair with one occurance of pants should be statistically significant


        end

        function TwoWayBF_Test(tc)

            % Number of Samples within a group and Number of groups
            rng(9001, "simdTwister")
            n_groups_A = 2;
            n_groups_B = 4;

            n_sub_reals = 15;

            bundle = cell(1, n_groups_A);
            Indicator_B = [];
            for K = 1:n_groups_A
                sub_bundle = [];
                COV = tc.COV_MATRIX{K};
                for KK = 1 :n_groups_B
                    temp_bundle = testFANOVA.Realizations(tc.eta_t, tc.t_grid, COV, n_sub_reals);
                    sub_bundle = [sub_bundle, temp_bundle];

                    temp_indicator = repmat(KK, n_sub_reals, 1);
                    Indicator_B = [Indicator_B; temp_indicator];
                end
                bundle{K} =  sub_bundle;
            end

            test2Way_BF = functionalANOVA(bundle, [-inf, inf], d_grid=tc.t_grid, ...
                SubgroupIndicator=Indicator_B,...
                PrimaryLabels=["Bladder", "Bottle"],...
                SecondaryLabels=["Pants", "Capri", "Skirt", "Shorts"]);

            test2Way_BF.verbose = false;

            test2Way_BF.TwoWayANOVA_BF(hypothesis="FAMILY");
            tc.verifyGreaterThan(test2Way_BF.TwoWay_BF_P_Table.("P-Value"), 0.05) 

            test2Way_BF.TwoWayANOVA_BF(hypothesis="PAIRWISE");
            tc.verifyGreaterThan(test2Way_BF.TwoWay_BF_P_Table{:, 2:end}, 0.05)

            % Test for Significance with Pants 
            rng(9001, "simdTwister")
            n_groups_A = 2;
            n_groups_B = 4;

            n_sub_reals = 15;

            bundle = cell(1, n_groups_A);
            Indicator_B = [];
            for K = 1:n_groups_A
                sub_bundle = [];
                COV = tc.COV_MATRIX{K};
                for KK = 1 :n_groups_B
                    temp_bundle = testFANOVA.Realizations(tc.eta_t, tc.t_grid, COV, n_sub_reals);

                    if KK == 1
                        temp_bundle = temp_bundle + ones(size(temp_bundle));
                    end

                    sub_bundle = [sub_bundle, temp_bundle];

                    temp_indicator = repmat(KK, n_sub_reals, 1);
                    Indicator_B = [Indicator_B; temp_indicator];
                end
                bundle{K} =  sub_bundle;
            end

            test2WayPants_BF = functionalANOVA(bundle, [-inf, inf], d_grid=tc.t_grid, ...
                SubgroupIndicator=Indicator_B,...
                PrimaryLabels=["Bladder", "Bottle"],...
                SecondaryLabels=["Pants", "Capri", "Skirt", "Shorts"]);

            test2WayPants_BF.verbose = false;

            test2WayPants_BF.TwoWayANOVA_BF(hypothesis="FAMILY");
            tc.verifyLessThan(test2WayPants_BF.TwoWay_BF_P_Table.("P-Value"), 0.05) 

            test2WayPants_BF.TwoWayANOVA_BF(hypothesis="PAIRWISE");
            pantsCountMask = count(test2WayPants_BF.TwoWay_BF_P_Table.Hypothesis, 'Pants') == 1;  % Only check where pants occurs once in pairwise
            subTable = test2WayPants_BF.TwoWay_BF_P_Table(pantsCountMask, :);
            tc.verifyLessThan(subTable{:, 2:end}, 0.05)  % Any pair with one occurance of pants should be statistically significant

            test2WayPants_BF.TwoWayANOVA_BF(hypothesis="primary");     % should return as not-significant. Constant was added to first level of secondary factor
            tc.verifyGreaterThan(test2WayPants_BF.TwoWay_BF_P_Table.("P-Value"),  0.05)  % Not signficant

            test2WayPants_BF.TwoWayANOVA_BF(hypothesis="secondary");     % should return assignificant. Constant was added to first level of secondary factor
            tc.verifyLessThan(test2WayPants_BF.TwoWay_BF_P_Table.("P-Value"),  0.05)  % signficant

            test2WayPants_BF.TwoWayANOVA_BF(hypothesis="interaction");
            tc.verifyGreaterThan(test2WayPants_BF.TwoWay_BF_P_Table.("P-Value"),  0.05)  % Not signficant

        end
    
        function PlotNullDist_Test(tc)
            % Find current Figures before class test
            startingChildren = get(0);
            startingChildren = startingChildren.Children;

           
            tc.genericANOVA.showSimulPlot = true;
            tc.genericANOVA.OneWayANOVA(ANOVA_Methods=["L2-Simul", "F-Simul"])


            % Find all Figures after class test
            endingChildren = get(0);
            endingChildren = endingChildren.Children;

            % Close all figures generated during test
            close(setdiff(endingChildren, startingChildren))

        end

        function PlotMean_Test(tc)
            % Find current Figures before class test
            startingChildren = get(0);
            startingChildren = startingChildren.Children;

           
            genericOne = tc.genericANOVA.copy();
            genericOne.PlotMeans()
            genericOne.PlotMeans(GroupLabels = ["Apple", "Pear", "Orange"])


            genericTwo = tc.genericTwoWayANOVA.copy();

            genericTwo.PlotMeans()
            genericTwo.PlotMeans(plotType='Primary')
            genericTwo.PlotMeans(plotType='Secondary')
            genericTwo.PlotMeans(plotType='Interaction')


            genericTwo.PlotMeans(PrimaryLabels=["Bladder", "Bottle"],...
            SecondaryLabels=["Pants", "Capri", "Skirt", "Shorts"])
            genericTwo.PlotMeans(plotType='Primary')
            genericTwo.PlotMeans(plotType='Secondary')
            genericTwo.PlotMeans(plotType='Interaction')

            % Find all Figures after class test
            endingChildren = get(0);
            endingChildren = endingChildren.Children;

            % Close all figures generated during test
            close(setdiff(endingChildren, startingChildren))
        end

        function PlotCovariances_Test(tc)
            % Find current Figures before class test
            startingChildren = get(0);
            startingChildren = startingChildren.Children;

            genericOne = tc.genericANOVA.copy();
            genericOne.PlotCovariances()
            genericOne.PlotCovariances(GroupLabels = ["Apple", "Pear", "Orange"])


            genericTwo = tc.genericTwoWayANOVA.copy();
            
            genericTwo.PlotCovariances()
            genericTwo.PlotCovariances(plotType='Primary')
            genericTwo.PlotCovariances(plotType='Secondary')
            genericTwo.PlotCovariances(plotType='Interaction')


            genericTwo.PlotCovariances(PrimaryLabels=["Bladder", "Bottle"],...
            SecondaryLabels=["Pants", "Capri", "Skirt", "Shorts"])
            genericTwo.PlotCovariances(plotType='Primary')
            genericTwo.PlotCovariances(plotType='Secondary')
            genericTwo.PlotCovariances(plotType='Interaction')

            % Find all Figures after class test
            endingChildren = get(0);
            endingChildren = endingChildren.Children;

            % Close all figures generated during test
            close(setdiff(endingChildren, startingChildren))
        end

        function DistributionCheck_Test(tc)

            % Find current Figures before class test
            startingChildren = get(0);
            startingChildren = startingChildren.Children;


            genericOne = tc.genericANOVA.copy();
            genericOne.PointwiseDistributionCheck(4, 'boxplot')
            genericOne.PointwiseDistributionCheck(4, 'qq')
            genericOne.PointwiseDistributionCheck(4, 'shapiro-wilkes')

            % Find all Figures after class test
            endingChildren = get(0);
            endingChildren = endingChildren.Children;

            % Close all figures generated during test
            close(setdiff(endingChildren, startingChildren))

        end
     
        function TwoSampleCovariance_Test(tc)
  
            DataStruct = testFANOVA.GenerateExampleData();
            DomainData = DataStruct(1).domainData;

            DataTwoWay = DataStruct(3).Data(1:2);

            test2SampleCov_different = functionalANOVA(DataTwoWay, [-inf, inf], d_grid=DomainData);
            test2SampleCov_different.N_boot=500;
            test2SampleCov_different.N_permutations=500;
            test2SampleCov_different.CovarianceTest_TwoSample()

            tc.verifyLessThan(test2SampleCov_different.COVAR_P_Table.("P-Value"), 0.05)


            DataTwoWay = DataStruct(1).Data(1:2);

            test2SampleCov_Same = functionalANOVA(DataTwoWay, [-inf, inf], d_grid=DomainData);
            test2SampleCov_Same.N_boot=500;
            test2SampleCov_Same.N_permutations=500;
            test2SampleCov_Same.CovarianceTest_TwoSample()

            tc.verifyGreaterThan(test2SampleCov_Same.COVAR_P_Table.("P-Value"), 0.05)


        end

        function KSampleCovariance_Test(tc)

            rng(9001, "simdTwister")
            n_reals = 30;
            n_groups = 3;

            bundle = cell(1, n_groups);
            for K =1:n_groups
                COV = tc.COV_MATRIX{K};
                bundle{K} =  testFANOVA.Realizations(tc.eta_t, tc.t_grid, COV, n_reals);
            end

            testKSampleCov_different = functionalANOVA(bundle, [-inf, inf], d_grid=tc.t_grid);
            testKSampleCov_different.verbose = false;
            testKSampleCov_different.N_boot=500;
            testKSampleCov_different.N_permutations=500;
            testKSampleCov_different.CovarianceTest()

            tc.verifyLessThan(testKSampleCov_different.COVAR_P_Table.("P-Value"), 0.05)

            rng(9100, "simdTwister")
            n_reals = 30;
            n_groups = 3;

            bundle = cell(1, n_groups);
            for K =1:n_groups
                COV = tc.COV_MATRIX{2};
                bundle{K} =  testFANOVA.Realizations(tc.eta_t, tc.t_grid, COV, n_reals);
            end

            testKSampleCov_Same = functionalANOVA(bundle, [-inf, inf], d_grid=tc.t_grid);
            testKSampleCov_Same.verbose = false;
            testKSampleCov_Same.N_boot=1000;
            testKSampleCov_Same.N_permutations=1000;
            testKSampleCov_Same.CovarianceTest()

            tc.verifyGreaterThan(testKSampleCov_Same.COVAR_P_Table.("P-Value"), 0.05)


        end

    end

    methods(Static)
        
        function f_x = LogisticFunction(L, k, x_0, x)
            f_x =  L ./ (1 + exp(-k.*(x-x_0)));
        end

        function K = radial_basis_kernel(s,t, l)
            K = exp(-(norm(s-t, 2)) ./ (2.*l^2));
        end

        function COV_MATRIX = CovarianceMatrix(Kernel_Handle, t)
            N = numel(t);
            COV_MATRIX = zeros(N, N);

            for K = 1 : N
                for KK = 1: N
                    if KK <= K
                        COV_MATRIX(K, KK) = Kernel_Handle(t(K), t(KK));
                    else
                        continue
                    end

                end
            end

            COV_MATRIX = tril(COV_MATRIX) + tril(COV_MATRIX, -1)';   % Utilize Symmetric Properties
        end

        function R = Realizations(eta_t, t_grid, COV_MATRIX, n_reals)
            R = eta_t + mvnrnd(zeros(n_reals, numel(t_grid)), COV_MATRIX);
            R = R';
        end

        function DataStruct = GenerateExampleData()
            
            DataStruct =  struct('RunType', [], 'covarianceType', [] ,'Data', [], 'Notes', [], ...
                'samplesPerGroup', [], 'levelsPerPrimaryFactor', ...
                [], 'levelsPerSecondaryFactor', [], 'domainData', [], 'IndicatorCell', []);

            baseline = 75;

            t_end = 25;
            t_grid = 0:0.01:t_end;

            N = length(t_grid);
            % Logistic Params

            L = 75;
            x_0 = t_end/2;
            k = 0.75;

            eta_t = testFANOVA.LogisticFunction(L, k, x_0, t_grid);

            % Kernel Params
            l = [0.1, 10, 100];
            sigma = [0.01, 0.01, 0.01];
            COV_MATRIX = cell(1, 3);

            for ii = 1:3
                Kernel_Handle = @(s,t) testFANOVA.radial_basis_kernel(s, t, l(ii));
                COV_MATRIX{ii} = testFANOVA.CovarianceMatrix(Kernel_Handle, t_grid);
            end

            %% OneWAY 
            % Generic Test 3 Groups with Same Mean Response
            rng(9001, "simdTwister")
            n_reals = 25;
            n_groups = 3;

            COV = COV_MATRIX{2};
            bundle = cell(1, n_groups);
            for K =1:n_groups
                bundle{K} =  testFANOVA.Realizations(eta_t, t_grid, COV, n_reals) + baseline;
                bundle{K} =  bundle{K} +  sigma(2)*randn(N, 1);
            end

            DataStruct(1).RunType = 'OneWay';
            DataStruct(1).covarianceType = 'Equal';
            DataStruct(1).Data = bundle;
            DataStruct(1).Notes = 'OneWay data set that is statistically not significant between groups';
            DataStruct(1).samplesPerGroup = n_reals;
            DataStruct(1).levelsPerPrimaryFactor = n_groups;
            DataStruct(1).levelsPerSecondaryFactor = nan;
            DataStruct(1).domainData = t_grid;

            % Make Echo Records
            if 8==exist('OneD', 'class')
                Rlabels = ["Cold", "Hot", "Ambient"];
                R_bundle = [];
                for K = 1 : numel(bundle)
                    for KK = 1:size(bundle{K}, 2)
                        r = TimeRecord(t_grid, bundle{K}(:,KK));
                        r.label('Temperature', Rlabels(K), 'SampleNumber', KK, 'RunType', DataStruct(1).RunType, 'covarianceType',  DataStruct(1).covarianceType, 'notes', DataStruct(1).Notes)
                        R_bundle = [R_bundle; r];
                    end
                end

                R_bundle = R_bundle.assignUnits('t', 'minute');
                R_bundle.setQuantity('t', 'Time', 'data', 'Heart Rate')
                DataStruct(1).EchoTimeRecords = R_bundle;
            end

            DataStruct(2).RunType = 'OneWay';
            DataStruct(2).covarianceType = 'Equal';
            modifiedBundle = bundle;
            modifiedBundle{1} = bundle{1} - 2;
            DataStruct(2).Data = modifiedBundle;
            DataStruct(2).Notes = 'OneWay data set that is statistically significant. level/Group 1 is different';
            DataStruct(2).samplesPerGroup = n_reals;
            DataStruct(2).levelsPerPrimaryFactor = n_groups;
            DataStruct(2).levelsPerSecondaryFactor = nan;

            % Make Echo Records
            if 8==exist('OneD', 'class')
                R_bundleModified = R_bundle.copy();
                R_bundleModified.removeLabels('notes');
                R_bundleModified.label('notes', 'OneWay data set that is statistically significant. level/Group 1 is different')
                R_ColdModified = R_bundleModified.pull('temperature', 'cold') - 2;
                DataStruct(2).EchoTimeRecords = [R_bundleModified.purge('temperature', 'cold'), R_ColdModified];
            end

            %% OneWay BF Data
            % Test 3 Groups with Same Mean Response
            rng(9000, "simdTwister")
            n_reals = 15;
            n_groups = 3;

            bundle = cell(1, n_groups);
            for K =1:n_groups
                COV = COV_MATRIX{K};
                bundle{K} =  testFANOVA.Realizations(eta_t, t_grid, COV, n_reals) + baseline;
                bundle{K} = bundle{K} +  sigma(K)*randn(N, 1);

            end

            DataStruct(3).RunType = 'OneWay';
            DataStruct(3).covarianceType = 'Unequal';
            DataStruct(3).Data = bundle;
            DataStruct(3).Notes = 'OneWay data set that is not statistically significant but different covariances';
            DataStruct(3).samplesPerGroup = n_reals;
            DataStruct(3).levelsPerPrimaryFactor = n_groups;
            DataStruct(3).levelsPerSecondaryFactor = nan;

            % Make Echo Records
            if 8==exist('OneD', 'class')
                Rlabels = ["Cold", "Hot", "Ambient"];
                R_bundle = [];
                for K = 1 : numel(bundle)
                    for KK = 1:size(bundle{K}, 2)
                        r = TimeRecord(t_grid, bundle{K}(:,KK));
                        r.label('Temperature', Rlabels(K), 'SampleNumber', KK, 'RunType', DataStruct(3).RunType, 'covarianceType',  DataStruct(3).covarianceType, 'notes', DataStruct(3).Notes)
                        R_bundle = [R_bundle; r];
                    end
                end
                R_bundle = R_bundle.assignUnits('t', 'minute');
                R_bundle.setQuantity('t', 'Time', 'data', 'Heart Rate')
                DataStruct(3).EchoTimeRecords = R_bundle;
            end
            %% TwoWay Data
            rng(9001, "simdTwister")
            n_groups_A = 2;
            n_groups_B = 4;
            n_sub_reals = 15;

            bundle = cell(1, n_groups_A);
            COV = COV_MATRIX{2};
            Indicator_B = [];
            Indicator_B_cell = {};

            PrimaryLabels=["Bladder", "Bottle"];
            SecondaryLabels=["Pants", "Shorts", "Skirt", "Dress"];
            R_bundle = [];

            DataStruct(4).RunType = 'TwoWay';
            DataStruct(4).covarianceType = 'equal';
            DataStruct(4).Data = bundle;
            DataStruct(4).Notes = 'TwoWay data set that is statistically significant. First level within secondary factor';


            for K = 1:n_groups_A
                sub_bundle = [];
                sub_indicator = [];

                for KK = 1 :n_groups_B
                    temp_bundle = testFANOVA.Realizations(eta_t, t_grid, COV, n_sub_reals) + baseline;
                    temp_bundle =  temp_bundle +  sigma(2)*randn(N, 1);


                    if KK == 1
                        temp_bundle = temp_bundle + 3;
                    end

                    sub_bundle = [sub_bundle, temp_bundle];
                    temp_indicator = repmat(KK, n_sub_reals, 1);
                    Indicator_B = [Indicator_B; temp_indicator];
                    sub_indicator = [sub_indicator; temp_indicator];


                    % Make Echo Records
                    if 8==exist('OneD', 'class')
                        for KKK = 1 : size(temp_bundle, 2)
                            r = TimeRecord(t_grid, temp_bundle(:,KKK));

                            r.label('HydrationType', PrimaryLabels(K), 'ClothingType',SecondaryLabels(KK), ...
                                'SampleNumber', KKK, ...
                                'RunType', DataStruct(4).RunType, ...
                                'covarianceType',  DataStruct(4).covarianceType, ...
                                'notes', DataStruct(4).Notes)

                            R_bundle = [R_bundle; r];
                        end
                    end
                end

                bundle{K} =  sub_bundle;
                Indicator_B_cell{K} = sub_indicator;

            end

            % Make Echo Records
            if 8==exist('OneD', 'class')
                R_bundle = R_bundle.assignUnits('t', 'minute');
                R_bundle.setQuantity('t', 'Time', 'data', 'Heart Rate')
            end
            
            DataStruct(4).Data = bundle;
            DataStruct(4).EchoTimeRecords = R_bundle;
            DataStruct(4).samplesPerGroup = n_sub_reals;
            DataStruct(4).levelsPerPrimaryFactor = n_groups_A;
            DataStruct(4).levelsPerSecondaryFactor = n_groups_B;
            DataStruct(4).IndicatorCell = Indicator_B_cell;


            %%


        end
    end

end