function PlotTestStats(p_value, alpha, NullDist, testStat, testName, Hypothesis, HypothesisLabel)
%PLOTTESTSTATS Plots Null Distribution and Test Statistic
% PLOTTESTSTATS(P_VALUE, ALPHA, NULLDIST, TESTSTAT, TESTNAME, HYPOTHESIS, HYPOTHESISLABEL)
%
% Plotting tool for funtionalANOVA class. Default behavior is to not plot the
% null distribution and test statistics for simulation based F-ANOVA
% methods. Can be turned on by setting the class property to "showSimulPlot=true". 
% PLOTTESTSTATS currently only works with the "OneWayANOVA" class method
% for homogeneous OneWay ANOVA analyses.
%
% <strong>Required Input</strong>
%
%         p_value ([1x1] Numeric)
%                 - P value from a F-ANOVA simulation method
%        NullDist ([Qx1] Numeric, default=10000)
%                 - N: Simulated Null Distribution
%                 - Q=self.N_simul
%        testStat ([1x1] Numeric)
%                 - Test statistic generated from either L2 or F-type test method
%        testName ([1x1] String)
%                 - Name of the ANOVA Method
%      Hypothesis ([1x1] String)
%                 - Name of the Hypothesis type.
%                 - Options: "FAMILY" or "PAIRWISE"
% HypothesisLabel ([1x1] String)
%                 - Label for the specific hypothesis
%                 - When Hypothesis="Family", HypothesisLabel="Family"
%                 - When Hypothesis="PAIRWISE", HypothesisLabel is the name of the
%                   levels used in the pairwise comparison
%
% See also FUNCTIONALANOVA, TWOWAYANOVA

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

if p_value <= alpha
    line_label = sprintf('%s Statistic P-Value: p=%0.4f <= %0.3f', testName, p_value, alpha);
    verdict_label = 'Verdict: Reject H_0';
else
    line_label = sprintf('%s Statistic P-Value: p=%0.4f > %0.3f', testName, p_value, alpha);
    verdict_label = 'Verdict: Fail to Reject H_0';
end

figure1 = figure('Name', sprintf('Null Distribution Plot for %s', testName), 'units', 'points', 'position', [90, 90, 1000, 800]);
axes1 = axes('Parent',figure1);
hold(axes1,'on');

hold on
H1 = histogram(NullDist, 'Parent',axes1, 100, "Normalization","pdf");
H2 = xline(testStat, '-r', testName, 'LineWidth',1.5, 'Color', 'r', 'LineStyle','--', 'LabelHorizontalAlignment','left');



if any(testStat > max(NullDist) * 1E2) || any(testStat < min (NullDist) * 1E2)
    set(axes1,'XMinorTick','on','XScale','log');
end

crit_value = quantile(NullDist, 1 - alpha);
ylim_vals = max(ylim);
y =[0, ylim_vals; 0, ylim_vals];
x = [crit_value, max(xlim)];
area(x, y, 'FaceAlpha', 0.3)
H3 = xline(crit_value, 'k' ,'Beginning of Critical Value Region', 'LineWidth',1.5, 'Color', 'k', 'LineStyle','-', 'LabelHorizontalAlignment','right');

if contains(testName, 'F','IgnoreCase', true)
    nullDistLabel = 'Simulated F-type Mixture Null Distribution';
    switch Hypothesis
        case "FAMILY"
            title_label = 'One-Way, Family, Functional ANOVA: F-type test';
        case "PAIRWISE"
            title_label = sprintf('One-Way, Pairwise (%s), Functional ANOVA: F-type test', HypothesisLabel);
    end
else
    nullDistLabel = 'Simulated \chi^2_{1}-type Mixture Null Distribution';
    switch Hypothesis
        case "FAMILY"
        title_label = 'One-Way, Family, Functional ANOVA: Squared L-2 Norm test';
        case "PAIRWISE"
        title_label = sprintf('One-Way, Pairwise (%s), Functional ANOVA: Squared L-2 Norm test', HypothesisLabel);       
    end
end


legend( [H1, H2, H3],{nullDistLabel, line_label, 'Region for Statistical Significance (Reject H_0)'}, Location="southoutside")

title({title_label ,verdict_label})
ylabel('PDF')
xlabel('Null Distribution')
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 18)
end
