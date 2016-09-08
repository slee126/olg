% clear; clc; close all;
load cps;
load popMat;
load cpsTable_Med;
gdp = xlsread('gdp.xlsx');
wageMat = tableAll(1:40, 5:4:end);
alphaMat = table2array(wageMat(:, 2:7))./table2array(repmat(wageMat(:, 1), 1, 6));

age_end = 75;
n_age = age_end - 20 + 1;
N_period = 2013 - 1976 + 1;
A = zeros(N_period*(age_end-20+1), 11);
for j = 1:N_period
    j_year = j + 1975;
    t_Period = tbl(tbl.YEAR == j_year, :);   
    gdpDiv = gdp(j);
    for i = 20:age_end
        temp2 = t_Period(t_Period.AGE == i & t_Period.INCWAGE > 0 & t_Period.WKSWORK1 > 0, :);
        popS = sum(temp2.WTSUPP);
        wageA = sum(temp2.INCWAGE.*temp2.WTSUPP)/popS;
        temp=[temp2.INCWAGE temp2.WTSUPP];
        temp = sortrows(temp, 1);
        temp = [temp, cumsum(temp(:, 2))];    
        med = find(temp(:, 3) > temp(end, 3)/2);
        wageM = temp(med(1), 1);
        
        hrA = sum(temp2.WKSWORK1.*temp2.UHRSWORKLY.*temp2.WTSUPP)/popS/52;
        temp = [temp2.WKSWORK1.*temp2.UHRSWORKLY temp2.WTSUPP];
        temp = sortrows(temp, 1);
        temp = [temp, cumsum(temp(:, 2))];    
        med = find(temp(:, 3) > temp(end, 3)/2);
        hrM = temp(med(1), 1)/52;
       
        hrWageA = sum(temp2.INCWAGE./(temp2.UHRSWORKLY.*temp2.WKSWORK1).*temp2.WTSUPP)/popS;
        temp = [temp2.INCWAGE./(temp2.UHRSWORKLY.*temp2.WKSWORK1) temp2.WTSUPP];
        temp = sortrows(temp, 1);
        temp = [temp, cumsum(temp(:, 2))];    
        med = find(temp(:, 3) > temp(end, 3)/2);
        hrWageM = temp(med(1), 1);
        
        temp2 = t_Period(t_Period.AGE == i &  t_Period.EDUC >= 110, :);
        popC = sum(temp2.WTSUPP);        
   
        A(n_age*(j-1)+i-19, :) = [popS wageA/gdpDiv*100 wageM/gdpDiv*100 hrA hrM ...
            hrWageA/gdpDiv*100 hrWageM/gdpDiv*100 hrWageM/gdpDiv*100/table2array(wageMat(j, 1)) popC j_year i];     
    end
end
 
popVec = simul_pop20(7:44, 3:n_age+2)';
popVec = popVec(:);
educVec = A(:, 9)./popVec;
educMat = reshape(educVec, n_age, 38)';
educMat_HP = hpfilter(educMat, 100);
educMat_SG = sgolayfilt(educMat, 3, 21, [], 2);

workingVec = A(:, 1)./popVec;
workingMat = reshape(workingVec, n_age, 38)';
workingMat_HP = hpfilter(workingMat, 100);
workingMat_SG = sgolayfilt(workingMat, 3, 21, [], 2);

effVec = A(:, 8);
effMat = reshape(effVec, n_age, 38)';
effMat_HP = hpfilter(effMat, 100);
effMat_SG = sgolayfilt(effMat, 3, 21, [], 2);

wageVec = A(:, 7);
wageMat = reshape(wageVec, n_age, 38)';
wageMat_HP = hpfilter(wageMat, 100);
wageMat_SG = sgolayfilt(wageMat, 3, 21, [], 2);

t_Year = 1976:2013;
year_label = cell(38, 1);
Y = zeros(38, n_age);
for i = 1:N_period
    year_label{i} = num2str(t_Year(i));
end

%===========================================================================
figure(1); clf;
set(1, 'defaulttextinterpreter', 'latex');
plot(20:age_end, effMat_HP(1, :), 'LineWidth', 2); hold on;
for i = 2:3:N_period-1
   plot(20:age_end, effMat_HP(i, :)); hold on; 
end
plot(20:age_end, effMat_HP(end, :), 'LineWidth', 2);
% lgnd = legend(year_label{1}, year_label{2:3:end}, year_label{end}, 'FontSize', 18, 'Location', 'best');
% set(lgnd, 'Color', 'none', 'box', 'off');
clickableLegend(year_label{1}, year_label{2:3:end-1}, year_label{end}, 'Location', 'best');
grid on; axis tight; title('Efficiency Weights HP Filtered', 'FontSize', 18);
set(gcf, 'Color', 'w'); set(gcf, 'Position', [0 0 800 600]);
export_fig figure2.pdf;
%--------------------------------------------------------------------------

figure(2); clf;
set(2, 'defaulttextinterpreter', 'latex')
plot(20:age_end, educMat_HP(1, :), 'LineWidth', 2); hold on;
for i = 2:3:N_period-1
    plot(20:age_end, educMat_HP(i, :)); hold on;
end
plot(20:age_end, educMat_HP(N_period, :), 'LineWidth', 2); hold on;
% lgnd = legend(year_label{1}, year_label{2:3:end}, year_label{end}, 'FontSize', 18, 'Location', 'best');
% set(lgnd, 'Color', 'none', 'box', 'off');
grid on; axis tight; title('\% of Cohorts with Bachelor''s or more', 'FontSize', 18);
clickableLegend(year_label{1}, year_label{2:3:end-1}, year_label{end}, 'Location', 'best');
set(gcf, 'Color', 'w'); set(gcf, 'Position', [0 0 800 600]);
export_fig figure3.png;

figure(3); clf;
set(3, 'defaulttextinterpreter', 'latex')
plot(20:age_end, workingMat_HP(1, :), 'LineWidth', 2); hold on;
for i = 2:3:N_period-1
    plot(20:age_end, workingMat_HP(i, :)); hold on;
end
plot(20:age_end, workingMat_HP(N_period, :), 'LineWidth', 2); hold on;
% lgnd = legend(year_label{1}, year_label{2:3:end}, year_label{end}, 'FontSize', 18, 'Location', 'best');
% set(lgnd, 'Color', 'none', 'box', 'off');
grid on; axis tight; title('\% of Cohorts Employed', 'FontSize', 18);
clickableLegend(year_label{1}, year_label{2:3:end-1}, year_label{end}, 'Location', 'best');
set(gcf, 'Color', 'w'); set(gcf, 'Position', [0 0 800 600]);
export_fig working.png;


figure(4); clf;
set(4, 'defaulttextinterpreter', 'latex')
plot(20:age_end, wageMat_HP(1, :), 'LineWidth', 2); hold on;
for i = 2:3:N_period-1
    plot(20:age_end, wageMat_HP(i, :)); hold on;
end
plot(20:age_end, wageMat_HP(N_period, :), 'LineWidth', 2); hold on;
% lgnd = legend(year_label{1}, year_label{2:3:end}, year_label{end}, 'FontSize', 18, 'Location', 'best');
% set(lgnd, 'Color', 'none', 'box', 'off');
grid on; axis tight; title('Wages Per Hour', 'FontSize', 18);
clickableLegend(year_label{1}, year_label{2:3:end-1}, year_label{end}, 'Location', 'best');
set(gcf, 'Color', 'w'); set(gcf, 'Position', [0 0 800 600]);
export_fig wageEff.png;
%--------------------------------------------------------------------------

save partA;