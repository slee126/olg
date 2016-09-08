% Get 1940 - 2014 Data and All  Projecting the population
clear;
close all;

xls_data = xlsread('../data/population.xlsx', 'census');
xls_data(isnan(xls_data)) = 0;
temp = xls_data(3:5:end, 2);

%1940-2014
all_=xls_data(3:5:end, 3:end);
all_ = [transpose(1940:2014)  sum(all_, 2) all_];

% Get birth rates
birth_1940_2014 = [transpose(1940:2014) xls_data(3:5:end, 3)./temp];
temp = [transpose(2015:2090), ones(76, 1)*birth_1940_2014(end, 2)];
birth_all = [birth_1940_2014; temp];

% Pseudo birth rate of 21 year olds
temp = sum(all_(:, 23:end), 2);
birth20_1940_2014 = [transpose(1940:2014) xls_data(3:5:end, 23)./temp];
temp = [transpose(2015:2090), ones(76, 1)*birth20_1940_2014(end, 2)];
birth_all20 = [birth20_1940_2014; temp];
birth_all20 = [birth_all20; [transpose(2091:2120) repmat(birth_all20(end, 2), 30, 1)]];

% Get death rates
death_2012_2090 = xlsread('../data/population.xlsx', 'death_prob_future');
death_1900_2011 = xlsread('../data/population.xlsx', 'death_prob_historical');
death_2012_2090 = death_2012_2090(2:end, :);
death_1900_2011 = death_1900_2011(2:end, :);
death_all =[death_1900_2011(41:end, :) ; death_2012_2090];
death_2014 = death_all(1:75, :);
death_all =[death_1900_2011(41:end, :) ; death_2012_2090];

end_year = 2170;
surv_prob_2014 = [transpose(1940:2014) 1 - death_all(1:75, 2:101)];
surv_all = [transpose(1940:2090) 1 - death_all(:, 2:101)];
surv_all = [surv_all; [transpose(2091:end_year) repmat(surv_all(end, 2:end), end_year-2090, 1)]];

%surv_all
surv_all = [surv_all(1:150, :); repmat(surv_all(151, :), 31, 1)];
% Simulation 1940 - 2014
% Get initial Population
bias = 1.0035;
year = 1970;
start_idx2 = year - 1940 + 1;
add_year = end_year-2120+1;
end_year_idx1 = end_year - year + 1;
end_year_idx2 = end_year -1940+ 1;
surv_simul = bias*surv_all(:, 2:end);
len = end_year-year+1;

% actual data--------------------------------------------------------------
actPop = [all_(31:end, 1) sum(all_(31:end, 23:end-1), 2) all_(31:end, 23:end-1)];
actSurv = log(actPop(2:end, 4:end)./actPop(1:end-1, 3:end-1));
actSimul = actPop;

temp = find(actPop(1, :) == 0);
if(~isempty(temp) || ~isempty(temp2))
    if(~isempty(temp))
        start = temp(1);
    elseif(~isempty(temp2))
        start = temp2(1);
    end
    for i = start-1:82
        actSimul(1, i) = actSimul(1, i-1)*surv_all(1, i);
    end
end

simul_pop20 = zeros(len, 80);
simul_pop20(1, :) = actSimul(1, 3:end);
%
% simul_pop20(:, 1) = 1970:1970+150;
% simul_pop20(:, 2) = sum(simul_pop20(:, 3:end), 2)

birth_all20 = [birth_all20; [transpose(2120:end_year), repmat(birth_all20(end, 2), add_year, 1)]];
surv_simul = [surv_simul; repmat(surv_simul(end, :), add_year, 1)];

for i = 2:(end_year-year+1)
    simul_pop20(i, 1) = sum(simul_pop20(i-1, :),2)*birth_all20(start_idx2+i-1, 2);
    simul_pop20(i, 2:end) = simul_pop20(i-1, 1:end-1).*surv_simul(start_idx2+i-1, 20:end-2);
end
simul_pop20 = [transpose(1970:end_year) sum(simul_pop20, 2) simul_pop20];

for h=1:44
    temp = find(actPop(h, :) == 0);
    if(~isempty(temp))
        if(~isempty(temp))
            start = temp(1);
        end
        for i = start-1:82
            actSimul(h, i) = actSimul(h, i-1)*surv_all(h, i+18);
        end
    end    
end


w1 = .8;
w2 = .1;
figure(1)
set(1, 'defaulttextinterpreter', 'latex');
bar(21:100, simul_pop20(1, 3:end)/1000, w1,'FaceColor',[1 1 1]*.75, 'EdgeColor', 'none'); hold on;
bar(21:100, simul_pop20(45, 3:end)/1000, .5,'FaceColor',[1 1 1]*.51, 'EdgeColor', 'none'); hold on;
bar(21:100, simul_pop20(end_year_idx1, 3:end)/1000, .1,'FaceColor',[1 1 1]*.1, 'EdgeColor', 'none'); hold on;
% axis([21 100 0 7e6/1000]);
ylabel('Thousands', 'Interpreter', 'latex');
lgnd = legend('1970', '2014', '2120', 'Location', 'best');
set(lgnd, 'box', 'off', 'color', 'none'); grid on;
title('Simulated Projection to 2120', 'FontSize', 18);
set(gcf, 'Color', 'w'); set(gcf, 'Position', [0 0 800 500]);
export_fig ../paper/popDynamicProj.png;

close all;
ind_Vec = round(linspace(2, 45, 10));
for i = 1:10
    figure(i)
    subplot(121)
    set(i, 'defaulttextinterpreter', 'latex');
    ind_ = ind_Vec(i);
    stem(21:100, actSimul(ind_, 3:end)/1000); hold on;
    bar(21:100, simul_pop20(ind_, 3:end)/1000, .5,'FaceColor',[1 1 1]*.8, 'EdgeColor', 'none'); hold on;
    ylabel('Thousands', 'Interpreter', 'latex');
    lgnd = legend('1970', '2014', 'Location', 'best');
    set(lgnd, 'box', 'off', 'color', 'none'); grid on;
    title(sprintf('Actual Population vs Simulated (%d)', actSimul(ind_, 1)), 'FontSize', 18);
    subplot(122)
    plot(21:100, actSimul(ind_, 3:end)/1000 - simul_pop20(ind_, 3:end)/1000);
    grid on;
    temp = sprintf('Actual-Simulated %d', actSimul(ind_, 1));
    title(temp);
    ylabel('Thousands', 'Interpreter', 'latex');
    set(gcf, 'Color', 'w'); set(gcf, 'Position', [0 0 1400 800]);
    temp = sprintf('../paper/ActualSimulated_%d.png', actSimul(ind_, 1));
    export_fig(temp);
end

workrat = zeros(end_year-year+1, 4);
workrat(:, 1) = 1970:end_year;
workrat(1:45, 2) = sum(actSimul(1:45, 47:end), 2);
workrat(46:end, 2) = sum(simul_pop20(46:end, 47:end), 2);

workrat(1:45, 3) = sum(actSimul(1:45, 3:46), 2);
workrat(46:end, 3) = sum(simul_pop20(46:end, 3:46), 2);

figure(11); clf;
set(11, 'defaulttextinterpreter', 'latex');
workrat(1:45, 4) = workrat(1:45, 2)./workrat(1:45, 3);
plot(1970:2011, workrat(1:42, 4)*100, '-k'); hold on;
workrat(46:end, 4) = sum(simul_pop20(46:end, 47:end), 2)./sum(simul_pop20(46:end, 3:46), 2);
plot(2011:2029, workrat(42:60, 4)*100, '-r'); hold on;
plot(2029:end_year, workrat(60:end, 4)*100, '-k'); 
set(gcf, 'Color', 'w'); set(gcf, 'Position', [0 0 800 600]);
grid on; axis tight;
ylabel('%');
title('Age Dependency Ratio $\frac{65+}{20-64}$', 'FontSize', 18);
export_fig ../paper/age_depend.png

