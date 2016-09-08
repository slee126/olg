%-----------------------------------------------------------------------
%survival probability
%source: https://www.ssa.gov/OACT/HistEst/DeathProbabilities2015.html
%produces survival.pdf
%-----------------------------------------------------------------------
clear;
close all;
%load historical
surv_proj = xlsread('../data/population.xlsx', 'death_prob_historical');
surv_proj = [surv_proj(61:end, 1), (1 -surv_proj(61:end, 2:end))];
prob85 = zeros(11, 1);
figure(1)
set(figure(1), 'defaulttextinterpreter', 'latex');
temp_ = '';
for i = 1:5:52
    if (i == 51)
        subplot(121)
        h(i)= plot(21:65, surv_proj(i, 23:67), '-r'); hold on;
        subplot(122)
        plot(65:100, surv_proj(i, 67:102), '-r'); hold on;
    else
        subplot(121)
        h(i) = plot(21:65, surv_proj(i, 23:67), 'color', [0 0 0] + 1 - .07*ceil(i/5)); hold on;
        subplot(122)
        plot(65:100, surv_proj(i, 67:102), 'color', [0 0 0] + 1 - .07*ceil(i/5)); hold on;
    end
    prob85(ceil(i/5)) = prod(surv_proj(i, 67:87));
    temp_ = sprintf('%s %g = %.3f\n', temp_, surv_proj(i, 1), prob85(ceil(i/5)));
end;
subplot(121)
grid on; axis tight
title('Survival Prob 21-65', 'Interpreter', 'latex', 'FontSize', 14);
xlabel('age');
lgnd = legend([h(1) h(6) h(11) h(16) h(21) h(26) h(31) h(36) h(41) h(46) h(51)], {'1960', '1965', '1970', ...
    '1975', '1980', '1985', '1990', '1995', '2000', '2005', '2010'}, 'Location', 'best');
set(lgnd, 'box', 'off', 'color', 'none');
subplot(122)
grid on; axis tight; xlabel('age');
set(gcf, 'Color', 'w'); set(gcf, 'Position', [0 0 800 500]);
title('Survival Prob 26-85', 'Interpreter', 'latex', 'FontSize', 14);
% Create textbox

temp_ = temp_(2: end);
annotation(figure(1),'textbox',...
    [0.573 .055410107705054 0.189769230769231 0.62121789560893],...
    'String',{'Probability to live up to 85', temp_},...
    'Interpreter','latex',...
    'FitBoxToText','off',...
    'EdgeColor','none');
figtitle('Historical', 'Interpreter', 'latex', 'FontSize', 20);
export_fig ../paper/survival_1.png;



temp_ = '';
surv_proj = xlsread('../data/population.xlsx', 'death_prob_future');
surv_proj = [surv_proj(:, 1), (1 -surv_proj(:, 2:end))];
prob85 = zeros(11, 1);
figure(1)
set(figure(2), 'defaulttextinterpreter', 'latex');
for i = 5:7:79
    if (i == 75)
        subplot(121)
        h(i)= plot(21:65, surv_proj(i, 23:67), '-r'); hold on;
        subplot(122)
        plot(65:100, surv_proj(i, 67:102), '-r'); hold on;
    else
        subplot(121)
        h(i) = plot(21:65, surv_proj(i, 23:67), 'color', [0 0 0] + 1 - .07*ceil(i/7)); hold on;
        subplot(122)
        plot(65:100, surv_proj(i, 67:102), 'color', [0 0 0] + 1 - .07*ceil(i/7)); hold on;
    end
    prob85(ceil(i/5)) = prod(surv_proj(i, 67:87));
    temp_ = sprintf('%s %g = %.3f\n', temp_, surv_proj(i, 1), prob85(ceil(i/5)));
end;
subplot(121)
grid on; axis tight
title('Survival Prob 21-65', 'Interpreter', 'latex', 'FontSize', 14);
xlabel('age');
lgnd = legend([h(5) h(12) h(19) h(26) h(33) h(40) h(47) h(54) h(61) h(68) h(75)], {'2016', '2023', '2030', ...
    '2037', '2044', '2051', '2058', '2065', '2072', '2079', '2086'}, 'Location', 'best');
set(lgnd, 'box', 'off', 'color', 'none');
subplot(122)
grid on; axis tight; xlabel('age');
set(gcf, 'Color', 'w'); set(gcf, 'Position', [0 0 800 500]);
title('Survival Prob 26-85', 'Interpreter', 'latex', 'FontSize', 14);
% Create textbox

temp_ = temp_(2: end);
annotation(figure(2),'textbox',...
    [0.573 .055410107705054 0.189769230769231 0.62121789560893],...
    'String',{'Probability to live up to 85', temp_},...
    'Interpreter','latex',...
    'FitBoxToText','off',...
    'EdgeColor','none');
figtitle('Projection', 'Interpreter', 'latex', 'FontSize', 20)
export_fig ../paper/survival_2.png;

