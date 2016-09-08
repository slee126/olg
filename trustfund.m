% Trust Fund -------------------------------------------------------------

T = xlsread('../data/trust_fund.csv');
close all;
figure(1)
set(1, 'defaulttextinterpreter', 'latex');
p1 = plot(T(1:42, 1), T(1:42,2)*100, '-k'); hold on;
plot(T(60:end, 1), T(60:end, 2)*100, '-k'); hold on;
p2 = plot(T(1:13, 5), T(1:13,7), '--k'); hold on;
p3 = plot(T(1:13, 5), T(1:13,6), ':k'); hold on;
plot([2015 2015], [-25 450]); hold on;
p4 = plot(T(1:47, 1), T(1:47,3)*100, ':b'); hold on;
p5 = plot(T(1:62, 1), T(1:62,4)*100, '-b'); hold on;
p6 = plot(T(42:60, 1), T(42:60, 2)*100, '-r'); hold on;
axis([1970 2040 -25 450]);
grid on; axis tight;
title('Asset reserves as \% of annual cost', 'FontSize', 18);
ylabel('%');
set(gcf, 'Color', 'w'); set(gcf, 'Position', [0 0 800 600]);

% Create textbox
annotation('textbox',...
    [0.594656471183015 0.854166666666667 0.066111981799797 0.0538333333333435],...
    'String',{'HIstorical'},...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation('textbox',...
    [0.688500000000001 0.854166666666667 0.066111981799797 0.0538333333333435],...
    'String','Estimated 2015+',...
    'FitBoxToText','off',...
    'EdgeColor','none');

lgnd = legend([p1 p2 p3 p4 p5 p6], 'OASI Intermediate', 'OASI High Cost', 'OASI Low Cost', 'DI', 'HI', 'Baby Boomers', 'Location', 'best');
set(lgnd, 'box', 'off', 'color', 'none');
export_fig '../figures/trustfund.png';

%-------------------------------------------------------------------------


%-------------------------------------------------------------------------
% Immigration
%-----------------------------------------------------------------------
T = xlsread('../data/population.xlsx', 'immigration');
figure(2)
set(2, 'defaulttextinterpreter', 'latex');
plot(T(:, 1), T(:,2)./T(:, 3)*100);
grid on; axis tight;
ylabel('%');
title('immigration \%', 'FontSize', 18);
set(gcf, 'Color', 'w'); set(gcf, 'Position', [0 0 800 600]);
export_fig '../figures/immigration.png';
