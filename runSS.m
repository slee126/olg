clear


% 
% % Changes to steady States
% load projectionMats.mat;
% load rates;
% 
% [aoptMat, coptMat, noptMat, pentMat, kbart, nbart, gov_sur] = ...
%     getSS_func2(partProj, wageProj, surv_, pop_, rates);
% 
% 
% year_label = cell(150, 1);
% for i = 1:150
%     year_label{i} = num2str(i+1975);
% end
% 
% figure(1); clf;
% set(1, 'defaulttextinterpreter', 'latex');
% for i = 1:10:150
%    plot(20:100, aoptMat(:, i)); hold on;
% end
% grid on; axis tight; title('Steady State Capital Profile', 'FontSize', 18);
% clickableLegend(year_label{1:10:150}, 'Location', 'best');
% 
% figure(2); clf;
% set(2, 'defaulttextinterpreter', 'latex');
% for i = 1:10:150
%    plot(20:100, coptMat(:, i)); hold on;
% end
% grid on; axis tight; title('Steady State Consumption Profile', 'FontSize', 18);
% clickableLegend(year_label{1:10:150}, 'Location', 'best');
% 
% figure(3); clf;
% set(3, 'defaulttextinterpreter', 'latex');
% for i = 1:10:150
%    plot(20:74, noptMat(1:55, i)); hold on;
% end
% grid on; axis tight; title('Steady State Labor Profile', 'FontSize', 18);
% clickableLegend(year_label{1:10:150}, 'Location', 'best');