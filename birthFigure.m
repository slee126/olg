close all;
clear;

%Birth Rates inferred from number of 0 year olds to the entire population

pop1939 = xlsread('../data/population.xlsx', '1939_pop');
educ = xlsread('../data/educWB2.xlsx', 'summary');
xls_data = xlsread('../data/population.xlsx', 'census');
xls_data(isnan(xls_data)) = 0;
temp = xls_data(3:5:end, 2);

%1940-2014 Want number of births/total population/25 and over
all_=xls_data(3:5:end, 3:end);
all_ = [transpose(1940:2014)  sum(all_, 2) all_];
all_ = [[1939 pop1939' zeros(1, 25)]; all_];
educ(1, 2:3) = all_(1, 2:3);
% Conferment vs Population growth
birth = educ(1:76, 3)./educ(1:76, 2);

figure(1);
plot(1939:2014, birth);
grid on; axis tight;
title('Birth Rate 1969-2014', 'Interpreter', 'latex', 'FontSize', 18);
set(gcf, 'Color', 'w'); set(gcf, 'Position', [0 0 800 600]);
export_fig ../paper/birthrate.png

