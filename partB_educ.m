clear;
load partA;
N_period = 2013 - 1976 + 1;
% Time series--------------------------------------------------------------
t_Year = 1977:2013;
year_labelh = cell(37, 1);
for i = 1:N_period-1
    year_labelh{i} = num2str(t_Year(i));
end

age_end = 75;
n_age2 = age_end - 20;

Y_t = educMat_HP(2:end, 2:end);
X_t = educMat_HP(1:end-1, 1:end-1);
Y_t_vec = reshape(Y_t', 37*n_age2, 1);

xAge = eye(n_age2);
xAge = repmat(xAge, 37, 1);
x4 = xAge;
xAge(:, 1) = [];

educHP_mat2 = zeros(37*n_age2, n_age2);

for i = 1:37
    educHP_mat2(n_age2*(i-1)+1:n_age2*i, :) = repmat(educMat_HP(i, 1:n_age2), n_age2, 1);
end

X_ts = [xAge x4.*educHP_mat2];
ts_mdl = LinearModel.fit(X_ts, Y_t_vec);
Y_t_t = reshape(Y_t_vec, n_age2, 37)';


%=========================================================================
% simulated Series
N_fut = 2170-2013;
lastSeries = educMat_HP(end, 1:n_age2);
yProj = zeros(N_fut, n_age2);

lastSeries2 = educMat_HP(end, 1:n_age2);
yProj2 = zeros(N_fut, n_age2);

temp2 = eye(n_age2);
temp3 = temp2(:, 2:end);
for i = 1:N_fut
    temp = repmat(lastSeries, n_age2, 1);
    X_proj = [temp3 temp2.*temp];
    yProj(i, :) = predict(ts_mdl, X_proj);
    lastSeries = [educMat_HP(end, 1) yProj(i, 1:end-1)];
end

t_Year = 2014:(2014+N_fut-1);
year_labelhh = cell(37, 1);
for i = 1:N_fut
    year_labelhh{i} = num2str(t_Year(i));
end

figure(5); clf;
set(5, 'defaulttextinterpreter', 'latex');
plot(21:age_end, yProj(1, 1:end), 'LineWidth', 2); hold on;
for i =2:3:N_fut-1
   plot(21:age_end, yProj(i, 1:end)); hold on;
end
plot(21:age_end, yProj(end, 1:end), '-r', 'LineWidth', 2);
grid on; axis tight; title('Projected \% of Bachelor''s or more', 'FontSize', 18);
clickableLegend(year_labelhh{[1 2:3:N_fut-1 N_fut]}, 'Location', 'best');
set(gcf, 'Color', 'w'); set(gcf, 'Position', [0 0 800 600]);
export_fig educatedFig.png;

%mechanical projection
educProjection2 = zeros(N_fut, n_age2+1);

last_ = educMat_HP(end, :);
enter_ = educMat_HP(end, 1:13);
for i = 1:N_fut
   educProjection2(i, 1:13) = enter_;
   educProjection2(i, 14:end) = last_(13:end-1);
   last_ = educProjection2(i, :);
end


figure(6); clf;
set(6, 'defaulttextinterpreter', 'latex');
plot(20:age_end, educProjection2(1, 1:end), 'LineWidth', 2); hold on;
for i =2:3:N_fut-1
   plot(20:age_end, educProjection2(i, 1:end)); hold on;
end
plot(20:age_end, educProjection2(end, 1:end), '-r', 'LineWidth', 2);
grid on; title('Projected \% of Bachelor''s or more', 'FontSize', 18);
clickableLegend(year_labelhh{[1 2:3:N_fut-1 N_fut]}, 'Location', 'best');
set(gcf, 'Color', 'w'); set(gcf, 'Position', [0 0 800 600]);
export_fig educatedFig2.png;

educProjection1 = [ones(157, 1)*educMat_HP(end, 1) yProj];
save educProjection yProj educProjection1 educProjection2;
%========================================================================