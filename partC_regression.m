% % fits===================================================================
close all; clc;
load partA;
load popMat;
load educProjection;

x = eye(n_age);
x = repmat(x, 38, 1);
x2 = x;
x(:, 1) = [];
surv75 = prod(surv_all(37:74, 65:75), 2);
educHP_mat = zeros(n_age*38, n_age);
surv75_mat = zeros(n_age*38, 1);
for i = 1:38
    educHP_mat(n_age*(i-1)+1:n_age*i, :) = repmat(educMat_HP(i, :), n_age, 1);
    surv75_mat(n_age*(i-1)+1:n_age*i) = repmat(surv75(i), n_age, 1);
end
X = [x x2.*educHP_mat surv75_mat];
popVec = simul_pop20(7:44, 3:n_age+2)';
popVec = popVec(:);

%  participation rate ====================================================
w20plus = A(:, 1)./popVec;
w20plus2 = reshape(transpose(workingMat_SG), 2128, 1);
mdl_part = LinearModel.fit(X, w20plus2);

N_fut = 112;

t_Year = 2014:(2014+N_fut-1);
year_labelhh = cell(37, 1);
for i = 1:N_fut
    year_labelhh{i} = num2str(t_Year(i));
end

educSG_mat_proj = zeros(n_age*N_fut, n_age);
for i = 1:N_fut
    educSG_mat_proj(n_age*(i-1)+1:n_age*i, :) = repmat(educProjection1(i, 1:n_age), n_age, 1);
end

x3 = eye(n_age);
x3 = repmat(x3, N_fut, 1);
x4 = x3;
x3(:, 1) = [];

surv75_proj = prod(surv_all(75:end, 65:75), 2);
ind_ = N_fut - length(surv75_proj);
surv75_proj = [surv75_proj; repmat(surv75_proj(end), ind_, 1)];
temp = repmat(surv75_proj', n_age, 1);
LE = temp(:);

X1 = [x3 x4.*educSG_mat_proj LE];
ypredpart = predict(mdl_part, X1);
ypredpartMat = reshape(ypredpart, n_age, N_fut)';
ypredpartMat =  sgolayfilt(ypredpartMat, 3, 21, [], 2);

figure(7); clf;
set(7, 'defaulttextinterpreter', 'latex');
for i = 1:4:40
    plot(20:(20+n_age-1), ypredpartMat(i, :)); hold on;
end
plot(20:(20+n_age-1), ypredpartMat(end, :));
grid on;
axis([20 75 0 .9]);
clickableLegend(year_labelhh{[1:4:40, N_fut]}, 'Location', 'best');
title('Projected Employment \%', 'FontSize', 18);
set(gcf, 'Color', 'w'); set(gcf, 'Position', [0 0 800 600]);
export_fig projectedEmploy.png;

%1976 to 2150
partProj = [workingMat_SG; ypredpartMat];
partProj = partProj(:, 1:end-1);
%=========================================================================
% 
% % shorter age regression
short_ = 11;
n_age2 = n_age - short_;
x = eye(n_age2);
x = repmat(x, 38, 1);
x2 = x;
x(:, 1) = [];
surv75 = prod(surv_all(37:74, 65:75), 2);
educHP_mat = zeros(n_age2*38, n_age2);
surv75_mat = zeros(n_age2*38, 1);
for i = 1:38
    educHP_mat(n_age2*(i-1)+1:n_age2*i, :) = repmat(educMat_HP(i, 1:n_age2), n_age2, 1);
    surv75_mat(n_age2*(i-1)+1:n_age2*i) = repmat(surv75(i), n_age2, 1);
end
X_short = [x x2.*educHP_mat surv75_mat];
popVec_short = simul_pop20(7:44, 3:n_age2+2)';
popVec_short = popVec_short(:);

x_old = eye(short_);
x_old = repmat(x_old, 38, 1);
x_old2 = x_old;
x_old(:, 1) = [];
surv75_2 = prod(surv_all(37:74, 65:75), 2);
educSG_mat_2 = zeros(short_*38, short_);
surv75_mat_2 = zeros(short_*38, 1);
for i = 1:38
    educSG_mat_2(short_*(i-1)+1:short_*i, :) = repmat(educMat_SG(i, end-(short_-1):end), short_, 1);
    surv75_mat_2(short_*(i-1)+1:short_*i) = repmat(surv75(i), short_, 1);
end
X_short_old = [x_old x_old2.*educSG_mat_2 surv75_mat_2];

%piecewise reduction in wages for ages 70 and 75
logwageMat = log(wageMat_HP);
X_rat70 = [logwageMat(:, 46) log(prod(surv_all(37:74, 65:75), 2))];
X_rat75 = [logwageMat(:, 46) log(prod(surv_all(37:74, 65:75), 2))];
Y_rat70 = log(logwageMat(:, 51)./logwageMat(:, 46));
Y_rat75 = log(logwageMat(:, 56)./logwageMat(:, 46));
mdl_rat1 = LinearModel.fit(X_rat70, Y_rat70);
mdl_rat2 = LinearModel.fit(X_rat75, Y_rat75);

%wages
wage20plus2 = reshape(transpose(wageMat_HP), 2128, 1);
mdl_wage = LinearModel.fit(X, log(wage20plus2));
fittedAll = reshape(mdl_wage.Fitted, n_age, 38)';

wage20short = reshape(transpose(wageMat_HP(:, 1:n_age2)), 38*n_age2, 1);
mdl_wage_short = LinearModel.fit(X_short, log(wage20short));
fittedShort = reshape(mdl_wage_short.Fitted, n_age2, 38)';
%========================================================================
% short
N_fut = 112;
educSG_mat_proj = zeros(n_age2*N_fut, n_age2);
for i = 1:N_fut
    educSG_mat_proj(n_age2*(i-1)+1:n_age2*i, :) = repmat(educProjection1(i, 1:n_age2), n_age2, 1);
end

x3 = eye(n_age2);
x3 = repmat(x3, N_fut, 1);
x4 = x3;
x3(:, 1) = [];

surv75_proj = prod(surv_all(75:end, 65:75), 2);
ind_ = N_fut - length(surv75_proj);
surv75_proj = [surv75_proj; repmat(surv75_proj(end), ind_, 1)];
temp = repmat(surv75_proj', n_age2, 1);
LE = temp(:);

X1 = [x3 x4.*educSG_mat_proj LE];
ypredwage = predict(mdl_wage_short, X1);
ypredwageMat = reshape(ypredwage, n_age2, N_fut)';
ypredwageMat =  sgolayfilt(ypredwageMat, 3, 21, [], 2);

% figure(8); clf;
% set(8, 'defaulttextinterpreter', 'latex');
% for i = 1:4:40
%     plot(20:(20+n_age2-1), ypredwageMat(i, :)); hold on;
% end
% plot(20:(20+n_age2-1), ypredwageMat(end, :));
% grid on; axis tight;
% clickableLegend(year_labelhh{[1:4:40, N_fut]}, 'Location', 'best');
% title('Projected Log Wages Per Hour', 'FontSize', 18);
% set(gcf, 'Color', 'w'); set(gcf, 'Position', [0 0 800 600]);


Xrat1 = [ypredwageMat(:, end) log(surv75_proj)];
y70proj = exp(predict(mdl_rat1, Xrat1)).*ypredwageMat(:, end);
Xrat2 = [ypredwageMat(:, end) y70proj log(surv75_proj)];
y75proj = exp(predict(mdl_rat2, Xrat1)).*ypredwageMat(:, end);

ages_ = 65:75;
y_s = [ypredwageMat(:, end) y70proj y75proj];
x_s = [64 70 75];
y_q = zeros(112, 11);
for i = 1:112
    y_q(i, :) = interp1(x_s, y_s(i, :), ages_);
end

y_proj_fin = [ypredwageMat(:, :) y_q(:, :)];
y_proj_fin = sgolayfilt(y_proj_fin, 3, 21, [], 2);
figure(9); clf;
set(9, 'defaulttextinterpreter', 'latex');
for i = 1:10:100
    plot(20:(20+n_age-1), y_proj_fin(i, :)); hold on;
end
plot(20:(20+n_age-1), y_proj_fin(end, :));
grid on; 
axis([20 75 2 3.4]);
clickableLegend(year_labelhh{[1:4:40, N_fut]}, 'Location', 'best');
title('Projected Log Wages Per Hour', 'FontSize', 18);
set(gcf, 'Color', 'w'); set(gcf, 'Position', [0 0 800 600]);
export_fig projectedWages;

wageProj = [log(wageMat_SG); y_proj_fin];
wageProj = wageProj(:, 1:end-1);

% % final matrices needed.....Survival, Productivity and Population 
%1976-2125
surv_all(:, 1) = 1970:2150;
surv_ = surv_all(7:156, 22:end-1);
pop_ = simul_pop20(7:156, 3:end);

save projectionMats wageProj partProj surv_ pop_;

