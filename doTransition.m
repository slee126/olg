clear;
close all;
clc;
tic;

load projectionMats;

% wageProj = exp(wageProj);
% wageProj = wageProj./repmat(mean(wageProj, 2), 1, 55);

surv_ = doDiags(surv_');
prod = doDiags(wageProj');
partProj = doDiags(partProj');

% tax rates and replacement rates
rates = csvread('../data/rates.csv');

%fixed parameters---------------------------------------------------------
G = 80;                         % Death Certainty
gama = 2;                       % risk aversion
eta = 2;                        % elasticity of labor supply
phi = 4.3;                      % disutility labor
alfa = .40;                     % share of capital
growth = 1.018;
bita = .99*growth^(-gama);      % discount factor
R = 55;
R_ss = 45;

% inline marginal utility
u_c = @(c) c.^(-gama);
r = 1/bita - 1;
delta = .07;

% minimum labor required
min_labor = .1;

% first year and last year of loop
year1 = 1976;
year2 = 2125;
N = 150;
NN = 150 + 79*2;

% 308 factor prices...this will have to be staggered for r and w

aopt = zeros(G + 1, 1);             
copt = zeros(G + 1, 1);            
nopt = .2*ones(G + 1, 1);
nopt(R+1:end) = 0;      


kbart = zeros(NN, 1);
nbart = zeros(NN, 1);
critt = zeros(N, 1);
gov_sur = zeros(NN, 1);

aoptMat = zeros(G + 1, NN);             
coptMat = zeros(G + 1, NN);            
noptMat = .2*ones(G + 1, NN);
noptMat(R+1:end, :) = 0; 

tax = [rates(1, 2)*ones(79, 1); rates(:, 2); rates(end, 2)*ones(190, 1)];
replacement = [rates(1, 3)*ones(79, 1); rates(:, 3); rates(end, 3)*ones(190, 1)];
%------------------------------------------------------------------------

% Loop parameters
crit_outer = 1; 
crit_main = 1;
crit_old = 1;
crit_diff = 1;
tol_outer = 1e-2;
tol_inner = 5e-2;
tol_main = 1e-2;
nq1 = 10000;
damp2 = .1;
damp1 = .1;
damp3 = .1;
iter_main = 0;

% 3 Loops
while(crit_main > tol_main)
    % needs to be wrapped in broylin update
    
    for h = 1:len_
        
        ind1 = h;
        ind2 = h+79;
        wt = w_mat(ind1:ind2);
        rt = r_mat(ind1:ind2);
        wt = (1 - alfa)*(kbart./nbart).^alfa;
        rt = alfa*(kbart./nbart).^(alfa - 1) - delta;
        if(h < G + 1)
            w = wt(1);
            r = rt(1);
            pen = pent(1);
        else
            w = wt(h - G + 1);
            r = rt(h - G + 1);     
            annual_income = w*(1-tau)*nopt(1:R_ss);
            pent(h - G + 1) = replacement*mean(annual_income);
            pen = pent(h - G + 1);
        end
        
        mass = mass_mat(:, h);
        surv_prob = survD(:, h);
        prod_vec = prodD(:, h);

        k_last = zeros(nq1, 1);
        k_first = zeros(nq1, 1);

        damp1 = .1;
        q1 = 0;
        crit_inner = 1;
        while (crit_inner > tol_inner);
            % iteration over q1: 
            if (q1 > 100)
                damp1 = .01;
            end;
            if (crit_main < .03 || iter_main > 8)
                tol_inner = 1e-3;
            end

            q1 = q1 + 1;
            if q1 == 1;
                kR = 0.1;
            elseif q1 == 2;
                kR = 1;
            else
            kR = k_last(q1 - 1) - damp1*((k_last(q1 - 1) - k_last(q1 - 2))/...
                (k_first(q1 - 1)- k_first(q1 - 2))*k_first(q1 - 1));
            end;

            aopt(G) = kR;
            k_last(q1) = kR;

            % decision rule for the retired
            for i = G-1 :-1:(R + 1);                   
                k2 = aopt(i + 2);
                k1 = aopt(i + 1);
                r0 = rt(i);
                r1 = rt(i+1);
                r2 = rt(i+2);
                % From Budget Constraint
                c1 = (1 + r2)*k1 - growth*k2 + pen;     
                % From Euler
                c0 = (u_c(c1)*surv_prob(i + 1)*bita*(1 + r1))^(1/-gama);
                k0 = (c0 + growth*k1 - pen)/(r0 + 1);
                aopt(i) = k0;
                copt(i) = c0;
                copt(i + 1) = c1;
            end;    

            % decision rule for the working and SS
            if R > R_ss
                for i = R:-1:R_ss;              
                  k2 = aopt(i + 2);
                  k1 = aopt(i + 1);
                  n1 = nopt(i + 1);
                    
                  if i == R
                    c1 = (1 + r)*k1  - growth*k2 + pen;
                  else
                    c1 = (1 + r)*k1 + (1 - tau)*w*n1*prod_vec(i + 1) - growth*k2 + pen;
                  end
                  %solve c0 and l0 ---> Euler/MRS
                  c0 = (u_c(c1)*surv_prob(i + 1)*bita*(1 + r))^(1/-gama);
                  % leisure/cons tradeoff
                  l0 = ((1 - tau)*w*prod_vec(i)/phi*c0^-gama)^(-1/eta);
                  n0 = 1 - l0;
                  if n0 < min_labor
                      n0 = 0;
                  end
                  k0 = (c0 + growth*k1 + n0*w*prod_vec(i)*(tau - 1) - pen)/(r + 1);
                  if c1 < 0
                    c1 = 0;
                  end;
                  aopt(i) = k0;
                  nopt(i) = n0;
                  copt(i) = c0;
                  copt(i + 1) = c1;
                end;
            end

            % decision rule for the working
            for i = R_ss-1:-1:1;              
                k2 = aopt(i + 2);
                k1 = aopt(i + 1);
                if i == R_ss - 1
                    n1 = nopt(i + 1);
                    c1 = (1 + r)*k1 + (1 - tau)*w*prod_vec(i + 1)*n1 - growth*k2 + pen; 
                else
                    n1 = nopt(i + 1);
                    c1 = (1 + r)*k1 + (1 - tau)*w*n1*prod_vec(i + 1) - growth*k2;
                end

                %solve c0 and l0 ---> Euler/MRS
                c0 = (u_c(c1)*surv_prob(i + 1)*bita*(1 + r))^(1/-gama);
                % leisure/cons tradeoff
                l0 = ((1 - tau)*w*prod_vec(i)/phi*c0^-gama)^(-1/eta);
                n0 = 1 - l0;
                if n0 < min_labor
                   n0 = 0; 
                end
                k0 = (c0 + growth*k1 + n0*w*prod_vec(i)*(tau - 1))/(r + 1);
                
                aopt(i) = k0;
                nopt(i) = n0;
                copt(i) = c0;
                copt(i + 1) = c1;
            end;

            k_first(q1) = aopt(1);     % stores all the first period k
            crit_inner = abs(aopt(1));
            fprintf('crit_inner: %.8f iter: %g  year: %g mainIter: %g main_crit: %.10f diff: %.10f\n', ...
                crit_inner, q1, h, iter_main, crit_main, crit_diff);
            if q1 > 2000
                break;
            end
        end;    
        aopt_all(:, h) = aopt;
        nopt_all(:, h) = nopt;
        copt_all(:, h) = copt; 
    end

    aopt_fullMat = [repmat(aoptMat(1:80, 1), 1, 79), aoptMat(1:end-1, :), ...
        repmat(aoptMat(1:80, end), 1, 79)];
    nopt_fullMat = [repmat(noptMat(1:55, 1), 1, 54), noptMat(1:55, :), ...
        repmat(noptMat(1:55, end), 1, 54)];
    
    aopt_diagMat = zeros(80, 229);
    nopt_diagMat = zeros(55, 229);
    
    start_ = 79;
    for i = 1:229
        ind1 = i;
        ind2 = 79+i;
        ind2_lab = 54+i;
        aopt_diagMat(:, i) = diag(aopt_fullMat(:, ind1:ind2));
        if i < 205
            nopt_diagMat(:, i) = diag(nopt_fullMat(:, ind1:ind2_lab));
        end
    end
    nopt_diagMat(:, 205:end) = repmat(nopt_diagMat(:, 204), 1, 25);
    

    for i = 0:len_- 1
        aoptD(:, i + 1) = diag(aoptMat, i);
        noptD(:, i + 1) = dcnbciag(noptMat, i);
        coptD(:, i + 1) = diag(coptMat, i);
    end

    for i = G:len_
        kbart(i - G + 1) =  mass_mat(1:G-1, i)'*(aoptD(2:G, i));
        nbart(i - G + 1) = mass_mat(1:R, i)'*(noptD(1:R, i).*prodD(:, i));
%         ret_pop(i - G +1) = sum(mass_mat(R_ss:end, i));
%         ann_income(i - G + 1) = mean(wt(i - G + 1)*(1 - tau)*noptD(1:R, i).*prodD(:, i) ...
%             + rt(i - G + 1)*aoptMat(1:R, i)); 
%         pent(i - G + 1) = .1*ann_income(i - G + 1);
        
        gov_income = noptD(1:R, i)'*(pop_matD(1:R, i).* prodD(:, i))*wt(i - G + 1)*tau/1e6;
        if (i-G+1 > 1)
            gov_surplus(i - G + 1) = gov_income - pent(i - G + 1)*sum(pop_matD(R+1:end, i)) + gov_surplus(i-G);
        else
            gov_surplus(i - G + 1) = gov_income - pent(i - G + 1)*sum(pop_matD(R+1:end, i));
        end
    end
    
    kbart = kbart*damp2 + kbart_old*(1 - damp2);
    nbart = nbart*damp2 + nbart_old*(1 - damp2);

    crit_main = max(abs(kbart - kbart_old));
    crit_diff = crit_main - crit_old;
    crit_old = crit_main;
    iter_main = iter_main + 1;
    
    kbart_old = kbart;
    if (crit_diff > 0)
        damp2 = damp2/2;
    end;
%     if iter_main > 10
%         damp2 = .02;
%     elseif iter_main > 15
%         damp2 = .01;
%     elseif iter_main > 20
%         damp2 = .005;
%     end
end