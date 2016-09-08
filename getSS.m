clear;
close all;
clc;
tic;

%fixed parameters---------------------------------------------------------
G = 80;                     % Death Certainty
gama = 2;                   % Risk Aversion
eta = 2;                    % elasticity of labor supply
phi = 4.3;                  % Disutility Labor
alfa = .40;                 % share of capital
growth = 1.018;
bita = .99*growth^(-gama); % discount factor
R = 55;                     % Mandatory Retirement Age 75
R_ss = 45;                  % SS age

% inline marginal utility
u_c = @(c) c.^(-gama);
r = 1/bita - 1;
delta = .07;

% minimum labor required
min_labor = .1;

year1 = 1976;
year2 = 2125;
N = 150;

aopt = zeros(G + 1, 1);             
copt = zeros(G + 1, 1);            
nopt = .2*ones(G + 1, 1);
nopt(R+1:end) = 0;                            

% load 4 projections: Participation, productivity, survival, and mass
load projectionMats;





% tax rates and replacement rates
rates = csvread('../data/rates.csv');

kbart = zeros(N, 1);
nbart = zeros(N, 1);
critt = zeros(N, 1);
gov_sur = zeros(N, 1);
workers = zeros(N, 1);
work_pop = zeros(N, 1);

aoptMat = zeros(G + 1, N);             
coptMat = zeros(G + 1, N);            
noptMat = .2*ones(G + 1, N);
noptMat(R+1:end, :) = 0; 
pentMat = zeros(2, 1);

yearMat = [1976; 2125];
tax = [rates(:, 2); rates(end, 2)*ones(111, 1)];
replacement = [rates(:, 3); rates(end, 3)*ones(111, 1)];
retpop = zeros(2, 1);

for h = 1:1
    tau = tax(h);
    replace = replacement(h);
    
    surv_prob = surv_(h, :);
    %computation of mass 
    mass =  pop_(h, :)./sum(pop_(h, :));
   
    retpop(h) = sum(mass(R_ss:end));    
    work_pop(h) = sum(mass(1:R));

    participation = partProj(h, :)';
    prod_vec = wageProj(h, :)'./mean(wageProj(h, :));
    
    iter = 0;
    crit_outer = 1; 
    tol_outer = 5e-3;
    tol_inner = 5e-2;
    nq1 = 10000;
    damp = .1;
    crit_old = 2;
    iter_outer = 0;

    % first guess
    nbar = .22;
    kbar = (alfa/(r + delta))^(1/(1 - alfa))*nbar;

    while (crit_outer > tol_outer)
        iter_outer = iter_outer + 1;
        w = (1 - alfa)*(kbar/nbar)^alfa;
        r = alfa*(kbar/nbar)^(alfa - 1) - delta;
        ret_pop = sum(mass(R_ss:end));
        annual_income = w*(1-tau)*nopt(1:R_ss);
        pen = replace*mean(annual_income);

        k_last = zeros(nq1, 1);
        k_first = zeros(nq1, 1);

        damp1 = .1;
        q1 = 0;
        crit_inner = 1;
        while (crit_inner > tol_inner);
            % iteration over q1: 
            if (q1 > 100)
                damp1 = .005;
            end;
            if (crit_outer < .01)
                tol_inner = 5e-3;
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

                % From Budget Constraint
                c1 = (1 + r)*k1 - growth*k2 + pen;     
                % From Euler
                c0 = (u_c(c1)*surv_prob(i)*bita*(1 + r))^(1/-gama);
                k0 = (c0 + growth*k1 - pen)/(r + 1);
                
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
                    c1 = (1 + r)*k1 + (1 - tau)*w*n1*prod_vec(i + 1)*participation(i + 1) - growth*k2 + pen;
                  end
                  %solve c0 and l0 ---> Euler/MRS
                  
                  %non linear solver for c and n if k0<0
                  c0 = (u_c(c1)*surv_prob(i)*bita*(1 + r))^(1/-gama);
                  
                  
                  % leisure/cons tradeoff
                  l0 = ((1 - tau)*w*prod_vec(i)/phi*c0^-gama)^(-1/eta);
                  lopt(i) = l0;
                  n0 = 1 - l0;
                  if n0 < min_labor
                      n0 = 0;
                  end
                  k0 = (c0 + growth*k1 + participation(i)*n0*w*prod_vec(i)*(tau - 1) - pen)/(r + 1);

                  if c0 < 0
                      c0 = 0;
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
                    c1 = (1 + r)*k1 + (1 - tau)*w*prod_vec(i + 1)*n1*participation(i+1) - growth*k2 + pen; 
                else
                    n1 = nopt(i + 1);
                    c1 = (1 + r)*k1 + (1 - tau)*w*n1*prod_vec(i + 1)*participation(i+1) - growth*k2;
                end

                %solve c0 and l0 ---> Euler/MRS
                c0 = (u_c(c1)*surv_prob(i)*bita*(1 + r))^(1/-gama);
                % leisure/cons tradeoff
                l0 = ((1 - tau)*w*prod_vec(i)/phi*c0^-gama)^(-1/eta);
                lopt(i) = l0;
                n0 = 1 - l0;
                if n0 < 0
                   n0 = 0; 
                end
                k0 = (c0 + growth*k1 + participation(i)*n0*w*prod_vec(i)*(tau - 1))/(r + 1);

                aopt(i) = k0;
                nopt(i) = n0;
                copt(i) = c0;
                copt(i + 1) = c1;
            end;
           k_first(q1) = aopt(1);     % stores all the first period k
           crit_inner = abs(aopt(1));
           fprintf('crit_outer: %.8f crit_inner: %.8f \t iter: %d outer: %d year: %g\n', crit_outer, abs(aopt(1)), q1, iter_outer, h);
           if q1 > 2000
               break;
           end
        end;    
        if(crit_outer < .01)
            damp = .01;
        end;

        knew = mass(1:G-1)*aopt(2:G); 
        crit_outer = abs((knew - kbar)/kbar);
        kbar = (1 - damp)*kbar + damp*knew;
        temp = mass(1:R).*participation';
        
        nnew = temp*(nopt(1:R).*prod_vec);
        nbar = (1 - damp)*nbar + damp*nnew;
        
        if(iter_outer > 100)
            break;
        end
    end;
        aoptMat(:, h) = aopt;
        coptMat(:, h) = copt;
        noptMat(:, h) = nopt;
        pentMat(h) = pen;
        critt(h) = crit_outer;
        kbart(h) = kbar;
        nbart(h ) = nbar;
        
        temp = mass(1:R).*participation';
        gov_income = (nopt(1:R).*prod_vec)'*temp'*w*tau;
        gov_out = pen*(1 - sum(mass(1:R_ss)));
        gov_sur(h) = gov_income - gov_out;
        temp = find(nopt ==  0);
        workers(h) = temp(1);
end

toc

% for i = 1:10:150
%     figure(1)
%     if h == i
%         plot(21:75, noptMat(1:55, i), '-r'); hold on;
%     else
%         plot(21:75, noptMat(1:55, i)); hold on;
%     end
%     grid on;
% end
% % toc;
% save guess1 kbart nbart gov_sur workers aoptMat coptMat noptMat pentMat