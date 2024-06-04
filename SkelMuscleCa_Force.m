
load PSO_25-Apr-2024.mat p pSol yinit

load Resistance_noSOCE6_5-9.mat Ca5 Ca_noSOCE6
load HIIT_5-18.mat Ca7 Ca_noSOCE8

%%
load Force_expt2.mat Ca_expt2 Force_expt2
%% 
Ca_expt = (10 .^(-1 .* Ca_expt2)) ./ (10^ (-6));  %pCa values
Ca_expt(end) = 0.393;
F_expt = Force_expt2 ; %Normalized
Ca_model = Ca5{1}; %pCa conversion
index = find(F_expt > 0.5);
Ca50_index = max(index);
Ca50_expt = Ca_expt(Ca50_index) + 1 ;
%Ca_model = -log(Ca5{1});
InterpExpt = interp1(Ca_expt,F_expt,min(Ca_expt):0.0001:max(Ca_expt));
%% 
initialGuess = [median(Ca_expt), 5,8,0]; % Initial guess for Kd, n and C

lb = [0, 1,0,0];   % Lower bounds 
ub = [max(Ca_expt), 10, min(F_expt), max(Ca_expt)]; % Upper bounds 
F_model = hillEquation(Ca_model, initialGuess);

options = optimoptions('lsqnonlin', 'Display', 'iter', 'Algorithm', 'trust-region-reflective','FunctionTolerance',1e-12);
[best_params, resnorm] = lsqnonlin(@(params) Est(params, Ca_expt, F_expt), initialGuess, lb, ub, options);

Kd_opt = best_params(1) ;
n_opt = best_params(2);
C_opt = best_params(3);
pCa50_opt = -log10(Kd_opt);
x0_opt = best_params(4);

Ca_interp = linspace(min(Ca_model), max(Ca_model), 5000);
F_fit = hillEquation(Ca_interp, [Kd_opt, n_opt, C_opt,x0_opt]);

figure;
semilogx(Ca_expt, F_expt, 'bo', 'DisplayName', 'Experimental Data'); 
hold on;
semilogx(Ca_interp, F_fit, 'r-', 'DisplayName', 'Fitted Model');
xlabel('Calcium Concentration');
ylabel('Force');
title('Fit of Hill Equation to Experimental Data');
legend show;

function F_predicted = hillEquation(Ca, params)
    Fmax = 1; % Assume Fmax is known or estimate it as the maximum observed force
    Kd = params(1);
    n = params(2);
    C = params(3);
    F_predicted = C + (Fmax * (Ca.^n) ./ ((Kd ).^n + Ca.^n));
end

function err = Est(params, Ca, F_expt)
    F_pred = hillEquation(Ca, params);
    err = sum((F_pred - F_expt) .^ 2);
end
