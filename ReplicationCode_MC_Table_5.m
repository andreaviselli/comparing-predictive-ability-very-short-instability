%==========================================================================
% MONTE CARLO SIMULATION (TABLE 5)
%==========================================================================
%
%   This script conducts a Monte Carlo simulation to evaluate the empirical
%   power of several predictive ability tests when the loss differential
%   exhibits serial correlation. The simulation is designed to replicate
%   the specific power analysis (Table 5) from the paper: "Comparing 
%   predictive ability in presence of instability over a very short time," 
%   by F. Iacone, L. Rossini, and A. Viselli (2024).
%
%   The Data Generating Process (DGP) under the alternative hypothesis is
%   an AR(1) process with a moderate degree of serial correlation (phi=0.5).
%   Instability is introduced through a constant mean shift ('mu') throughout
%   the sample, combined with an additional, large mean shift ('delta') in
%   the final 'm=3' observations. 
%
%   1. Diebold-Mariano (DM) test.
%   2. Giacomini-Rossi (GR) Fluctuation test.
%   3. Andrews End-of-Sample (S) Instability test.
%   4. The MAX test for detecting short-lived instabilities.
%
%   The script calculates the empirical rejection frequency (power) for each
%   test under this specific DGP. The results are compiled and displayed in a
%   table in the command window, corresponding to Table 5 of the paper.
%
%   This script requires the following custom functions in the MATLAB path:
%   - dmtest.m
%   - grtest.m
%   - esitest.m
%   - maxtest.m
%   - longrunvariance.m
%
% -------------------------------------------------------------------------

clear; clc
rng(0)

% --- MC Settings ---

sigma2 = 1.8;                               % Innovation variance
phi = 0.5;                                  % Coefficient value

T = 80;                                     % Sample size
S = 20;                                     % Post-break sample size (for Fl test)
nsim = 10000;                               % MC simulations

mu = -1.3;                                  % Differential for DM and Fl
delta = (-53.3+1.3);                        % End-of-sample differential
m = 3;                                      % Length of the instability (S)

bw = [3, 8];                                % Bandwidth for the long run variance
k = [0.1, 0.3];                             % Size of the rolling window (Fl)
Flcv = [3.393, 3.012];                      % Critical value for the Fl test

% ---

scenarios = "end";
res = splitvars(table(zeros(1, 9)));

res.Properties.VariableNames = ["-mu", "-delta", "DM8", ...
    "Fl_k01", "Fl_k03", "S_eye", "S_tilde", "S_hat", "MAX"];

% ---

% Pre-allocate simulation shocks in advance
eps = randn(T, nsim) * sqrt(sigma2);

row_index = 1;
for i_m = 1:length(m)
    for i_scenario = 1:length(scenarios)
        for i_mu = 1:length(mu)
            for i_delta = 1:length(delta)
                
                % Initialize mu_current, delta_current, and scenario
                mu_current = mu(i_mu);
                delta_current = delta(i_delta);
                scenario = scenarios(i_scenario);
                
                % Cases for simulations
                cond1 = scenario == "end";
                
                % Initialize vectors for storage
                pass12 = zeros(nsim,1);                         % DM_bw8
                pass21 = zeros(nsim,1);                         % Fl_k01_bw3
                pass22 = zeros(nsim,1);                         % Fl_k03_bw3
                pass3 = zeros(nsim,1);                          % S with eye(m)
                pass4 = zeros(nsim,1);                          % S with Omega_tilde
                pass5 = zeros(nsim,1);                          % S with Omega_hat
                pass6 = zeros(nsim,1);                          % MAX
                
                % Check simulation case and proceed
                if cond1
                    for i_sim = 1:nsim
    
                        % Generate the loss from an AR(1)
                        d = zeros(T, 1);
                        d(1) = sqrt( 1/(1-phi^2) ) * eps(1, i_sim);
        
                        for t = 2:T
                            d(t, :) = phi*d(t-1, :) + eps(t, i_sim);
                        end

                        switch scenario
                            case "const"
                            % CASE 1: delta_s = mu_current for all s
                                d = mu_current + d;
                            case "end"
                            % CASE 3: delta_s = mu_current for s = T-m+1
                                d = mu_current + d;
                                d(T-m(i_m)+1:T)= delta_current + d(T-m(i_m)+1:T);
                        end
    
                        % Diebold-Mariano test (bw = 8):
                        [DM_2,~] = dmtest(d,bw(2));
    
                        % Critical value for the DM test (fixed smoothing)
                        DMcv_2 = 1.96+2.9694*(bw(2)/T)+0.4160*((bw(2)/T)^2)-0.5324*((bw(2)/T)^3);
                        
                        if abs(DM_2) > DMcv_2
                            pass12(i_sim) = 1;
                        end
                        
                        % Giacomini-Rossi test (k = 0.1; bw = 3):
                        [GR_1,~,~,~] = grtest(d,0,fix(k(1)*T),[],bw(1),[]);
    
                        if any(abs(GR_1) > Flcv(1))
                            pass21(i_sim) = 1;
                        end
    
                        % Giacomini-Rossi test (k = 0.3; bw = 3):
                        [GR_2,~,~,~] = grtest(d,0,fix(k(2)*T),[],bw(1),[]);
    
                        if any(abs(GR_2) > Flcv(2))
                            pass22(i_sim) = 1;
                        end
                        
                        % Andrews test (using sigma_hat: stability part of the sample):
                        [S_n,q_n] = esitest(d,[],T-m(i_m)+1,"nplusone");
                    
                        if S_n > q_n
                            pass5(i_sim) = 1;
                        end
                    
                        % Andrews test (using sigma_tilde, entire sample):
                        [S_nplusm,q_nplusm] = esitest(d,[],T-m(i_m)+1,"nplusm");
                    
                        if S_nplusm > q_nplusm
                            pass4(i_sim) = 1;
                        end
                        
                        % Andrews test (using the identity matrix):
                        [S_eye,q_eye] = esitest(d,[],T-m(i_m)+1,"eye");
                    
                        if S_eye > q_eye
                            pass3(i_sim) = 1;
                        end
                        
                        % MAX test:
                        if floor(0.95*T) > T-m(i_m)+1
                            error('foo:bar', ['Error (MAX test): "m" is set such that "taustar" ' ...
                                'occurs later than "T-m+1".\nDecrease "m" or increase "taustar";' ...
                                'but note that the FPR of the test changes too.'])
                        end
                        [MAX,cv,~] = maxtest(d,floor(0.95*T),"T1",T);
                    
                        if MAX > cv
                            pass6(i_sim)=1;
                        end
                    end
    
                    % Store the results for each scenario and parameter choice
                    res{row_index,:} = [abs(mu_current), abs(delta_current), mean(pass12), ...
                        mean(pass21), mean(pass22),mean(pass3), mean(pass4), ...
                        mean(pass5), mean(pass6)];
                    row_index = row_index+1;
    
                end
            end
        end
    end
end

% ------------------------------------------------------------------

% Create Table 5
format short g
res.Variables = round(res.Variables, 3);
fprintf("Table 5: Results for the MC with phi=0.5, sigma2=1.8, T=80, m=3.\n")
disp(res)