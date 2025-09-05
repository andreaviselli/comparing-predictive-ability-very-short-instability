%==========================================================================
% MONTE CARLO SIMULATION (TABLE 1)
%==========================================================================
%
%   This script conducts a Monte Carlo simulation to evaluate and compare
%   the empirical size of several tests for equal predictive ability. The
%   simulation is designed to replicate the size analysis (Table 1) from
%   the paper: "Comparing predictive ability in presence of instability
%   over a very short time," by F. Iacone, L. Rossini, and A. Viselli (2024).
%
%   The data generating process (DGP) under the null hypothesis is a mean-zero
%   AR(1) process for the loss differential series.
%
%   The following tests are implemented and evaluated:
%   1. Diebold-Mariano (DM) test with different bandwidths.
%   2. Giacomini-Rossi (GR) Fluctuation test with different window sizes.
%   3. Andrews End-of-Sample (S) Instability test with various settings.
%   4. The MAX test for detecting short-lived instabilities.
%
%   The script calculates the empirical rejection frequency (size) for each
%   test and configuration. The results are compiled and displayed in a
%   table in the command window.
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
phi_vec = [0, 0.25, 0.5];                   % Coefficient value

T = 80;                                     % Sample size
nsim = 10000;                               % MC simulations

m = [1, 3];                                 % Length of the instability (S, MAX)
bw = [3, 8];                                % Bandwidth for the long run variance
k = [0.1, 0.3];                             % Size of the rolling window (Fl)
Flcv = [3.393, 3.012];                      % Critical value for the Fl test

% ---

res = splitvars(table(zeros(3, 10)));
res.Properties.VariableNames = ["phi", "DM3", "DM8", "Fl_k01", ...
    "Fl_k03", "S_m1", "S_eye_m3", "S_tilde_m3", "S_hat_m3", "MAX_m3"];

% ---

% Pre-allocate simulation shocks in advance
eps = randn(T, nsim) * sqrt(sigma2);

row_index = 1;
for i_phi = 1:3

    % Initialize vectors for storage
    pass11 = zeros(nsim,1);                         % DM_bw3
    pass12 = zeros(nsim,1);                         % DM_bw8
    pass21 = zeros(nsim,1);                         % Fl_k01_bw3
    pass22 = zeros(nsim,1);                         % Fl_k03_bw3
    pass3 = zeros(nsim,1);                          % S with m=1
    pass4 = zeros(nsim,1);                          % S with eye(m), m=3
    pass5 = zeros(nsim,1);                          % S with Omega_tilde, m=3
    pass6 = zeros(nsim,1);                          % S with Omega_hat, m=3
    pass7 = zeros(nsim,1);                          % MAX, m=3
    
    % Check simulation case and proceed
    for i_sim = 1:nsim

        % Generate the loss from an AR(1)
        d = zeros(T, 1);
        d(1) = sqrt( 1/(1-phi_vec(i_phi)^2) ) * eps(1, i_sim);

        for t = 2:T
            d(t, :) = phi_vec(i_phi)*d(t-1, :) + eps(t, i_sim);
        end

        % Diebold-Mariano test (bw = 2):
        [DM_1,~] = dmtest(d,bw(1));

        % Critical value for the DM test (fixed smoothing)
        DMcv_1 = 1.96+2.9694*(bw(1)/T)+0.4160*((bw(1)/T)^2)-0.5324*((bw(1)/T)^3);
        
        if abs(DM_1) > DMcv_1
            pass11(i_sim) = 1;
        end

        % Diebold-Mariano test (bw = 8):
        [DM_2,~] = dmtest(d,bw(2));

        % Critical value for the DM test (fixed smoothing)
        DMcv_2 = 1.96+2.9694*(bw(2)/T)+0.4160*((bw(2)/T)^2)-0.5324*((bw(2)/T)^3);
        
        if abs(DM_2) > DMcv_2
            pass12(i_sim) = 1;
        end
        
        % Giacomini-Rossi test (k = 0.1; bw = 2):
        [GR_1,~,~,~] = grtest(d,0,fix(k(1)*T),[],bw(1),[]);

        if any(abs(GR_1) > Flcv(1))
            pass21(i_sim) = 1;
        end

        % Giacomini-Rossi test (k = 0.3; bw = 2):
        [GR_2,~,~,~] = grtest(d,0,fix(k(2)*T),[],bw(1),[]);

        if any(abs(GR_2) > Flcv(2))
            pass22(i_sim) = 1;
        end
        
        % Andrews test (m=1):
        [S_eye,q_eye] = esitest(d,[],T-m(1)+1);
    
        if S_eye > q_eye
            pass3(i_sim) = 1;
        end
        
        % Andrews test (using sigma_hat: stability part of the sample, m=3):
        [S_n,q_n] = esitest(d,[],T-m(2)+1,"nplusone");
    
        if S_n > q_n
            pass6(i_sim) = 1;
        end
    
        % Andrews test (using sigma_tilde, entire sample, m=3):
        [S_nplusm,q_nplusm] = esitest(d,[],T-m(2)+1,"nplusm");
    
        if S_nplusm > q_nplusm
            pass5(i_sim) = 1;
        end
        
        % Andrews test (using the identity matrix, m=3):
        [S_eye,q_eye] = esitest(d,[],T-m(2)+1,"eye");
    
        if S_eye > q_eye
            pass4(i_sim) = 1;
        end
        
        % MAX test:
        if floor(0.95*T) > T-m(2)+1
            error('foo:bar', ['Error (MAX test): "m" is set such that "taustar" ' ...
                'occurs later than "T-m+1".\nDecrease "m" or increase "taustar";' ...
                'but note that the FPR of the test changes too.'])
        end
        [MAX,cv,~] = maxtest(d,floor(0.95*T),"T1",T);
    
        if MAX > cv
            pass7(i_sim)=1;
        end
    end

    % Store the results for each scenario and parameter choice
    res{row_index,:} = [phi_vec(i_phi), mean(pass11), mean(pass12), mean(pass21), ...
        mean(pass22), mean(pass3), mean(pass4), mean(pass5), mean(pass6), mean(pass7)];
    row_index = row_index+1;

end


% ------------------------------------------------------------------

% Create Table 1 (subject to the choices of phi)
format short g
res.Variables = round(res.Variables, 3);
fprintf("Table 1: Results for the MC with sigma2=1.8, T=80; size study.\n")
disp(res)