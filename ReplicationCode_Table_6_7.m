%==========================================================================
% EMPIRICAL APPLICATION (TABLE 6 & 7)
%==========================================================================
%
%   This script replicates the empirical application from the paper "Comparing
%   predictive ability in presence of instability over a very short time,"
%   by F. Iacone, L. Rossini, and A. Viselli (2024).
%
%   The application compares the median nowcast from the Survey of
%   Professional Forecasters (SPF) against a naive, zero-growth benchmark.
%   The analysis is performed on two distinct sample periods to highlight the
%   effect of the instability:
%   1. The pre-COVID period: 2000:Q1 - 2019:Q4.
%   2. The full sample including the shock: 2000:Q1 - 2020:Q3.
%
%   A suite of predictive ability tests is applied to the squared forecast
%   errors of the two models. These include standard full-sample tests
%   (Diebold-Mariano) as well as tests designed to detect instability
%   (Giacomini-Rossi Fluctuation, Andrews End-of-Sample, and the MAX procedure).
%
%   The script produces two main outputs:
%   1. Two figures visualizing the forecast errors and the nowcasts against
%      the actual data series (corresponding to Figures 1 & 2 in the paper).
%   2. Two tables (Tables 6 & 7) printed to the command window, summarizing
%      the results of the predictive ability tests for the different sample
%      periods and test configurations.
%
%   This script requires the `SPF.mat` data file and the following custom
%   functions to be in the MATLAB path:
%   - dmtest.m
%   - grtest.m
%   - esitest.m
%   - maxtest.m
%
% -------------------------------------------------------------------------

% DATA INITIALIZATION

clc; clear

% Loading the SPF data:
load('SPF.mat');

%   Description of the dataset:

% • NGDP1 is the actual value of GDP for the last period;
% • NGDP2 is the nowcast, i.e. the current-quarter forecast;
% • NGDP3 is the one-step ahead forecast, and so on...

yLag = NGDP1;                               % One-period lagged GDP growth;                                      
ySPF = NGDP2;                               % SPF nowcast;

y = [NGDP1(2:size(NGDP1,1))' NaN(1)]';      % Actual GDP growth;

% We nowcast the growth rate with the SPF median forecast
% with a 0 growth rate using the naive benchmark.

e1 = 100*(y-ySPF)./yLag;            % Error associated to SPF;
e2 = 100*(y-yLag)./yLag;            % Error associated to the naive benchmark;

tStart = 126;                       % 126 corresponds to Q1:2000;

e1 = e1(tStart:size(e1,1)-1);       % Trim pre-2000 observations
e2 = e2(tStart:size(e2,1)-1);       % and drops the last one (NaN);


% ------------------------------------------------------------------

% FIGURES

% Visualization settings:

Time = datetime(YEAR,1,1)+calquarters(QUARTER-1);
Time = Time(tStart:end-1);

legendItem_1 = ["$e_{1t}$, SPF","$e_{2t}$, naive"];
legendItem_2 = ["$\hat{y}_{1t}$, SPF","$\hat{y}_{2t}$, naive","$y_t$, observed"];

linewidth = 3;
fontSizeAxis = 30;
fontSizeLegend = 35;

% Computing the forecasts:

yHat_1 = 100*(ySPF-yLag)./yLag; 
yHat_2 = zeros(size(yHat_1,1),1); 
yAct = 100*(y-yLag)./yLag;

yHat_1 = yHat_1(tStart:size(yHat_1,1)-1); 
yHat_2 = yHat_2(tStart:size(yHat_2,1)-1); 
yAct = yAct(tStart:size(yAct,1)-1);

% Plot:

figure(1)

plot(Time(1:83),e1(1:83),"LineWidth",linewidth);
hold on
plot(Time(1:83),e2(1:83),"LineWidth",linewidth);
grid on
axis padded
set(gca,'FontSize',fontSizeAxis);
legend(legendItem_1,'Location','northwest','Fontsize',fontSizeLegend, ...
    'Interpreter','latex'); 
legend boxoff

figure(2)

plot(Time(1:83),yHat_1(1:83),"LineWidth",linewidth);
hold on
plot(Time(1:83),yHat_2(1:83),"LineWidth",linewidth);
hold on
plot(Time(1:83),yAct(1:83),"k","LineWidth",linewidth,"LineStyle","-.");
grid on
axis padded
set(gca,'FontSize',fontSizeAxis);
legend(legendItem_2,'Location','northwest','Fontsize',fontSizeLegend, ...
    'Interpreter','latex'); 
legend boxoff


% ------------------------------------------------------------------

% NOWCAST EVALUATION 

% Compute the squared residuals:

e1sq = e1.*e1;
e2sq = e2.*e2;

% Function inputs:

n = 80;             % Number of pre-break observations;
nplusone = 81;      % First post-break observation;
nplusm = 83;        % Total number of observations;

% RMSE Ratios:

RMSE_SPF_n = sqrt( mean(e1sq(1:n)) );
RMSE_SPF_nplusm = sqrt( mean(e1sq(1:nplusm)) );

RMSE_NAIVE_n = sqrt(  mean(e2sq(1:n)) );
RMSE_NAIVE_nplusm = sqrt( mean(e2sq(1:nplusm)) );

Ratio_n = RMSE_SPF_n / RMSE_NAIVE_n;
Ratio_nplusm = RMSE_SPF_nplusm / RMSE_NAIVE_nplusm;

% Diebold-Mariano test:

bw_DM = 8;

[DM_n,~] = dmtest(e1sq(1:n)-e2sq(1:n),bw_DM);
[DM_nplusm,~] = dmtest(e1sq(1:nplusm)-e2sq(1:nplusm),bw_DM);

disp("Results of the empirical application:")
disp("DM Pre  DM Post")
disp([DM_n,DM_nplusm])

% Giacomini-Rossi test (k = 0.1):

bw_GR = 3;

m_n = fix(0.1*n);
[GR_n_01,~,~,~] = grtest(e1sq(1:n)-e2sq(1:n),0,m_n,[],bw_GR,"two");

m_nplusm = fix(0.1*nplusm);
[GR_nplusm_01,~,~,~] = grtest(e1sq(1:nplusm)-e2sq(1:nplusm),0,m_nplusm,[],bw_GR,"two");

disp("  GR_l Pre  GR_u Pre  GR_l Post  GR_u Post")
disp([GR_n_01,GR_nplusm_01])

% Giacomini-Rossi test (k = 0.3):

m_n = fix(0.3*n);
[GR_n_03,~,~,~] = grtest(e1sq(1:n)-e2sq(1:n),0,m_n,[],bw_GR,"two");

m_nplusm = fix(0.3*nplusm);
[GR_nplusm_03,~,~,~] = grtest(e1sq(1:nplusm)-e2sq(1:nplusm),0,m_nplusm,[],bw_GR,"two");

disp("  GR_l Pre  GR_u Pre  GR_l Post  GR_u Post")
disp([GR_n_03,GR_nplusm_03])

% Andrews test:

[S_n,q_n,~] = esitest(e1sq(1:nplusm)-e2sq(1:nplusm),[],nplusone,"nplusone");
[S_nplusm,q_nplusm,~] = esitest(e1sq(1:nplusm)-e2sq(1:nplusm),[],nplusone,"nplusm");
[S_eye,q_eye,~] = esitest(e1sq(1:nplusm)-e2sq(1:nplusm),[],nplusone,"eye");

disp("      S_n     q_n     S_nplusm    q_nplusm    S_eye   q_eye   (Q2-Q4:2020)")
disp([S_n/10^3,q_n/10^3,S_nplusm,q_nplusm,S_eye/10^3,q_eye/10^3])

% MAX procedure:

[MAX,MAX_cv,alpha] = maxtest(e1sq(1:nplusm)-e2sq(1:nplusm),n,"T1",nplusm);

disp("      MAX         q   (Q2-Q4:2020)")
disp([MAX^(1/2),MAX_cv^(1/2)])


% ------------------------------------------------------------------

disp("Results of the empirical application.")

% Create Table 6
varNames1 = ["Sample", "RMSE_SPF", "RMSE_NAIVE", "RATIO", "DM", "Fl_L_01", "Fl_U_01", "Fl_L_03", "Fl_U_03"];
Sample = ["Q1:2000 - Q4:2019"; "Q1:2000 - Q3:2020"];
RMSE_SPF = [RMSE_SPF_n; RMSE_SPF_nplusm];
RMSE_NAIVE = [RMSE_NAIVE_n; RMSE_NAIVE_nplusm];
Ratio = [Ratio_n; Ratio_nplusm];
DM = [DM_n; DM_nplusm];
Fl_L_01 = [GR_n_01(1); GR_nplusm_01(1)];
Fl_U_01 = [GR_n_01(2); GR_nplusm_01(2)];
Fl_L_03 = [GR_n_03(1); GR_nplusm_03(1)];
Fl_U_03 = [GR_n_03(2); GR_nplusm_03(2)];

res1 = table(Sample, RMSE_SPF, RMSE_NAIVE, Ratio, DM, Fl_L_01, Fl_U_01, Fl_L_03, Fl_U_03, 'VariableNames', varNames1);
disp(res1)
fprintf("Table 6: Columns RMSE_SPF and RMSE_NAIVE are the average RMSE for the SPF and Naive " + ...
    "benchmark,\nrespectively. Column RATIO refers to the ratio RMSE_SPF/RMSE_NAIVE. Columns " + ...
    "DM , Fl_L, and Fl_U\nare the DM and lower and upper Fl test statistics. The 5%% critical " + ...
    "values for two-sided tests \nare 2.261 for the DM test and 3.393 and 3.01 for the two Fl " + ...
    "tests.\n")

% Create Table 7
varNames2 = ["S(I)", "q_S(I)", "S(Sigma_Tilde)", "q_S(Sigma_Tilde)", "S(Sigma_Hat)", ...
    "q_S(Sigma_Hat)", "MAX", "q_MAX"];

res2 = table(S_eye, q_eye, S_nplusm, q_nplusm, S_n, q_n, MAX, MAX_cv, 'VariableNames', varNames2);
disp(res2)
fprintf("Table 7: Columns S(I), S(eΣ), and S(bΣ), denote the S test statistics when the identity " + ...
    "matrix,\nthe restricted residuals, and the unrestricted residuals are used to weight the errors, " + ...
    "respecti-\nvely. Columns qS(I), qS( eΣ), and qS( eΣ), are the respective critical values. The " + ...
    "theoretical\nsize is 5%%. Column MAX denotes the maximum procedure over the period 2020:Q1 to " + ...
    "2020:Q3, while\ncolumn qMAX denotes the MAX over the period 2000:Q1 to 2019:Q4. The false " + ...
    "positive rate of the\nprocedure is 3.6%%.\n")