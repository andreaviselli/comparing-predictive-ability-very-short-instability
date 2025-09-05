function [S,q,p,S_vec] = esitest(y,X,nplusone,cov,a,ssc)
%ESITEST End-of-sample instability test by Andrews (2003).
%
%   Syntax:
%   [S,q,p,S_vec] = esitest(y,X,nplusone)
%   [S,q,p,S_vec] = esitest(y,X,nplusone,cov,a,ssc)
%
%   Description:
%   This function performs the end-of-sample instability test proposed by
%   Andrews (2003) for a linear regression model. The test is designed to
%   detect parameter instability or outliers occurring at the very end of
%   the sample, accommodating post-change sample sizes that can be
%   arbitrarily small (as small as one).
%
%   The null hypothesis is that the model parameters are stable over the
%   entire sample. The test statistic is constructed from post-change
%   residuals and compared against an empirical distribution of similar
%   statistics generated from the pre-change data. The test is robust to
%   serial correlation, provided the data are strongly stationary and ergodic.
%
%   Input Arguments:
%   • y: A T-by-1 numeric vector of the dependent variable, where T is the
%        total sample size (n+m).
%   • X: A T-by-K numeric matrix of regressors. Do not include a column
%        of ones for the intercept, as it is added automatically. K is the
%        number of regressors.
%   • nplusone: A scalar integer indicating the starting index of the
%               post-change period. The pre-change sample size is thus
%               n = nplusone - 1.
%   • cov: (Optional) A string specifying the covariance matrix estimation
%          method. Default is "nplusone". Options are:
%          - "nplusone": Uses pre-change observations (1 to n) to estimate
%                        the covariance matrix of residuals.
%          - "nplusm": Uses the full sample (1 to n+m).
%          - "eye": Uses the identity matrix (assumes homoskedastic and
%                   serially uncorrelated errors).
%   • a: (Optional) A scalar numeric value for the quantile level used to
%        determine the critical value. Default is 0.95.
%   • ssc: (Optional) A flag (1 or 0) to apply a small-sample correction
%          as suggested by Andrews (2003) to improve the test's size.
%          Default is 1 (correction applied).
%
%   Output Arguments:
%   • S: The scalar value of the end-of-sample test statistic, S_{n+1}.
%   • q: The sample quantile (critical value) from the empirical
%        distribution of subsample statistics, calculated at level 'a'. The
%        null hypothesis is rejected if S > q.
%   • p: The p-value of the test, calculated as the proportion of
%        subsample statistics that are greater than or equal to S.
%   • S_vec: A vector containing the subsample statistics {S_j} for
%            j=1, ..., n-m+1. This vector forms the empirical distribution
%            used for inference.
%
%   Reference:
%   • Andrews, Donald WK. "End‐of‐sample instability tests." Econometrica
%     71.6 (2003): 1661-1694.


%% Initialization
if nargin < 3
    error(['Error: The dependent variable y, regressors X, and', ...
        ' post-change index nplusone are mandatory arguments.', ...
        ' Please provide them.'])
end
if nargin < 4; cov = "nplusone"; a = 0.95; ssc = 1; end
if nargin < 5; a = 0.95; ssc = 1; end
if nargin < 6; ssc = 1; end


%% Testing procedure

nplusm = size(y,1);                         % Full sample size;
n = nplusone - 1;                           % Pre-change sample size;
m = nplusm - n;                             % Post-change sample size;

X = [ones(nplusm,1),X];                     % Design matrix;

Xn = X(1:n,:);
yn = y(1:n);

beta_n_hat = (Xn'*Xn)\(Xn'*yn);             % OLS estimator (pre-change);     
e_u = y - X*beta_n_hat;                     % Residuals (unrestricted model);

beta_nplusm_hat = (X'*X)\(X'*y);            % OLS estimator (full sample);
e = y - X*beta_nplusm_hat;                  % Residuals (restricted model);

if cov == "nplusone"
    nplus = nplusone-m;
elseif cov == "nplusm"
    e_u = e;
    nplus = nplusone;
end

% Sample covariance matrix:
if cov == "nplusone" || cov == "nplusm"
    sigma_sum = zeros(m,m);
    for j = 1:nplus
        e_j = e_u(j:j+m-1);
        sigma_sum = sigma_sum + e_j*e_j';
    end
    sigma_hat = sigma_sum./nplus;
elseif cov == "eye"
    sigma_hat = eye(m);
end

% Subsample statistics S_j:
S_vec = zeros(n-m+1,1);         % Vector of subsample statistics S_j;
d = size(X,2);                  % N. of regressors (including the constant);
if d <= m

    if ssc == 0
        for j = 1:n-m+1
            
            X_j = X(j:j+m-1,:);
            e_j = e(j:j+m-1);
            
            A_j = X_j'/sigma_hat*e_j;
            V_j = X_j'/sigma_hat*X_j;
            S_vec(j) = A_j'/V_j*A_j;        
        end
   
    elseif ssc == 1
        for j = 1:n-m+1
            
            y_j = yn; y_j(j:j+ceil(m/2)-1) = [];
            X_j = Xn; X_j(j:j+ceil(m/2)-1,:) = [];
            beta_2j_hat = (X_j'*X_j)\(X_j'*y_j);
            
            y_j = y(j:j+m-1);
            X_j = X(j:j+m-1,:);
            e_j = y_j - X_j*beta_2j_hat;
            
            A_j = X_j'/sigma_hat*e_j;
            V_j = X_j'/sigma_hat*X_j;
            S_vec(j) = A_j'/V_j*A_j;
        end
    end

elseif d > m

    if ssc == 0
        for j = 1:n-m+1
    
            e_j = e(j:j+m-1);
            S_vec(j) = e_j'/sigma_hat*e_j;        
        end
    
    elseif ssc == 1
        for j = 1:n-m+1

            y_j = y(1:n); y_j(j:j+ceil(m/2)-1) = [];
            X_j = X(1:n); X_j(j:j+ceil(m/2)-1,:) = [];
            beta_2j_hat = (X_j'*X_j)\(X_j'*y_j);
            
            e_j = y - X*beta_2j_hat;
            e_j = e_j(j:j+m-1);

            S_vec(j) = e_j'/sigma_hat*e_j;
        end
    end
end


%% Results

F_S = sort(S_vec);              % Empirical distribution function of S_j;
q = F_S(ceil(a*(n-m+1)));       % Sample quantile (at level a);

e_nplusone = e(nplusone:nplusone+m-1);
X_nplusone = X(nplusone:nplusone+m-1);

% Test statistic, S = S_n+1;
if d <= m
    
    A = X_nplusone'/sigma_hat*e_nplusone;
    V = X_nplusone'/sigma_hat*X_nplusone;
    S = A'/V*A;

elseif d > m
    
    A = e_nplusone;
    V = sigma_hat;
    S = A'/V*A;

end

% Compute the p-value:
p = mean(S <= S_vec);
if p == 0; p = "<e-04"; end