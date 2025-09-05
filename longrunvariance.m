function Omega_hat = longrunvariance(X,h)
%LONGRUNVARIANCE Newey-West estimator of the long-run covariance matrix.
%
%   Syntax:
%   Omega_hat = longrunvariance(X,h)
%
%   Description:
%   This function computes the Newey-West (1987) Heteroskedasticity and
%   Autocorrelation Consistent (HAC) estimator of the long-run covariance
%   matrix for a multivariate time series. This estimator is robust to
%   both general forms of heteroskedasticity and autocorrelation.
%
%   The estimator uses the Bartlett kernel for weighting the autocovariance
%   terms and is guaranteed to be positive semi-definite by construction.
%   The bandwidth parameter 'h' determines the number of autocovariance
%   terms included in the estimation.
%
%   Input Arguments:
%   • X: A T-by-N numeric matrix representing a multivariate time series,
%        where T is the number of observations (time periods) and N is the
%        number of variables (series). The series are demeaned internally.
%   • h: A scalar integer specifying the bandwidth (lag truncation
%        parameter). This determines the number of autocovariances used in
%        the computation. This argument is mandatory. For guidance on
%        choosing h, refer to Newey and West (1994).
%
%   Output Arguments:
%   • Omega_hat: The N-by-N estimated long-run covariance matrix.
%
%   References:
%   • Newey, Whitney K., and Kenneth D. West. "A simple, positive
%     semidefinite, heteroskedasticity and autocorrelation consistent
%     covariance matrix." Econometrica 55.3 (1987): 703-708.
%   • Newey, Whitney K., and Kenneth D. West. "Automatic lag selection
%     in covariance matrix estimation." The Review of Economic Studies
%     61.4 (1994): 631-653.


X = X';                                 % Transpose for better
[N,T] = size(X);                        % handling time-series;

X = X - ones(N,T)*mean(X,2);            % De-meaned series;
Omega_hat = (1/T)*(X*X');               % Gamma_0 (Variance);

if h > 0
   for j = 1:h
        
       X_0 = X(:,1+j:T);
       X_j = X(:,1:T-j);

       Gamma_hat_j = (1/T)*(X_0*X_j');  % Gamma_j;
       w = 1 - (j/h);                   % Bartlett kernel;
       
       % Long-run variance:
       Omega_hat = Omega_hat + w*(Gamma_hat_j+Gamma_hat_j');

   end
end