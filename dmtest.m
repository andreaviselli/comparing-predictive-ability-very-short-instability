function [DM,pvalue] = dmtest(d,b)
%DMTEST Diebold-Mariano test for equal predictive accuracy.
%
%   Description:
%   This function performs the Diebold-Mariano (DM) test for the null
%   hypothesis of unconditional equal predictive accuracy between two
%   competing forecasts. The test is based on the sample mean of a loss
%   differential series.
%
%   The null hypothesis is H0: E[d(t)] = 0, where d(t) is the loss
%   differential. Under H0, the DM statistic is asymptotically N(0,1).
%   The null is rejected at a significance level 'alpha' if the returned
%   p-value is less than 'alpha'.
%
%   Input:
%   • d: A T-by-1 numeric vector of the loss differential series.
%        The loss differential is defined as d(t) = L(e1(t)) - L(e2(t)),
%        where L(.) is a loss function and e1, e2 are forecast errors.
%        NaN values are automatically removed.
%   • b: A scalar integer specifying the bandwidth, or lag truncation
%        parameter, for the long-run variance estimation. For an h-step
%        ahead forecast, a common choice is b = h. If b is provided as
%        empty ([]), a default rule is used: b = floor(max(1, T^0.3)) + 1,
%        where T is the effective sample size.
%
%   Output:
%   • DM: The scalar value of the Diebold-Mariano test statistic.
%   • pvalue: The two-sided p-value for the test, computed using the
%             standard normal distribution.
%
%   Reference:
%   • Diebold, Francis X., and Robert S. Mariano. "Comparing predictive
%     accuracy". Journal of Business & economic statistics 20.1 (2002): 134-144.


%% Run the test

% Remove Nan values;
d = rmmissing(d);

T = size(d,1);
if isempty(b)                       % Cut-off rule;
    b = floor(max(1,T^0.3))+1;      % Bandwidth;
end

% Compute the DM statistic:
aVar = longrunvariance(d,b);
DM = sqrt(T)*mean(d)/sqrt(aVar);

% Compute the p-value:
pvalue = 2*normcdf(-abs(DM));      % Two-sided test;