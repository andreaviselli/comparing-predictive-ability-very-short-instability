function [sup_F,F,k_a,rej,alt] = grtest(d,h,m,a,b,alt)
%GRTEST Giacomini-Rossi (2010) fluctuation test for predictive ability.
%
%   Syntax:
%   [sup_F,F,k_a,rej,alt] = grtest(d,h,m)
%   [sup_F,F,k_a,rej,alt] = grtest(d,h,m,a,b,alt)
%
%   Description:
%   This function performs the fluctuation test of Giacomini and Rossi
%   (2010) for the null hypothesis of equal "local" unconditional
%   predictive ability between two competing forecasts. The test evaluates
%   the stability of relative forecast performance over time by computing
%   test statistics on centered rolling windows of the data. The null is
%   rejected if the sequence of statistics fluctuates too much, exceeding
%   tabulated critical values.
%
%   Input Arguments:
%   • d: A T-by-1 numeric vector of the loss differential series, defined
%        as d(t) = L(e1(t)) - L(e2(t)), where L is a loss function and
%        e1, e2 are forecast errors from model 1 and model 2, respectively.
%        NaN values are automatically removed.
%   • h: A scalar integer for the forecast horizon (steps ahead). Nowcasts
%        (h=0) and backcasts (h=-1) are treated as one-step-ahead (h=1).
%   • m: A scalar integer specifying the size of the centered rolling
%        windows. Note: 'm' must be less than the effective number of
%        out-of-sample observations.
%   • a: (Optional) A string specifying the significance level of the
%        test. Must be either "0.05" or "0.10". Default is "0.05".
%   • b: (Optional) A scalar integer for the bandwidth of the HAC long-run
%        variance estimator. If empty ([]), a default rule is used:
%        b = floor(T^0.3) + 1.
%   • alt: (Optional) A string specifying the type of alternative
%          hypothesis. Default is "two". Options are:
%          - "one": One-sided test, checks for any significant deviation.
%          - "two": Two-sided test, checks for superiority of one model
%                   over the other.
%
%   Output Arguments:
%   • sup_F: The value of the fluctuation test statistic. For a two-sided
%            test, this is a 1-by-2 vector [min(F), max(F)]. For a
%            one-sided test, it is the scalar max(abs(F)).
%   • F: The full vector of fluctuation test statistics F_{t,m} computed
%        over each rolling window.
%   • k_a: The scalar critical value from Giacomini and Rossi (2010) for
%          the specified significance level 'a' and window size 'm'.
%   • rej: A string indicating the outcome of the test.
%          - "NO": The null hypothesis is not rejected at any point in time.
%          - "YES": For a one-sided test, the null is rejected.
%          - For a two-sided test, where d = L1 - L2:
%            - "YES_1": Model 1 is significantly worse than model 2 at
%                       some point (max(F) > k_a).
%            - "YES_2": Model 2 is significantly worse than model 1 at
%                       some point (min(F) < -k_a).
%            - "YES_12": Both cases occur at different points in time.
%   • alt: A string, "one-sided" or "two-sided", indicating the type of
%          test performed.
%
%   Reference:
%   • Giacomini, Raffaella, and Barbara Rossi. "Forecast comparisons in
%     unstable environments." Journal of Applied Econometrics 25.4 (2010):
%     595-620.

%% Initialization

if nargin < 5
    a = "0.05";
    b = [];
    alt = "two";
end
if nargin < 7
    alt = "two";
end

d = rmmissing(d);                       % Remove missing observations;

T = size(d,1);                         % Out-of-sample # of obs.;
mh = ceil(m/2);                         % Size of the rolling window;
mu = round(m/T,1);                      % Proportion of # of rolling window
                                        % obs. to the # of o.o.s. obs;
if mu > 1
    error("The number of observation in the rolling window is " + ...
        "greater than the number of out-of-sample observations");
end

if isempty(a)                           % Default size of the test;
    a = "0.05";
end

if isempty(b)                           % Rule for the bandwidth of the
    b = floor(max(1,T^0.3))+1;            % var-cov matrix estimator;
end

if isempty(alt)
    alt = "two";                        % Default for the alternative
end                                     % hypothesis;

if h == -1 || h == 0                    % Backcasts and nowcasts are treated
    h = 1;                              % like one-step-ahead forecasts;
end


%% Running the test:

aVar = longrunvariance(d,b);            % Asymptotic (long-run) variance;
F = nan(T-2*mh-h+2,1);                  % Fluctuation test statistics;

for t = h+mh:T-mh+1

    d_t = d((t-mh:t+mh-1));                         % (Rolling) loss differential;
    F(t-h-mh+1) = sum(d_t)/( sqrt(m)*sqrt(aVar) );  % Computing the test statistic;
    
end


%% Results

mu_seq = 0.1:0.1:0.9;                           % Values for mu (read below);
r = mu == mu_seq;                               % Matrix row (read below);
if a == "0.05"; c = 1; else; c = 2; end         % Matrix column;

% The following matrices display the asymptotic critical values, where the
% 1st column refers to the 5% confidence level, whereas the 2nd to the 10%.
% Each row is associated to a peculiar value of mu, see the reference.

if alt == "two"                     % Two-sided test;
    acv = [3.393, 3.170; ...        % Asymptotic critical values;
           3.179, 2.948; ...
           3.012, 2.766; ...
           2.890, 2.626; ...
           2.779, 2.500; ...
           2.634, 2.356; ...
           2.560, 2.252; ...
           2.433, 2.130; ...
           2.248, 1.950];
elseif alt == "one"                 % One-sided test;
    acv = [3.176, 2.928; ...        % Asymptotic critical values;
           2.938, 2.676; ...
           2.770, 2.482; ...
           2.624, 2.334; ...
           2.475, 2.168; ...
           2.352, 2.030; ...
           2.248, 1.904; ...
           2.080, 1.740; ...
           1.975, 1.600];
end

k_a = acv(r,c);                     % Critical value;

res = "";
m1 = "1";                           % Print if model 1 is locally superior;
m2 = "2";                           % Print if model 2 is locally superior;

if alt == "two"
    if max(F) > k_a; res = strcat(res,m1); end
    if min(F) < -k_a; res = strcat(res,m2); end
    if res == ""; rej = "NO"; else; rej = strcat("YES_",res); end
    sup_F = [min(F),max(F)];
elseif alt == "one"
    if max(abs(F)) > k_a; rej = "YES"; else; "NO"; end
    sup_F = [max(abs(F))];
end

alt = strcat(alt,"-sided");         % Type of test (output);