function [MAX,cv,alpha] = maxtest(d,T0,varargin)
%MAXTEST MAX procedure for detecting instabilities in forecast accuracy.
%
%   Syntax:
%   [MAX,cv,alpha] = maxtest(d,T0)
%   [MAX,cv,alpha] = maxtest(d,T0,'Name',Value)
%
%   Description:
%   This function performs the MAX procedure proposed by Harvey et al.
%   (2021) to detect instabilities in predictive ability that occur over a
%   very short sample. This implementation is a specific adaptation
%   (as in Iacone, Rossini, and Viselli, 2024).
%
%   The procedure splits the data into a training period (1 to T0) and a
%   monitoring period (T0+1 to T1). The test statistic (MAX) is the maximum
%   squared loss differential observed during the monitoring period. This is
%   compared against a critical value (cv) derived from the maximum
%   squared loss differential in the training period. A break is signaled
%   if MAX > cv.
%
%   Input Arguments:
%   • d: A T-by-1 numeric vector of the forecast loss differential series.
%   • T0: A scalar integer specifying the end of the training period.
%         Observations from 1 to T0 are used to compute the critical value.
%
%   Optional Name-Value Pair Arguments:
%   • 'alpha': A scalar numeric value for the desired False Positive Rate
%              (FPR) of the procedure. The end of the monitoring period,
%              T1, will be calculated based on this value. This option is
%              mutually exclusive with 'T1'. Default: 0.05.
%   • 'T1': A scalar integer specifying the end of the monitoring period.
%           If provided, the FPR 'alpha' will be calculated from T1. This
%           option is mutually exclusive with 'alpha'.
%
%   Output Arguments:
%   • MAX: The scalar value of the MAX test statistic, calculated as the
%          maximum of the squared loss differential over the monitoring
%          period (T0+1 to T1).
%   • cv: The scalar critical value for the test, calculated as the maximum
%         of the squared loss differential over the training period (1 to T0).
%         The null hypothesis of no break is rejected if MAX > cv.
%   • alpha: The actual False Positive Rate (FPR) of the procedure, which is
%            either the specified input or the value calculated from T1.
%
%   Reference:
%   • Harvey, David I., et al. "Real‐time detection of regimes of
%     predictability in the US equity premium." Journal of Applied
%     Econometrics 36.1 (2021): 45-70.

% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ %

% Optional parameters:
for n = 1:2:size(varargin,2)
    if strcmp('T1',varargin{n})
        if ~exist('alpha','var')
            T1 = varargin{n+1};
            alpha = round( (T1-T0)/T1 ,2);
        else
            error('foo:bar',"ERROR: T1 and alpha cannot be input at the same " + ...
            "time.\n Please either input T1 or alpha.")
        end
    end
    if strcmp('alpha',varargin{n})
        if ~exist('T1','var')
            alpha = varargin{n+1};
            T1 = fix( (T0-alpha)/(1-alpha) );
        else
            error('foo:bar',"ERROR: T1 and alpha cannot be input at the same " + ...
            "time.\n Please either input T1 or alpha.")
        end
    end
end

% Default parameters:
if ~exist('alpha','var') && ~exist('T1','var')
    alpha = 0.05;
    T1 = fix( (T0-alpha)/(1-alpha) );
end

% Warning check:
if T1 > size(d,1)
    T1 = size(d,1);
    alpha = round( (T1-T0-m+1)/(T1-2*m+1) ,2);
    warning('foo:bar',"T1 exceeds the length of l_1:\n T1 is set equal to " + ...
        "the size of l_1 and the correct alpha is calculated. \n If you have " + ...
        "input alpha, it must be smaller or equal than %.2f.", alpha)
end

% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ %


% Training statistics:
dTrain = (d(1:T0)).^2;

% Monitoring statistics:
dMonitor = (d(T0+1:T1)).^2;

% MAX test statistic:
MAX = max(dMonitor);

% Critical value:
cv = max(dTrain);