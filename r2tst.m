function r = r2tst(tf, x, y, df, s)
% r = r2tst(test_flag, x, y, df, s)
%
% Goodness of fit tests. The type of test is specified by test_flag:
%   'sse': sum of squared error
%   'chi2': (reduced) Chi^2 statistic
%   'r2': R^2 coefficient of determination
%   'normr': norm of the residuals
%   'stderr': standard deviation of the residuals, normalized by # of pts
%   'me': mean absolute value of error
% The above can also be specified by integer values 0 to 5. Other arguments:
%    x = expected outcomes or frequencies
%    y = observations
%    df = degrees of freedom
%    s = standard deviation of observations
% df and s are only required for the reduced Chi^2, but df can be used as the
% normalizer for tests 2, 4, or 5 (eg. if one wants the sample standard 
% deviation etc). Note: x and y can be any size / dimension as long as they
% agree; only one number is returned in any case. As well, if x or y contain
% NaN's these will be ignored.

% SSE = sum{ (x_i - y_i)^2    note: SS_tot ~ variance eg sum{ (y_i - mean(y))^2
% Chi^2 = weighted SSE sum{ (x_i - y_i)^2/sigma_i    reduced: divide this by df
% R^2 = 1 - SSE/SS_tot    adj: 1 - (1 - R^2)*(n - 1)/(n - df - 1)
% Norm of residuals: sqrt(SSE)
% Standard error: sqrt(SSE/n)    sample adjustment: n => n - df + 1
% Mean abs dev:

if ischar(tf)
    if strncmpi(tf, 'sse', min(length(tf),3))
        tf = 0;
    elseif strncmpi(tf, 'chi2', min(length(tf),4))
        tf = 1;
    elseif strncmpi(tf, 'r2', min(length(tf),2))
        tf = 2;
    elseif strncmpi(tf, 'normr', min(length(tf),5))
        tf = 3;
    elseif strncmpi(tf, 'stderr', min(length(tf),6))
        tf = 4;
    elseif strncmpi(tf, 'me', min(length(tf),2))
        tf = 5;
    else
        error(['Unknown test flag: ' tf]);
    end
elseif tf < 0 || tf > 5
    error(['Unknown test flag: ' num2str(tf)]);
end

if nargin < 4 || isempty(df)
    df = 0;
end
if nargin < 5
    s = [];
end

ii = ~isnan(x) & ~isnan(y); n = nnz(ii);
if tf < 5
    c1 = (x - y).^2;
else
    c1 = abs(x - y);
end

if tf == 1
    if ~isempty(s)
        jj = ii & (~isnan(s) & (s ~= 0));
        c2 = sum(c1(jj) ./ s(jj).^2);
    else
        c2 = sum(c1(ii) ./ abs(x(ii)));
    end
else
    c2 = sum(c1(ii));
end

if tf == 2
    c3 = c2/sum((y(ii) - mean(y(ii))).^2);
elseif tf > 3
    c3 = c2 / (n - df);
else
    c3 = c2;
end

if tf == 1 && df > 0
    r = c3 / df;
elseif tf == 2
    if df > 0
        r = 1 - c3*(n - 1)/(n - df - 1);
    else
        r = 1 - c3;
    end
elseif tf == 3 || tf == 4
    r = sqrt(c3);
else
    r = c3;
end
