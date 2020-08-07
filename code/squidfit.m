function [pset, opstr] = squidfit(ydata, p1, y1, fixvals, usefun, varargin)
% [pset, erval] = squidfit(ydata, p1, y1, fixvals, usefun, ...)
%
% Wrapper for SQUID model fitting. ydata is either country or US state
% data (normalized); p1 gives initial parameter guess, and y1 gives initial
% conditions; fixvals is an indicator vector (logical or index) specifying
% parameters to hold fixed (empty if none are). Arg usefun specifies which
% fitting function to use (0 = lsqcurvefit, default; 1 = fminsearch);
% subsequent arguments are sent to the fitting routine.

[m, n] = size(ydata);
if n == 3
    retvals = {4:6,5,6};
elseif n == 2
    retvals = {4:6,6};
else
    retvals = {4:6};
end
xx = 1:m;
resc = (m*n)/norm(ydata);
ytarg = resc*ydata;

if nargin < 3 || isempty(y1)
    error('SQUIDFIT: Require initial condition');
end

k = length(p1); pp = [p1(:); y1(:)]';
for i = 9:2:k    % check if pulse dates need to be normalized
    if pp(i) > 1
        pp(i) = (pp(i) - 1)/(m - 1);
    end
end
lb = zeros(size(pp)); lb(8:2:k) = -1;
ub = ones(size(pp));
if nargin > 3 && ~isempty(fixvals)
    lb(fixvals) = pp(fixvals);
    ub(fixvals) = pp(fixvals);
end

if nargin < 6 || usefun == 0
    useopts = optimoptions(@lsqcurvefit, varargin{:});
    [pset, erval, ~, ~, opstr] = lsqcurvefit(@(x,xdata) resc*squid2sol(x(1:k),xdata,x(k+1:end),retvals), pp,xx,ytarg,lb,ub,useopts);
    opstr.resnorm = sqrt(erval);
elseif usefun == 1
    useopts = optimset(varargin{:});
    [pset, erval, eflag, opstr] = fminsearch(@(x) norm(resc*squid2sol(min(max(x(1:k),lb(1:k)),ub(1:k)),xx,min(max(x(k+1:end),lb(k+1:end)),ub(k+1:end)),retvals) - ytarg),pp,useopts);
    pset = min(max(pset, lb), ub);
    opstr.exitflag = eflag; opstr.resnorm = erval;
end
for i = 9:2:k    % un-normalize pulse dates in result
    pset(i) = pset(i)*(m - 1) + 1;
end

end

