function [pset, opstr] = squidfit2d(ydata, p1, y1, fixvals, varargin)
% [pset, erval] = squidfit2d(ydata, p1, y1, fixvals, ...)
%
% Wrapper for 2D version of SQUID model fitting. ydata is either country or
% US state data (normalized); p1 gives initial parameter guess, and y1 gives
% initial conditions guess; fixvals is an indicator vector (logical or index)
% specifying parameters to hold fixed (empty if none are). Changes from 1D
% version: no specification of which compartments are fit to (it is assumed
% cumulative cases and deaths); no specification of fitting function (only
% lsqcurvfit is used); data rescaling is done separately for each
% compartment (so that death data is given equal weight even though the
% number of deaths is small in comparison to the number of cases).

[m, n] = size(ydata); assert(n == 2);
xx = 1:m; retvals = {4:6,6};
mresc = min(m/norm(ydata(:,1)), 1e4);                % rescale data so that error values are in
mresc = [mresc, min(m/norm(ydata(:,1)), 20*mresc)];  % ok range for fitter and deaths are not
mresc = diag(mresc); ytarg = ydata*mresc;            % ignored due to low counts

if nargin < 3 || isempty(y1)
    error('SQUIDFIT: Require initial condition');
end

k = length(p1); pp = [p1(:); y1(:)]';
for i = [7, 11:2:k]    % check if pulse and death rate change dates need to be normalized
    if pp(i) > 1
        pp(i) = (pp(i) - 1)/(m - 1);
    end
end
lb = zeros(size(pp)); lb(12:2:k) = -1;    % q_i's after 1st could all release from Q
ub = ones(size(pp)); ub(2) = Inf;   % power law exponent allowed to go above 1
if nargin > 3 && ~isempty(fixvals)
    lb(fixvals) = pp(fixvals); ub(fixvals) = pp(fixvals);
end

useopts = optimoptions(@lsqcurvefit, varargin{:});
[pset, erval, ~, ~, opstr] = lsqcurvefit(@(x,xdata) squid2sol2d(x(1:k),xdata,x(k+1:end),retvals)*mresc, pp,xx,ytarg,lb,ub,useopts);
opstr.resnorm = sqrt(erval);

for i = [7,11:2:k]    % un-normalize pulse etc dates in result
    pset(i) = pset(i)*(m - 1) + 1;
end

end

