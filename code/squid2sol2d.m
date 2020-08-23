function [yy, ry] = squid2sol2d(pars, xdays, INIT, opfrm)
% [yy, ry] = squid2sol2d(pars, xdays, INIT_COND, opfrm)
%
% 2D version of squid2sol function, which merely substitutes 2D appropriate
% functions or values in relevant places in code. For general info see help
% for squid2sol. Extra output ry is a log of detected recovereds who are
% recycled into the susceptible population, which can be used as a corrective
% for long term cumulative case counts.

% Wrapper for ode45 solution to SQID+ model using parameter set pars,
% for days xdays. For a description of the model and expected paramater set
% see the help for function squid2dt.m . Other arguments:
%     xdays: vector if integers giving days to return values for. If given
%           as a scalar then days 1:xdays will be used.
%     INIT: vector of initial population levels to use for solution.
%           If given as scalar > 1 is assumed to be total population, and
%           initial condition will be set to 1 person infected in an
%           otherwise susceptible population.
%
% Return value yy is governed by the last argument opfrm. By default the
% values for all non-susceptibles are returned in a matrix (S can easily be
% recovered by subtracting the sum of these from 1). If opfrm is a scalar
% or vector then the specified column(s) will be returned (where 1
% indicates the S column). IF opfrm is a cell array each column of yy will
% be the sum of the columns specified in the respective cell.

if nargin < 4 || isempty(opfrm)
    opfrm = 2:7;
end
if length(INIT) == 1
    if INIT > 1
        INIT = 1/INIT;
    end
    INIT = [1 - INIT; INIT; zeros(5,1)];
elseif length(INIT) < 7
    INIT = [1 - sum(INIT); INIT(:); zeros(6 - length(INIT),1)];
else
    INIT = INIT(:);
end
if nargout > 1
    INIT = [INIT; 0];
end
if length(xdays) == 1
    xdays = 1:xdays; d1 = 1; d2 = xdays;
else
    d1 = min(xdays); d2 = max(xdays);
end

% Note: our own spec fit fun uses normalized time params, so it can never be used
% for xdays going below 1
if d1 >= 1
    for i = [7,11:2:length(pars)]  % to use wrapper for times <= first xday use time <= 0 normalized by run length
        if pars(i) <= 1
            pars(i) = d1 + pars(i)*(d2 - d1);
        end
    end
end

[~, yall] = ode45(@(t,x) squid2dt2d(t,pars,x), xdays, INIT, odeset('MaxStep',0.5));

if ~iscell(opfrm)
    yy = yall(:,opfrm);
else
    yy = zeros(length(xdays), length(opfrm));
    for i = 1:length(opfrm)
        yy(:,i) = sum(yall(:,opfrm{i}),2);
    end
end
if nargout > 1
    ry = yall(:,end);
end

% to fit
% lsqcurvefit(@(x,xdata) rescale*squid2sol(x, xdata, init_conds, op_fmt), init_parms, #days, rescale*y_data, LB, UB)

end

