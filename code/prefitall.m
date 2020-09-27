function Tout = prefitall(T, n1, nq, doproj, p1, init1, rchk)
% Tout = prefitall(T, n1, nq, doproj, p1, init1, rchk)
%
% Fit empirical covid data in table T using the original 1D versions of
% SQUIDER code over a specified range of the data. Main purpose of function
% is to do prefits of data in T for subsequent extension and refitting with
% a bifurcated death rate. Additional arguments:
%     n1 = number of days to fit data to (default = 101);
%     nq = number of q interventions to generate (default = 1);
%     doproj = if non-negative add projected compartment values to results
%         table. By default these are not included at this stage.
%     p1 = starting guess for parameters to use. See comments for defaults.
%     init1 = starting guess for initial conditions. See comments for defaults.
%     rchk = R^2 test value for accepting fit. Def = 0.95 ...

%{
 EXTRA INFO:
rchk - If R^2 value for a fit is not high enough then fit is redone using
    different starting guesses until success, or all alternative starting
    guesses are exhausted; in the latter case the fit with the best R^2 is
    returned.
doproj - If a positive value is given than an ODE run of doproj days is
    generated, returning all compartments except S (= 1 - total) and added
    to the table. If logical true is given as an argument then a run of the
    same period as the fit data is generated.
p1 - For general default values, see code (ll. 109..). If nq is given and 
    p1 has length 2*nq then it assumed to specify default interventions
    only and general defaults are used for the rates.
   - To do a general 1D refit p1 can be given as a table, cell array, or 
    matrix with the same number of rows as T. In this case the default n1 
    will match the number of days of data in T, and any given nq values will 
    be ignored (use extfitall function to adjust number of interventions).
   - Initial conditions included in tables, m x 2 cell arrays, or even width 
    matrices will override defaults or given values.
   - To convert 2D results use the input {'convert', p1}, which will refit
    with the mean death rate as starting guess.
init1 - By default the starting U_0 values (which is also a fit parameter) 
    is set to 1 over the population of the state/country in question, and 
    all other compartments (except S) are 0. Use this if you need to set Q
    or other bins to nonzero values at the start.
%}

doconv = false;
if nargin < 2
    n1 = []; % def = 101
end
if nargin < 3
    nq = []; % def = 1
end
if nargin < 4
    doproj = [];
end
if nargin < 5
    p1 = []; % default: see below
end
if nargin < 6
    init1 = []; % default: see below
end
if nargin < 7 || isempty(rchk)
    rchk = 0.95;
end
cn = T.Properties.RowNames; m1 = length(cn); cnn = T{:,1};   % 2nd used for consol output (could be State or Country)
n2 = length(T.confirmed{1}); cpop = num2cell(T.Population);

% set up starting guess(es)
cp0 = cell(size(cn)); ci0 = cell(size(cn)); % starting param & initial condition
defp1 = []; defq1 = [];
if ~isempty(p1)
    if iscell(p1) && (numel(p1) == 2 && isequal(p1{1},'convert'))
        doconv = true; p1 = p1{2};  % conversion from 2D
    end
    [m0, k0] = size(p1);
    if m0 > 1           % refit of previously fitted data. We assume that rows in p1 correspond
        if m0 ~= m1     % to rows of T. As well, any nq value given will be ignored ...
            error('Starting parameter guesses: size mismatch');
        end
        if istable(p1)    % fit results table such as returned by this or refitall
            init1 = p1.init; p1 = p1.params;
            if isempty(doproj)    % if refitting from table, do projection if prev did and not specified otherwise
                kk = find(contains(p1.Properties.VariableNames, 'proj'));
                if length(kk) == 1
                    doproj = cellfun('size',p1{1,kk},1);
                end
            end
        end
        if iscell(p1)    % cell array of starting guesses. First col = params only
            cp0 = p1(:,1);
            if k0 >= 2  % ow use default initial conditions
                init1 = p1(:,2);
            end
        else   % matrix given, may include init conds
            if mod(k0, 2) == 0  % matrix includes init conds, assume last col = U_0
                init1 = p1(:,k0); p1 = p1(:,1:k0-1);
            end
            cp0 = num2cell(p1, 2);
        end
        if isempty(n1)    % if n1 not specifically given assume we want full refit
            n1 = n2;
        end
        nq = [];    % since now presumably irrelevant
    else   % single tuple, using one starting guess
        if ~isempty(nq) && k0/2 == nq
            defq1 = p1;
        elseif k0 < 7
            error('Incorrect format for starting param guess');
        elseif k0 == 7
            defp1 = p1;
        else
            if mod(k0, 2) == 0
                init1 = p1(k0); p1 = p1(1:k0-1);
            end
            defp1 = p1(1:7); defq1 = p1(8:end); nq = length(defq1)/2;
        end
    end
end
if ~isempty(init1)
    if iscell(init1)
        ci0 = init1;
    elseif size(init1,1) == m1
        ci0 = num2cell(init1, 2);
    else
        ci0(:) = {init1};
    end
end
if isempty(n1)
    n1 = 101;
end
if isempty(doproj)
    doproj = 0;
elseif isequal(doproj, true)
    doproj = n1;
end
ll = cellfun(@isempty, cp0);
if any(ll)
    if isempty(nq)
        nq = 1;
    end
    if isempty(defp1)
        defp1 = [3/4, 0.99, 1/8, 1/2, 1/2, 0.02, 0.01];
    end
    if isempty(defq1)
        defq1 = zeros(1,2*nq); defq1(1:4:end) = 1/5; defq1(3:4:end) = -1/3;
        defq1(2:2:end) = n1*(2:(nq + 1))/(nq + 2);
    end
    cp0(ll) = {[defp1,defq1]};
end
ll = cellfun(@isempty, ci0);
if any(ll)
    ci0(ll) = num2cell(1./cell2mat(cpop(ll)));
end
if doconv
    for i = 1:m1
        cp0{i} = [cp0{i}(1:5),mean(cp0{i}(6:2:8)),cp0{i}(9:end)];
    end
end

disp(['Fitting to day ', num2str(n1)]);
if ~isempty(nq)
    disp(['Generating ', num2str(nq), ' intervention(s)']);
elseif doconv
    disp('Converting from 2D fits');
else
    disp('Refitting on prior results');
end
disp(['Fit R^2 target: ', num2str(rchk)]);
if doproj > 0
    disp('Including projections');
end

cyy = cell(m1,1); cp1 = cell(m1,1); ci1 = cell(m1,1); co1 = cell(m1,1); cr1 = cell(m1,1);
for i = 1:m1
    disp([num2str(i), ': ', cnn{i}]);
    cyy{i} = [T.confirmed{i}(1:n1); T.deaths{i}(1:n1)]'/cpop{i};
    [ctmp, co1{i}] = squidfit(cyy{i}, cp0{i}, ci0{i});
    k = length(cp0{i}); cp1{i} = ctmp(1:k); ci1{i} = ctmp(k+1:end);
    ytmp = squid2sol(cp1{i}, n1, ci1{i}, {4:6,6});
    cr1{i} = r2tst('r2', cyy{i}, ytmp, length(ctmp));
    disp({'resnorm',co1{i}.resnorm;'r2',cr1{i}})
end

disp('Checking R^2 vals...');
ll = c2m(cr1) < rchk; kk = find(ll);
if isempty(kk)
    disp('... all good.');
elseif all(ll)   % if refit version of this works, use same response here ...
    disp('Something went wrong, exiting.');
    Tout = cat(2, cn, cpop, cyy, cp1, ci1, co1, cr1);
    return;
else
    disp('...');
    ttmp = sortrows(cat(2,cn,cf(@getfield,co1,'resnorm'),cr1,num2cell(1:length(cn))'),[3,2,1,4]);
    disp(flipud(ttmp(1:length(kk)+1,:)));

    tmpx = complex(T.Long, T.Lat);
    tdist = abs(bsxfun(@minus, tmpx, tmpx(ll).'));
    tried = false(size(tdist));
    while any(ll)
        for i = 1:length(kk)
            if ~ll(kk(i))
                continue
            end
            jj = find(~ll & ~tried(:,i));
            if isempty(jj)
                disp('Problem with refits')
                Tout = {cn, cpop, cyy, cp1, co1, cr1, ll, kk, i, tdist, tried};
                return;
            end
            [~, jk] = min(tdist(jj,i)); tried(jj(jk),i) = true;
            disp(['Refitting ' cn{kk(i)} ' based on ' cn{jj(jk)}]);
            [ptmp, otmp] = squidfit(cyy{kk(i)}, cp1{jj(jk)}, ci1{jj(jk)});
            ytmp = squid2sol(ptmp(1:end-1), n1, ptmp(end), {4:6,6});
            rtmp = r2tst('r2', cyy{kk(i)}, ytmp, length(ptmp));
            disp({'resnorm',otmp.resnorm;'r2',rtmp})
            if rtmp > cr1{kk(i)}
                k = length(cp1{jj(jk)});
                cp1{kk(i)} = ptmp(1:k); ci1{kk(i)} = ptmp(k+1:end);
                co1{kk(i)} = otmp; cr1{kk(i)} = rtmp;
                if rtmp >= rchk
                    disp('success')
                    ll(kk(i)) = false;
                end
            end
        end
    end
end
Tout = cat(2,cpop,cyy,cp1,ci1,co1,cr1);
vnams = {'pop','ndata','params','init','ofit','R2'};
if doproj > 0
    cy1 = cell(m1,1); n3 = doproj;
    disp(['Calculating projections for ' num2str(n3), ' days']);
    for i = 1:m1
        cy1{i} = squid2sol(cp1{i}, n3, ci1{i});
    end
    Tout = cat(2, Tout, cy1);
    if doproj > 284
        vnams = cat(2, vnams, ['proj',num2str(round(n3/365)), 'y']);
    elseif doproj > 14
        vnams = cat(2, vnams, ['proj',num2str(round(n3/30)), 'm']);
    else
        vnams = cat(2, vnams, ['proj',num2str(n3), 'd']);
    end
end
Tout = cell2table(Tout,'VariableNames',vnams,'RowNames',cn);

end
