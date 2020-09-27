function Tout = extfitall(T, Tprv, qnew, doconv)
% Tout = extfitall(T_data, T_prevfit, q_new, doconv)
%
% Extends fits for all jurisdictions in the table T_prevfit which are present
% also in table T_data over the entire extent of T_data time series data,
% optionally adding (or removing) q interventions and converting from 1D to
% 2D fits. Arg qnew = number of new interventions to add (if a tuple then
% values are used in starting guess, while arg doconv, if set to false, tells
% function that Tprv param sets already have 2 death rates; default = true,
% since main purpose of function is to be middle stage of 3 pass
% multi-jurisdiction fit. Fits are not redone if R^2 value is low, but
% interventions in param set results are checked for consistency and
% reordered if need be. Projections are not generated, as they optionally
% are for prefitall and refitall functions.

if nargin < 4 || isempty(doconv)
    doconv = true;
end
tmplate = [];
if nargin < 3 || isempty(qnew)
    qnew = 0;    % we expect this to be set, but if not don't assume anything
elseif length(qnew) > 1
    tmplate = qnew; qnew = length(qnew)/2;
end
if qnew == 0 && ~doconv    % nothing to do, pass the input back
    Tout = Tprv; return;
end
if istable(Tprv)
    cn = Tprv.Properties.RowNames; ntmp = T.Properties.RowNames;
    if ~isequal(cn, ntmp)   % following allows one of tables to be subset of other wrt jurisdictions included
        ll = ismember(cn, ntmp);
        if ~any(ll)
            error('Table mismatch between new data and prev results');
        elseif ~all(ll)
            cn = cn(ll);
        end
        T = T(cn,:);
    end
    Tprv = table2cell(Tprv(cn,{'params','init'}));
elseif size(T,1) ~= size(Tprv,1)    % if matrix or cell array rows must correspond to T
    error('Size mismatch for new data and prev results');
else
    cn = T.Properties.RowNames;
end

% set up initial guesses
[m1, k1] = size(Tprv); ci0 = [];
if iscell(Tprv) % extend to previous fits to use as new starting guesses
    cp0 = Tprv(:,1);
    if k1 > 1
        ci0 = Tprv(:,2);
    end
else    % Tprv is matrix
    if mod(k1,2) == 1
        cp0 = num2cell(Tprv, 2);
    else   % even param length means initial vals are included, assume only U_0 value
        cp0 = num2cell(Tprv(:,1:end-1), 2); ci0 = num2cell(Tprv(:,end));
    end
end
if isempty(ci0)
    ci0 = num2cell(1./T.population);
end
n2 = length(T.confirmed{1});
if doconv    % add space for 2nd death rate; starting guess puts time at 1/2 way pt
    for i = 1:m1
        cp0{i} = [cp0{i}(1:6), n2/2, cp0{i}(6:end)];
    end
end
resc = false;
if qnew > 0
    if isempty(tmplate)
        resc = true; tmplate = NaN(1,2*qnew);
    end
    for i = 1:m1
        if resc    % construct default extension: alt q vals 1/5 and -1/3, eq space from last intervention to end of run
            k = (length(cp0{i}) - 9)/2;    % current # of interventions, we don't assume same for each case
            tmplate(1:2:end) = 8*mod(k+1:k+qnew, 2)/15 - 1/3;    % q values
            tmplate(2:2:end) = (n2 - cp0{i}(end))*(1:qnew)/(qnew+1) + cp0{i}(end);
        end
        cp0{i} = [cp0{i}, tmplate];
    end
elseif qnew < 0    % removing interventions (presumably too many generated previously)
    cp0{i} = cp0{i}(1,max(end+2*qnew,9));   % probably safest to just remove from end
end
cpop = num2cell(T.Population); cnn = T{:,1};    % cnn = name of state/country
cyy = cf(@rdivide,cf(@transpose,cf(@cat,1,T.confirmed,T.deaths)),cpop);
cp2 = cell(m1,1); ci2 = cell(m1,1); co2 = cell(m1,1); cr2 = cell(m1,1);

disp(['Extending to day ', num2str(n2)]);
if qnew > 0
    disp(['Adding ', num2str(qnew), ' intervention(s)']);
elseif qnew < 0
    disp(['Removing ', num2str(abs(qnew)), ' intervention(s)']);
else
    disp('Nmber of interventions will be unchanged');
end
if doconv
    disp('Converting 1D to 2D fits');
end
kk = [1:5,9];    % hold these params fixed for now; allow death rates, interventions, and initial conditions to float
for i = 1:m1
    disp([num2str(i), ': ', cnn{i}]); k = length(cp0{i});
    [ctmp, co2{i}] = squidfit2d(cyy{i}, cp0{i}, ci0{i}, kk);
    cp2{i} = ctmp(1:k); ci2{i} = ctmp(k+1:end);
    ytmp = squid2sol2d(cp2{i}, n2, ci2{i}, {4:6,6});
    cr2{i} = r2tst('r2', cyy{i}, ytmp, length(ctmp)-length(kk));
    disp({'resnorm',co2{i}.resnorm;'r2',cr2{i}})
end

% Check interventions and put in proper order, removing -ve or very low 1st
% interventions that would have no effect on results anyways.
for i = 1:m1
    tmp1 = cp2{i}(10:end); k = length(tmp1)/2;
    if any(diff(tmp1(2:2:end)) < 0) || tmp1(1) < 5e-5
        tmp1 = sortrows(reshape(tmp1, 2, k)', 2);
        for j = 1:k
            if tmp1(1,1) < 5e-5    % place null interventions at end of fit period where subsequent refits might have use for them
                tmp1 = [tmp1(2:k,:); 0, (tmp1(k,k) + n2)/2];
            else
                break    % note: if loop goes whole distance there are effectively no interventions left
            end
        end
        cp2{i}(10:end) = reshape(tmp1', 1, 2*k);
    end
end

disp('Jurisdictions with low R^2:')
cat(2, cn(c2m(cr2)<.95), cr2(c2m(cr2)<.95))    % just show for now, will deal with on 3rd pass if need be
Tout = cat(2, cpop, cyy, cp2, ci2, co2, cr2);
vnams = {'pop','ndata','params','init','ofit','R2'};
Tout = cell2table(Tout, 'VariableNames',vnams, 'RowNames',cn);

end
