function Tout = refitall(T, Tprv, doproj, rchk)
% Tout = refitall(T, Tprv, doproj, rchk)
%
% Refits all jurisdictions in the table of fit results T_prv which are also
% present in table T_data, over the entire extent of T_data time series. No
% interventions are added or removed (use extfit for this), and Tprv fit
% data is assumed to be already in 2D format. Main use is for the last leg
% of a 3 pass multijurisdiction fit calculation. Additional arguments:
% doproj, if false, suppresses the 2 year (from 1st day of available data)
% projection that is added to the result table by default; for other time
% periods set doproj to the desired number of days. rchk gives the desired
% value for the measured coefficient of determination R^2 used as a
% criterion for keeping a fit; if test fails, fit is retried using
% alternate starting guesses. By default fit is kept if mean of separate
% R^2's for confirmed cases and deaths is 0.95 (see comments for more).

% Set up criteria for keeping fit. Options are mean (def), min (ie both
% >=), max (at least one), either confirmed or deaths alone, or compare
% entire matrix of results in single R^2 op (like in prefitall function  etc)
rfunc = @mean; fstr = 'mean';
if nargin < 4 || isempty(rchk)
    rchk = 0.95;
elseif iscell(rchk)    % use alternative to mean for R^2 test combo
    if isequal(rchk{1}, 'min')
        rfunc = @min; fstr = 'min';
    elseif isequal(rchk{1}, 'max')
        rfunc = @max; fstr = 'max';
    elseif isequal(rchk{1},1) || isequal(rchk{1}(1:min(end,4)),'conf')
        rfunc = @(x) x(1); fstr = 'confirmed';    % just check R^2 for confirmed cases
    elseif isequal(rchk{1},2) || isequal(rchk{1},'deaths')
        rfunc = @(x) x(2); fstr = 'deaths';    % just check R^2 for deaths
    elseif isequal(rchk{1}(1:min(end,4)),'comb')
        rfunc = []; fstr = 'combined';    % revert to combined check used in prefitall etc.
    elseif ~isequal(rchk{1},'mean')
        error(['Unknown argument: ', rchk{1}]);
    end
    if length(rchk) > 1
        rchk = rchk{2};
    else
        rchk = 0.95;
    end
end
if (nargin < 3 || isempty(doproj)) || isequal(doproj, true)
    doproj = 730;
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
elseif size(T,1) ~= size(Tprv,1)    % if matrix or cell array rows must correspond
    error('Size mismatch for new data and prev results');
else
    cn = T.Properties.RowNames;
end

% set up initial guesses; if init conds missing give error
[m1, k1] = size(Tprv); ci0 = [];
if iscell(Tprv) % extend to previous fits to use as new starting guesses
    if k1 > 1
        cp0 = Tprv(:,1); ci0 = Tprv(:,2);
    end
elseif mod(k1,2) == 0   % Tprv is matrix, initial conds must be last col
    cp0 = num2cell(Tprv(:,1:end-1), 2); ci0 = num2cell(Tprv(:,end));
end
if isempty(ci0)
    error('Require initial values for refits');
end
n2 = length(T.confirmed{1});

% cn = T.Properties.RowNames; v1 = T.Properties.VariableNames{1};
cpop = num2cell(T.Population); cnn = T{:,1};
cyy = cf(@rdivide,cf(@transpose,cf(@cat,1,T.confirmed,T.deaths)),cpop);
cp3 = cell(m1,1); ci3 = cell(m1,1); co3 = cell(m1,1); cr3 = cell(m1,1);
rcmp = cell(m1,1);

disp(['Fitting to day ', num2str(n2)]);
disp(['Criterion for keeping fits: ', fstr, ' R^2 at least ', num2str(rchk)]);
if ~doproj || doproj < 0
    disp('Projections will not be added to table');
else
    disp(['Adding projections for ', num2str(doproj), ' days from beginning of data']);
end
for i = 1:m1
    disp([num2str(i), ': ', cnn{i}]);
    [ctmp, co3{i}] = squidfit2d(cyy{i}, cp0{i}, ci0{i});
    k = length(cp0{i}); cp3{i} = ctmp(1:k); ci3{i} = ctmp(k+1:end);
    ytmp = squid2sol2d(cp3{i}, n2, ci3{i}, {4:6,6});
    if ~isempty(rfunc)
        k = length(ctmp);
        cr3{i} = [r2tst('r2',cyy{i}(:,1),ytmp(:,1),k),r2tst('r2',cyy{i}(:,2),ytmp(:,2),k)];
        rcmp{i} = rfunc(cr3{i});
    else
        cr3{i} = r2tst('r2', cyy{i}, ytmp, length(ctmp));
        rcmp{i} = cr3{i};
    end
    disp({'resnorm',co3{i}.resnorm;'r2',rcmp{i}})
end

disp('Checking R^2 vals...');
ll = c2m(rcmp) < rchk; kk = find(ll);
if isempty(kk)
    disp('... all good.');
elseif all(ll)
    disp('Something might have gone wrong, all fits have low R^2, exiting.');
else   % Idea: refit using as initial guess params based on successful fits of other jurisdictions.
    disp('...');
    ttmp = sortrows(cat(2,cn,cf(@getfield,co3,'resnorm'),rcmp,num2cell(1:length(cn))'),[3,2,1,4]);
    disp(flipud(ttmp(1:length(kk)+1,:)));

    tmpx = complex(T.Long, T.Lat);   % Starting guesses selected in order by distance from target jurisdiction
    tdist = abs(bsxfun(@minus, tmpx, tmpx(ll).'));
    tried = false(size(tdist));
    while any(ll)
        for i = 1:length(kk)
            if ~ll(kk(i))
                continue
            end
            jj = find(~ll & ~tried(:,i));
            if isempty(jj)
                if nnz(~tried(~ll,ll(kk))) > 0
                    continue;   % still a chance for other jurisdictions, so carry on
                else
                    disp('Some fits still have low R^2, returning best results')
                    disp(cat(2, cn(ll), rcmp(ll))); ll(ll) = false;    % clear ll so we can exit outer loop
                    break;
                end
            end
            [~, jk] = min(tdist(jj,i)); tried(jj(jk),i) = true;
            disp(['Refitting ' cn{kk(i)} ' based on ' cn{jj(jk)}]);
            [ptmp, otmp] = squidfit2d(cyy{kk(i)}, cp3{jj(jk)}, ci3{jj(jk)});
            k = length(cp3{jj(jk)});
            ytmp = squid2sol2d(ptmp(1:k), n2, ptmp(k+1:end), {4:6,6});
            if ~isempty(rfunc)
                r2mp = [r2tst('r2',cyy{kk(i)}(:,1),ytmp(:,1),length(ctmp)),r2tst('r2',cyy{kk(i)}(:,2),ytmp(:,2),length(ctmp))];
                rtmp = rfunc(r2mp);
            else
                r2mp = r2tst('r2', cyy{kk(i)}, ytmp, length(ctmp));
                rtmp = r2mp;
            end
            disp({'resnorm',otmp.resnorm;'r2',rtmp})
            if rtmp > rcmp{kk(i)}   % save best so far, in any case
                cp3{kk(i)} = ptmp(1:k); ci3{kk(i)} = ptmp(k+1:end);
                co3{kk(i)} = otmp; cr3{kk(i)} = r2mp; rcmp{kk(i)} = rtmp;
                if rtmp >= rchk
                    disp('success')
                    ll(kk(i)) = false;
                end
            end
        end
    end
end

Tout = cat(2,cpop, cyy, cp3, ci3, co3, cr3);
vnams = {'pop','ndata','params','init','ofit','R2'};
if doproj > 0
    cy3 = cell(m1,1); n3 = doproj;
    for i = 1:m1
        cy3{i} = squid2sol2d(cp3{i}, n3, ci3{i});
    end
    Tout = cat(2, Tout, cy3);
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
