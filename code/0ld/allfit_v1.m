function Tout = allfit(T, nq, n1, p1)
% Tout = allfit(T, nq, n1, p1)
%
% Combines all fit ops that are done and redone in the notes..txt files
% into a single function for simplicity. T = covid data table, nq = number
% of q lockdowns/relaxations to apply, n1 = length of initial fit period
% (default = 101 days), p1 = init guess for parameters; if a duple, this is
% just a guess for the 1st q pulse, and default vals are used for the other
% params.

if nargin < 3
    n1 = 101;
end
if nargin < 4
    p1 = [3/4, 0.99, 1/8, 1/2, 1/2, 0.02, 0.01, 1/5, round(2*n1/3)];
    q1 = p1(8:9);
elseif length(p1) == 2
    q1 = p1; p1 = [3/4, 0.99, 1/8, 1/2, 1/2, 0.02, 0.01, q1];
end
n2 = length(T.confirmed{1});
cn = T.Properties.RowNames;
cpop = num2cell(T.Population);
cyy = cf(@rdivide,cf(@transpose,cf(@cat,1,T.confirmed,T.deaths)),cpops);

disp(['First pass: to day ', num2str(n1)]);
cp1 = cell(size(cn)); co1 = cell(size(cn)); cr1 = cell(size(cn));
for i = 1:length(cn)
    disp([num2str(i), ': ', T.State{cn{i}}]);
    [cp1{i},co1{i}] = squidfit(cyy{i}(1:n1,:), p1, 1/cpop{i});
    ytmp = squid2sol(cp1{i}(1:end-1), n1, cp1{i}(end), {4:6,6});
    cr1{i} = r2tst('r2', cyy{i}(1:n1,:), ytmp, length(cp1{i}));
    disp({'resnorm',co1{i}.resnorm;'r2',cr1{i}})
end

disp('Checking R^2 vals...');
ll = c2m(cr1) < .95; kk = find(ll);
if isempty(kk)
    disp('... all good.');
elseif all(ll)
    disp('Something went wrong, exiting.');
    Tout = cat(2, cn, cpop, cyy, cp1, co1, cr1);
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
            [ptmp, otmp] = squidfit(cyy{kk(i)}(1:n1,:), cp1{jj(jk)}(1:end-1), cp1{jj(jk)}(end));
            ytmp = squid2sol(ptmp(1:end-1), n1, ptmp(end), {4:6,6});
            rtmp = r2tst('r2', cyy{kk(i)}(1:n1,:), ytmp, length(ptmp));
            disp({'resnorm',otmp.resnorm;'r2',rtmp})
            if rtmp > cr1{kk(i)}
                cp1{kk(i)} = ptmp; co1{kk(i)} = otmp; cr1{kk(i)} = rtmp;
                if rtmp >= 0.95
                    disp('success')
                    ll(kk(i)) = false;
                end
            end
        end
    end
end
for i = 1:length(cn)
    cp1{i}(8) = max(cp1{i}(8), 0);    % -ve initial q rate has no effect
end
% if nq == 1
%     % wrap up current results
%     Tout = cat(2, cn, cpop, cyy, cp1, co1, cr1);
%     return;
% end

disp(['Second pass: extend to day ', num2str(n2)]);
cp2 = cell(size(cn)); co2 = cell(size(cn)); cr2 = cell(size(cn));
tmplate = -.5*ones(1,2*(nq-1)); tmplate(3:4:end) = q1(1);
tmplate(2:2:end) = (2:2:(2*nq-1))/(2*nq-1);
for i = 1:length(cn)
    disp([num2str(i), ': ', T.State{cn{i}}]);
    qtmp = tmplate; qtmp(2:2:end) = (n2 - cp1{i}(9))*qtmp(2:2:end) + cp1{i}(9);
    [cp2{i},co2{i}] = squidfit(cyy{i}, [cp1{i}(1:end-1), qtmp], cp1{i}(end), 1:7);
    ytmp = squid2sol(cp2{i}(1:end-1), n2, cp2{i}(end), {4:6,6});
    cr2{i} = r2tst('r2', cyy{i}, ytmp, length(cp2{i})-8);
    disp({'resnorm',co2{i}.resnorm;'r2',cr2{i}})
end

for i = 1:length(cn)    % -ve initial q rate again, now we have to check all q's
    for j = 1:nq
        if cp2{i}(8 + 2*(j-1)) <= 0
            cp2{i}(8 + 2*(j-1)) = 0;
        else
            break;
        end
    end
end
disp('States with low R^2:')
cat(2,cn(c2m(cr2)<.95),cr2(c2m(cr2)<.95))    % just show for now, will deal with on 3rd pass if need be
disp('Relative change in initial condition:')
feval(@(x) [min(x),median(x),max(x)], c2m(cf(@(x1,x2) 2*abs(x1-x2)/(x1+x2), cf(@id,cp1,'end'),cf(@id,cp2,'end')))')

disp('Third pass: full refit + 2yr projection');
cp3 = cell(size(cn)); co3 = cell(size(cn)); cr3 = cell(size(cn)); cy3 = cell(size(cn));
for i = 1:length(cn)
    disp([num2str(i), ': ', T.State{cn{i}}]);
    [cp3{i},co3{i}] = squidfit(cyy{i}, cp2{i}(1:end-1), cp2{i}(end));
    cy3{i} = squid2sol(cp3{i}(1:end-1), 730, cp3{i}(end));
    cr3{i} = r2tst('r2', cyy{i}, [sum(cy3{i}(1:n2,3:5),2),cy3{i}(1:n2,5)], length(cp3{i}));
    disp({'resnorm',co3{i}.resnorm;'r2',cr3{i}})
end

disp('Checking R^2 vals...');
ll = c2m(cr3) < .95; kk = find(ll);
if isempty(kk)
    disp('... all good.');
elseif all(ll)
    disp('Something went wrong, exiting.');
    Tout = {cn, cpop, cyy, cp1, co1, cr1, cp2, co2, cr2, cp3, co3, cr3, cy3, ll};
    return;
else
    disp('...');
    ttmp = sortrows(cat(2,cn,cf(@getfield,co3,'resnorm'),cr3,num2cell(1:length(cn))'),[3,2,1,4]);
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
                Tout = {cn, cpop, cyy, cp1, co1, cr1, cp2, co2, cr2, cp3, co3, cr3, cy3, ll, kk, i, tdist, tried};
                return;
            end
            [~, jk] = min(tdist(jj,i)); tried(jj(jk),i) = true;
            disp(['Refitting ' cn{kk(i)} ' based on ' cn{jj(jk)}]);
            [ptmp, otmp] = squidfit(cyy{kk(i)}, cp3{jj(jk)}(1:end-1), cp3{jj(jk)}(end));
            ytmp = squid2sol(ptmp(1:end-1), 730, ptmp(end));
            rtmp = r2tst('r2', cyy{kk(i)}, [sum(ytmp(1:n2,3:5),2),ytmp(1:n2,5)], length(ptmp));
            disp({'resnorm',otmp.resnorm;'r2',rtmp})
            if rtmp > cr3{kk(i)}
                cp3{kk(i)} = ptmp; co3{kk(i)} = otmp; cr3{kk(i)} = rtmp; cy3{kk(i)} = ytmp;
                if rtmp >= 0.95
                    disp('success')
                    ll(kk(i)) = false;
                end
            end
        end
    end
end

Tout = cell2table(cat(2,cpop,cyy,cp1,cp3,co3,cr3,cy3),'VariableNames',{'pop','ndata','porig','pfinal','ofit','R2','proj2y'},'RowNames',cn);

end
