function outable = genhubfile(rtable, dvec, loctable, wks, fnamstub)
% outable = genhubfile(intable, dvec, loctable, wks, fnamstub)
%
% Generate a Matlab table of predictions from fit results and projections
% in rtable which can uploaded to the COVID-19 ForecastHub. dvec gives dates
% over which the projections in rtable run, loctable gives a listing of US
% state FIPS codes (which need to be included), and wks is an integer
% giving the number of weeks for which predictions are desired (by default, 1).
% If the last row of loctable refers to the U.S. as a whole (with location
% given as 'US') then aggregate results for the entire country will also be
% included.
%   Returns the table of predictions (row-sorted); if string fnamstub is
% provided function also saves a csv file using the ForecastHub naming rules.

vnames = {'forecast_date','target','target_end_date','location','type','quantile','value'};
cdata = cell(size(vnames));
outable = table(cdata{:}, 'VariableNames', vnames);
itrg = 2; ited = 3; iloc = 4; ival = 7;
cdata{5} = 'point'; cdata{6} = 'NA';

targs = {' wk ahead cum death',' wk ahead inc death',' wk ahead inc case'};
if nargin < 4 || isempty(wks)
    wks = 1;
end
aggus = zeros(wks,3);

% do some datetime math to get the values we want
curd = datetime('today','Format','yyyy-MM-dd'); cdata{1} = curd;
jc = datenum(curd) - datenum(dvec(1)) + 1;    % current date wrt projection period
jw = weekday(curd); jj = (jc - jw) + 7*(0:wks);
tards = datetime(dvec(jj(2:end)), 'Format','yyyy-MM-dd');

for i = 1:size(loctable,1)
    cdata{iloc} = loctable.location{i};
    curst = loctable.abbreviation{i};
    if isequal(curst, 'us')    % create aggregate predictions from aggus
        vals = aggus;
    else
        tmpv = rtable.pop(curst) * rtable.proj2y{curst}(jj,3:5);
        vals = [tmpv(2:end,end), diff(tmpv(:,end)), diff(sum(tmpv,2))];
        aggus = aggus + vals;
    end
    for j = 1:wks
        cdata{ited} = tards(j);
        for k = 1:3
            cdata{itrg} = [num2str(j), targs{k}];
            cdata{ival} = round(1000*vals(j,k))/1000;
            outable = [outable; cdata];
        end
    end
end
outable = sortrows(outable);

if nargin > 4 && ~isempty(fnamstub)
    writetable(outable, [datestr(curd,'yyyy-mm-dd'), '-', fnamstub, '.csv']);
end

end

