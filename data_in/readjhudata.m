function T = readjhudata(f_jhu)
% T = readjhudata(f_jhu)
%
% Read the John Hopkins University data in the file f_jhu and put into a
% normal form for combination with other datasets. This involves mainly
% aggregating the data that appears by province for some countries (eg
% Australia) rather than nationally, and promoting entries from some of the
% French etc colonial possessions into the County/Region category, as well
% as collecting the time series into cells so multiple stats can be held in
% one table.

warning('OFF', 'MATLAB:table:ModifiedAndSavedVarnames');
T = readtable(f_jhu, 'ReadVariableNames',true);
warning('ON', 'MATLAB:table:ModifiedAndSavedVarnames');
if isequal(T.Properties.VariableNames{1}, 'UID')   % indicates US data broken out by municipality
    T = T(:,contains(T.Properties.VariableNames,{'Admin2','Province','Lat','Long'}) | ~cellfun(@isempty,T.Properties.VariableDescriptions));
    T.Properties.VariableNames(1:2) = {'Province_State','Country_Region'};
    desig = 'State';
else
    desig = 'Country';
end

% Stage 1: collect out non-empty Province/State data to either aggregate or
% read as its own Country/Region
kk = ~cellfun(@isempty, T.Province_State);
Tpr = T(kk,:); T(kk,:) = [];

% Get names of countries to which province/states belong, and check if any
% have their own entry in main table; these are overseas possessions that
% can be simply promoted to the country column for our purposes
N0 = unique(sort(Tpr.Country_Region));
kk = ismember(N0, T.Country_Region); N1 = N0(~kk);    % + save aside names with no independent row
ll = ismember(Tpr.Country_Region, N0(kk));
Tpr.Country_Region(ll) = Tpr.Province_State(ll);
T = [T; Tpr(ll,:)]; Tpr(ll,:) = [];

% Remaining countries in N1 have to have province data combined into
% single row, which will as well be readded to T
for i = 1:length(N1)
    jj = find(ismember(Tpr.Country_Region, N1{i}));
    tmp = Tpr{jj,3:4}; Tpr{jj(1),3:4} = mean(tmp(any(tmp ~= 0, 2),:));    % use mean of non-null longitude and latitude over provinces (null denotes eg ships at sea)
    Tpr{jj(1),5:end} = sum(Tpr{jj,5:end});  % sum statistic of interest over all provinces
    Tpr(jj(2:end),:) = [];     % clear other provinces
end
T = [T; Tpr];

% Remove the province/state column and resort
T.Province_State = []; T = sortrows(T);

% Stage 2: collect time series data into arrays (+ collect dates themselves
% from variable names)
D = datetime(strrep(strrep(T.Properties.VariableNames(4:end),'x',''),'_','/'), 'inputformat','MM/dd/yy', 'pivotyear',2000);
if isequal(D, D(1):D(end))
    D = D([1 end]);
end
tmp = num2cell(T{:,4:end}, 2);
T.Properties.VariableNames(1:5) = {desig,'Lat','Long','Range','Data'};
T(:,6:end) = [];
T.Range = cell(height(T),1); T{:,4} = {D};
T.Data = tmp;

end

