function T = collatejhudata(f_pop, f_conf, f_dead, f_recv)
% T = collatejhudata(f_conf, f_dead, f_recv)
%
% Collate the John Hopkins University epidemiological data in files f_conf,
% f_dead, and f_recv into a single table that also includes population and
% density data (in file f_pop).

T = readtable(f_pop, 'ReadVariableNames',true);
if contains('tag', T.Properties.VariableNames)
    Nr = T.tag; T.tag = [];    % if there is a tag column these will be used as row names for easier access
else
    Nr = [];
end

T1 = readjhudata(f_conf); T = join(T, T1);
T.Properties.VariableNames{end} = 'confirmed';

if nargin > 3 && ~isempty(f_recv)
    T1 = readjhudata(f_recv); T = join(T, T1(:,[1 end]));
    T.Properties.VariableNames{end} = 'recovered';
end

T1 = readjhudata(f_dead); T = join(T, T1(:,[1 end]));
T.Properties.VariableNames{end} = 'deaths';

if ~isempty(Nr)
    [Nr, ii] = sort(Nr);
    T = T(ii,:); T.Properties.RowNames = Nr;
    if isequal(T.Properties.VariableNames{1},'Country')
        T.Country{'us'} = 'United States';
    end
end

end

