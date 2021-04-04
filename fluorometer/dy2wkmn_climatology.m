function [wk_avg, wk_std, yd_wk, mn_avg, mn_std, yd_mn] = dy2wkmn_climatology(dy_data, yearlist)
%function [wk_avg, wk_std, yd_wk, mn_avg, mn_std, yd_mn] = dy2wkmn_climatology(dy_data, yearlist)
%given a matrix of daily data and and a yearlist, produced from Heidi's timeseries2ydmat
%for each week and month, calculates the average and std dev from all the
%points in the timeseries. This results in a slightly different number than
%if one were to construct the climatologies with precomputed weekly or
%monthly averages for each year.

%weekly climatology:

yd_wk = (1:7:364)';
wk_avg = nan(52,1);
wk_std = nan(52,1);

for week = 1:52
    
    if week < 52
        ii = week*7-6:week*7;
        temp_data=dy_data(ii,:); temp_data=temp_data(:);
        wk_avg(week) = nanmean(temp_data);
        wk_std(week) = nanstd(temp_data); 
    else
        ii= week*7-6:366;
        temp_data=dy_data(ii,:); temp_data=temp_data(:);
        wk_avg(week) = nanmean(temp_data);
        wk_std(week) = nanstd(temp_data); 
    end;
end

%
%monthly climatology:

%from Heidi's code:
%find which month each day belongs to for the data:
numyrs = length(yearlist);
mdate_year = datenum(yearlist,0,0);
yd = (1:366)';
mdate = repmat(yd,1,numyrs)+repmat(mdate_year,length(yd),1);
[~,month,~,~,~,~] = datevec(mdate);

mn_avg=nan(12,1);
mn_std=nan(12,1);

%yd_mn = datenum([zeros(12,1) (1:12)' ones(12,1)]); %corresponds to first of month for a leap year
yd_mn = datenum([repmat(2003,12,1) (1:12)' ones(12,1)])-datenum('1-0-2003'); %corresponds to first of month for a non-leap year

if any(size(dy_data)~=size(month))
    keyboard
end

for mn = 1:12,
    ii=find(month==mn);
    mn_avg(mn) = nanmean(dy_data(ii));
    mn_std(mn) = nanstd(dy_data(ii)); 
end;




