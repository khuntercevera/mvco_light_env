load /Users/kristenhunter-cevera/MVCO_light_at_depth/fluorometer/ecochl_all.mat %load all in situ fluorometer results (ecofl series sensors)
plot(dat(:,1), dat(:,7), 'c-', 'linewidth', 2)
ylim([0 15])
datetick('x')
[y,m,d,h,mi,s] = datevec(dat(:,1));
ind_night = find(h < 9); %UTC hours for middle of local night
ind_midday = find(h >= 14 & h <= 18); %UTC hours for mid-day local
hold on
plot(dat(ind_night,1), dat(ind_night,7), 'b.')
plot(dat(ind_midday,1), dat(ind_midday,7), 'g.')
ylabel('Chl (mg m^{-3})')

%calculate the mean and std dev for values each night and each mid-day
day = floor(dat(:,1));
unqdays = unique(day);
ecochl_mean  = NaN(length(unqdays),2);
ecochl_std  = ecochl_mean;
for count = 1:length(unqdays),
    ind = find(day == unqdays(count) & h < 9);
    ecochl_mean(count,1) = nanmean(dat(ind,7)); %night
    ecochl_std(count,1) = nanstd(dat(ind,7),0,1);    
    ind = find(day == unqdays(count) & h <= 16 & h >= 14);
    ecochl_mean(count,2) = nanmean(dat(ind,7)); %mid-day
    ecochl_std(count,2) = nanstd(dat(ind,7),0,1);
end;

plot(unqdays+1/6, ecochl_mean(:,1), 'b^', 'markerfacecolor', 'b')
plot(unqdays+2/3, ecochl_mean(:,2), 'g^', 'markerfacecolor', 'g')


%% Now, load in the fluorometer deployments to see if any went awry:

filename = '/Users/kristenhunter-cevera/MVCO_light_at_depth/fluorometer/Fluor_MVCO_deploy_matlab_readable.txt';
delimiter = '\t';
startRow = 2;
formatSpec = '%s%s%s%s%s%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);
fclose(fileID);

%% reorganize into matrix:
fluor_deployment_titles={'Deployment Number';'Serial number';'Date In';'Date Out';'Calibration'};

fluor_deployment=cell(33,5);

temp=str2num(char(dataArray{:,1}));
for j=1:33 %predetermined :)
    
    jj=find(temp==j);
    
    fluor_deployment{j,1}=j;
    fluor_deployment{j,2}=char(dataArray{2}(jj(1)));

    if ~strcmp(char(dataArray{2}(jj(1))),char(dataArray{2}(jj(2))))
        keyboard
    end
    
    fluor_deployment{j,3}=datenum(char(dataArray{3}(jj(1))));
    if j~=33
    fluor_deployment{j,4}=datenum(char(dataArray{3}(jj(2))));
    else
        fluor_deployment{j,4}=datenum('sept-18-2017');
    end
    
    fluor_deployment{j,5}=datenum(char(dataArray{5}(jj(1))));
end

%% and now color coded figure!
%can organize by deployment, by calibration or by serial number
cc=parula(33);
for j=1:33
    jj=find(cellfun(@(x) x==j,fluor_deployment(:,1))==1);
    
    ii=find(dat(:,1) >= fluor_deployment{jj,3} & dat(:,1) <= fluor_deployment{jj,4});
    plot(dat(ii,1), dat(ii,7),'-', 'color',cc(j,:), 'linewidth', 2)
end


%%
load /Users/kristenhunter-cevera/MVCO_light_at_depth/fluorometer/CHLatASIT.mat %load discrete sample extracted chl results
hold on
plot(FL_matdate, FL_chl(:,1), 'r*') %fluorometric analysis of extracts
ylabel('Chl (mg m^{-3})')

legend('eco, in situ', 'eco night', 'eco day','eco avg night', 'eco avg day', 'extract-fl')
set(gcf, 'position', [29 378 1388 420])

%find the set of in situ eco results that "match up" with discrete samples
%consider both night before, night after (and day time too, but won't use these)
FLday = floor(FL_matdate);
ecochl_match = NaN(length(FLday),4);
ecochl_match_std = ecochl_match;
for count = 1:length(FLday),
    ind = find(unqdays == FLday(count));
    if ~isempty(ind),
        ecochl_match(count,1) = ecochl_mean(ind,1); %night before
        ecochl_match_std(count,1) = ecochl_std(ind,1); %night before
        ecochl_match(count,3) = ecochl_mean(ind,2); %day before
    end;
    ind = find(unqdays == FLday(count)+1);
    if ~isempty(ind),
        ecochl_match(count,2) = ecochl_mean(ind,1); %night after
        ecochl_match_std(count,2) = ecochl_std(ind,1); %night after
        ecochl_match(count,4) = ecochl_mean(ind,2); %day after
    end;
end;

[y,m,d,h,mi,s] = datevec(FL_matdate);
figure
plot(FL_matdate-datenum(y,1,0), ecochl_match(:,1)./FL_chl(:,1), '.') %night before
hold on
plot(FL_matdate-datenum(y,1,0), ecochl_match(:,2)./FL_chl(:,1), '.c') %night after
ylim([0 4]), xlim([0 365])
line([0 365], [1 1])
ylabel('Chl eco/discrete-fl')
xlabel('Year day')
legend('night before', 'night after')

figure
plot(FL_matdate, FL_chl(:,1), '.-')
hold on
errorbar(FL_matdate, ecochl_match(:,1), ecochl_match_std(:,1), 'r*')
errorbar(FL_matdate, ecochl_match(:,2), ecochl_match_std(:,2), 'm*')
ylabel('Chl (mg m^{-3})')
legend('discrete-fl', 'eco, night before', 'eco, night after')
datetick('x')
set(gcf, 'position', [29 378 1388 420])

figure
plot(FL_matdate, ecochl_match(:,1)./FL_chl(:,1), '.') %night before
hold on
plot(FL_matdate, ecochl_match(:,2)./FL_chl(:,1), '.c') %night after
datetick('x')
ylim([0 4])
set(gca, 'xgrid', 'on')
line(xlim, [1 1])
ylabel('Chl eco/discrete-fl')

%load DailySolar %MJ m^-2
%yd = (1:365)';
%Solar_matday = [yd+datenum(2003,1,0); yd+datenum(2004,1,0); yd+datenum(2005,1,0); yd+datenum(2006,1,0); yd+datenum(2007,1,0); yd+datenum(2008,1,0); yd+datenum(2009,1,0); yd+datenum(2010,1,0)];

%Solar_match = NaN(1,length(unqdays));
%for count = 1:length(unqdays),
%    ind = find(unqdays(count) == Solar_matday);
%    if ~isempty(ind),
%        Solar_match(count) = DailySolar(ind);
%    end;
%end;

%Solar_matchFL = NaN(1,length(FLday));
%for count = 1:length(FLday),
%    ind = find(FLday(count) == Solar_matday);
%    if ~isempty(ind),
%        Solar_matchFL(count) = DailySolar(ind);
%    end;
%end;

%find the set of in situ eco results that "match up" with discrete samples
%this time for HPLC analysis
%consider both night before, night after (and day time too, but won't use these)
HPLCday = floor(HPLC_matdate);
ecochl_match_hplc = NaN(length(HPLCday),4);
ecochl_match_hplc_std = ecochl_match_hplc;
for count = 1:length(HPLCday),
    ind = find(unqdays == HPLCday(count));
    if ~isempty(ind),
        ecochl_match_hplc(count,1) = ecochl_mean(ind,1); %night before
        ecochl_match_hplc_std(count,1) = ecochl_std(ind,1); %night before
        ecochl_match_hplc(count,3) = ecochl_mean(ind,2); %day before
    end;
    ind = find(unqdays == HPLCday(count)+1);
    if ~isempty(ind),
        ecochl_match_hplc(count,2) = ecochl_mean(ind,1); %night after
        ecochl_match_hplc_std(count,2) = ecochl_std(ind,1); %night after
        ecochl_match_hplc(count,4) = ecochl_mean(ind,2); %day after
    end;
end;


figure
plot(FL_chl(:,1), nanmean(ecochl_match(:,1:2),2), '.')
hold on
plot(HPLC_chl, nanmean(ecochl_match_hplc(:,1:2),2), 'r.')
ylabel('Chl (mg m^{-3}), est. in situ')
xlabel('Chl (mg m^{-3}), extract')
legend('FL', 'HPLC')
line([0 12], [0 12])

y = 0:.1:4.5;
x = .5*y.^2;
hold on
plot(x, y, 'g-')
