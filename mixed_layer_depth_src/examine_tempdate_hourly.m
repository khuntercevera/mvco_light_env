%% investigate whether or not see diel changes in temperature from 4m to 12m:
% uses scripts from Heidi that aggregate Tday products


% temperature from the beam: every hour:

mdate = [];
T = [];

for ii = 4:6
    eval(['load /Volumes/Lab_data/MVCO/MVCO_Moxanode/code/CTDsmoothhr20' num2str(ii,'%02d')]);
    mdate = [mdate; CTDsmooth(:,1)];
    T = [T; CTDsmooth(:,2)];
end;

%early part of 2007 (from mininode period)
load /Volumes/Lab_data/MVCO/MVCO_Moxanode/code/CTDsmoothhr2007a
mdate = [mdate; CTDsmooth(:,1)];
T = [T; CTDsmooth(:,2)];

% transition to moxanode
for ii = 7:17
    eval(['load /Volumes/Lab_data/MVCO/MVCO_Moxanode/code/CTDsmoothhr20' num2str(ii,'%02d')]);
    mdate = [mdate; CTDsmooth(:,1)];
    T = [T; CTDsmooth(:,2)];
end;

Tbeam=T;
mdate_beam=mdate;

%% load in the node data: every 20 min:

datadir = '/Volumes/Lab_data/MVCO/EnvironmentalData/';
load([datadir 'other2001'])
load([datadir 'other2002'])
load([datadir 'other2003_04'])
load([datadir 'other2005'])
load([datadir 'other2006'])
load([datadir 'other2007'])
load([datadir 'other2008'])
load([datadir 'other2009'])
load([datadir 'other2010'])
load([datadir 'other2011'])
load([datadir 'other2012'])
load([datadir 'other2013'])
load([datadir 'other2014'])
load([datadir 'other2015'])
load([datadir 'other2016'])
load([datadir 'other2017'])

%this will lead to averaging of seacat and node temps for overlap days
yd_ocn2003 = [yd_ocn2003; yd_seacat2003];
Temp2003 = [Temp2003; temp_seacat2003];
yd_ocn2004 = [yd_ocn2004; yd_seacat2004];
Temp2004 = [Temp2004; temp_seacat2004];


%% okay, so now that have raw data loaded, ask:
%difference within a day?
%temperature difference within a day?
year=2001:2017;

mdate_node=[]; Tnode=[];
for count = 1:length(year),
    eval(['yd_ocn = yd_ocn' num2str(year(count)) ' + datenum(''1-0-' num2str(year(count)) ''');'])
    eval(['Temp = Temp' num2str(year(count)) ';'])
    mdate_node=[mdate_node; yd_ocn];
    Tnode=[Tnode; Temp];
end;

%% find the hourly value for each year:

%this takes awhile! ~30 min!!!
yearlist=2003:2017;

Tbeam_hour=nan(366*24*length(yearlist),1);
Tnode_hour=nan(366*24*length(yearlist),1);
time_hour=nan(366*24*length(yearlist),1);

vec_mdate_Tbeam=datevec(mdate_beam);
vec_mdate_Tnode=datevec(mdate_node);

count=0;
tic
for yr=1:length(yearlist)
    for d=1:366
        for hr=0:23
            count=count+1;
            
            tempvec=datevec(datenum(['1-0-' num2str(yearlist(yr))]) + d);
            tempvec(4)=hr;
            time_hour(count)=datenum(tempvec);
            
            bb=find(vec_mdate_Tbeam(:,1)==tempvec(1) & vec_mdate_Tbeam(:,2)==tempvec(2) ...
                & vec_mdate_Tbeam(:,3)==tempvec(3) & vec_mdate_Tbeam(:,4)==tempvec(4)); %locate year-month-day-hour
            Tbeam_hour(count)=nanmean(Tbeam(bb));
            
            nn=find(vec_mdate_Tnode(:,1)==tempvec(1) & vec_mdate_Tnode(:,2)==tempvec(2) ...
                & vec_mdate_Tnode(:,3)==tempvec(3) & vec_mdate_Tnode(:,4)==tempvec(4));
            Tnode_hour(count)=nanmean(Tnode(nn));
        end
    end
    
end
toc

%this take half an hour! :(
%%

% remove non-leap year day repeats:

for j=[2004 2006 2007 2008 2010 2011 2012 2014 2015 2016]
    
    ii=find(floor(time_hour)==datenum(['1-1-' num2str(j)]));
    if length(ii) ~= 48
        keyboard
    end
    time_hour(ii(1):ii(24)) = [];
    Tbeam_hour(ii(1):ii(24)) =[];
    Tnode_hour(ii(1):ii(24)) =[];
    
end


%% okay, now, let's plot and look at difference!

figure, plot(mdate_beam,Tbeam,'.-')
hold on
plot(mdate_node,Tnode,'.-')
plot(time_hour,Tnode_hour,'k.-')
plot(time_hour,Tbeam_hour,'b.-')

%looking good except for a few stray points...

%% remove known bad points in Tbeam:

%manually:
jj=find(diff(Tbeam_hour) > 3);
Tbeam_hour(jj+1) %these should be it!
Tbeam_hour(jj+1)=NaN;


%% and now difference! Yay!

%ask, within a day, what is the maximum and minimum difference? How many
%hours within a day are above a 0.6 degree difference?

unqdays=unique(floor(time_hour-4/24));
unqdays=unqdays(2:end-1); %remove a dec 31 2002 and Jan 1 2018
mindiff_day=nan(length(unqdays),1);
maxdiff_day=nan(length(unqdays),2);
min_day=nan(length(unqdays),2);
max_day=nan(length(unqdays),2);

for q=1:length(unqdays)
    ii=find(floor(time_hour-4/24)==unqdays(q));
    
    [~,im]=min(abs(Tbeam_hour(ii)-Tnode_hour(ii)));
    mindiff_day(q)=Tbeam_hour(ii(im))-Tnode_hour(ii(im));  %min abs temp diff temp within a day
    
    [~,ma]=max(abs(Tbeam_hour(ii)-Tnode_hour(ii)));
    maxdiff_day(q)=Tbeam_hour(ii(ma))-Tnode_hour(ii(ma)); %max abs temp diff temp within a day
    maxdiff_day(q,2)=max(Tbeam_hour(ii)-Tnode_hour(ii)); %max temp diff temp within a day
    
    min_day(q,:)=[min(Tbeam_hour(ii)) min(Tnode_hour(ii))]; %min temp within a day
    max_day(q,:)=[max(Tbeam_hour(ii)) max(Tnode_hour(ii))]; %max temp within a day
    
end


%% more plots...

subplot(2,3,1,'replace')
plot(find_yearday(time_hour),Tbeam_hour-Tnode_hour,'.')
xlabel('year day')
ylabel('Tbeam (4 m) - Tnode (12 m) (\circC)')
title('hourly differences')
xlim([1 366])

subplot(2,3,2,'replace')
plot(max_day(:,2),max_day(:,1),'.')
line([0 22],[0 22],'color','r')
xlabel('maximum temp within a day at 12 m node')
ylabel('maximum temp within a day at 4 m node')
title('within a day')

subplot(2,3,3,'replace')
plot(min_day(:,2),min_day(:,1),'.')
line([0 22],[0 22],'color','r')
xlabel('minimum temp wihtin a day at 12 m node')
ylabel('minimum temp within a day at 4 m node')
title('within a day')

subplot(2,3,4,'replace')
plot(find_yearday(unqdays),mindiff_day,'.')
xlabel('year day')
ylabel('minimum difference in day between 4 m beam and 12 m node')
xlim([1 366])

subplot(2,3,5,'replace')
plot(find_yearday(unqdays),maxdiff_day(:,1),'.')
xlabel('year day')
ylabel('maximum abs difference in a day between 4 m beam and 12 m node')
xlim([1 366])

subplot(2,3,6,'replace')
plot(find_yearday(unqdays),maxdiff_day(:,2),'.')
xlabel('year day')
ylabel('maximum difference in a day between 4 m beam and 12 m node')
xlim([1 366])

%% save this figure for tex document:

%shows spread of values, but also that distinct seasonality:

figure
plot(find_yearday(time_hour),Tbeam_hour-Tnode_hour,'.')
xlabel('year day')
ylabel('Tbeam (4 m) - Tnode (12 m) (\circC)')
title('hourly \DeltaT by yearday')
xlim([1 366])
set(gca,'fontsize',20)
%%
set(gcf,'color','w')
addpath /Users/kristenhunter-cevera/Documents/MATLAB/matlab_tools/export_fig_2016/

export_fig /Users/kristenhunter-cevera/MVCO_light_at_depth/figures_for_tex_doc/hourly_deltaT_by_yearday.pdf

%% compare with light data!

total_Solar=[];
solar_time=[];
total_dawn=[];

for year=2003:2017
    
    switch year
        case 2003
            filelabel='May';
        case 2004
            filelabel='Apr';
        case 2005
            filelabel='Apr';
        case 2006
            filelabel='May';
        case 2007
            filelabel='Mar';
        otherwise
            filelabel='Jan';
    end
    
    eval(['load /Volumes/Lab_data/MVCO/FCB/MVCO_' filelabel num2str(year) '/model/solar' num2str(year) '.mat'])
    Solar(Solar <0)=0;
    
    total_Solar=[total_Solar; Solar];
    solar_time=[solar_time; date_met];
    total_dawn=[total_dawn; dawn];
    
end

% record rough dawn and dusks:
dawn_dusk=nan(length(unqdays),2);
for q=1:length(unqdays)
    jj=find(floor(solar_time-4/24)==unqdays(q));
    ee=find(total_Solar(jj) > 5); %find range that is above 5 W/2
    if ~isempty(ee)
        [~,~,~,h1,~,~]=datevec(solar_time(jj(ee(1)))-4/24);
        [~,~,~,h2,~,~]=datevec(solar_time(jj(ee(end)))-4/24);
        dawn_dusk(q,:)=[h1 h2];
    end
end

% just use median values to find light window:
temp=find_yearday(unqdays);
med_dd=nan(366,2);
for j=1:366
    ii=find(temp==j); %current yearday
    med_dd(j,:)=[nanmedian(dawn_dusk(ii,1)) nanmedian(dawn_dusk(ii,2))];
end

%a bit rough, but okay for now...

figure, plot(1:366,med_dd(:,1),'.')
hold on
plot(1:366,med_dd(:,2),'.')

%% okay, now go through and check what time of day stratification happens:

% a slighlty more complicated question: how many hours in daylight would
% the water column be considered 'stratified'? Continuously?

above_06=nan(length(unqdays),3);

for q=1:length(unqdays)
    
    tempyrdy=find_yearday(unqdays(q)); %unqdays is in local time
    
    %find daylight portion:
    ii=find((time_hour-4/24) >= (unqdays(q) + med_dd(tempyrdy,1)/24) &  (time_hour-4/24) <= (unqdays(q) + med_dd(tempyrdy,2)/24)  );
    
    if isempty(ii) %should always be able to find this!
        keyboard
    end
    
    above_06(q,3)=length(ii); %record length of daylight as sanity check
    
    %how many hours would be considered stratified?
    
    delta=Tbeam_hour(ii) - Tnode_hour(ii);   
    delta=delta(~isnan(delta));
        
    if ~isempty(delta) %meaning some real values:
        jj=find(delta >= 0.6); %over 0.6 looks like stratification
        if ~isempty(delta)
            above_06(q,1:2)=[length(jj) length(delta)];
        else
            above_06(q,1:2)=[0 length(delta)];
        end
    end
    
end

%% sort so can put in a bar chart:

%bin by week:
unq_yrdy=find_yearday(unqdays);
wk_days=[1:7:365 367];
rec=nan(length(wk_days),12);

for w=1:length(wk_days)-1
    jj=find(unq_yrdy >= wk_days(w) &  unq_yrdy < wk_days(w+1) & ~isnan(above_06(:,2)) & (above_06(:,2)./above_06(:,3) > 0.6));
    %find yearday within week and only consider days where at least 0.6 of the hours of the day were available
    rec(w,1:11)=histc(above_06(jj,1)./above_06(jj,2),0:0.1:1)';
    rec(w,12)=sum(rec(w,1:11));
end

%%
clf
bar(wk_days,rec(:,1:11),2.5,'stacked') %hooray!
xlim([-2 368]) %breathing space
colormap jet
hbar=colorbar;
set(hbar,'yticklabel',cellstr(num2str((0:0.1:1)')))
xlabel('Year Day')
ylabel('Frequency')
title('Percentage of day light hours "stratified" for available days in dataset')

%%
set(gcf,'color','w')
addpath /Users/kristenhunter-cevera/Documents/MATLAB/matlab_tools/export_fig_2016/

export_fig /Users/kristenhunter-cevera/MVCO_light_at_depth/figures_for_tex_doc/percentage_hours_stratified.pdf

%% see how many days would have some type of stratification in the
clf
plot(find_yearday(unqdays),above_06(:,2),'.')
hold on
plot(find_yearday(unqdays),above_06(:,1),'.')

%% at any rate, it's a low percentage of days that have more than 2/3 possibly stratified:
kk=find(above_06(:,2)./above_06(:,3) > 0.6); %2/3 days
jj=find(above_06(kk,1)./above_06(kk,2) > 0.5); %of those days, how many hours are stratified?
length(jj)./length(kk)
