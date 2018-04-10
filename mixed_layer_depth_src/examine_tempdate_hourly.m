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

mdate_node=[]; Tnode=[];
for count = 1:length(year),    
    eval(['yd_ocn = yd_ocn' num2str(year(count)) ' + datenum(''1-0-' num2str(year(count)) ''');'])
    eval(['Temp = Temp' num2str(year(count)) ';'])
    mdate_node=[mdate_node; yd_ocn];
    Tnode=[Tnode; Temp];
end;

%% find the hourly value for each year:

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

%% remove bad points:

%manually:
%datestr(unqdays(find(maxdiff_day > 4)))
%jj=find(floor(time_hour)==datenum('20-Jul-2017'));
%Tbeam_hour(jj(6))=NaN:


%% and now difference! Yay!

%ask, within a day, what is the maximum and minimum difference? How many
%hours within a day are above a 0.6 degree difference?

unqdays=unique(floor(time_hour-4/24));
mindiff_day=nan(length(unqdays),1);
maxdiff_day=nan(length(unqdays),2);
min_day=nan(length(unqdays),2);
max_day=nan(length(unqdays),2);
above_06=nan(length(unqdays),2);

for q=1:length(unqdays)
    ii=find(floor(time_hour-4/24)==unqdays(q));
    
    [~,im]=min(abs(Tbeam_hour(ii)-Tnode_hour(ii)));
    mindiff_day(q)=Tbeam_hour(ii(im))-Tnode_hour(ii(im));   
    
    [~,ma]=max(abs(Tbeam_hour(ii)-Tnode_hour(ii)));
    maxdiff_day(q)=Tbeam_hour(ii(ma))-Tnode_hour(ii(ma));
    maxdiff_day(q,2)=max(Tbeam_hour(ii)-Tnode_hour(ii));
    
    min_day(q,:)=[min(Tbeam_hour(ii)) min(Tnode_hour(ii))];
    max_day(q,:)=[max(Tbeam_hour(ii)) max(Tnode_hour(ii))];
    
end

%% more plots...

subplot(2,3,1,'replace')
plot(max_day(:,2),max_day(:,1),'.')
line([0 22],[0 22],'color','r')
xlabel('max temp 12 m node')
ylabel('max temp 4 m node')

subplot(2,3,2,'replace')
plot(min_day(:,2),min_day(:,1),'.')
line([0 22],[0 22],'color','r')
xlabel('min temp 12 m node')
ylabel('min temp 4 m node')

subplot(2,3,3,'replace')
plot(find_yearday(time_hour),Tbeam_hour-Tnode_hour,'.')
xlabel('year day')
ylabel('Tbeam - Tnode')

subplot(2,3,4,'replace')
plot(find_yearday(unqdays),maxdiff_day(:,1),'.')
xlabel('year day')
ylabel('max abs difference between 4 m beam and 12 m node')

subplot(2,3,5,'replace')
plot(find_yearday(unqdays),maxdiff_day(:,2),'.')
xlabel('year day')
ylabel('max difference between 4 m beam and 12 m node')

subplot(2,3,6,'replace')
plot(find_yearday(unqdays),mindiff_day,'.')
xlabel('year day')
ylabel('min difference between 4 m beam and 12 m node')


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


for q=1:length(unqdays)
    jj=find(floor(solar_time-4/24)==unqdays(q));
    ii=find(floor(time_hour-4/24)==unqdays(q)); %unqdays is local time


end
% a slighlty more complicated question: how many hours in daylight would
% the water column be considered 'stratified'? Continuously?

%above_06(q,:)= [length(find( (Tbeam_hour(ii)-Tnode_hour(ii)) > 0.6))   length(ii)];

