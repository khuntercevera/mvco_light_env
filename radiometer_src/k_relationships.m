%% Collect and organize k values:

%gather k values from each casts
% take various means and aggregates to look at seasonality


%% A time series plot of K_d values!

%sourcepath=fullfile('\\sosiknas1\lab_data\mvco\HyperPro_Radiometer\processed_radiometer_files\'); % path to folders with raw data...
sourcepath=fullfile('/Volumes/Lab_data/MVCO/HyperPro_Radiometer/');
processed_path=fullfile(sourcepath,'/processed_radiometer_files/');

%How many of these do we have?
d = dir(processed_path);
isub = [d(:).isdir]; %# returns logical vector
foldernames = {d(isub).name}'; %names of folders
temp=regexp(foldernames,'\d{1,2}\w{3,4}20\d{2}'); %find all folders with this naming scheme
datafolders = foldernames(cellfun('isempty',temp)==0);

load(fullfile(sourcepath,'good_data_folders.mat'))

%% Gather K values from all the data to plot by time and by lat & lon...

k_values=[];
k_record={};
lambdas=[];
k_lambdas=[];

for foldernum=good_data' %folders with viable casts in them
    
    %load in all the data!
    matsource=fullfile(processed_path,datafolders{foldernum},'/mat_outfiles/');
    
    eval(['load ' matsource 'data_' datafolders{foldernum} '.mat'])
    eval(['tempdata=data_' datafolders{foldernum} ';'])
    
    eval(['load ' matsource 'K_PAR_' datafolders{foldernum} '.mat'])
    eval(['K_PAR=K_PAR_' datafolders{foldernum} ';'])
    
    eval(['load ' matsource 'location_' datafolders{foldernum} '.mat'])
    eval(['location=location_' datafolders{foldernum} ';'])
    
    eval(['load ' matsource 'k_lambda_' datafolders{foldernum} '.mat'])
    eval(['k_lambda=k_lambda_' datafolders{foldernum} ';'])
    
    %because some of the casts were split, need to account for this:
    for filenum=1:length(K_PAR);
        
        if K_PAR(filenum).flag==0;
            
            k_values=[k_values; datenum(datestr(datafolders{foldernum})) location(filenum).lat location(filenum).lon K_PAR(filenum).K(2) K_PAR(filenum).stats(1) K_PAR(filenum).flag];
            k_record=[k_record; {datafolders{foldernum}} {foldernum} {K_PAR(filenum).file} {location(filenum).file} location(filenum).eventnum {K_PAR(filenum).NOTES} {location(filenum).notes}];
            k_lambdas=[k_lambdas; -k_lambda(filenum).k_wv(:,3)'];  %attenuation coefficient for each wavelength
            lambdas=[lambdas; k_lambda(filenum).k_wv(:,1)']; %wavelength used
            
        elseif K_PAR(filenum).flag==3;
            
            k_values=[k_values; datenum(datestr(datafolders{foldernum})) location(filenum).lat location(filenum).lon K_PAR(filenum).K1(2) NaN K_PAR(filenum).flag];
            k_record=[k_record; {datafolders{foldernum}} {foldernum} {K_PAR(filenum).file} {location(filenum).file} location(filenum).eventnum {K_PAR(filenum).NOTES} {location(filenum).notes}];
            k_lambdas=[k_lambdas; -k_lambda(filenum).k_wv1(:,3)']; %attenuation coefficient for each wavelength
            lambdas=[lambdas; k_lambda(filenum).k_wv1(:,1)']; %wavelength used
        end
    end
    
end

if sum(sum(lambdas + repmat(-lambdas(1,:),size(lambdas,1),1)))==0
    disp('all recorded wavelengths are the same!')
    lambdas=lambdas(1,:);
else
    disp('Ruh oh - not all wavelengths are the same?')
end

%add in yearday:
k_values=[find_yearday(k_values(:,1)) k_values];

k_record_titles={'date';'folder number';'K_PAR file name';'location file name';'event number';'K_PAR notes';'location notes'};
k_value_titles={'year day';'matlab date';'lat';'lon';'k';'r2';'flag';'station'};
%% organize these points based on position:

%a quick plot of where these points lie:
subplot(2,3,1,'replace'), hold on
%plot station points:
plot(-70.567,41.325,'o','markersize',16,'color',[0.5 0.5 0.5]) %tower
plot(-70.555,41.335,'o','markersize',16,'color',[0.5 0.5 0.5]) %node
plot(-70.6275,41.325,'o','markersize',16,'color',[0.5 0.5 0.5]) %station
plot(-70.505,41.325,'o','markersize',16,'color',[0.5 0.5 0.5]) %station
plot(-70.45,41.3275,'o','markersize',16,'color',[0.5 0.5 0.5]) %station
plot(-70.567,41.255,'o','markersize',16,'color',[0.5 0.5 0.5]) %station
plot(-70.567,41.2,'o','markersize',16,'color',[0.5 0.5 0.5]) %station
plot(-70.567,41.145,'o','markersize',16,'color',[0.5 0.5 0.5]) %station

plot(k_values(:,4),k_values(:,3),'r.')

%% Now code each cast more or less to a station number:

approx_station_loc=[4 -70.567 41.325; %tower
    3 -70.5564 41.3366; %node
    1 -70.45 41.3275; %station 1
    2 -70.505 41.325; %station 2
    5 -70.6275 41.325;  %station 5
    6 -70.567 41.255; %station 6
    7 -70.567 41.2; %station 7
    8 -70.567 41.145]; %station 8

shore_pos=[41.3499,-70.5267];

% set up the indexes:

for j=1:length(approx_station_loc)
    temp_ind=find(k_values(:,3) > approx_station_loc(j,3)-0.01 & k_values(:,3) < approx_station_loc(j,3)+0.01 ... %lat
        & k_values(:,4) > approx_station_loc(j,2)-0.01 & k_values(:,4) < approx_station_loc(j,2)+0.01); %lon
    k_values(temp_ind,8)=approx_station_loc(j,1);
end

%and color coordinated plot!
scatter(k_values(:,4),k_values(:,3),30,k_values(:,8),'filled')
set(gca,'box','on','fontsize',14)
xlabel('Longitude')
ylabel('Latitude')
%excellent - just the one outlier!

%% hmmm...should double check and take average of mulitple casts per event number...maybe this helps the relationships?
%still leaves them separate on the day though...

%average within an event - this way matches to chl...
eventlist=unique(k_record(:,5));
k_avg=[]; %matched by event number...
k_low=[]; %lowest k recorded for event...
k_high=[]; %highest k_recorded for event...

for q=1:length(eventlist)%length(daylist)
    
    qq=find(cellfun('isempty',strfind(k_record(:,5),eventlist{q}))==0);     %matching by event number
    k_avg=[k_avg; k_values(qq(1),1:4) mean(k_values(qq,5)) length(qq) qq(1) k_values(qq(1),8)]; %average within event num
    k_low=[k_low; k_values(qq(1),1:4) min(k_values(qq,5)) length(qq) qq(1) k_values(qq(1),8)]; %min within event num
    k_high=[k_high; k_values(qq(1),1:4) max(k_values(qq,5)) length(qq) qq(1) k_values(qq(1),8)]; %max within event num

end

k_avg_titles={'year day';'matdate';'lat';'lon';'avg k within event';'number of obs';'index into k_values';'station'};

%% K averages taken over the day, merging events:

daylist=unique(floor(k_values(:,2)));
k_avg_byday=[]; %matched by day

for q=1:length(daylist)
    
    qq=find(k_values(:,2)==daylist(q));
    st=unique(k_values(qq,8)); %average within station
    for j=1:length(st)
        jj=find(k_values(qq,8)==st(j));
        k_avg_byday=[k_avg_byday; k_values(qq(jj(1)),1:4) mean(k_values(qq(jj),5)) length(jj) qq(jj(1)) st(j) ];
        if length(unique(k_record(qq(jj),5))) > 1 % a discrepancy between station number and event number
            disp(st(j))
            disp(unique(k_record(qq(jj),5)))
        end
       
    end
    
end

k_avg_byday_titles={'year day';'matdate';'lat';'lon';'avg k over day';'number of obs';'index into k_values';'station'};


%%  Average casts that were less than 10 yeardays apart:

%excludes outershelf stations 7 and 8
%k_avg is average over event num, keeping them separate on the day

unq_yrdy=unique(k_avg(:,1));
kaggregate=nan(8,4);
kaggregate_hilow=nan(8,4);
count=0;

for q=1:length(unq_yrdy)
    
    if ismember(unq_yrdy(q),[66 72 73 82 255 267])
        switch unq_yrdy(q)
            case 66
                count=count+1;
                qq=find(k_avg(:,1) >= 66 & k_avg(:,1) <= 82 & k_avg(:,8)~=8 & k_avg(:,8)~=7 & k_avg(:,8)~=0); %& k_avg(:,8)~=1); %exclude outer shelf stations
                kaggregate(count,1)=mean([66 72 73 82]);
                kaggregate_hilow(count,1)=mean([66 72 73 82]);
            case 267 %255 only has casts at 7 and 8...
                count=count+1;
                qq=find(k_avg(:,1) >= 255 & k_avg(:,1) <= 267 & k_avg(:,8)~=8 & k_avg(:,8)~=7 & k_avg(:,8)~=0); % & k_avg(:,8)~=1); %exclude outer shelf stations
                kaggregate(count,1)=267;
                kaggregate_hilow(count,1)=267;
        end
        
        kaggregate(count,2)=nanmean(k_avg(qq,5));
        kaggregate(count,3)=nanstd(k_avg(qq,5));
        kaggregate(count,4)=length(qq); 
        
        kaggregate_hilow(count,2)=min(k_low(qq,5));
        kaggregate_hilow(count,3)=max(k_high(qq,5));
        kaggregate_hilow(count,4)=length(qq); 
        
    else
    count=count+1;    
    qq=find(k_avg(:,1)==unq_yrdy(q) & k_avg(:,8)~=8 & k_avg(:,8)~=7 & k_avg(:,8)~=0); % & k_avg(:,8)~=1); %exclude outer shelf stations
    kaggregate(count,1)=unq_yrdy(q);
    kaggregate(count,2)=nanmean(k_avg(qq,5));
    kaggregate(count,3)=nanstd(k_avg(qq,5));
    kaggregate(count,4)=length(qq); 
    
    kaggregate_hilow(count,1)=unq_yrdy(q);
    kaggregate_hilow(count,2)=min(k_low(qq,5));
    kaggregate_hilow(count,3)=max(k_high(qq,5));
    kaggregate_hilow(count,4)=length(qq); 
    
    end

end
%average by month:
% yrdy_list=nan(12,1);
% for mn=1:12
%     yrdy_list(mn)=find_yearday(datenum([num2str(mn) '-1-03']));
% end
% yrdy_list=[yrdy_list; 367];
% 
% %% use k_avg as this is the average within an event:
% 
% mn_avg=nan(12,3);
% for q=1:12
%     qq=find(k_avg(:,1) >= yrdy_list(q) & k_avg(:,1) < yrdy_list(q+1) & k_avg(:,8)~=8 & k_avg(:,8)~=7 & k_avg(:,8)~=0 & k_avg(:,8)~=1); %exclude outer shelf casts
%     mn_avg(q,1)=nanmean(k_avg(qq,5));
%     mn_avg(q,2)=nanstd(k_avg(qq,5));
%     mn_avg(q,3)=length(qq);
% end

%% save vars:


save k_lite k_avg* k_record k_values k_low k_high kaggregate*



%% now examine k, by yearday, location, etc...

subplot(2,3,2,'replace'), hold on

scatter(k_values(:,1),-k_values(:,5),30,k_values(:,8),'filled')
colorbar
set(gca,'fontsize',14,'box','on')
ylabel('Attenuation coefficient, K')
xlabel('Year Day')
xlim([1 366])
title('k for all stations')
% unqdays=unique(k_values(:,1));
% for q=1:length(unqdays)
%     line([unqdays(q) unqdays(q)],ylim,'color',[0.6 0.6 0.6])
% end

subplot(2,3,3,'replace'), hold on
%light gray vertical lines for each sampling day:
xlim([1 366]); ylim([0.15 0.5])
unqdays=unique(k_values(:,1));
for q=1:length(unqdays)
    line([unqdays(q) unqdays(q)],ylim,'color',[0.6 0.6 0.6])
end
ii=find(k_values(:,8)==3 | k_values(:,8)==4 | k_values(:,8)==1 | k_values(:,8)==2);
scatter(k_values(ii,1),-k_values(ii,5),30,k_values(ii,8),'filled')

ii=find(k_avg(:,8)==3 | k_avg(:,8)==4);
plot(k_avg(ii,1),-k_avg(ii,5),'ko')
ii=find(k_avg_byday(:,8)==3 | k_avg_byday(:,8)==4);
plot(k_avg_byday(ii,1),-k_avg_byday(ii,5),'ro')

%scatter(k_avg(ii,1),-k_avg(ii,5),30,k_avg(ii,7),'filled')


caxis([0 8]), colorbar
set(gca,'fontsize',14,'box','on')
ylabel('Attenuation coefficient, K')
xlabel('Year Day')
xlim([1 366])
title('k for just tower and node locations')

%%

subplot(1,2,2,'replace'), hold on


% tt=find(k_avg(:,8)~=7 & k_avg(:,8)~=8 & k_avg(:,8)~=0 & k_avg(:,8)~=1);
% scatter(k_avg(tt,1),-k_avg(tt,5),30,k_values(tt,8),'filled')

%plot the averages: with connections:
%calc slope...
% m1=kaggregate(end,2);
% m2=kaggregate(1,2);
% m=(m2-m1)/(365-347 + 16); 
% x=[1; kaggregate(:,1); 365];
% y=[-(m*(366-347)+m1); -kaggregate(:,2); -(m*(365-347)+m1)];

%or just pad:
x=[kaggregate(end,1)-365; kaggregate(:,1); kaggregate(1,1)+365];
y=[-kaggregate(end,2); -kaggregate(:,2); -kaggregate(1,2)];
plot(x, y,':','color',[0.3 0.3 0.3],'linewidth',2)

x=[kaggregate_hilow(end,1)-365; kaggregate_hilow(:,1); kaggregate_hilow(1,1)+365];
y=[-kaggregate_hilow(end,2); -kaggregate_hilow(:,2); -kaggregate_hilow(1,2)];
plot(x, y,':','color',[0.6 0.6 0.6],'linewidth',2)

x=[kaggregate_hilow(end,1)-365; kaggregate_hilow(:,1); kaggregate_hilow(1,1)+365];
y=[-kaggregate_hilow(end,3); -kaggregate_hilow(:,3); -kaggregate_hilow(1,3)];
plot(x, y,':','color',[0.6 0.6 0.6],'linewidth',2)

%plot(kaggregate(:,1), -kaggregate(:,2),'.','markersize',12,'color',[0.5 0.5 0.5])

tt=find(k_values(:,8)~=7 & k_values(:,8)~=8 & k_values(:,8)~=0 & k_values(:,8)~=1);
scatter(k_values(tt,1),-k_values(tt,5),30,k_values(tt,8),'filled')

colorbar
caxis([1 8])
set(gca,'fontsize',14,'box','on')
ylabel('Attenuation coefficient, K')
xlabel('Year Day')
xlim([1 366])
title('k for all stations')

colormap jet