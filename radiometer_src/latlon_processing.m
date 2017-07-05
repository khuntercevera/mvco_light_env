%Figure out the latitudes and longitudes for each cast

%find the data folders:
sourcepath=fullfile('/Volumes/Lab_data/MVCO/HyperPro_Radiometer/processed_radiometer_files/');

d = dir(sourcepath);
isub = [d(:).isdir]; %# returns logical vector
foldernames = {d(isub).name}'; %names of folders
temp=regexp(foldernames,'\d{1,2}\w{3,4}20\d{2}'); %find all folders with this naming scheme
datafolders = foldernames(cellfun('isempty',temp)==0);


%to place lat/lon and station that goes along with each cast, use info in nutrient file:
load /Users/kristenhunter-cevera/Documents/MATLAB/MVCO_Syn_analysis/MVCO_environmental_data/nut_data_repsB.mat

%Search this log for dates that match the folder dates:
MVCO_eventdays=cellfun(@(x) datenum(x),MVCO_nut_reps(:,3));

mvco_match={};
for q=1:length(datafolders)
    tempday=datenum(datafolders{q});
    ii=find(MVCO_eventdays==tempday);
    
    if ~isempty(ii)
        mvco_match=[mvco_match; cellstr(repmat(datestr(tempday),length(ii),1)) MVCO_nut_reps(ii,[1 5:6]) num2cell(station4nut(ii))];
    else
        mvco_match=[mvco_match; {datestr(tempday)} cell(1,4)];
    end
    
end


%excellent - the ones that don't match are empty anyways!
good_data=[2:3 5:10 12 14:16 18:20];
bad_data=setxor(good_data,1:20);


%% Okay, now walk through each day of data - see if lat/lons match to MVCO event data
%If so, record station number
%If not, try to piece that back together!


for foldernum=[7:10 12 14:16 18:20] %[2:3 5:10 12 14:16 18:20]; %files that have data in them!
    
    matsource=fullfile(sourcepath,datafolders{foldernum},'/mat_outfiles/');
    
    eval(['load ' matsource 'data_' datafolders{foldernum} '.mat'])
    eval(['tempdata=data_' datafolders{foldernum} ';'])
    
    ii=find(strcmp(datestr(datafolders{foldernum}),mvco_match(:,1))==1);
    mvco_latlon=cell2mat(mvco_match(ii,3:5));
    
    for q=1:length(tempdata)
        
        location(q).file = tempdata(q).file;
        
        templon=regexp(tempdata(q).lon,'(?<degree>\d{2})\s(?<decmin>\d*\.\d*)','names');
        templat=regexp(tempdata(q).lat,'(?<degree>\d{2})\s(?<decmin>\d*\.\d*)','names');
        
        if ~isempty(templon) & ~isempty(templat)
            templon=-[str2num(templon.degree) + str2num(templon.decmin)/60];
            templat=[str2num(templat.degree) + str2num(templat.decmin)/60];
            
            location(q).lat=templat;
            location(q).lon=templon;
            
            %find the station: look at the closet distances in the record?
            [~, ia]=min(abs(templat-mvco_latlon(:,1)));
            [~, io]=min(abs(templon-mvco_latlon(:,2)));
            
            if mvco_latlon(io,3) ~= mvco_latlon(ia,3) %if stations do not match, then have a problem...
                location(q).station=[];
            else
                location(q).station=mvco_latlon(io,3);
            end
                
        else %need to look at mvco records...
            disp(['Need to look at the MVCO records for this one! ' tempdata(q).file '  ' datafolders{foldernum} '  ' num2str(foldernum)])   
        end
        
    end %files within tempdata
    
    eval(['location_' datafolders{foldernum} '=location;'])
    eval(['save ' matsource 'location_' datafolders{foldernum} ' location_' datafolders{foldernum}])
    clear location
    
end %foldernum


%% for the ones that don't have lat/lon - entry manually... 

%23Sept2008  18 CANNOT FIND WHICH MVCO EVENTS WENT WITH THESE CASTS! - CHECK
%THE SHIP LOGS! string of single casts, likely at just one station, but day was multi-station

foldernum=19;

matsource=fullfile(sourcepath,datafolders{foldernum},'/mat_outfiles/');
eval(['load ' matsource 'data_' datafolders{foldernum} '.mat'])
eval(['tempdata=data_' datafolders{foldernum} ';'])
     
clf, hold on
for filenum=1:length(tempdata)
    plot(tempdata(filenum).mprtime, tempdata(filenum).depth,'.-')
end
set(gca,'Ydir','reverse')
datetick

%%
rr=14;
for j=rr
plot(tempdata(j).mprtime, tempdata(j).depth,'k.-')
end

%%
event='MVCO_206';
ii=find(strcmp(event,MVCO_nut_reps(:,1))==1);
station4nut(ii)

%%
for q=rr
    location(q).station = station4nut(ii(1));
    location(q).lon = MVCO_nut_reps{ii(1),6};
    location(q).lat = MVCO_nut_reps{ii(1),5};
    location(q).file = tempdata(q).file;
    location(q).notes= ['from event info ' event '; matched based on depth and event# order'];
end

%%
eval(['location_' datafolders{foldernum} '=location;'])
eval(['save ' matsource 'location_' datafolders{foldernum} ' location_' datafolders{foldernum}])
clear location