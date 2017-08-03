%Figure out the latitudes and longitudes for each cast
%Sometimes this is recorded in the header of each file, but sometimes is not
%Also check what the MVCO event log says of lat/lon (which should hopefully match any header information)

%find the data folders:
sourcepath=fullfile('/Volumes/Lab_data/MVCO/HyperPro_Radiometer/processed_radiometer_files/');

d = dir(sourcepath);
isub = [d(:).isdir]; %# returns logical vector
foldernames = {d(isub).name}'; %names of folders
temp=regexp(foldernames,'\d{1,2}\w{3,4}20\d{2}'); %find all folders with this naming scheme
datafolders = foldernames(cellfun('isempty',temp)==0);

%To avoid chasing down data for empty casts, identify the folders that have good data in them:
good_data=[];
bad_data=[];
for foldernum=1:length(datafolders)
    
    eval(['load ' sourcepath datafolders{foldernum} '/mat_outfiles/data_' datafolders{foldernum} '.mat'])
    eval(['tempdata=data_' datafolders{foldernum} ';'])
    
    if any(cellfun(@(x) any(x==0),{tempdata(:).emptyflag})==1) % meaning at least one good cast
        disp('Found good casts...')
        good_data=[good_data; foldernum];
    else
        disp('Only empty files for this folder...')
        bad_data=[bad_data; foldernum];
    end
end

%% Okay, the plan is to ...
%loop through good folders and see if lat/lon was recorded.
% Also see if lat/lon was recorded in MVCO event log and hopefully these
% compare with the entered values!

%to place lat/lon and station that goes along with each cast, use info in nutrient file:
%load /Users/kristenhunter-cevera/Documents/MATLAB/MVCO_Syn_analysis/MVCO_environmental_data/nut_data_repsB.mat

%trying out 'readtable' - can also be accomplished by textscan if using
%older version of MATLAB:
mvco_events=readtable('/Users/kristenhunter-cevera/MVCO_light_at_depth/MVCO_event_log.txt'); %results in a table structure
%see list of variables: mvco_events.Properties.VariableNames
mvco_event_num=mvco_events{:,{'Event_Number'}};
mvco_time=datenum(mvco_events{:,{'Start_Date'}});
mvco_lon=mvco_events{:,{'Longitude'}};
mvco_lat=mvco_events{:,{'Latitude'}};
mvco_start_time=mvco_events{:,{'Start_Time_UTC'}};
mvco_depth=mvco_events{:,{'Water_Depth_m'}};


%% This is a very manual process....

for foldernum=good_data(5)' %files that have data in them!
    %%
    matsource=fullfile(sourcepath,datafolders{foldernum},'/mat_outfiles/');
    
    eval(['load ' matsource 'data_' datafolders{foldernum} '.mat'])
    eval(['tempdata=data_' datafolders{foldernum} ';'])
    
    %Okay, would like to take a look at the casts all together, with info from MVCO data...
    if foldernum==7
        disp('adjusting mprtime for these files')
        for filenum=1:length(tempdata)
            tempdata(filenum).mprtime = tempdata(filenum).mprtime-(1/24); %suspected hour offset with UTC time
        end
    end
    
    figure(21), clf, hold on
    for filenum=1:length(tempdata)
        plot(tempdata(filenum).mprtime, tempdata(filenum).depth,'.-')
        %text(mean(tempdata(filenum).mprtime),-2,{(templat);(templon)})
    end
    set(gca,'Ydir','reverse')
    datetick
    
    ii=find(mvco_time==floor(tempdata(1).timestamp)); %find the dates that correspond
    localt=datenum([repmat('1-0-03 ',length(ii),1) char(mvco_start_time(ii))])-datenum('1-0-03')+floor(tempdata(1).timestamp);
    plot(localt,mvco_depth(ii),'p','markersize',12)
    set(gca,'xlim',[min(localt(1),tempdata(1).timestamp) max(localt(end),tempdata(end).timestamp)],'xgrid','on','ygrid','on')
    %dsiplay for reference:
    mvco_events(ii,{'Start_Date','Start_Time_UTC','Latitude','Longitude','Water_Depth_m'})
    datetick
    
    for q=1:length(tempdata)
        %%
        location(q).file = tempdata(q).file;
        
        if tempdata(q).emptyflag~=1
            
            %any recorded values:
            templon=regexp(tempdata(q).lon,'(?<degree>\d{2})\s(?<decmin>\d*\.\d*)','names');
            templat=regexp(tempdata(q).lat,'(?<degree>\d{2})\s(?<decmin>\d*\.\d*)','names');
            
            %find closest mvco values with regard to time and cast:
            [temp, itime]=sort(localt-tempdata(q).timestamp);
            
            %maybe look for the cast right before and after radiometer cast:
            temp(temp<0)=0; temp(temp>0)=1;
            if ~any(temp==1) %meaning all zeros -> end of the day casts
                itime=itime([length(temp) length(temp)-1]); 
            elseif ~any(temp==0) %meaing all ones -> beginning of day casts
                itime=itime(1:2);
            else
                ww=find(diff(temp==1));
                itime=itime(max(ww-1,1):min(ww+1,length(temp))); %for the cases where first and last cast
            end
            
            %[~, itime]=sort(abs(localt-tempdata(q).timestamp));
            
            max_depth=max(tempdata(q).depth);
            idep=find(mvco_depth(ii)+0.25-max_depth > 0);
            
            %only evaluate indices with depths that are valid:
            itime=intersect(itime,idep,'rows','stable');
            
            if ~isempty(templon) && ~isempty(templat) %then there are few set of logic statements that we can undertake:
                
                templon=-[str2num(templon.degree) + str2num(templon.decmin)/60]; %convert to decimal
                templat=[str2num(templat.degree) + str2num(templat.decmin)/60];
                
                %also find the closet values with respect to recorded event lat/lon:
                [~, ilat]=sort(abs(templat-mvco_lat(ii)));
                [~, ilon]=sort(abs(templon-mvco_lon(ii)));
                
                ilat=intersect(itime,ilat,'rows','stable');
                ilon=intersect(itime,ilon,'rows','stable');
                
                %so, then we could just choose the mode out of these...
                %Most likely cast:
                if ilat(1) == ilon(1)
                    im=ilat(1);
                else
                    disp(['hmmm...for ' num2str(q) ' out of: ' num2str(length(tempdata)) ' seem to have a problem finding the corresponding cast'])
                    keyboard
                end
                
                %THE MOST LIKELY CORRESPONDING CAST:
                mvcotemplat=mvco_lat(ii(im));
                mvcotemplon=mvco_lon(ii(im));
                plot(localt(im),mvco_depth(ii(im)),'p','markersize',12,'markerfacecolor','b')
                
                if abs(mvcotemplat-templat) < 0.01 && abs(mvcotemplon-templon) < 0.01 %close enough!
                    location(q).lat=templat;
                    location(q).lon=templon;
                    location(q).notes=['lat/lon from hdr used; matches event record values with closest record: ' mvco_event_num{ii(im)}];
                    disp(['for file: ' num2str(q) ' out of: ' num2str(length(tempdata)) ' lat/lon entered with good agreement between record values, with closest record: ' mvco_event_num{ii(im)}])
                else
                    disp(['for file: ' num2str(q) ' out of: ' num2str(length(tempdata)) ', not a good match between entered and record data'])
                    keyboard
                    decision=input('Do you trust entered log lat/lon or recorded mvco lat/lon? Enter either: "MVCO" or "log"\n');
                                        
                    if strcmp(decision,'MVCO')
                        location(q).lat=mvcotemplat;
                        location(q).lon=mvcotemplon;
                        location(q).notes=['lat/lon entered but seems incorrect compared to log; using mvco event: ' mvco_event_num{ii(im)} ' lat/lon for position'];
                    elseif strcmp(decision,'log')
                        
                        location(q).lat=templat;
                        location(q).lon=templon;
                        location(q).notes=['lat/lon entered does not match closest entry compared to mvco event records ' mvco_event_num{ii(im)} ', but using entered lat/lon for position'];
                        
                    end
                     
                end
                
            else %no entered lat/lon -> use data from MVCO event number:
                
                %THE MOST LIKELY CORRESPONDING CAST:
                mvcotemplat=mvco_lat(ii(itime(1)));
                mvcotemplon=mvco_lon(ii(itime(1)));
                
                in=itime(1);
                plot(localt(itime(1)),mvco_depth(ii(itime(1))),'p','markersize',12,'markerfacecolor','b')
                
                if    1
                else
                    disp(['hmmm...for ' num2str(q) ' out of: ' num2str(length(tempdata)) ' seem to have a problem finding the corresponding cast'])
                    keyboard
                end
                
                location(q).lat=mvcotemplat;
                location(q).lon=mvcotemplon;
                location(q).notes=['no lat/lon entered; using mvco event: ' mvco_event_num{ii(in)} ' lat/lon for position'];
                disp(['for file: ' num2str(q) ' out of ' num2str(length(tempdata)) ' no lat/lon entered, using mvco record values'])
                
            end
            
        else %empty file
            location(q).lat=[];
            location(q).lon=[];
            location(q).notes='empty_file';
        end
        
        keyboard
        
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