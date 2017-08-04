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

%trying out 'readtable' - can also be accomplished by textscan if using older version of MATLAB:
mvco_events=readtable('/Users/kristenhunter-cevera/MVCO_light_at_depth/MVCO_event_log.txt'); %results in a table structure
%see list of variables: mvco_events.Properties.VariableNames
mvco_event_num=mvco_events{:,{'Event_Number'}};
mvco_time=datenum(mvco_events{:,{'Start_Date'}});
mvco_lon=mvco_events{:,{'Longitude'}};
mvco_lat=mvco_events{:,{'Latitude'}};
mvco_start_time=mvco_events{:,{'Start_Time_UTC'}};
mvco_depth=mvco_events{:,{'Water_Depth_m'}};


%% This is a very manual process....

for foldernum=good_data(15)' %files that have data in them! Go through one by one....
    
    matsource=fullfile(sourcepath,datafolders{foldernum},'/mat_outfiles/');
    
    eval(['load ' matsource 'data_' datafolders{foldernum} '.mat'])
    eval(['tempdata=data_' datafolders{foldernum} ';'])
    
    %It turns out there are some discrepancies in UTC offset, timestamps
    %and records. So far, these have been identified in:
    
    switch datafolders{foldernum}
        case '12Aug2009'
            hour_offset=-(1/24);
            timestamp_offset=-(1/24);
        case '21Oct2009'
            hour_offset=-(1/24);
            timestamp_offset=-(1/24);
        case '23Mar2011'
            hour_offset=-(1/24);
            timestamp_offset=-(1/24);
        case '13Mar2008'
            hour_offset=(1/24);
            timestamp_offset=(1/24);
        case '16Jan2008'
            hour_offset=(1/24);
            timestamp_offset=(1/24);
        case '13Dec2007'
            hour_offset=(1/24);
            timestamp_offset=1+(1/24);
        case '27Apr2009'
            timestamp_offset=-(1/24);
            hour_offset=-(1/24);
        case '4June2009'
            timestamp_offset=-(1/24);
            hour_offset=-(1/24);
        otherwise
            hour_offset=0;
            timestamp_offset=0;
    end
    
    if hour_offset~=0 || timestamp_offset~=0
        for filenum=1:length(tempdata)
            tempdata(filenum).mprtime = tempdata(filenum).mprtime+hour_offset; %suspected hour offset with UTC time
            tempdata(filenum).timestamp = tempdata(filenum).timestamp+timestamp_offset;
        end
    end
    
    ii=find(mvco_time==floor(tempdata(end).timestamp)); %find the dates that correspond
    localt=datenum([repmat('1-0-03 ',length(ii),1) char(mvco_start_time(ii))])-datenum('1-0-03')+floor(tempdata(end).timestamp);
    
    figure(21), clf, hold on
    for filenum=1:length(tempdata)
        plot(tempdata(filenum).mprtime, tempdata(filenum).depth,'.-')
        %text(mean(tempdata(filenum).mprtime),-2,{(templat);(templon)})
    end
    
    plot(localt,mvco_depth(ii),'p','markersize',12)
    title(datafolders{foldernum})
    set(gca,'Ydir','reverse','xlim',[min(localt(1),tempdata(1).timestamp) max(localt(end),tempdata(end).timestamp)],'xgrid','on','ygrid','on')
    datetick
    
    %dsiplay for reference:
    mvco_events(ii,{'Start_Date','Event_Number','Start_Time_UTC','Latitude','Longitude','Water_Depth_m'})
    
    for q=1:length(tempdata)
        
        location(q).file = tempdata(q).file;
        
        if tempdata(q).emptyflag~=1
            
            %any recorded values:
            templon=regexp(tempdata(q).lon,'(?<degree>\d{2})\s(?<decmin>\d*\.\d*)','names');
            templat=regexp(tempdata(q).lat,'(?<degree>\d{2})\s(?<decmin>\d*\.\d*)','names');
            
            %find closest mvco values with regard to time and cast:
            [~, itime_abs]=sort(abs(localt-tempdata(q).timestamp)); %absolute in time closest values
            [temp, ~]=sort(localt-tempdata(q).timestamp); %sorted to find inbetween timepoints
            
            %maybe look for the cast right before and after radiometer cast:
            temp(temp<0)=0; temp(temp>0)=1;
            if ~any(temp==1) %meaning all zeros -> end of the day casts
                itime=itime_abs(1:2);
                %itime=itime([length(temp) length(temp)-1]);
            elseif ~any(temp==0) %meaing all ones -> beginning of day casts
                itime=itime_abs(1:2);
            else
                ww=find(diff(temp==1));
                rr=(max(ww-1,1):min(ww+1,length(temp))); %range
                %itime=itime(max(ww-1,1):min(ww+1,length(temp))); %for the cases where first and last cast
                itime=intersect(itime_abs,rr','rows','stable');
            end
            
            %[~, itime]=sort(abs(localt-tempdata(q).timestamp));
            
            max_depth=max(tempdata(q).depth);
            idep=find(mvco_depth(ii)+0.8-max_depth > 0);
            
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
                    disp(['for file: ' num2str(q) ' out of: ' num2str(length(tempdata)) ', not a good match between entered and record data: ' mvco_event_num{ii(im)}])
                    keyboard
                    decision=input('Do you trust entered log lat/lon or recorded mvco lat/lon? Enter either: "MVCO" or "log" or "other" to return to keyboard\n');
                    
                    if strcmp(decision,'MVCO')
                        location(q).lat=mvcotemplat;
                        location(q).lon=mvcotemplon;
                        location(q).notes=['lat/lon entered but seems incorrect compared to log; using mvco event: ' mvco_event_num{ii(im)} ' lat/lon for position'];
                    elseif strcmp(decision,'log')
                        
                        location(q).lat=templat;
                        location(q).lon=templon;
                        location(q).notes=['lat/lon entered does not match closest entry compared to mvco event records ' mvco_event_num{ii(im)} ', but using entered lat/lon for position'];
                    elseif strcmp(decision,'other')
                        keyboard
                    end
                    
                end
                
            else %no entered lat/lon -> use data from MVCO event number:
                
                %THE MOST LIKELY CORRESPONDING CAST:
                mvcotemplat=mvco_lat(ii(itime(1)));
                mvcotemplon=mvco_lon(ii(itime(1)));
                
                plot(localt(itime(1)),mvco_depth(ii(itime(1))),'p','markersize',12,'markerfacecolor','b')
                
                location(q).lat=mvcotemplat;
                location(q).lon=mvcotemplon;
                
                location(q).notes=['no lat/lon entered; using mvco event: ' mvco_event_num{ii(itime(1))} ' lat/lon for position'];
                disp(['for file: ' num2str(q) ' out of ' num2str(length(tempdata)) ' no lat/lon entered, using mvco record value: ' mvco_event_num{ii(itime(1))}])
                
            end
            
        else %empty file
            disp('empty file')
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


%% And now a quick plot to show where the data is!


%If need to reload data:
location_mat=[]; location_rec={};

for foldernum=good_data'
    
    matsource=fullfile(sourcepath,datafolders{foldernum},'/mat_outfiles/');
    eval(['load ' matsource 'location_' datafolders{foldernum} '.mat'])
    eval(['location=location_' datafolders{foldernum} ';'])
    
    location_rec=[location_rec; repmat({datafolders{foldernum}},length(location),1) repmat({foldernum},length(location),1) {location(:).file}'  {location(:).notes}'];
    location_mat=[location_mat; [cell2mat({location(:).lat}')  cell2mat({location(:).lon}')]];
    
end

jj=find(strcmp(location_rec(:,4),'empty_file')~=1);
location_rec=location_rec(jj,:);

%% and plot!

clf
plot(location_mat(:,2),location_mat(:,1),'.'), hold on

%% hmmm, if some do look suspicious:
[x y]=ginput;
hold on
plot(x,y,'r.')

%% now let's find those points and see if want to adjust:

q=6; %how many points you need to correct:
plot(x(q),y(q),'kp','markersize',12)

qq=find((1e-2)*floor(100*location_mat(:,2))==((1e-2)*floor(100*x(q))) & (1e-2)*floor(100*location_mat(:,1))==((1e-2)*floor(100*y(q))));

disp([location_rec{qq,1:2} location_rec{qq,4}])

foldername=location_rec{qq,1};
eval(['templist={location_' foldername '(:).file}'';'])

w=find(strcmp(templist,location_rec{qq,3})==1);

mvconum=regexp(location_rec{qq,4},'MVCO_\d{3}','match');
jj=find(strcmp(mvco_event_num,mvconum)==1);

%%
eval(['location_' foldername  '(w).lat=mvco_lat(jj);'])
eval(['location_' foldername '(w).lat=mvco_lat(jj);'])
eval(['location_' foldername '(w).notes=''lat/lon entered but seems incorrect compared to log; using mvco event: ' mvconum{:} ' lat/lon for position'';'])

plot(mvco_lon(jj),mvco_lat(jj),'rp','markersize',12)

%% If that looks good, can save:

eval(['save ' matsource 'location_' datafolders{foldernum} ' location_' datafolders{foldernum}])

