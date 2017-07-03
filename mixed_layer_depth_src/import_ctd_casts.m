%script to import all the Tioga processed data from the SurveyCruises
%folder:

%cast_type='Dbin'; %Ubin %files can be processed seqparately as just
%downward or upward casts,
%BUT these seem to have strange time jump errors and are 0.5m binned...not
%sure what processing is taking place!
cast_type='both'; %if want both the downcast and upcast, more raw data

pathname = '/Volumes/TaylorF/from_Samwise/data/MVCO/SurveyCruises/'; %once maddie is mounted
%pathname = '/Volumes/J_data/MVCO/SurveyCruises/';

folder_list=dir(pathname);
folder_names=extractfield(folder_list,'name'); %cell array of folder names

%find the Tioga folders:
ii=find(cellfun('isempty',strfind(folder_names','Tioga'))==0);

TEMP=struct('cruise_name',{},'cruise_data',{});

for q=1:length(ii)
    
    if strcmp(cast_type,'both')
            filelist=dir([pathname folder_names{ii(q)} '/CTD/converted/*.cnv']); %Ubin.cnv, have a choice of extracting downward or upward cast (if processed)
        filepath=[pathname folder_names{ii(q)} '/CTD/converted/'];
    else
    filelist=dir([pathname folder_names{ii(q)} '/CTD/binned/*' cast_type '.cnv']); %Ubin.cnv, have a choice of extracting downward or upward cast (if processed)
    filepath=[pathname folder_names{ii(q)} '/CTD/binned/'];
    end
    disp(folder_names{ii(q)})
    
    TEMP(q).cruise_name=folder_names{ii(q)};
    if ~isempty(filelist)
        temp=struct('filename',{},'lat',{},'lon',{},'UTC',{},'upload_time',{},'col_headers',{},'data',{});
        
        for j=1:length(filelist)
            filename=filelist(j).name;
            [lat,lon, UTC_time,upload_time,header,data]=import_cnvfile([filepath filename]);
            
            temp(j).filename = filename;
            temp(j).lat=lat;
            temp(j).lon=lon;
            temp(j).UTC=UTC_time;
            temp(j).upload_time=upload_time;
            temp(j).col_headers={header};
            temp(j).data={data};
            
        end
    else
        keyboard
    end
    
    TEMP(q).cruise_data=temp;
end


%% Extract data from structure file and calculat potential density:
%Slightly clunky:
addpath /Users/kristenhunter-cevera/Dropbox/MVCO_mixed_layer_depth/seawater_ver3_2/

position=[];
file_time={};

pos_titles={'lat_deg','lat_min','lat_sec','lon_deg','lon_min','lon_sec'};

data={}; hdr={};
for q=1:length(TEMP)
    
    temp=TEMP(q).cruise_data; %temp is now the the structure that contains the data for each cast on that cruise
    for j=1:length(temp)
        
        file_time=[file_time; {temp(j).filename} {temp(j).UTC} {temp(j).upload_time}];
        lat=temp(j).lat;
        lon=temp(j).lon;
        
        if ~isempty(lat) && ~isempty(lon)
            position=[position;  str2num(lat.deg) str2num(lat.min) str2num(lat.sec) str2num(lon.deg) str2num(lon.min) str2num(lon.sec)];
        elseif isempty(lat) && ~isempty(lon)
            position=[position;  nan(1,3) str2num(lon.deg) str2num(lon.min) str2num(lon.sec)];
        elseif isempty(lon) && ~isempty(lat)
            position=[position; str2num(lat.deg) str2num(lat.min) str2num(lat.sec) nan(1,3) ];
        else
            position=[position;  nan(1,6)];
        end
        
        %calculate and potential density:
        temp_data=temp(j).data{:};
        temp_hdr=temp(j).col_headers{:};
        if any(~cellfun('isempty',temp_data)) %if non-empty cell exists:
            %depth=temp_data{1}; press=temp_data{2}; temperature=temp_data{3}; sal=temp_data{5};
            pdens=sw_pden(temp_data{5},temp_data{3},temp_data{2},0);
            temp_data{end+1}=pdens;
            temp_hdr{end+1}='Potential Density';
        end
        %now extract the data:
        data=[data; {temp_data}];
        hdr=[hdr; {temp_hdr}];
        
    end
end

if strcmp(cast_type,'Ubin') %default is downcasts
    TEMP_upcast=TEMP;
    data_upcast=data;
    hdr_upcast=hdr;
    position_upcast=position;
end

%% find only trips to tower or node:

tower_ind=find((position(:,2) == 19 | position(:,2) == 20) & (position(:,5) == 33 | position(:,5) == 34));

% in general how to access data out of the structure array:
%TEMP(index).cruise_data.field;


