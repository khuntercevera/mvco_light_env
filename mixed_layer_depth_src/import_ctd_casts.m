%script to import all the Tioga processed CTD data (.cnv files) from the
%processed_CTD_data folder:

addpath /Users/kristenhunter-cevera/MVCO_light_at_depth/seawater_ver3_2/
pathname='/Volumes/Lab_data/MVCO/processed_CTD_casts/';
%pathname = '/Volumes/TaylorF/from_Samwise/data/MVCO/SurveyCruises/'; %once maddie is mounted
%pathname = '/Volumes/J_data/MVCO/SurveyCruises/';

filelist=dir(pathname);
filenames=extractfield(filelist,'name'); %cell array of folder names
temp=regexp(filenames,'\.cnv'); %find only the .cnv files
ii=find(cellfun('isempty',temp)==0);
filenames=filenames(ii)';

%find the Tioga folders:
%ii=find(cellfun('isempty',strfind(folder_names','Tioga'))==0);
load(fullfile(pathname,'list_and_location_of_raw_ctd_files.mat'))

CTD=struct('cast_name',{},'file_location',{},'lat',{},'lon',{},'UTC',{},'upload_time',{},'col_headers',{},'data',{});

for q=1:length(filenames)
    
    disp(filenames{q})
    
    CTD(q).cast_name=filenames{q};
    
    jj=find(cellfun('isempty',regexp(datafiles,filenames{q}(1:end-4)))==0);
    CTD(q).file_location=datafiles{jj};
    
    [lat,lon, UTC_time,upload_time,header,data]=import_cnvfile([pathname filenames{q}]);
    
    if ~isempty(lat) || ~isempty(lon)            
        CTD(q).lat=[str2num(lat.deg)+str2num(lat.min)/60+str2num(lat.sec)/3600];
        CTD(q).lon=-[str2num(lon.deg)+str2num(lon.min)/60+str2num(lon.sec)/3600];
    else
        CTD(q).lat=NaN;
        CTD(q).lon=NaN;
    end
    
    CTD(q).UTC=UTC_time;
    CTD(q).upload_time=upload_time;
    
    %calculate and potential density:
    if any(~cellfun('isempty',data)) %if non-empty cell exists:
        
        %make sure you've got the right headers:
        s=find(cellfun('isempty',regexp(header,'Salinity'))==0);
        t=find(cellfun('isempty',regexp(header,'Temperature'))==0);
        p=find(cellfun('isempty',regexp(header,'Pressure'))==0);
        
        pdens=sw_pden(data{s},data{t},data{p},0); %0 refers to reference pressure
        data{end+1}=pdens;
        header{end+1}='Potential Density';
    end
    
    CTD(q).data_hdr=header;
    CTD(q).data=data;
    
end


%% How many files don't have lat/lon?

ii=find(cellfun('isempty',{CTD(:).lat}')==1);
{CTD(ii).cast_name}'
%Okay, only about 8 or so actual cruises where we are missing this data...


%% find only trips to tower or node:

templat=cell2mat({CTD(:).lat}');
templon=cell2mat({CTD(:).lon}');

box_tower=[-70.58 -70.53 41.315 41.33];
box_node=[-70.58 -70.53 41.33 41.345];

mvco_ind=find(templat > 41.3 & templat < 41.35 & templon < -70.53 & templon > -70.60);
%this should be about ~99 casts...

%% and if curious, see where all the casts come from...

plot(templon,templat,'.')
patch([-70.53 -70.53 -70.60 -70.60],[41.3 41.35 41.35 41.3],'k','facecolor','none')