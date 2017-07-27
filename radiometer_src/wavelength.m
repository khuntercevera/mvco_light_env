% Script to look at how the spectral quality of light changes with time:

%% A time series plot of K_d values!

%sourcepath=fullfile('\\sosiknas1\lab_data\mvco\HyperPro_Radiometer\processed_radiometer_files\'); % path to folders with raw data...
sourcepath=fullfile('/Volumes/Lab_data/MVCO/HyperPro_Radiometer/processed_radiometer_files/');

%How many of these do we have?
d = dir(sourcepath);
isub = [d(:).isdir]; %# returns logical vector
foldernames = {d(isub).name}'; %names of folders
temp=regexp(foldernames,'\d{1,2}\w{3,4}20\d{2}'); %find all folders with this naming scheme
datafolders = foldernames(cellfun('isempty',temp)==0);

%%
cast_record=[];

for foldernum=[2:3 5:10 12 14:16 18:20]; %are these the ones with lat/lon?

    %find dat data!
    matsource=fullfile(sourcepath,datafolders{foldernum},'/mat_outfiles/');
    
    eval(['load ' matsource 'data_' datafolders{foldernum} '.mat'])
    eval(['tempdata=data_' datafolders{foldernum} ';'])    
    
    eval(['load ' matsource 'K_PAR_' datafolders{foldernum} '.mat'])
    eval(['K_PAR=K_PAR_' datafolders{foldernum} ';'])

    eval(['load ' matsource 'location_' datafolders{foldernum} '.mat'])
    eval(['location=location_' datafolders{foldernum} ';'])

    %timestamp={tempdata.timestamp}'; timestamp=cell2mat(timestamp);
    %INCORRECT DATES LISTED HERE!!!
    
    timestamp=repmat(datenum(datestr(datafolders{foldernum})),1,length(tempdata))';
    if foldernum ~= 18
        lat={location.lat}'; 
        lon={location.lon}';
        station={location.station}';  
    else
        lat=cell(length(timestamp),1);
        lon=lat;
        station=lat;
    end
    
    %now that we have the data in, let's look at those wavelenghts!
    
    %find the depth for 4m, get the wave distribution, then keep and find
    %max...
    
    wv=data_06Jul2011(1).wavelen_downwell;
    edl=data_06Jul2011(1).adj_edl;
    depth=data_06Jul2011(1).depth;
    edl_ind=data_06Jul2011(1).edl_ind;
    
    figure, hold on
    depthcolor=hsv(111);
    for j=1:size(edl,1)
        scatter(wv,edl(j,:),40,depth(edl_ind(j))*ones(size(wv)),'filled')
    end
    
    %need downcast colors!!!
end


