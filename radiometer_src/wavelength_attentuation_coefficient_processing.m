% %% Calculation of attentuation coefficients for each k for each lambda :)


%sourcepath=fullfile('\\sosiknas1\lab_data\mvco\HyperPro_Radiometer\processed_radiometer_files\'); % path to folders with raw data...
sourcepath=fullfile('/Volumes/Lab_data/MVCO/HyperPro_Radiometer/');
processed_path=fullfile(sourcepath,'/processed_radiometer_files/');

d = dir(processed_path);
isub = [d(:).isdir]; %# returns logical vector
foldernames = {d(isub).name}'; %names of folders
temp=regexp(foldernames,'\d{1,2}\w{3,4}20\d{2}'); %find all folders with this naming scheme
datafolders = foldernames(cellfun('isempty',temp)==0);
%there should be (as of summer 2017) only 20 folders...

%load in the list of good data folders:
load(fullfile(sourcepath,'good_data_folders.mat'))

%% right now, manually entered, but can put into a for loop

foldernum=good_data(1);
load(fullfile(processed_path,datafolders{foldernum},['/mat_outfiles/data_' datafolders{foldernum} '.mat']))
eval(['tempdata=data_' datafolders{foldernum} ';'])
load(fullfile(processed_path,datafolders{foldernum},['/mat_outfiles/K_PAR_' datafolders{foldernum}]));
eval(['K_PAR=K_PAR_' datafolders{foldernum} ';'])


for filenum=1:length(K_PAR);
    
    if K_PAR(filenum).flag==0;  
        
        filename=tempdata(filenum).file;
        depth=tempdata(filenum).depth;
        adj_edl=tempdata.adj_edl; %dark corrected downwelling light amounts for each wavelength with depth
        wavelen_downwell=tempdata_wavelen_downwell; %corresponding wavelengths
        edl_PAR=tempdata(filenum).edl_PAR; %calculated PAR
        edl_ind=tempdata(filenum).edl_ind; %index that matches to MPR sensor (depth)

        K=K_PAR(filenum).K;
        impr=K_PAR(filenum).depth_index; %index that matches to downcast and used to do PAR regression
        ipar=K_PAR(filenum).par_index; %index that matches to downcast and used to do PAR regression

        %calculate k for each wavelength and record:
        k_wv=nan(size(wavelen_downwell,1),4);
        for wv=1:length(wavelen_downwell)
            [k,~,~,~,stats] = regress(log(adj_edl(ipar,wv)),[ones(size(depth(impr))) depth(impr)]);
            k_wv(wv,:)=[wavelen_downwell(wv) k(1) k(2) stats(1)];   
        end
        
        
        if plotflag
              subplot(1,2,2,'replace')
              plot(log(adj_edl(:,wv)),depth(edl_ind),'.'), hold on
              plot(log(adj_edl(ipar,wv)),depth(impr),'.')
              plot(k(1)+k(2)*depth(impr),depth(impr),'-')
            set(gca,'YDir','reverse')
            title(['Folder: ' num2str(foldernum) ', File: ' num2str(filenum) ' out of ' num2str(length(K_PAR)) '; ' filename])
%             xl=get(gca,'xlim');  yl=get(gca,'ylim');
%         text(0.80*diff(xl)+xl(1),0.95*diff(yl)+yl(1),{['K: ' num2str(K(2))];['R2: ' num2str(STATS(1))] })
         
            
        end
        
    elseif K_PAR(filenum).flag==3; %split casts! Can look at both :)                
    
    
    end
end
