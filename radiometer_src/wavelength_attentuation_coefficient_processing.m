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

foldernum=good_data(6);
load(fullfile(processed_path,datafolders{foldernum},['/mat_outfiles/data_' datafolders{foldernum} '.mat']))
eval(['tempdata=data_' datafolders{foldernum} ';'])
load(fullfile(processed_path,datafolders{foldernum},['/mat_outfiles/K_PAR_' datafolders{foldernum}]));
eval(['K_PAR=K_PAR_' datafolders{foldernum} ';'])
matsource=fullfile(processed_path,datafolders{foldernum},'/mat_outfiles/');

%%
for filenum=1:length(K_PAR); 
    %
    k_lambda(filenum).file=K_PAR(filenum).file;
    k_lambda(filenum).NOTES=K_PAR(filenum).NOTES;
    k_lambda(filenum).flag=K_PAR(filenum).flag;
    
    if K_PAR(filenum).flag==0;
        
        filename=tempdata(filenum).file;
        depth=tempdata(filenum).depth;
        adj_edl=tempdata(filenum).adj_edl; %dark corrected downwelling light amounts for each wavelength with depth
        wavelen_downwell=tempdata(filenum).wavelen_downwell; %corresponding wavelengths
        edl_PAR=tempdata(filenum).edl_PAR; %calculated PAR
        edl_ind=tempdata(filenum).edl_ind; %index that matches to MPR sensor (depth)
        
        K=K_PAR(filenum).K;
        impr=K_PAR(filenum).depth_index; %index that matches to downcast and used to do PAR regression
        ipar=K_PAR(filenum).par_index; %index that matches to downcast and used to do PAR regression
        
        %[~, i400]=min(abs(wavelen_downwell-400)); %indexes of closest wv to 400
        [~, i720]=min(abs(wavelen_downwell-720)); %indexes of closest wv to 700
        
        % calculate k for each wavelength and record:
        k_wv=nan(i720,4);
        for wv=1:i720 %don't go beyond this as water absorbs very quickly...
            
            if any(adj_edl(ipar,wv) < 0)
                k_wv(wv,:)=[wavelen_downwell(wv) NaN NaN NaN];
            else
                [k,~,~,~,stats] = regress(log(adj_edl(ipar,wv)),[ones(size(depth(impr))) depth(impr)]);
                
                 % k_wv(wv,:)=[wavelen_downwell(wv) k(1) k(2) stats(1)];
                if stats(1) < 0.70
                    k_wv(wv,:)=[wavelen_downwell(wv) NaN NaN NaN]; %don't keep bad fits...
                else
                    k_wv(wv,:)=[wavelen_downwell(wv) k(1) k(2) stats(1)];
                end
%                 
                if plotflag1 && floor(wv/5)==wv/5 %plot every 10
                    figure(1), clf, hold on
                    plot(log(adj_edl(:,wv)),depth(edl_ind),'.')
                    plot(log(adj_edl(ipar,wv)),depth(impr),'.')
                    plot(k(1)+k(2)*depth(impr),depth(impr),'-')
                    set(gca,'YDir','reverse')
                    title({['Folder: ' num2str(foldernum) ', File: ' num2str(filenum) ' out of ' num2str(length(K_PAR)) '; ' filename]; ['wv=' num2str(wavelen_downwell(wv)) ' R2:' num2str(stats(1))]})
                    pause
                end
            end
            
        end
        %
        k_lambda(filenum).k_wv=k_wv;
        k_lambda(filenum).depth_index=impr;
        k_lambda(filenum).par_index=ipar;
        
        figure(6), hold on
        plot(k_wv(:,1),-k_wv(:,3),'.')
        
    elseif K_PAR(filenum).flag==3; %split casts! Can look at both :)
        
        filename=tempdata(filenum).file;
        depth=tempdata(filenum).depth;
        adj_edl=tempdata(filenum).adj_edl; %dark corrected downwelling light amounts for each wavelength with depth
        wavelen_downwell=tempdata(filenum).wavelen_downwell; %corresponding wavelengths
        edl_PAR=tempdata(filenum).edl_PAR; %calculated PAR
        edl_ind=tempdata(filenum).edl_ind; %index that matches to MPR sensor (depth)
        
        impr1=K_PAR(filenum).depth_index1; %index that matches to downcast and used to do PAR regression
        ipar1=K_PAR(filenum).par_index1; %index that matches to downcast and used to do PAR regression
        impr2=K_PAR(filenum).depth_index2; %index that matches to downcast and used to do PAR regression
        ipar2=K_PAR(filenum).par_index2; %index that matches to downcast and used to do PAR regression
        
        [~, i720]=min(abs(wavelen_downwell-720)); %indexes of closest wv to 720
        k_wv1=nan(i720,4);
        for wv=1:i720
            if any(adj_edl(ipar1,wv) < 0)
                k_wv1(wv,:)=[wavelen_downwell(wv) NaN NaN NaN];
            else
                [k1,~,~,~,stats1] = regress(log(adj_edl(ipar1,wv)),[ones(size(depth(impr1))) depth(impr1)]);
                if stats1(1) < 0.95
                    k_wv1(wv,:)=[wavelen_downwell(wv) NaN NaN NaN]; %don't keep bad fits...
                else
                    k_wv1(wv,:)=[wavelen_downwell(wv) k1(1) k1(2) stats1(1)];
                end
            end
        end
        
        k_wv2=nan(i720,4);
        for wv=1:i720
            if any(adj_edl(ipar2,wv) < 0)
                k_wv2(wv,:)=[wavelen_downwell(wv) NaN NaN NaN];
            else
                [k2,~,~,~,stats2] = regress(log(adj_edl(ipar2,wv)),[ones(size(depth(impr2))) depth(impr2)]);
                if stats2(1) < 0.95
                    k_wv2(wv,:)=[wavelen_downwell(wv) NaN NaN NaN]; %don't keep bad fits...
                else
                    k_wv2(wv,:)=[wavelen_downwell(wv) k2(1) k2(2) stats2(1)];
                end
            end
        end
        
        k_lambda(filenum).k_wv1=k_wv1;
        k_lambda(filenum).depth_index1=impr1;
        k_lambda(filenum).par_index1=ipar1;
        k_lambda(filenum).k_wv2=k_wv2;
        k_lambda(filenum).depth_index2=impr2;
        k_lambda(filenum).par_index1=ipar2; 
        
        figure(6), hold on
        plot(k_wv1(:,1),-k_wv1(:,3),'.')
        plot(k_wv2(:,1),-k_wv2(:,3),'.')
    end
    
    pause 
    
end

%%
%save the attenuation fits:
eval(['k_lambda_' datafolders{foldernum} '=k_lambda;'])
eval(['save ' matsource 'k_lambda_' datafolders{foldernum} '.mat k_lambda_' datafolders{foldernum}])
figure(6), clf
clear k_lambda tempdata K_PAR mprtime edl_ind edl_PAR depth impr* ipar*



