%% to be incorporated in anotherscript...

%what we want to do in this script is for each "set" of data,
%see what the casts look like over the day
%calculate attentuation coefficient for each .raw file
%see what k looks like over the casts

%other fancy plots at the end :)

% first find all the folders with the processed matlab files...

%sourcepath=fullfile('\\sosiknas1\lab_data\mvco\HyperPro_Radiometer\processed_radiometer_files\'); % path to folders with raw data...
sourcepath=fullfile('/Volumes/Lab_data/MVCO/HyperPro_Radiometer/processed_radiometer_files/');

%How many of these do we have?
d = dir(sourcepath);
isub = [d(:).isdir]; %# returns logical vector
foldernames = {d(isub).name}'; %names of folders
temp=regexp(foldernames,'\d{1,2}\w{3,4}20\d{2}'); %find all folders with this naming scheme
datafolders = foldernames(cellfun('isempty',temp)==0);
%there should be (as of summer 2017) only 20 folders...

%Load in a log about the data:
load(fullfile('\\sosiknas1\lab_data\mvco\HyperPro_Radiometer\','initial_data_notes.mat'))
 
%load in the list of good data folders:
load(fullfile('\\sosiknas1\lab_data\mvco\HyperPro_Radiometer\','good_data_folders.mat')) 
%%
%can just manually enter in a foldernum for now .... put into a for loop later....
foldernum=good_data(1);

%find the data!
matsource=fullfile(sourcepath,datafolders{foldernum},'/mat_outfiles/');

%where to save:
savepath=fullfile(matsource,'mat_outfiles/'); %this directory should already exist!

eval(['load ' matsource 'data_' datafolders{foldernum} '.mat'])
eval(['tempdata=data_' datafolders{foldernum} ';'])
eval(['sourcefiles = {data_' datafolders{foldernum} '.file};']) %the files in the structure
numfiles=length(sourcefiles);

% first an overview plot:
figure(1), clf, hold on

eval(['load ' matsource 'location_' datafolders{foldernum} '.mat'])
eval(['location=location_' datafolders{foldernum} ';'])
for filenum=1:numfiles
    plot(tempdata(filenum).mprtime, tempdata(filenum).depth,'.-')
end
set(gca,'Ydir','reverse')
datetick


%% Okay, if that looks reasonable, then move onto attenuation calculations!

%disp(['#files: ' num2str(numfiles)])
for filenum=1:numfiles
    
    filename=sourcefiles{filenum};
    ii=find(cellfun('isempty',strfind(data_notes(:,2),filename(1:end-4)))==0); %find the file
    
    %Is it a good cast?
    if isempty(data_notes{ii,4}) || ~cellfun('isempty',(regexpi(data_notes(ii,4),'not a cast')))
        
        disp(['File ' num2str(filenum) ' is an empty file or has no valid casts - skipping: ' filename])
        K_PAR(filenum).file=filename;
        K_PAR(filenum).NOTES='empty file / not a cast';
    
    else
        disp(['Processing: file number: ' num2str(filenum) ' out of ' num2str(numfiles) ' | ' filename ' for attenuation coefficient'])
        
        %a reasonable plan may be to look at the sign of the change in
        %depth - a positive, continuous value implies a falling profiler -
        %look for these to find the downcasts.
        %Remove top and bottom 2/3 as well for calculation
        
        %for easier handling:
        depth=tempdata(filenum).depth;
        edl_ind=tempdata(filenum).edl_ind;
        edl_PAR=tempdata(filenum).edl_PAR;
        
        dz=smooth2(diff(depth),15); %change in depth, smoothed to reduce noise
        dc=find(dz > 0.025); %dc for downcast
        
        jj=find(depth(dc) > 2 & depth(dc) < max(depth)-1.5); %exclude top and bottom measurements
        %jj=find(depth(dc) > 2 & depth(dc) < max(depth)-0.5); %special cases - closer to bottom, for short casts      
        %jj=find(depth(dc) > 2 & depth(dc) < max(depth)-20); %special case for deep casts, change in K with deeper depths
        
        %matching indexes for PAR:
        [~,impr,ipar]=intersect(dc(jj),edl_ind);
        impr=dc(jj(impr)); %ipar already indexes into edl_PAR...
        
        %now, calculate k:
        [K,~,~,~,STATS] = regress(log(edl_PAR(ipar)),[ones(size(depth(impr))) depth(impr)]);
        
        if STATS(1) < 0.9
            disp('Hmmm...maybe not such a great regression fit to find k?')
        end
        
        %plot those results!
        figure(filenum), set(gcf,'position',[33         468        1218         510])
        subplot(1,2,1,'replace'), hold on
        plot(tempdata(filenum).mprtime,depth,'.-')
        plot(tempdata(filenum).mprtime(impr),depth(impr),'.-')
        datetick
        
        plot(tempdata(filenum).mprtime(edl_ind),0.1*edl_PAR-20,'.-')
        plot(tempdata(filenum).mprtime(impr),0.1*edl_PAR(ipar)-20,'.-')
        %plot(tempdata(filenum).mprtime(impr),0.1*tempdata(filenum).esl_PAR(ipar)+20,'.-')
        %set(gca,'Ydir','reverse')
        
        subplot(1,2,2,'replace')
        plot(log(edl_PAR),depth(edl_ind),'.'), hold on
        plot(log(edl_PAR(ipar)),depth(impr),'.')
        plot(K(1)+K(2)*depth(impr),depth(impr),'-')
        set(gca,'YDir','reverse')
        title(['File: ' num2str(filenum) ' ; ' filename])
        text(2,4,['K: ' num2str(K(2))])
        
        keyboard
        
        %if need to exclude a cast:
        % casts=find(diff(impr)>100);
        % cc=7;
        % [K,~,~,~,STATS] = regress(log(edl_PAR(ipar(casts(cc):end))),[ones(size(depth(impr(casts(cc):end)))) depth(impr(casts(cc):end))]);
        
        notes='';
        K_PAR(filenum).file=filename;
        K_PAR(filenum).stats=STATS;
        K_PAR(filenum).K=K(2);
        K_PAR(filenum).NOTES=notes;
    end
    
end

clear tempdata mprtime edl_ind edl_PAR depth

%save the attenuation fits:
eval(['K_PAR_' datafolders{foldernum} '=K_PAR;'])
eval(['save ' matsource 'K_PAR_' datafolders{foldernum} '.mat K_PAR_' datafolders{foldernum}])

close all
clear K_PAR


%% Okay, and if you want to check the fits for a certain set of data:

foldernum=good_data(1);
load(fullfile(sourcepath,datafolders{foldernum},['/mat_outfiles/data_' datafolders{foldernum}]))
eval(['tempdata=data_' datafolders{foldernum} ';'])
load(fullfile(sourcepath,datafolders{foldernum},['/mat_outfiles/K_PAR_' datafolders{foldernum}]));
eval(['K_PAR=K_PAR_' datafolders{foldernum} ';'])
%%
filenum=1; 

depth=tempdata(filenum).depth;
edl_ind=tempdata(filenum).edl_ind;
edl_PAR=tempdata(filenum).edl_PAR;
K=K_PAR(filenum).K;

figure(filenum), set(gcf,'position',[33         468        1218         510])
        subplot(1,2,1,'replace'), hold on
        [ax h1 h2]=plotyy(tempdata(filenum).mprtime,depth,tempdata(filenum).mprtime(edl_ind),edl_PAR)
        %plot(tempdata(filenum).mprtime(impr),depth(impr),'.-')
        %plot(tempdata(filenum).mprtime(impr),0.1*edl_PAR(ipar)-20,'.-')
        datetick(ax(2))
        datetick(ax(1))
        set(ax(1),'YDir','reverse')     
        ylabel(ax(1),'Depth')
        ylabel(ax(2),'PAR')

        subplot(1,2,2,'replace')
        plot(log(edl_PAR),depth(edl_ind),'.'), hold on
        %plot(log(edl_PAR(ipar)),depth(impr),'.')
        %plot(K(1)+K(2)*depth(impr),depth(impr),'-')
         plot(K(1)+K(2)*depth,depth,'-')
        set(gca,'YDir','reverse')
        title(['File: ' num2str(filenum) ' ; ' filename])
        text(2,4,['K: ' num2str(K(2))])
