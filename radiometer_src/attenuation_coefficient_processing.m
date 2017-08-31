%% Calculation of attentuation coefficient K:

%what we want to do in this script is for each "set" of data,
%see what the casts look like over the day
%calculate attentuation coefficient for each .raw file
%see what k looks like over the casts

%other fancy plots at the end :)

% first find all the folders with the processed matlab files...

%sourcepath=fullfile('\\sosiknas1\lab_data\mvco\HyperPro_Radiometer\processed_radiometer_files\'); % path to folders with raw data...
sourcepath=fullfile('/Volumes/Lab_data/MVCO/HyperPro_Radiometer/');
processed_path=fullfile(sourcepath,'/processed_radiometer_files/');
%How many of these do we have?
d = dir(processed_path);
isub = [d(:).isdir]; %# returns logical vector
foldernames = {d(isub).name}'; %names of folders
temp=regexp(foldernames,'\d{1,2}\w{3,4}20\d{2}'); %find all folders with this naming scheme
datafolders = foldernames(cellfun('isempty',temp)==0);
%there should be (as of summer 2017) only 20 folders...

%Load in a log about the data:
load(fullfile(sourcepath,'initial_data_notes.mat'))

%load in the list of good data folders:
load(fullfile(sourcepath,'good_data_folders.mat'))


%%
%can just manually enter in a foldernum for now .... put into a for loop later....
foldernum=good_data(15);

%find the data!
matsource=fullfile(processed_path,datafolders{foldernum},'/mat_outfiles/');

%where to save:
savepath=fullfile(matsource,'mat_outfiles/'); %this directory should already exist!

eval(['load ' matsource 'data_' datafolders{foldernum} '.mat'])
eval(['tempdata=data_' datafolders{foldernum} ';'])
eval(['sourcefiles = {data_' datafolders{foldernum} '.file};']) %the files in the structure
numfiles=length(sourcefiles);
eval(['load ' matsource 'location_' datafolders{foldernum} '.mat'])
eval(['location=location_' datafolders{foldernum} ';'])

% first an overview plot:
figure(21), clf
subplot(1,2,1,'replace'), hold on
for filenum=1:numfiles
    plot(tempdata(filenum).mprtime, tempdata(filenum).depth,'.-')
end
set(gca,'Ydir','reverse','xgrid','on','ygrid','on')
datetick

subplot(1,2,2,'replace')
hold on
%plot station points:
plot(-70.567,41.325,'o','markersize',16,'color',[0.5 0.5 0.5]) %tower
plot(-70.555,41.335,'o','markersize',16,'color',[0.5 0.5 0.5]) %node
plot(-70.6275,41.325,'o','markersize',16,'color',[0.5 0.5 0.5]) %station
plot(-70.505,41.325,'o','markersize',16,'color',[0.5 0.5 0.5]) %station
plot(-70.45,41.3275,'o','markersize',16,'color',[0.5 0.5 0.5]) %station
plot(-70.567,41.255,'o','markersize',16,'color',[0.5 0.5 0.5]) %station
plot(-70.567,41.2,'o','markersize',16,'color',[0.5 0.5 0.5]) %station
plot(-70.567,41.145,'o','markersize',16,'color',[0.5 0.5 0.5]) %station

%and the lat lon:
pos=[cell2mat({location(:).lon}') cell2mat({location(:).lat}')];
plot(pos(:,1),pos(:,2),'rp')

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
        K_PAR(filenum).flag=1;
        
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
        
        %         impr=impr(32:end);
        %         ipar=ipar(32:end);
        
        %now, calculate k:
        [K,~,~,~,STATS] = regress(log(edl_PAR(ipar)),[ones(size(depth(impr))) depth(impr)]);
        
        %for a split k:
        %          d1=find(depth(impr) < 16);
        %          d2=find(depth(impr) > 17);
        %         [K1,~,~,~,STATS1] = regress(log(edl_PAR(ipar(d1))),[ones(size(depth(impr(d1)))) depth(impr(d1))]);
        %         [K2,~,~,~,STATS2] = regress(log(edl_PAR(ipar(d2))),[ones(size(depth(impr(d2)))) depth(impr(d2))]);
        
        if STATS(1) < 0.9
            disp('Hmmm...maybe not such a great regression fit to find k?')
        end
        
        %plot those results!
        figure(filenum), set(gcf,'position',[33         468        1218         510])
        
        subplot(1,2,1,'replace')
        [ax h1 h2]=plotyy(tempdata(filenum).mprtime,depth,tempdata(filenum).mprtime(edl_ind),edl_PAR);
        hold(ax(1)); hold(ax(2));
        plot(ax(1),tempdata(filenum).mprtime(impr),depth(impr),'.','markersize',12)
        plot(ax(2),tempdata(filenum).mprtime(impr),edl_PAR(ipar),'.','markersize',12)
        set(ax(1),'YDir','reverse','xgrid','on','ygrid','on')
        ylabel(ax(1),'Depth')
        ylabel(ax(2),'PAR')
        datetick(ax(2))
        datetick(ax(1))
        
        subplot(1,2,2,'replace')
        plot(log(edl_PAR),depth(edl_ind),'.'), hold on
        plot(log(edl_PAR(ipar)),depth(impr),'.')
        plot(K(1)+K(2)*depth(impr),depth(impr),'-')
        set(gca,'YDir','reverse')
        title(['File: ' num2str(filenum) ' ; ' filename])
        xl=get(gca,'xlim');  yl=get(gca,'ylim');
        text(0.80*diff(xl)+xl(1),0.95*diff(yl)+yl(1),{['K: ' num2str(K(2))];['R2: ' num2str(STATS(1))] })
        
        if max(depth) < 7
            disp('too short of a cast')
            K_PAR(filenum).file=filename;
            K_PAR(filenum).NOTES='depth too short for reliable cast, excluding...';
            K_PAR(filenum).flag=2;
            
        else
            
            K_PAR(filenum).file=filename;
            K_PAR(filenum).stats=STATS;
            K_PAR(filenum).K=K;
            K_PAR(filenum).depth_index=impr;
            K_PAR(filenum).par_index=ipar;
            K_PAR(filenum).NOTES=notes;
            K_PAR(filenum).flag=0;
            
            %split cast:
            %             K_PAR(filenum).file=filename;
            %             K_PAR(filenum).stats=STATS1;
            %             K_PAR(filenum).K1=K1;
            %             K_PAR(filenum).depth_index1=impr(d1);
            %             K_PAR(filenum).par_index1=ipar(d1);
            %             K_PAR(filenum).stats=STATS2;
            %             K_PAR(filenum).K2=K2;
            %             K_PAR(filenum).depth_index2=impr(d2);
            %             K_PAR(filenum).par_index2=ipar(d2);
            %             K_PAR(filenum).NOTES='split cast from 0-15 and 16-max';
            %             K_PAR(filenum).flag=3;
            
        end
        
        keyboard
        
        %if need to exclude a cast:
        % casts=find(diff(impr)>100);
        % cc=7;
        % [K,~,~,~,STATS] = regress(log(edl_PAR(ipar(casts(cc):end))),[ones(size(depth(impr(casts(cc):end)))) depth(impr(casts(cc):end))]);
        
    end
    
end

clear tempdata mprtime edl_ind edl_PAR depth

%save the attenuation fits:
eval(['K_PAR_' datafolders{foldernum} '=K_PAR;'])
eval(['save ' matsource 'K_PAR_' datafolders{foldernum} '.mat K_PAR_' datafolders{foldernum}])

%disp K's:
k_rec=[];
for filenum=1:numfiles
    if K_PAR(filenum).flag ==0
        k_rec=[k_rec; K_PAR(filenum).K(2)];
    end
end

disp(k_rec)
disp({K_PAR(:).flag}')

close all
clear K_PAR



%% IF WANT TO CHECK K AND REGRESSION FITS FOR THE DATA:

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

%% put into a for loop or just go one by one....

clearvars K_PAR tempdata

foldernum=good_data(15);
load(fullfile(processed_path,datafolders{foldernum},['/mat_outfiles/data_' datafolders{foldernum} '.mat']))
eval(['tempdata=data_' datafolders{foldernum} ';'])
load(fullfile(processed_path,datafolders{foldernum},['/mat_outfiles/K_PAR_' datafolders{foldernum}]));
eval(['K_PAR=K_PAR_' datafolders{foldernum} ';'])

set(gcf,'position',[33         468        1218         510])

%%
for filenum=1:length(K_PAR);
    
    if K_PAR(filenum).flag==0;
        
        filename=tempdata(filenum).file;
        depth=tempdata(filenum).depth;
        edl_ind=tempdata(filenum).edl_ind;
        edl_PAR=tempdata(filenum).edl_PAR;
        K=K_PAR(filenum).K;
        impr=K_PAR(filenum).depth_index;
        ipar=K_PAR(filenum).par_index;
        STATS=K_PAR(filenum).stats;
        
        subplot(1,2,1,'replac e')
        [ax h1 h2]=plotyy(tempdata(filenum).mprtime,depth,tempdata(filenum).mprtime(edl_ind),edl_PAR);
        hold(ax(1)); hold(ax(2));
        plot(ax(1),tempdata(filenum).mprtime(impr),depth(impr),'.','markersize',12)
        plot(ax(2),tempdata(filenum).mprtime(impr),edl_PAR(ipar),'.','markersize',12)
        set(ax(1),'YDir','reverse','xgrid','on','ygrid','on')
        title(datafolders{foldernum})
        ylabel(ax(1),'Depth')
        ylabel(ax(2),'PAR')
        datetick(ax(2))
        datetick(ax(1))
        
        subplot(1,2,2,'replace')
        plot(log(edl_PAR),depth(edl_ind),'.'), hold on
        plot(log(edl_PAR(ipar)),depth(impr),'.')
        plot(K(1)+K(2)*depth(impr),depth(impr),'-')
        set(gca,'YDir','reverse')
        title(['Folder: ' num2str(foldernum) ', File: ' num2str(filenum) ' out of ' num2str(length(K_PAR)) '; ' filename])
        xl=get(gca,'xlim');  yl=get(gca,'ylim');
        text(0.80*diff(xl)+xl(1),0.95*diff(yl)+yl(1),{['K: ' num2str(K(2))];['R2: ' num2str(STATS(1))] })
        
        pause
        clf
        
    elseif K_PAR(filenum).flag==3; %split casts! Can look at both :)
        
        filename=tempdata(filenum).file;
        depth=tempdata(filenum).depth;
        edl_ind=tempdata(filenum).edl_ind;
        edl_PAR=tempdata(filenum).edl_PAR;

        STATS1=K_PAR(filenum).stats;
        K1=K_PAR(filenum).K1;
        impr1=K_PAR(filenum).depth_index1;
        ipar1=K_PAR(filenum).par_index1;
        STATS2=K_PAR(filenum).stats;
        K2=K_PAR(filenum).K2;
        impr2=K_PAR(filenum).depth_index2;
        ipar2=K_PAR(filenum).par_index2;

        subplot(1,2,1,'replace')
        [ax h1 h2]=plotyy(tempdata(filenum).mprtime,depth,tempdata(filenum).mprtime(edl_ind),edl_PAR);
        hold(ax(1)); hold(ax(2));
        plot(ax(1),tempdata(filenum).mprtime(impr1),depth(impr1),'.','markersize',12,'color',[0.2081    0.1663    0.5292])
        plot(ax(2),tempdata(filenum).mprtime(impr1),edl_PAR(ipar1),'.','markersize',12,'color',[0.2081    0.1663    0.5292])
        plot(ax(1),tempdata(filenum).mprtime(impr2),depth(impr2),'.','markersize',12,'color',[0.9763    0.9831    0.0538])
        plot(ax(2),tempdata(filenum).mprtime(impr2),edl_PAR(ipar2),'.','markersize',12,'color',[0.9763    0.9831    0.0538])
        set(ax(1),'YDir','reverse','xgrid','on','ygrid','on')
        title(datafolders{foldernum})
        ylabel(ax(1),'Depth')
        ylabel(ax(2),'PAR')
        datetick(ax(2))
        datetick(ax(1))
        
        subplot(1,2,2,'replace')
        plot(log(edl_PAR),depth(edl_ind),'.'), hold on
        plot(log(edl_PAR(ipar1)),depth(impr1),'.','color',[0.2081    0.1663    0.5292])
        plot(K1(1)+K1(2)*depth(impr1),depth(impr1),'-')
        plot(log(edl_PAR(ipar2)),depth(impr2),'.','color',[0.9763    0.9831    0.0538])
        plot(K2(1)+K2(2)*depth(impr2),depth(impr2),'-')
        set(gca,'YDir','reverse')
        title(['Folder: ' num2str(foldernum) ', File: ' num2str(filenum) ' out of ' num2str(length(K_PAR)) '; ' filename ' SPLIT CAST!'])
        xl=get(gca,'xlim');  yl=get(gca,'ylim');
        text(0.80*diff(xl)+xl(1),0.85*diff(yl)+yl(1),{['K1: ' num2str(K1(2))];['R2: ' num2str(STATS1(1))] })
        text(0.80*diff(xl)+xl(1),0.95*diff(yl)+yl(1),{['K2: ' num2str(K2(2))];['R2: ' num2str(STATS2(1))] })
        
        pause
        clf
    end
    
end
