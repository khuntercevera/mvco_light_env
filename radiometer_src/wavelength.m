% Script to look at how the spectral quality of light changes with time:

%% A time series plot of K_d values!

%sourcepath=fullfile('\\sosiknas1\lab_data\mvco\HyperPro_Radiometer\processed_radiometer_files\'); % path to folders with raw data...
sourcepath=fullfile('/Volumes/Lab_data/MVCO/HyperPro_Radiometer/');
processed_path=fullfile(sourcepath,'/processed_radiometer_files/');

%How many of these do we have?
d = dir(processed_path);
isub = [d(:).isdir]; %# returns logical vector
foldernames = {d(isub).name}'; %names of folders
temp=regexp(foldernames,'\d{1,2}\w{3,4}20\d{2}'); %find all folders with this naming scheme
datafolders = foldernames(cellfun('isempty',temp)==0);

load(fullfile(sourcepath,'good_data_folders.mat'))

%%
wv_record=[];
wv_log={};
wv_spectra={};
plotflag=0;
extraplotflag=0;

for foldernum=good_data' %are these the ones with lat/lon?
    
    %load in all the data!
    matsource=fullfile(processed_path,datafolders{foldernum},'/mat_outfiles/');
    
    eval(['load ' matsource 'data_' datafolders{foldernum} '.mat'])
    eval(['tempdata=data_' datafolders{foldernum} ';'])
    
    eval(['load ' matsource 'K_PAR_' datafolders{foldernum} '.mat'])
    eval(['K_PAR=K_PAR_' datafolders{foldernum} ';'])
    
    eval(['load ' matsource 'location_' datafolders{foldernum} '.mat'])
    eval(['location=location_' datafolders{foldernum} ';'])
    
    eval(['load ' matsource 'k_lambda_' datafolders{foldernum} '.mat'])
    eval(['k_lambda=k_lambda_' datafolders{foldernum} ';'])
    
    % because some of the casts were split, need to account for this:
    for filenum=1:length(K_PAR);
        
        
        if K_PAR(filenum).flag ~=1 && K_PAR(filenum).flag ~=2 %not a cast or two short a cast
            
            wv=tempdata(filenum).wavelen_downwell;
            edl=tempdata(filenum).adj_edl;
            depth=tempdata(filenum).depth;
            edl_ind=tempdata(filenum).edl_ind;
            esl_ind=tempdata(filenum).esl_ind;
            
            wv_sol=tempdata(filenum).wavelen_solarst;
            esl=tempdata(filenum).adj_esl;
            
            if K_PAR(filenum).flag==0 %a legit cast, 1 is empty file, 2 is data, but not a cast
                
                impr=K_PAR(filenum).depth_index;
                ipar=K_PAR(filenum).par_index;
                
            elseif K_PAR(filenum).flag==3;
                
                impr=K_PAR(filenum).depth_index1;
                ipar=K_PAR(filenum).par_index1;
                
            end
            
            %record solar std spectrum:
            [~, im]=max(esl,[],2);
            max_wv_sol=[0 mode(wv_sol(im))];
            
            %temp matrix of max wv at each depth
            [~, ii]=max(edl(ipar,:),[],2); %find max light level at each depth record
            max_wv=[depth(impr) wv(ii)]; %temp record
            
            %indexes that correspond to different depth windows:
            i4=find(max_wv(:,1)>3.5 & max_wv(:,1)<4.5);
            i8=find(max_wv(:,1)>7.5 & max_wv(:,1)<8.5);
            i12=find(max_wv(:,1)>11.5 & max_wv(:,1)<12.5);
            
            %find corresponding solar std spectra...a bit of a mission!
            %Note, this may not be an exact match up in time...will need to think on this more...
            if ~isempty(i4)
                i4_sol=nan(size(i4));
                for j=1:length(i4)
                    [~, jj]=min(abs(esl_ind-impr(i4(j)))); %find closest index
                    i4_sol(j)=jj;
                end
            end
            
            %if need to check:
            %time=tempdata(filenum).mprtime;
            %[time(esl_ind(i4_sol))-time((impr(i4)))]
            
             if ~isempty(i8)
                i8_sol=nan(size(i8));
                for j=1:length(i8)
                    [~, jj]=min(abs(esl_ind-impr(i8(j)))); %find closest index
                    i8_sol(j)=jj;
                end
             end
            
              if ~isempty(i12)
                i12_sol=nan(size(i12));
                for j=1:length(i12)
                    [~, jj]=min(abs(esl_ind-impr(i12(j)))); %find closest index
                    i12_sol(j)=jj;
                end
              end
            
            if isempty(i4), r4=[NaN NaN]; else r4=[nanmean(max_wv(i4,1)) mode(max_wv(i4,2))]; end
            if isempty(i8), r8=[NaN NaN]; else r8=[nanmean(max_wv(i8,1)) mode(max_wv(i8,2))]; end
            if isempty(i12), r12=[NaN NaN]; else r12=[nanmean(max_wv(i12,1)) mode(max_wv(i12,2))]; end
            
            wv_record=[wv_record; datenum(datestr(datafolders{foldernum})) location(filenum).lat location(filenum).lon max_wv_sol r4 r8 r12];
            wv_log=[wv_log; {datafolders{foldernum}} {foldernum} {K_PAR(filenum).file} {location(filenum).file} location(filenum).eventnum {K_PAR(filenum).NOTES} {location(filenum).notes}];
            
            if isempty(i4), sp4={''}; sol4={''}; else sp4={edl(ipar(i4),:)}; sol4={esl(i4_sol,:)};end
            if isempty(i8), sp8={''}; sol8={''}; else sp8={edl(ipar(i8),:)}; sol8={esl(i8_sol,:)};end
            if isempty(i12), sp12={''}; sol12={''};  else sp12={edl(ipar(i12),:)}; sol12={esl(i12_sol,:)}; end
            
            wv_spectra=[wv_spectra; sp4 sol4 sp8 sol8 sp12 sol12];
            
            if plotflag
                
                subplot(1,4,1,'replace')
                %a look at the cast:
                plot(time,depth,'-'), hold on
                plot(time(impr),depth(impr),'.','markersize',12)
                set(gca,'YDir','reverse')
                datetick
                
                subplot(1,4,2,'replace'), hold on
                plot(wv_sol,esl,'-','color',[0.5 0.5 0.5])
                title(tempdata(filenum).file,'Interpreter','none')
                
                for j=1:length(ipar)
                    scatter(wv,edl(ipar(j),:),40,depth(impr(j))*ones(size(wv)),'filled')
                end
                
                subplot(1,4,3,'replace'),hold on
                if ~isempty(i4), plot(wv,edl(ipar(i4),:),'-','color',[0 0.5 1]); end
                if ~isempty(i8),plot(wv,edl(ipar(i8),:),'-','color',[0 0 1]); end
                if ~isempty(i12),plot(wv,edl(ipar(i12),:),'-','color',[0 0 0.7]); end
                
                subplot(1,4,4,'replace'), hold on
                plot(max_wv(:,1),max_wv(:,2),'.','markersize',12)
                if ~isempty(i4), plot(wv_record(end,6),wv_record(end,7),'p','color',[0 0.5 1]); end
                if ~isempty(i8),plot(wv_record(end,8),wv_record(end,9),'p','color',[0 0 1]); end
                if ~isempty(i12),plot(wv_record(end,10),wv_record(end,11),'p','color',[0 0 0.7]); end
                
                xlabel('Depth (m)')
                ylabel('Wavelength of maximum energy (nm)')
                
                keyboard
            end
        end
    end
end

wv_record=[find_yearday(wv_record(:,1)) wv_record];

%% recode for easier station finding:
approx_station_loc=[4 -70.567 41.325; %tower
    3 -70.5564 41.3366; %node
    1 -70.45 41.3275; %station 1
    2 -70.505 41.325; %station 2
    5 -70.6275 41.325;  %station 5
    6 -70.567 41.255; %station 6
    7 -70.567 41.2; %station 7
    8 -70.567 41.145]; %station 8

shore_pos=[41.3499,-70.5267];

% set up the indexes:

for j=1:length(approx_station_loc)
    temp_ind=find(wv_record(:,3) > approx_station_loc(j,3)-0.01 & wv_record(:,3) < approx_station_loc(j,3)+0.01 ... %lat
        & wv_record(:,4) > approx_station_loc(j,2)-0.01 & wv_record(:,4) < approx_station_loc(j,2)+0.01); %lon
    wv_record(temp_ind,13)=approx_station_loc(j,1);
end

% find the tower and node casts:
tn=find(wv_record(:,13)==3 | wv_record(:,13)==4);

%%
figure, hold on
plot(wv_record(:,1),wv_record(:,8),'o','color',[0 0.7 0]) %4m max wv
plot(wv_record(:,1),wv_record(:,10),'p','color',[0 0.5 1]) %8m max wv
plot(wv_record(:,1),wv_record(:,12),'s','color',[0 0 0.8]) %12m max wv

%hmmm - cool, looks like mainly steady at 540-560, but some interesting
%features at 500? to explore more!

%% fancier plot for time of year:

figure, ylim([0 1.1]), hold on
cc=jet(366);

[~,is]=sort(wv_record(:,1));
for q=1:size(wv_record,1)
    plot(wv,wv_spectra{is(q),1}./repmat(max(wv_spectra{is(q),1},[],2),1,size(wv_spectra{is(q),1},2)),'.-','color',cc(wv_record(is(q),1),:))
    pause(0.2)
end


%What I'd also like to do is record a set of casts, and then plot these
%overlayed, color coded by time of year to identify any shifts??





