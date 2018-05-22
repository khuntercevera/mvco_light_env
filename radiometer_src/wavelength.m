% Script to look at how the spectral quality of light changes with time:

%below is a bit convoluted...lots of different thigns recorded, etc...not
%the easiest to deal with or revisit...

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
wv_record=[]; %information about files
wv_log={};
wv_spectra={}; %holds spectra
plotflag=0;
%extraplotflag=1;
wv_spectra_titles={'dw spectra at 4m'; 'corresponding solal spectra for 4m';
   'dw spectra at 8m'; 'corresponding solal spectra for 8m';
   'dw spectra at 12 m'; 'corresponding solal spectra for 12m';};


for foldernum=good_data' %are these the ones with lat/lon?
    
    %load in all the data!
    matsource=fullfile(processed_path,datafolders{foldernum},'/mat_outfiles/');
    
    eval(['load ' matsource 'data_' datafolders{foldernum} '.mat'])
    eval(['tempdata=data_' datafolders{foldernum} ';'])
    
    eval(['load ' matsource 'K_PAR_' datafolders{foldernum} '.mat'])
    eval(['K_PAR=K_PAR_' datafolders{foldernum} ';'])
    
    eval(['load ' matsource 'location_' datafolders{foldernum} '.mat'])
    eval(['location=location_' datafolders{foldernum} ';'])
    
   
    % because some of the casts were split, need to account for this:
    for filenum=1:length(K_PAR);
        
        
        if K_PAR(filenum).flag ~=1 && K_PAR(filenum).flag ~=2 %not a cast or two short a cast
            
            wv=tempdata(filenum).wavelen_downwell;
            edl=tempdata(filenum).adj_edl;
            depth=tempdata(filenum).depth;
            time=tempdata(filenum).mprtime;
            edl_ind=tempdata(filenum).edl_ind;
            esl_ind=tempdata(filenum).esl_ind;
            
            wv_sol=tempdata(filenum).wavelen_solarst;
            esl=tempdata(filenum).adj_esl;
            
            if K_PAR(filenum).flag==0 %a legit cast, 1 is empty file, 2 is data, but not a cast
                
                impr=K_PAR(filenum).depth_index;
                ipar=K_PAR(filenum).par_index;
                
            elseif K_PAR(filenum).flag==3; %k-par was piecewise linear, only investigate first fit
                
                impr=K_PAR(filenum).depth_index1; %downcast 
                ipar=K_PAR(filenum).par_index1;
                
            end
            
            %record solar std spectrum:
            [~, im]=max(esl,[],2);
            max_wv_sol=[0 mode(wv_sol(im))];
            
            %temp matrix of max wv at each depth
            [~, ii]=max(edl(ipar,:),[],2); %find max light level at each depth record
            max_wv=[depth(impr) wv(ii)]; %[depth    wavelength of max energy at that depth]
            
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
            
            if isempty(i4), sp4={[]}; sol4={[]}; else sp4={edl(ipar(i4),:)}; sol4={esl(i4_sol,:)};end
            if isempty(i8), sp8={[]}; sol8={[]}; else sp8={edl(ipar(i8),:)}; sol8={esl(i8_sol,:)};end
            if isempty(i12), sp12={[]}; sol12={[]};  else sp12={edl(ipar(i12),:)}; sol12={esl(i12_sol,:)}; end
            
            wv_spectra=[wv_spectra; sp4 sol4 sp8 sol8 sp12 sol12];
            
            if plotflag
                
                set(gcf,'color','w')
                subplot(1,4,1,'replace')
                %a look at the cast:
                plot(time,depth,'-'), hold on
                plot(time(impr),depth(impr),'.','markersize',12)
                set(gca,'YDir','reverse','fontsize',14)
                datetick
                title(['File: ' num2str(filenum) ' out of ' num2str(length(K_PAR)) ],'Interpreter','none')
                ylabel('Depth (m)')
                
                subplot(1,4,2,'replace'), hold on
                plot(wv_sol,esl,'-','color',[0.5 0.5 0.5])
                title(tempdata(filenum).file,'Interpreter','none')                
                colormap(flipud(parula(length(impr))))
                for j=1:length(ipar)
                    scatter(wv,edl(ipar(j),:),40,depth(impr(j))*ones(size(wv)),'filled')
                end
                caxis([0 max(depth(ipar))])
                hbar=colorbar; ylabel(hbar,'depth')
                set(gca,'fontsize',14,'box','on')
                ylabel('Irradiance (uW/cm^2/nm)')
                xlabel('Wavelength (nm)')
                xlim([345 805])
                ylim([-2 max(max(esl))+2])
                
                
                subplot(1,4,3,'replace'),hold on
                if ~isempty(i4), h1=plot(wv,edl(ipar(i4),:),'-','color',[0 0.5 1]); end
                if ~isempty(i8),h2=plot(wv,edl(ipar(i8),:),'-','color',[0 0 1]); end
                if ~isempty(i12),h3=plot(wv,edl(ipar(i12),:),'-','color',[0 0 0.7]); end
                
                if ~isempty(i4) && ~isempty(i8) && ~isempty(i12)
                legend([h1(1); h2(1); h3(1)],'4 m','8 m', '12 m')
                elseif ~isempty(i4) && ~isempty(i8) && isempty(i12)
                    legend([h1(1); h2(1)],'4 m','8 m')
                else ~isempty(i4) && isempty(i8) && isempty(i12)
                    legend(h1(1),'4 m')
                end
                
                set(gca,'fontsize',14,'box','on')
                ylabel('Irradiance (uW/cm^2/nm)')
                xlabel('Wavelength (nm)')
                xlim([345 805])
                ylim([-2 max(max(edl(ipar(i4),:)))+2])
                
                subplot(1,4,4,'replace'), hold on
                plot(max_wv(:,1),max_wv(:,2),'.','markersize',12,'color',[0.4 0.4 0.4])
                %lines:
                line([4 4],ylim,'color',[0.5 0.5 0.5])
                line([8 8],ylim,'color',[0.5 0.5 0.5])
                line([12 12],ylim,'color',[0.5 0.5 0.5])
                if ~isempty(i4), plot(wv_record(end,6),wv_record(end,7),'p','markerfacecolor',[0 0.5 1],'markeredgecolor',[0 0.5 1],'markersize',10); end
                if ~isempty(i8),plot(wv_record(end,8),wv_record(end,9),'p','markerfacecolor',[0 0 1],'markeredgecolor',[0 0 1],'markersize',10); end
                if ~isempty(i12),plot(wv_record(end,10),wv_record(end,11),'p','markerfacecolor',[0 0 0.6],'markeredgecolor',[0 0 0.6],'markersize',10); end
                set(gca,'fontsize',14,'box','on')
                xlabel('Depth (m)')
                ylabel('Wavelength of maximum irradiance (nm)')
                
                %keyboard
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
subplot(1,3,3,'replace'), hold on
plot(wv_record(:,1),wv_record(:,8),'.','color',[0 0.8 1],'markersize',12) %4m max wv
plot(wv_record(:,1),wv_record(:,10),'.','color',[0 0.3 1],'markersize',12) %8m max wv
plot(wv_record(:,1),wv_record(:,12),'.','color',[0 0 0.6],'markersize',12) %12m max wv

legend('4m','8m','12m')
%hmmm - cool, looks like mainly steady at 540-560, but some interesting
%features at 500? to explore more!
set(gca,'box','on')
xlabel('Year Day')
ylabel('Max wavelength (nm)')
%% fancier plot for time of year:

figure, ylim([0 1.1]), hold on
cc=jet(366);

[~,is]=sort(wv_record(:,1));
for q=1:size(wv_record,1)
    plot(wv,wv_spectra{is(q),1}./repmat(max(wv_spectra{is(q),1},[],2),1,size(wv_spectra{is(q),1},2)),'.-','color',cc(wv_record(is(q),1),:))
    pause
end

%What I'd also like to do is record a set of casts, and then plot these
%overlayed, color coded by time of year to identify any shifts??

%% Or could make a 3D plot!

%figure, 
subplot(1,3,1,'replace')
hold on, grid on
cc=jet(366);

[~,is]=sort(wv_record(:,1));
for q=1:size(wv_record,1)
    plot3(wv,wv_record(is(q),1)*ones(size(wv)),wv_spectra{is(q),1}./repmat(max(wv_spectra{is(q),1},[],2),1,size(wv_spectra{is(q),1},2)),'.-','color',cc(wv_record(is(q),1),:))

    %plot(wv,wv_spectra{is(q),1}./repmat(max(wv_spectra{is(q),1},[],2),1,size(wv_spectra{is(q),1},2)),'.-','color',cc(wv_record(is(q),1),:))
    %pause(0.2)
end

zlim([0 1.1])
ylim([1 366])
set(gca,'YDir','reverse')

xlabel('Wavelength (nm)')
ylabel('Year Day')
zlabel('Relative distrbution')
%view([32  34]) %for subplot2
view([-24 42.8]) %for subplot1

%% Or even better a heatmap!

 addpath /Users/kristenhunter-cevera/MVCO_V6V8/MVCO_Syn_V6V8_analysis/sequence_analysis/compositional_stats_analysis/
 %what is the best way to normalize spectra?
 
%let's only use tower casts for now...
%and spectra at 4 m ....

% find the tower casts:
tn=find(wv_record(:,13)==4);

unqdaylist=unique(wv_record(tn,2));

rel_sp4=nan(137,14); 
rel_sp8=nan(137,14);  
rel_sp12=nan(137,14); 
y=[];
 
for q=1:size(unqdaylist,1)
  
    qq=find(wv_record(tn,2)==unqdaylist(q));
    disp([num2str(length(qq)) ' tower casts: ' datestr(unqdaylist(q))])
    %spectra at 4m:
%     if wv_record(tn(is(q)),2) ~= tempday
%         rel_sp=[rel_sp NaN(length(wv),3)];
%         tempday=wv_record(tn(is(q)),2); %if not, replace with current date being evaluated
%     end
    
    all_spectra4=cell2mat(wv_spectra(tn(qq),1));    
    all_spectra4(all_spectra4<0)=0; %hmmm...some measurement noise and offsets?
    temp4=closure(all_spectra4);
    norm_spectra4=center(temp4);
    %take compositional mean of these? or have just max at 1?

    all_spectra8=cell2mat(wv_spectra(tn(qq),3));    
    all_spectra8(all_spectra8<0)=0; %hmmm...some measurement noise and offsets?
    temp8=closure(all_spectra8);
    norm_spectra8=center(temp8);
    
    all_spectra12=cell2mat(wv_spectra(tn(qq),5));      
    all_spectra12(all_spectra12<0)=0; %hmmm...some measurement noise and offsets?
    temp12=closure(all_spectra12);
    norm_spectra12=center(temp12);

    subplot(1,3,1,'replace')
    plot(wv, temp4,'.-','color',[0.5 0.5 0.5]), hold on
    plot(wv,norm_spectra4,'b.-')
    title('Spectra at 4m')
    
    subplot(1,3,2,'replace')
    title('Spectra at 8m')
    if ~isempty(all_spectra8)
    plot(wv, temp8,'.-','color',[0.5 0.5 0.5]), hold on
    plot(wv,norm_spectra8,'b.-')
    end
    
    subplot(1,3,3,'replace')
    title('Spectra at 12m')
    if ~isempty(all_spectra12)
    plot(wv, temp12,'.-','color',[0.5 0.5 0.5]), hold on
    plot(wv,norm_spectra12,'b.-')
    end
    
    pause
    %if this looks good, record for that day:

    rel_sp4(:,q)= norm_spectra4';
    if ~isempty(all_spectra8) rel_sp8(:,q)=norm_spectra8'; end
    if ~isempty(all_spectra12), rel_sp12(:,q)=norm_spectra12'; end
    
    y=[y unqdaylist(q)];
    
end

%%
save wv-lite wv* rel*

%% line figures?

cc=jet(366);
figure, hold on
for j=1:14
plot(wv, rel_sp4(:,j),'-','color',cc(find_yearday(unqdaylist(j)),:),'linewidth',3)
end

xlim([350 750])
ylim([0 0.025])

colormap jet
colorbar
%%

[~,is]=sort(find_yearday(y));
subplot(1,3,1,'replace')
imagesc(find_yearday(y(is)),wv,rel_sp4(:,is))

subplot(1,3,2,'replace')
imagesc(find_yearday(y(is)),wv,rel_sp8(:,is))

subplot(1,3,3,'replace')
imagesc(find_yearday(y(is)),wv,rel_sp12(:,is))

%% towards publication figure....

figure, 
imagesc(find_yearday(y(is)),wv,rel_sp12(:,is))
xlabel('Yearday')
ylabel('Wavelength (nm)')
hbar=colorbar;
line(xlim,[540 540],'color',[0.5 0.5 0.5])
line(xlim,[495 495],'color',[0.5 0.5 0.5])
line(xlim,[560 560],'color',[0.5 0.5 0.5])
set(hbar,'location','northoutside')

%Ooh! Put the syn emission spectra right next to them!
%Ooh! and then and then put the lambda specific wavelengths next to them!


%% An example plot of solar and wavelengths for 4m, 8m, 12m depths
q=30;
figure,hold on
%spectra at depth:
plot(wv,wv_spectra{q,1},'.-','color',[0 0.6 1])
plot(wv,wv_spectra{q,3},'.-','color',[0 0.3 1])
plot(wv,wv_spectra{q,5},'.-','color',[0 0 0.6])
%solar spectra:
plot(wv,wv_spectra{q,2},'.-','color',[0.6 0.6 0.6])
plot(wv,wv_spectra{q,4},'.-','color',[0.4 0.4 0.4])
plot(wv,wv_spectra{q,6},'.-','color',[0.2 0.2 0.2])

