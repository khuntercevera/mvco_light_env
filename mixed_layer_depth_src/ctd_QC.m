%script to try to see if the water column is stratified or not:
addpath /Users/kristenhunter-cevera/Dropbox/MVCO_mixed_layer_depth/seawater_ver3_2/
%addpath /Users/Kristen/Documents/MATLAB/12.808/seawater_ver3_2/

%a look at the potential density and Brunt Vaisala Frequency:
for j=1:length(tower_ind);

    col_hdr=hdr{tower_ind(j)};
    temp_data=data{tower_ind(j)};
    dens=temp_data{6};
    
    depth=temp_data{1};
    press=temp_data{2};
    temperature=temp_data{3};
    sal=temp_data{5};
  
    time=file_time{tower_ind(j),3};
    %find time
    ii=find(cellfun('isempty',regexp(col_hdr,'Time, Elapsed'))==0);
    cast_time=temp_data{ii};
    
    %calculate and potential density:
    pdens=sw_pden(sal,temperature,press,0);
    temp_data{end+1}=pdens;
    
    %just density:
    clf
    subplot(1,3,1,'replace')
    plot(dens,depth,'*k')
    hold on
    plot(pdens,depth,'r*')
    set(gca,'ydir','reverse','fontsize',14)
    xlabel('Density (kg/m^3)') %ylabel(col_hdr{6})
    title([file_time{tower_ind(j),2} ' - ' datestr(file_time{tower_ind(j),3})])
    legend('Obs \rho','Potential \rho','location','NorthWest')
    
    %calculate the Brunt-Vaisala Frequnecy:
    N2=sw_bfrq(sal,temperature,press);
    temp_data{end+1}=N2;
    subplot(1,3,2,'replace')
    plot(N2,depth(2:end),'*k-')
    set(gca,'ydir','reverse','fontsize',14)
    xlabel('Brunt-Vaisala Freq') %ylabel(col_hdr{6})
    title('stratified at N2 ~ 0.5-1 10e-4  s^{-2}')
    line([1e-4 1e-4], get(gca,'ylim'),'color','r')
    
    subplot(1,3,3,'replace')
    plot(cast_time,depth,'*k-')
    set(gca,'ydir','reverse','fontsize',14)
    xlabel('Time') %ylabel(col_hdr{6})
    
    pause
    col_hdr=[col_hdr; {'potential density'};{'Brunt-Vaisala Frequency'}];
    hdr{tower_ind(j)}=col_hdr;
    data{tower_ind(j)}=temp_data;
   
end


%% Edit the casts that have 0 salinity (residual fresh water that hasn't cleared the conductivity sensor?)

%cd /Users/Kristen/Documents/MATLAB/12.808/seawater_ver3_2

% tower_ind=tower_ind_DC;
% data=data_downcasts;
% hdr=hdr_downcasts;
% file_time=file_time_dc;

data2=cell(length(tower_ind),1);
data2_time=zeros(length(tower_ind),1);
for j=1:length(tower_ind);

    col_hdr=hdr{tower_ind(j)};
    temp_data=data{tower_ind(j)};
    temp_data2=cell(1,6);
    
    %extract the data from the cell array:
    depth=temp_data{1};
    press=temp_data{2};
    temperature=temp_data{3};
    sal=temp_data{5};
    
    if ~isempty(file_time{tower_ind(j),2});
        data2_time(j)=file_time{tower_ind(j),2};
    else
        data2_time(j)=file_time{tower_ind(j),3};
    end
    
    %remove first 2 meters of data on the down cast and then any points
    %that have anomalous salinity values:
    ii=find(depth > 3 & sal > 30); 
    
    depth=depth(ii);
    press=press(ii);
    temperature=temperature(ii);
    sal=sal(ii);
    
    %calculate and potential density and Brunt-Vaisala Frequnecy:
    pdens=sw_pden(sal,temperature,press,0);
    N2=sw_bfrq(sal,temperature,press);
    
    %put into new cell array:
    temp_data2{1}=depth;
    temp_data2{2}=press;
    temp_data2{3}=temperature;
    temp_data2{4}=sal;
    temp_data2{5}=pdens;
    temp_data2{6}=[nan; N2];
    
    data2{j}=temp_data2;
    
%     %plots:
%     clf
%     subplot(121)
%     plot(pdens,depth,'r*')
%     set(gca,'ydir','reverse')
%     ylabel('Depth')
%     title([file_time{tower_ind(j),2} ' : ' datestr(file_time{tower_ind(j),3})])
%     
%     subplot(122)
%     plot(N2,depth(2:end),'*k-')
%     set(gca,'ydir','reverse')
%     ylabel('Depth')
%     title('Brunt-Vaisala Freq')
%     
%     pause
   
end

%% Now find the mixed layer:

mixed_layer_rec=zeros(length(data2),10);
ml_titles={'MLD','ref_pdens','dens_at_MLD','temp at MLD','depth_12','dens_at_12','temp_at_12','depth_at_4','dens_at_4','temp_at_4'}

for j=1:length(data2)
    %
    temp_data=data2{j};
    
    N2=temp_data{6};
    depth=temp_data{1};
    temperature=temp_data{3};
    pdens=temp_data{5};
    
    ref_pdens=min(pdens(depth > 4)); %just to avoid noisy readings at surface
    test=pdens-ref_pdens;
    jj=find(test > 0.08); %must have an abs density difference
    
    if isempty(jj)
        mld=depth(end);
        k=length(depth);
    else
        ii=find(N2(jj) > 0.8e-3); %Is this too high?
        if ~isempty(ii)
            mld=depth(jj(ii((1))));
            k=jj(ii((1)));
        else
        mld=depth(jj(1));
        k=jj(1);
        end
    end
    
    
%     ii=find(depth > 5);
%     [maxN2 ia]=max(N2(ii));
    
   % mld(j)=depth(ii(ia));
    
    %     %plots:
%     clf
%     subplot(141)
%     plot(pdens,depth,'k*')
%     set(gca,'ydir','reverse')
%     ylabel('Depth')
%     hold on
%     plot(pdens(k),mld,'ro')
%     title(datestr(data2_time(j)))
%     
%     subplot(142)
%     plot(N2,depth,'*b-')
%     hold on
%      plot(N2(k),mld,'ro')
%     set(gca,'ydir','reverse')
%     ylabel('Depth')
%     title('Brunt-Vaisala Freq')
%     
%      subplot(143)
%   
%      plot(N2./maxN2,depth,'bo')
%     set(gca,'ydir','reverse')
%     ylabel('Depth')
%     title('Relative Brunt-Vaisala Freq')
%     
%      subplot(144)
%   
%      plot(temperature,depth,'ro')
%     set(gca,'ydir','reverse')
%     ylabel('Depth')
%     title('Temperature')
    %
    ind12=find(depth < 12.1);
    ind12=ind12(end); %so this could be 12m or just the bottom of the cast...
    ind4=find(depth >= 4);
    ind4=ind4(1); 
    mixed_layer_rec(j,:)=[mld ref_pdens pdens(k) temperature(k) depth(ind12) pdens(ind12) temperature(ind12) depth(ind4) pdens(ind4) temperature(ind4)];
    %
  %keyboard
    

end

%% Plots to look at the relationship between mixed layer depth and differences in temperature/density

%for some reason (should still resolve this as of July 5th, 2014!) the
%casts after 1/1/2011 seem a little strange, so for now, use only those
%from 2007-2010:

k=find(data2_time < datenum(['1-1-2011']));
sub_mld=mixed_layer_rec(k,:);
sub_time=data2_time(k);

figure %mld over time
plot(sub_time(:,1),sub_mld(:,1),'b*')
hold on
plot(sub_time(:,1),sub_mld(:,1),'g')

figure %mld against difference in temperature, density between 4m and 12m and relationship btwn diff in temp and diff in dens
subplot(131)
plot(sub_mld(:,1),sub_mld(:,7)-sub_mld(:,10),'bo')
xlabel('mixed layer depth')
ylabel('Temperature difference btwn 4 and 12m')

subplot(132)
plot(sub_mld(:,1),sub_mld(:,6)-sub_mld(:,9),'bo')
xlabel('mixed layer depth')
ylabel('Potential Density difference btwn 4 and 12m')

subplot(133)
plot(sub_mld(:,7)-sub_mld(:,10),sub_mld(:,6)-sub_mld(:,9),'bo')
xlabel('Temperature difference btwn 4 and 12m')
ylabel('Potential Density difference btwn 4 and 12m')