%script to explore the density structure around the tower over time:

addpath /Users/Kristen/Documents/MATLAB/tools/plotxx/
%construct density plots:
% a quick look the tower profiles with both the upcasts and downcasts
figure, hold on
for j=1:length(tower_ind_DC)
    
    %downcasts:
    col_hdr_dc=hdr_downcasts{tower_ind_DC(j)};
    temp_data_dc=data_downcasts{tower_ind_DC(j)};
    dens_dc=temp_data_dc{6};
    depth_dc=temp_data_dc{1};
    temperature_dc=temp_data_dc{3};
    sal_dc=temp_data_dc{5};
    
    time_dc=file_time_dc{tower_ind_DC(j),3};
    
    %plots:
    %just density:
%     plot(dens_dc,depth_dc,'*k')
%     set(gca,'ydir','reverse')
%     ylabel(col_hdr_dc{6})
%        title([file_time_dc{tower_ind_downcast(j),2} ' : ' datestr(file_time_dc{tower_ind_downcast(j),3})])
    figure(1)
    clf
    subplot(121)
    [ax, h1,h2]=plotxx(dens_dc,depth_dc,temperature_dc,depth_dc);
    set(ax(1),'ydir','reverse')
    set(ax(2),'ydir','reverse')
    set(h1,'marker','*','color','k')
    set(h2,'marker','*','color','r')
    ylabel(col_hdr{6})
    title(['Downcast: ' file_time{tower_ind(j),2} ' : ' datestr(file_time{tower_ind(j),3})])

    
    %Upcasts:
    col_hdr_uc=hdr_upcast{tower_ind_UC(j)};
    temp_data_uc=data_upcast{tower_ind_UC(j)};
    dens_uc=temp_data_uc{6};
    depth_uc=temp_data_uc{1};
    temperature_uc=temp_data_uc{3};
    sal_uc=temp_data_uc{5};
    
%      plot(dens_uc,depth_uc,'*b')
%      title(['DC: ' datestr(file_time_dc{tower_ind_DC(j),3}) ' : UC: ' datestr(file_time_uc{tower_ind_UC(j),3})])
% 
     
    
    figure(2)
    clf
    [ax2, h1b,h2b]=plotxx(dens_uc,depth_uc,temperature_uc,depth_uc);
    set(ax2(1),'ydir','reverse')
    set(ax2(2),'ydir','reverse')
    set(h1b,'marker','*','color','k')
    set(h2b,'marker','*','color','r')
    ylabel(col_hdr{6})
    title(['Upcast: ' file_time{tower_ind(j),2} ' : ' datestr(file_time{tower_ind(j),3})])

    %cool scatter plots:
       % scatter(time*ones(length(depth),1),depth,70,dens,'filled'), caxis([1021 1026])    
       % scatter(time*ones(length(temperature),1),depth,70,temperature,'filled'), caxis([0 21])

   % xlim([1021 1026])
    pause
    
end


%% A quick gather of the salinity values for each cast:

addpath /Users/Kristen/Documents/MATLAB/plotxx/
%construct density plots:
% a quick look the tower profiles with both the upcasts and downcasts
figure, hold on
time_rec=[];
for j=1:length(tower_ind_DC)
    
    %downcasts:
    col_hdr_dc=hdr_downcasts{tower_ind_DC(j)};
    temp_data_dc=data_downcasts{tower_ind_DC(j)};
    dens_dc=temp_data_dc{6};
    depth_dc=temp_data_dc{1};
    temperature_dc=temp_data_dc{3};
    sal_dc=temp_data_dc{5};
    
    time_dc=file_time_dc{tower_ind_DC(j),3};
    time_rec=[time_rec; time_dc];
    %plots:
    %just density:
%     plot(dens_dc,depth_dc,'*k')
%     set(gca,'ydir','reverse')
%     ylabel(col_hdr_dc{6})
%        title([file_time_dc{tower_ind_downcast(j),2} ' : ' datestr(file_time_dc{tower_ind_downcast(j),3})])
    figure(1)
    clf
    subplot(121)
    [ax, h1,h2]=plotxx(dens_dc,depth_dc,sal_dc,depth_dc);
    set(ax(1),'ydir','reverse')
    set(ax(2),'ydir','reverse','xlim',[29 34])
    set(h1,'marker','*','color','k')
    set(h2,'marker','*','color','r')
    ylabel(col_hdr{5})
    title(['Downcast: ' file_time{tower_ind(j),2} ' : ' datestr(file_time{tower_ind(j),3})])
    
    pause

end