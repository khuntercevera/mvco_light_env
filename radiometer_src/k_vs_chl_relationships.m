

%% RELATIONSHIP TO CHLOROPHYLLS at MVCO:

load /Users/kristenhunter-cevera/Documents/MATLAB/MVCO_Syn_analysis/mvco_chlrep.mat

%load /Volumes/Lab_data/Grazers/EnvironmentalData/grazer_environmental_data.mat

%chlavg - those titles are:
%col 1: average of all replicates for whole chl
%col 2: average of all replicates for < 10 um chl
%col 3: average of all replicates for < 80 um chl
%%
subplot(2,3,4,'replace'), hold on

chl_yearday=find_yearday(matdate);
%plot(chl_yearday,chlavg(:,1),'.','color',[0 0.7 0])
set(gca,'fontsize',14,'box','on')
xlim([1 366])
xlabel('Year day')
ylabel('Chlorophyll Avg')
%Also a nice seasonal pattern!

%and for plotting, find all chl that came from box around sites close to
%tower:
chl_tower=find(lat > approx_station_loc(1,3)-0.01 & lat < approx_station_loc(1,3)+0.01 & ...
    lon > approx_station_loc(1,2)-0.01 & lon < approx_station_loc(1,2)+0.01);
chl_node=find(lat > approx_station_loc(2,3)-0.01 & lat < approx_station_loc(2,3)+0.01 & ...
    lon > approx_station_loc(2,2)-0.01 & lon < approx_station_loc(2,2)+0.01);
plot(chl_yearday(chl_tower),chlavg(chl_tower,1),'.','color',[0 0.6 0])
plot(chl_yearday(chl_node),chlavg(chl_node,1),'.','color',[0 0.3 0])
title('Avg total Chl values at tower and node')

%% RELATIONSHIP? MATCH CHL TO CASTS

%use event number to match...
% chl_match=nan(length(k_record),2);
% for q=1:length(k_record)
%     
%     ii=find(cellfun('isempty',regexp(cellstr(event),k_record{q,5}))==0); %match by event number
%     
%     if ~isempty(ii)
%         chl_match(q,1)=k_values(q,2);
%         chl_match(q,2)=nanmean(chlavg(ii,1)); %for now, just use the average of all the measurements at depth
%     else
%         keyboard
%     end
% end

chl_match=nan(length(k_avg),2);
for q=1:length(k_avg)
    
    ii=find(cellfun('isempty',regexp(cellstr(event),eventlist{q}))==0); %match by event number
    
    if ~isempty(ii)
        chl_match(q,1)=k_avg(q,2);
        chl_match(q,2)=nanmean(chlavg(ii,1)); %for now, just use the average of all the measurements at depth
    else
        keyboard
    end
end


%%
subplot(2,3,5,'replace')
hold on
plot(chl_match(:,2),-k_avg(:,5),'o','markersize',4,'color',[0.5 0.5 0.5])
%scatter(chl_match(:,2),-k_avg(:,5),30,k_avg(:,1),'filled')
%plot(chl_match(:,3),-cast_record(:,2),'ro')
ylim([0 0.6])
%%
%how about the ones just at the tower?
%jj=find(k_avg(:,8)==3 | k_avg(:,8)==4);
%or the latter half of the year:
jj=find(k_avg(:,1) > 150);
scatter(chl_match(jj,2),-k_avg(jj,5),30,k_avg(jj,1),'filled')
% plot(chl_match(jj,2),-k_avg(jj,5),'.','markersize',14,'color',[0 0.5 1])

set(gca,'box','on','fontsize',14)
xlabel('Chlorophyll Avg')
ylabel('Attenuation coefficient, K')


%% and a plot of K_par vs Chl with Morel 1988 equation and then fitted equation of that form:

%CAREFUL - XSCALE IS IN LOG SCLAE IN MOREL 1988!!!! Dah!!!!

% The curve from Morel 1988:
%K_PAR= 0.121 * C(mg/m3)^).428
plot(0.1:0.1:10,0.121*(0.1:0.1:10).^0.428,'-','linewidth',2,'color',[0.5 0.5 0.5])
title('note log scale for x-axis')
set(gca,'xscale','log')
xlim([0 10])

% fit a power curve through the data:
%all the data:
x = chl_match(:,2);
y = -k_avg(:,5);
nn=find(~isnan(x) & ~isnan(y));

[X1,~,~,EXITFLAG] = lsqnonlin(@(theta) fit_powercurve(theta,x(nn),y(nn)),[0.2, 2.2]);
plot(sort(x),X1(1)*sort(x).^X1(2),'-','color',[0 0 0])

%just the tower and node:
%tn=find(k_avg(:,8)==3 | k_avg(:,8)==4); %tower or node values
tn=find(k_avg(:,1) > 150);
x = chl_match(tn,2);
y = -k_avg(tn,5);
nn=find(~isnan(x) & ~isnan(y));
[X2,~,~,EXITFLAG] = lsqnonlin(@(theta) fit_powercurve(theta,x(nn),y(nn)),[0.2, 2.2]);
plot(sort(x),X2(1)*sort(x).^X2(2),'-','color',[0 0.5 1])

%this also works!
% f = fit(x(jj),y(jj),'power')
% figure,
% plot(f,x,y)

text(0.2,0.45,['Tower data: ' num2str(X1(1)) 'x^{' num2str(X1(2)) '}'],'color',[0 0.5 1])
text(0.2,0.42,['All data: ' num2str(X2(1)) 'x^{' num2str(X2(2)) '}'])
text(0.2,0.39,'Morel 1988: 0.121x^{0.428}','color',[0.4 0.4 0.4])

%% what about a linear form??

tn=find(k_avg(:,8)==3 | k_avg(:,8)==4); %tower or node values
x = chl_match(tn,2);
y = -k_avg(tn,5);
nn=find(~isnan(x) & ~isnan(y));
[b,~,~,~,stats] = regress(y(nn),[ones(length(nn),1) x(nn)]);

%%
x = chl_match(:,2);
y = -k_avg(:,5);
nn=find(~isnan(x) & ~isnan(y));
[b,~,~,~,stats] = regress(sqrt(y(nn)),[ones(length(nn),1) x(nn)]);

%% Climatologies...

addpath /Users/kristenhunter-cevera/Documents/MATLAB/mvco_tools/

% mean weekly chl:
[time_chl, daily_chl, chl_years] = timeseries2ydmat(matdate, chlavg(:,1));
[chl_wk_avg, chl_wk_std, yd_wk, chl_mn_avg, chl_mn_std, yd_mn] = dy2wkmn_climatology(daily_chl, chl_years);

[~,is]=sort(k_avg(:,2));
[time_k, daily_k, k_years] = timeseries2ydmat(k_avg(is,2), -k_avg(is,5));
[k_wk_avg, k_wk_std, k_yd_wk, k_mn_avg, k_mn_std, k_yd_mn] = dy2wkmn_climatology(daily_k, k_years);

itn=union(chl_tower,chl_node); %indexes into chl for tower and node:
[time_chltn, daily_chltn, chltn_years] = timeseries2ydmat(matdate(itn), chlavg(itn,1));
[chltn_wk_avg, chltn_wk_std, yd_wk, chltn_mn_avg, chltn_mn_std, yd_mn] = dy2wkmn_climatology(daily_chltn, chltn_years);

[~,is]=sort(k_avg(tn,2));
[time_ktn, daily_ktn, ktn_years] = timeseries2ydmat(k_avg(tn(is),2), -k_avg(tn(is),5));
[ktn_wk_avg, ktn_wk_std, ktn_yd_wk, ktn_mn_avg, ktn_mn_std, ktn_yd_mn] = dy2wkmn_climatology(daily_ktn, ktn_years);

%%
subplot(2,3,6,'replace')
nn=find(~isnan(k_wk_avg));
[ax, h1, h2]=plotyy(yd_wk,chl_wk_avg,k_yd_wk(nn),k_wk_avg(nn));
hold(ax(2))
hold(ax(1))
set(ax(1),'xlim',[1 366],'ycolor',[0 0.7 0],'fontsize',14)
set(ax(2),'xlim',[1 366],'ycolor',[0 0.5 1],'fontsize',14,'visible','on')
ylabel(ax(1),'Chlorophyll Avg')
ylabel(ax(2),'K')
set(h1,'color',[0 0.7 0],'marker','.','linewidth',2,'markersize',10)
set(h2,'color',[0 0.5 1],'marker','o','linestyle','none','markersize',8)
xlabel('Year Day')

h3=plot(ax(1),yd_wk,chltn_wk_avg,'.-','linewidth',2,'color',[0 0.3 0]);
h4=plot(ax(2),ktn_yd_wk,ktn_wk_avg,'o','markersize',8,'color',[0 0.3 0]);

legend([h1(1); h2(1); h3(1); h4(1)],'all chl','all k','tower/node chl','tower/node k')
title('Weekly averages')

%% hmmm...may need to look at monthly averages...


[chl_month, time_chl_mn] = ydmat2monthlymat(daily_chl, chl_years);
[k_month, time_k_mn] = ydmat2monthlymat(daily_k, k_years);



%% save that figure!

addpath /Users/kristenhunter-cevera/Documents/MATLAB/matlab_tools/export_fig_2016/
set(gcf,'color','w')
export_fig /Users/kristenhunter-cevera/MVCO_light_at_depth/PAR_k_values.pdf

