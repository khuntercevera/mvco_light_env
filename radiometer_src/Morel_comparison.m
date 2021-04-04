%% exploration and comparison to Morel equations:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%K(LAMBDAS)...comparison with Morel...

%So, each k at each wavelength was plotted against chl value and then a
%power function was fit through it (after subtracting kw)...

% import the values for comparison:
filename = '/Users/kristenhunter-cevera/MVCO_light_at_depth/radiometer_src/Morel_Maritorena_2001_k_fits.txt';
delimiter = '\t';
formatSpec = '%f%f%f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);
fclose(fileID);

Morel2001=cell2mat(dataArray(1:4));
titles_Morel2001={'wavelength';'Kw';'e param';'chi param'};

%initial plot for attenuation of water:
figure
plot(Morel2001(:,1), Morel2001(:,2),'.-')

%interpolate these to match lambdas used on our radiometer:
k_water=interp1(Morel2001(:,1), Morel2001(:,2),lambdas);

%% Ok, so let's see what those the Kbio looks like after has substracted Kw plotted against chlorophyll:
%for first pass, just compare directly wiht Morel2001, but should try to calculate kw
%directly from aw+1/2(bw)...

% [~, im]=min(abs(Morel2001(:,1)-lambdas(w))); %if don't want to interpolate and just use closest wavelength....
%k_bio=k_lambdas(:,w)-Morel2001(im,2);
test=[0:0.2:10];
k_bio=k_lambdas-repmat(k_water,92,1);

%do a linear regression of log-log transformed data to match with
%Morel:
rec=nan(length(lambdas),6);
for w=1:length(lambdas)
    
    x=[ones(size(chl_match(:,2))) log(chl_match(:,2))];
    y=log(k_bio(:,w)); %wavelength by wavelength
    [b,~,~,~,stats]=regress(y,x);
    
    [~, im]=min(abs(Morel2001(:,1)-lambdas(w))); %find closest to Morel
    
    if ~any(imag(b(:))) && ~any(isinf(b(:))) && ~any(isnan(b(:)))
        
        rec(w,:)=[lambdas(w) b(1) b(2) stats(1) stats(3) im]; %Morel2001(im,1:4)
        
        subplot(1,2,1,'replace'), hold on %linear view
        plot(log(chl_match(:,2)),log(k_bio(:,w)),'.')
        line([-2 10], b(1) + b(2)*[-2 10])
        title(['\lambda ' num2str(lambdas(w))])
        xlabel('Chlorophyll mg/m^{3}')
        ylabel('k_{bio} [k(\lambda) - k_w(\lambda)]')
        set(gca,'box','on','fontsize',14)
        
        %power function view:
        subplot(1,2,2,'replace'), hold on %linear view
        plot(chl_match(:,2),k_bio(:,w),'.')
        plot(test,exp(b(1)).*test.^b(2),'.-')
        plot(test,Morel2001(im,4).*(test.^Morel2001(im,3)),'.-')
        set(gca,'xscale','log','box','on','fontsize',14)
        set(gca,'yscale','log')
        xlim([0.01 100])
        ylim([0.001 10])
        title(['Morel match: ' num2str(Morel2001(im,1))])
        xlabel('Chlorophyll mg/m^{3}')
        ylabel('k_{bio} [k(\lambda) - k_w(\lambda)]')
        %because non negative y's can't appear on log scale
        
        set(gcf,'color','w')
        
        pause
        
    end
end

rec_titles={'lambda' 'log chi' 'e' 'R2' 'p-value' 'index to closest lambda in Morel'};

%% a figure to show the consistent differences:

ll=[420 480 550 650]; %wavelengths to see
for i=1:4
    
    [~, w]=min(abs(lambdas-ll(i)));
    [~, im]=min(abs(Morel2001(:,1)-lambdas(w)));
    
    subplot(2,2,i,'replace'), hold on 
    plot(chl_match(:,2),k_bio(:,w),'.')
    plot(test,exp(rec(w,2)).*test.^rec(w,3),'-')
    plot(test,Morel2001(im,4).*(test.^Morel2001(im,3)),'-')
    set(gca,'xscale','log','box','on','fontsize',14)
    set(gca,'yscale','log')
    xlim([0.01 100])
    ylim([0.001 10])
    title(['\lambda: ' num2str(lambdas(w)) ' nm, Morel match: ' num2str(Morel2001(im,1)) ' nm'])
    xlabel('Chlorophyll mg/m^{3}')
    ylabel('k_{bio} [k(\lambda) - k_w(\lambda)]')
    %because non negative y's can't appear on log scale
end
set(gcf,'color','w')


%% hmmm...these plots look a bit different from Morel...

%figure 4 in Morel & Maritorena:
figure
subplot(1,2,1,'replace'), hold on
plot(rec(:,1),rec(:,3),'.-') %MVCO e
hold on
plot(Morel2001(:,1),Morel2001(:,3),'.-') %Morel e
title('e parameter')

subplot(1,2,2,'replace'), hold on
plot(rec(:,1),exp(rec(:,2)),'.-') %MVCO chi
plot(Morel2001(:,1),Morel2001(:,4),'.-') %Morel chi
title('Chi parameter')

%% let's construct their k_total vs wavelength plot and compare too:

%Figure 5 in Morel & Maritorena:

figure, hold on
k_MM = Morel2001(:,4).*(0.03).^Morel2001(:,3) + Morel2001(:,2);
plot(Morel2001(:,1),k_MM,'b.-')
k_MM = Morel2001(:,4).*(0.3).^Morel2001(:,3) + Morel2001(:,2);
plot(Morel2001(:,1),k_MM,'r.-')
k_MM = Morel2001(:,4).*(1).^Morel2001(:,3) + Morel2001(:,2);
plot(Morel2001(:,1),k_MM,'g.-')
k_MM = Morel2001(:,4).*(3).^Morel2001(:,3) + Morel2001(:,2);
plot(Morel2001(:,1),k_MM,'m.-')

k_MVCO = exp(rec(:,2)).*(0.03).^rec(:,3) + k_water';
plot(rec(:,1),k_MVCO,'b.:')
k_MVCO = exp(rec(:,2)).*(0.3).^rec(:,3) + k_water';
plot(rec(:,1),k_MVCO,'r.:')
k_MVCO = exp(rec(:,2)).*(1).^rec(:,3) + k_water';
plot(rec(:,1),k_MVCO,'g.:')
k_MVCO = exp(rec(:,2)).*(3).^rec(:,3) + k_water';
plot(rec(:,1),k_MVCO,'m.:')

xlim([400 700])


%% color coded k-lambda plot?

clf, hold on
cc=jet(5);
tn=find(k_values(:,8)==3 | k_values(:,8)==4); 

count1=0; count2=0; count3=0; count4=0;

seasons=[1 90;
91 180;
181 274;
275 366];

%marker_types={'.','o','^','p','s'};
dy=unique(k_values(tn,2));

for q=1:length(dy)
    
    qq=find(k_values(tn,2)==dy(q));
    
    temp=unique(k_values(tn(qq),1));
    if  temp >= seasons(1,1) && temp <= seasons(1,2)
        plotnum=1; count1=count1+1; count=count1;
    elseif temp >= seasons(2,1) && temp <= seasons(2,2)
        plotnum=2;  count2=count2+1; count=count2;
    elseif temp >= seasons(3,1) && temp <= seasons(3,2)
        plotnum=3; count3=count3+1; count=count3;
    elseif temp >= seasons(4,1) && temp <= seasons(4,2)
        plotnum=4; count4=count4+1; count=count4;
    end
    
    subplot(1,4,plotnum)
    plot(lambdas, k_lambdas(tn(qq),:),'.-','color',cc(count,:))
%     plot(lambdas, k_lambdas(tn(qq),:),'-','color',cc(temp,:),'marker',marker_types{count})
    
    hold on
    plot(lambdas,k_water,'k:')
end


