function [smu, new_mu_est, no_new_est]=replace_suspicious_mus(time_mu,daily_mu,plotflag)
%highlight suspicious mu's based on what we know about temperature response
%and see if can find a replacement in that day's model runs:

load ~/Documents/MATLAB/MVCO_Syn_analysis/mvco_envdata_26Jul2016.mat  %(or most current processed MVCO Enviornmental data)

smu=find(daily_mu > 0.20 & Tbeam_corr <= 4);
smu=[smu; find(daily_mu > 0.25 & (Tbeam_corr <= 5 & Tbeam_corr > 4))];
smu=[smu; find(daily_mu > 0.30 & (Tbeam_corr <= 6 & Tbeam_corr > 5))];
smu=[smu; find(daily_mu > 0.35 & (Tbeam_corr <= 7 & Tbeam_corr > 6))];
%smu=[smu; find(daily_mu > 0.45 & (Tbeam_corr <= 8 & Tbeam_corr > 7))];
smu=unique(smu);

% if want to see:
% figure, plot(Tbeam_corr,daily_mu,'b.')
% hold on
% plot(Tbeam_corr(smu),daily_mu(smu),'r.')
%
no_new_est=[];
new_mu_est=[];
%Now, load in model data and as - is there an alternative growth rate that
%is close in likelihood and lower?

%% load in those model runs:
for q=1:length(smu)
    
    day=time_mu(smu(q));
    [yearlabel,~,~]=datevec(day);
    disp(['Day: ' datestr(day) ' : temp:' num2str(Tbeam_corr(smu(q))) ' : mu: ' num2str(daily_mu(smu(q)))])
    
    switch yearlabel
        case 2003
            filelabel='May';
        case 2004
            filelabel='Apr';
        case 2005
            filelabel='Apr';
        case 2006
            filelabel='May';
        case 2007
            filelabel='Mar';
        otherwise
            filelabel='Jan';
    end
    
    eval(['load /Volumes/Lab_data/MVCO/FCB/MVCO_' filelabel num2str(yearlabel) '/model/output_July2016/mvco_14par_dmn_' num2str(yearlabel) '.mat'])
    
    qq=find(modelresults(:,1)==day);
    temp=allmodelruns{qq,1};
    
    %were there other modelfits that resulted in a lower growth rate?
    [ss is]=sort(temp(:,15));
    
    if Tbeam_corr(smu(q)) <= 4
        tt=find(temp(:,16) < .20);
    elseif  Tbeam_corr(smu(q)) > 4 & Tbeam_corr(smu(q)) <= 5
        tt=find(temp(:,16) < .25);
    elseif  Tbeam_corr(smu(q)) > 5 & Tbeam_corr(smu(q)) <= 6
        tt=find(temp(:,16) < .30);
    elseif  Tbeam_corr(smu(q)) > 6 & Tbeam_corr(smu(q)) <= 7
        tt=find(temp(:,16) < .35);
    end
    
    if ~isempty(tt)
        [s1 is1]=sort(temp(tt,15)); %sort log likelihoods
        fprintf('Difference in -logL from next reasonable value: %6.2f, %6.2f, and diff: %4.3f\n',ss(1),s1(1),ss(1)-s1(1))
        fprintf('with that division rate being: %1.3f\n',temp(tt(is1(1)),16))
        new_mu_est=[new_mu_est; smu(q) temp(tt(is1(1)),16)];
    else
        [s1 is1]=sort(temp(:,16)); %sort division rates -> choose the lowest one?
        fprintf('Difference in -logL from lowest mu found: %6.2f, %6.2f, and diff: %4.3f\n',ss(1),temp(is1(1),15),ss(1)-temp(is1(1),15))
        fprintf('with that division rate being: %1.3f\n',temp(is1(1),16))
        no_new_est=[no_new_est; smu(q) size(temp,1) temp(is1(1),16)];
    end
    
    %plot ofr sanity check
    if plotflag
        h1=subplot(1,3,1,'replace'); plot(1:length(temp),temp(is,15),'.-')
        h2=subplot(1,3,2,'replace'); plot(1:length(temp),temp(is,16),'.-')
        h3=subplot(1,3,3,'replace'); plot(Tbeam_corr,daily_mu,'b.'), hold on, plot(Tbeam_corr(smu(q)),daily_mu(smu(q)),'rp')
        if ~isempty(tt)
            ww=find(is==tt(is1(1))); %find lowest -logL that matches the new temp criteria
            axes(h1), hold on, plot(ww,temp(tt(is1(1)),15),'ro')
            axes(h2), hold on, plot(ww,temp(tt(is1(1)),16),'ro')
        else
            ww=find(is==is1(1));
            axes(h1), hold on, plot(ww,temp(is1(1),15),'ko')
            axes(h2), hold on, plot(ww,temp(is1(1),16),'ko')   
        end
        title(['Day: ' datestr(day) ' : temp:' num2str(Tbeam_corr(smu(q))) ' : mu: ' num2str(daily_mu(smu(q)))])
        pause
    end
    
    
end


