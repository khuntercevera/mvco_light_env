% processPROII_KRHC
% SRL
% May 2016
% from processSAS
% Additionally editted by KRHC, specific to MVCO data needs, saves intermediate txt files
% Script imports and processes data from 3 sensors on HyperPro


%addpath /Users/kristenhunter-cevera/Documents/MATLAB/mvco_tools/

clear all;
close all;
close all; %for any open files
plot_flag=1;

raw_data_log={};

mastersourcepath=fullfile('\\sosiknas1\lab_data\mvco\HyperPro_Radiometer\raw_data\'); % path to folders with raw data...

%mastersourcepath=fullfile(pwd);%'/Volumes/Lab_data/SummerStudents/2013_Marco/MVCO_hyperpro'; %top level directory
%mastersourcepath=fullfile('/Users/kristenhunter-cevera/Desktop/forKristin');
masteroutpath=fullfile('\\sosiknas1\lab_data\mvco\HyperPro_Radiometer\','processed_radiometer_files'); %where you'd like final folders to be stored
%calfiledir = fullfile(pwd, 'calibration_files'); %where the calibration files are located
calfiledir = fullfile('\\sosiknas1\lab_data\mvco\HyperPro_Radiometer\', 'calibration_files\CD_Nov2008\Calibration_Files\');

% Raw files for MVCO casts are stored in different folders labeled by date
%To find the correct files in each folder:

d = dir(mastersourcepath);
isub = [d(:).isdir]; %# returns logical vector
foldernames = {d(isub).name}'; %names of folders
temp=regexp(foldernames,'\d{1,2}\w{3,4}20\d{2}'); %find all folders with this naming scheme
%temp=regexp(foldernames,'testdata'); %test case
datafolders = foldernames(cellfun('isempty',temp)==0);

%HAVE NOT INCORPORATED THIS YET, BUT MAY BE NICE:
% a choice of processing options:
% save_txtfiles=1; %if want to save the intermediate txt files that are generated from satcon
% do_satcon=1; %if do not need to run satcon and want to only work wth preprocessed txt files

%%
%Now for each folder, go through and process the .raw files
for foldernum = 1:length(datafolders)
    
    %Keep the date-folder structure for now and store each processed text
    %and mat files in separate folders under the date:
    txtoutdir = fullfile(masteroutpath,datafolders{foldernum},'converted_txtfiles'); %where you want .txt files to be stored
    matoutdir = fullfile(masteroutpath,datafolders{foldernum},'mat_outfiles'); % where you'd like processed .mat files to be
    
    %If these directories don't exist, make them!
    if ~isequal(exist(fullfile(masteroutpath,datafolders{foldernum}),'dir'),7)
        mkdir(fullfile(masteroutpath,datafolders{foldernum}))
    end
    if ~isequal(exist(txtoutdir,'dir'),7)
        mkdir(txtoutdir)
    end
    if ~isequal(exist(matoutdir,'dir'),7)
        mkdir(matoutdir)
    end
    
    
    %find dat data!
    tempdatasource=fullfile(mastersourcepath,datafolders{foldernum});
    sourcefiles=dir([tempdatasource '/*.raw']);
    
    %if filelist comes back empty, there may be a separate folder that
    %contains the raw files:
    if isempty(sourcefiles)
        
        fprintf('No .raw files in top directory...searching for raws folder under %s....\n',datafolders{foldernum});
        subd = dir([tempdatasource '/']);
        isub = [subd(:).isdir]; %# returns logical vector
        subfoldernames = {subd(isub).name}'; %names of folders
        temp=regexp(subfoldernames,'raws');
        subfolder=subfoldernames(cellfun('isempty',temp)==0);
        
        tempdatasource=fullfile(tempdatasource,subfolder{:});
        sourcefiles=dir([tempdatasource '/*.raw']);
        
        if isempty(sourcefiles)
            fprintf('Uh oh - cannot find a raws folder either!\n')
        end
    end
    
    
    % the plan:
    % loop through all the fnum *.raw files in that directory
    % each time:
    % -> convert raw into .txt file using MPR cal file
    % -> read in that .txt file and convert to matlab variables
    % -> convert raw into .txt using 284 (downwelling irrad) cal file
    % -> read into .txt and convert into other matlab variables
    % -> repeat for 285 (solar reference)
    % make some informative plots and flag bad data
    % save all matlab vars that you care about into mat file names as input raw file
    
    
    % now using satcon 1.5.5_2 NEW SWTCHES
    % HELPFUL TO READ THE MANUAL!
    % xtra output info, convert phys vals + immersion, proc time tags & format,
    % precision 4
    % buffersize 32768, overwrite old files, kill end window, no conversion log
    
    %satconswitches = ' /x /ci /tf /p10 /OverWrite=yes /EW=y /cl=nul';
    
    %Conversion switches:
    %/x extra output heading (info about conversion, which file used, immersion-corrected, etc.)
    %/c converts measured voltages into physical/calibrated/useful values
    %/ci converts measured voltages into physical/calibrated/useful values and then applies an immersion coefficient for sensors operated under water
    %/tf append time for each frame with readable format
    %/p precision for floating point values
    
    %SatCON Override Switches:
    %/OverWrite=yes: overwrite existing files if present
    %/EW (EndWindow): Close conversion progress window after completed conversion
    %/CL (Conversion log)=nul: do not maintain a conversion log
    %
    numfiles = length(sourcefiles);
    
    for filenum = 1:numfiles
        
        infile = fullfile(tempdatasource, sourcefiles(filenum).name);
        outfile = fullfile(matoutdir, [sourcefiles(filenum).name(1:end-3) 'mat']);
        
        fprintf('Processing %s into text into matlab files\n', sourcefiles(filenum).name);
        empty_flag=0;
        
        for ftype = 1:5,
            
            switch ftype
                
                case 1 %MPR sensor (depth, temperature, etc.)
                    calfilestr = [calfiledir '\MPR106a.cal'];
                    txtfile = fullfile(txtoutdir,[sourcefiles(filenum).name(1:end-4) '-MPR.txt']);
                    satconswitches = ' /x /c /tf /p10 /OverWrite=yes /EW=y /cl=nul';
                    varstub ='mpr';
                    
                case 2 %Solar standard, light measurments
                    calfilestr = [calfiledir '\Hse285c.cal'];
                    txtfile = fullfile(txtoutdir,[sourcefiles(filenum).name(1:end-4) '-ESTDL.txt']);
                    satconswitches = ' /x /c /tf /p10 /OverWrite=yes /EW=y /cl=nul';
                    varstub = 'esl';
                    
                case 3 %Solar standard, dark measurments
                    calfilestr = [calfiledir '\HED285c.cal'];
                    txtfile = fullfile(txtoutdir,[sourcefiles(filenum).name(1:end-4) '-ESTDD.txt']);
                    satconswitches = ' /x /c /tf /p10 /OverWrite=yes /EW=y /cl=nul';
                    varstub = 'esd';
                    
                case 4 %Downwelling radiation, light measurments
                    calfilestr = [calfiledir '\HPE284c.cal'];
                    txtfile = fullfile(txtoutdir,[sourcefiles(filenum).name(1:end-4) '-EDOWNL.txt']);
                    satconswitches = ' /x /ci /tf /p10 /OverWrite=yes /EW=y /cl=nul'; %want immersion corrected
                    varstub = 'edl';
                    
                case 5 %Downwelling radiation, dark measurments
                    calfilestr = [calfiledir '\PED284c.cal'];
                    txtfile = fullfile(txtoutdir,[sourcefiles(filenum).name(1:end-4) '-EDOWND.txt']);
                    satconswitches = ' /x /ci /tf /p10 /OverWrite=yes /EW=y /cl=nul'; %want immersion corrected
                    varstub = 'edd';
                    
                otherwise
                    fprintf('Error in switch\n');
                    fclose all;
                    return;
                    
            end;
            
            
            % need to satcon these data into temp files -----------------------------------------------------------
            %fnamestr = ['"C:\Program Files (x86)\Satlantic\SatCon\SatCon.exe" ' infile ' ' txtfile ' ' calfilestr];
            fnamestr = ['"C:\Program Files (x86)\Satlantic\SatCon\SatCon.exe" ' infile ' ' txtfile ' ' calfilestr];
            cmdstr = [fnamestr satconswitches];
            [s,w] = system(cmdstr);
            %             %-------------------------------------------------------------------------------------------------------
            % As of 9/22, these strings were working from the sosiknas1
            % hyperpro directory....I'm wondering if it's just not a file
            % issue.....rather than a string issue?
            
            %             calfilestr ='\\sosiknas1\Lab_data\MVCO\HyperPro_Radiometer\calibration_files\CD_Nov2008\Calibration_Files\HPE284c.cal';
            %             infile='\\sosiknas1\Lab_data\MVCO\HyperPro_Radiometer\11Oct2007\2007-284-182304.raw';
            %             txtfile='\\sosiknas1\Lab_data\MVCO\HyperPro_Radiometer\text.txt';
            
            %Once have made the txt files, reimport to matlab to save as .mat files:
            fprintf('reading in %s data\n', varstub);
            
            %could in theory find the formatSpec once for .txt import .... but doesn't look too taxing at the moment...
            
            %find out how many columns you have in the file:
            fid = fopen(txtfile);
            line = fgetl(fid); k=1;
            while ~strcmp(line,'')   %look for this space between header lines and column labels
                line = fgetl(fid);
                k=k+1;
            end
            line = fgetl(fid);
            fclose(fid);
            
            tmp = textscan(line, '%s'); %returns cell array
            n = length(tmp{:}); %number of columns
            
            %textscan magic:
            %first gather column header info
            fileID = fopen(txtfile,'r');
            formatSpec=repmat('%s',1,n);
            textscan(fileID, '%[^\n\r]', k-1, 'ReturnOnError', false); %throw away until row 3
            column_headers = textscan(fileID, formatSpec, (k+1)-k+1, 'Delimiter', '\t', 'ReturnOnError', false);
            
            %reshuffle header info:
            column_headers=cat(2,column_headers{:})';
            
            %find known string fields:
            strind1=find(strcmp('DATETAG',column_headers(:,1))==1);
            strind2=find(strcmp('TIMETAG2',column_headers(:,1))==1);
            
            formatSpec=repmat('%f',1,n);
            formatSpec(2*strind1-1:2*strind1)='%s';
            formatSpec(2*strind2-1:2*strind2)='%s';
            
            dataArray = textscan(fileID, formatSpec, 'Delimiter', '\t', 'EmptyValue' ,NaN,'ReturnOnError', false,'TreatAsEmpty','Error 6');
            fclose(fileID);
            
            %combine and reformat time to something useful:
            %gather time pieces separately:
            time1=regexp(dataArray{strind1},'(?<year>\d{4})-(?<yrday>\d{1,3})','names');
            time1=cell2mat(time1);
            
            time2=regexp(dataArray{strind2},'(?<hour>\d{2}):(?<min>\d{2}):(?<sec>\d{2})\.(?<nano>\d{3})','names');
            time2=cell2mat(time2);
            
            %converted string to numeric values and calculate time:
            if ~isempty(time1) && ~isempty(time2)
                jl_time=arrayfun(@(x) str2num(x.yrday),time1) + ... %julian day
                    (1/24)*arrayfun(@(x) str2num(x.hour),time2) + ... %hours
                    (1/24/60)*arrayfun(@(x) str2num(x.min),time2) + ... %minutes
                    (1/24/60/60)*(arrayfun(@(x) str2num(x.sec),time2) + 0.001*arrayfun(@(x) str2num(x.nano),time2)); %seconds
                
                matdate = jl_time + datenum([repmat('1-0-',length(jl_time),1) cell2mat({time1(:).year}')]);
                
                dataArray(:,n+1)={jl_time}; column_headers{n+1,1}='JULIAN_TIME';
                dataArray(:,n+2)={matdate}; column_headers{n+2,1}='MATLAB_TIME';
                
            else
                empty_flag=1;
                fprintf('This file seems to have no data in it....\n')
            end
            
            eval([varstub '_data=dataArray;'])
            eval([varstub '_hdr=column_headers;'])
            
            clear dataArray column_headers txtfile
            
            if ftype==1 %just do this once :)
                %PRESSURE TARE - this is recorder in a header region of each .raw file:
                %--------------------------------------------------------------------------------
                delimiter = ' ';
                startRow = 11;
                endRow = 11;
                % Format string for each line of text:
                %   column1: text (%s)
                %	column2: text (%s)
                formatSpec = '%s%s%*s%*s%*s%*s%*s%*s%*s%*s%[^\n\r]';
                fileID = fopen(infile,'r');
                textscan(fileID, '%[^\n\r]', startRow-1, 'ReturnOnError', false);
                dataArray = textscan(fileID, formatSpec, endRow-startRow+1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'ReturnOnError', false);
                fclose(fileID);
                
                pressure_tare = [dataArray{1:end-1}];
                pressure_tare=str2num(pressure_tare{2}); %+ distance back up to downwelling sensor??
                
                clearvars filename delimiter startRow endRow formatSpec fileID dataArray
            end
            
            %%
            %GET LAT AND LONGITUDE AS WELL!!!!
            
            
            %If can't find a pressure tare?
            
        end; % of ftype
        
        %
        %PHASE 2 - CLEANING UP AND PROCESSING OF DATA:
        
        % Once we have all the variables in the workspace
        % parse out wavelength data into matrices
        % correct light measurements with shutter dark measurements
        % create a common time axis for the light measurements
        % find all the downcasts
        % calculate PAR, attentuation coefficient
        
        %Wavelength processing: pull out the light data into a matrix for
        %easier processing:
        
        varstubs={'mpr','edl','edd','esl','esd'};
        
        if empty_flag==0
            for q=2:5;
                eval(['hdr=' varstubs{q} '_hdr;'])
                temp=regexp(hdr,'\d{3}\.\d{2}','match');
                ind=find(cellfun('isempty',temp)==0);
                temp=temp(ind);
                temp=[temp{:}]';
                
                eval([varstubs{q} '_lambdas=str2num(char(temp));'])
                
                %extract just the wavelength measurments:
                eval([varstubs{q} '_wv_data=cell2mat(' varstubs{q} '_data(:,ind));'])
            end
            
            %Check that all wavelengths are the same!
            if sum(edl_lambdas-edd_lambdas) < 0.01;
                lambdas_downwell=edl_lambdas;
            end
            if sum(esl_lambdas-esd_lambdas) < 0.01;
                lambdas_solarstd=esl_lambdas;
            end
            
            clear *_lambdas
            
            % CORRECT LIGHT MEASUREMETNS WITH DARK MEASUREMENTS:
            %--------------------------------------------------------------------------------------------------
            % Find the nearest dark measurement and subtract this value:
            
            
            nd=length(edl_data{end-1}); %how many measurments, julian time column
            edd2edl=nan(nd,1); %index matching dark to light measurements
            for j=1:nd
                [mm, ii]=min(abs(edd_data{end-1}-edl_data{end-1}(j)));
                edd2edl(j)=ii;
            end
            
            adj_edl=edl_wv_data-edd_wv_data(edd2edl,:);
            
            %and the same for the solar standard:
            ns=length(esl_data{end-1}); %how many light measurments
            esd2esl=nan(ns,1);
            for j=1:ns
                [mm, ii]=min(abs(esd_data{end-1}-esl_data{end-1}(j)));
                esd2esl(j)=ii;
            end
            
            adj_esl=esl_wv_data-esd_wv_data(esd2esl,:);
            
            % CALCULATE PAR
            %--------------------------------------------------------------------------------------------------
            
            h=6.625e-34; %Planck's constant
            c=3e8; %speed of light
            
            %Equation from ProSoft manual:
            %Integral from 400 to 700 (nm) of (lambda/hc)*E(lambda)*dlambda
            %This is finding out how many photons were there for each wavelength:
            %So, per photon of a given wavelength, the energy (in Joules) will be:
            %(E = hc/lambda)
            % We divide our recorded energy by this to get number of
            % photons and then just integrate that!
            
            [~, i400]=min(abs(lambdas_solarstd-400)); %indexes of closest wv to 400
            [~, i700]=min(abs(lambdas_solarstd-700)); %indexes of closest wv to 700
            
            photons = repmat((1e-9/(h*c))*lambdas_solarstd(i400:i700)',size(adj_esl(:,i400:i700),1),1)*(1e-6).*adj_esl(:,i400:i700); % number of photons per cm2 per second
            esl_PAR = (6.02214e-23)*1e6*100*trapz(lambdas_solarstd(i400:i700),photons,2); %in micromol photons/m^2/s (converts to moles, to micromoles, then from cm2 to m2)
            
            clearvars i400 i700 photons
            
            [~, i400]=min(abs(lambdas_downwell-400)); %indexes of closest wv to 400
            [~, i700]=min(abs(lambdas_downwell-700)); %indexes of closest wv to 700
            
            photons = repmat((1e-9/(h*c))*lambdas_downwell(i400:i700)',size(adj_edl(:,i400:i700),1),1)*(1e-6).*adj_edl(:,i400:i700); % number of photons per cm2 per second
            edl_PAR = (6.02214e-23)*1e6*100*trapz(lambdas_downwell(i400:i700),photons,2); %in micromol photons/m^2/s (converts to moles, to micromoles, then from cm2 to m2)
            
            %             edl_PAR = trapz(edl_lambdas(i400:i700),adj_edl(:,i400:i700),2); %in uW/cm^2
            %             edl_PAR = (1e4/1e6)*edl_PAR;  %in W/m^2
            
            
            if (max(esl_PAR)-min(esl_PAR))/max(esl_PAR) > 0.10 %10% variation?
                solarstd_flag=1;
            else
                solarstd_flag=0;
            end
            
            % TIME SYNC WITH MPR SENSOR:
            %--------------------------------------------------------------------------------------------------
            %assumes more time measurements on mpr sensor than for light
            %sensors:
            mpr2esl=nan(ns,1);
            for j=1:ns
                ii=find(mpr_data{end-1} == esl_data{end-1}(j));
                if isempty(ii) %meaning no exact time match
                    [mm, ii]=min(abs(mpr_data{end-1}-esl_data{end-1}(j)));
                end
                mpr2esl(j)=ii(1);
            end
            
            mpr2edl=nan(nd,1);
            for j=1:nd
                ii=find(mpr_data{end-1} == edl_data{end-1}(j));
                if isempty(ii) %meaning no exact time match
                    [mm, ii]=min(abs(mpr_data{end-1}-edl_data{end-1}(j)));
                end
                mpr2edl(j)=ii(1); %just in case there are two equal options
            end
            
            mpr_mattime=mpr_data{end};
            mpr_jultime=mpr_data{end-1};
            
            % Some useful quantities from the MPR:
            %-----------------------------------------------------------------------
            
            jj= strcmp(mpr_hdr,'Pres')==1;
            depth=mpr_data{:,jj}-pressure_tare;
            
            %to make some of the calculations easier, trim off the first
            %few meters and above the bottom:
            tdepth=depth; %t for trim
            tdepth(tdepth<1)=NaN;
            tdepth(tdepth > max(tdepth)-1) = NaN;
            %this avoids potentially dealing with wave focusing or profile
            %bobbing up and down on the bottom!
            
            fprintf('Pressure tare: %2.2f         min deth: %2.2f max deth: %2.2f \n',pressure_tare,min(depth),max(depth));
                     
            
            if plot_flag==1
                
                clf
                % SANITY CHECK PLOTS:
                %dark and light adjusted data
                subplot(2,3,1,'replace'), hold on
                for j=20:10:130 %a few representative wavelengths
                    plot(esd_data{end-1},esd_data{j},'.-','color',[0.4 0.4 0.4])
                end
                plot(esd_data{end-1},min(esd_wv_data,[],2),'k-','linewidth',2)
                plot(esd_data{end-1},max(esd_wv_data,[],2),'k-','linewidth',2)
                datetick('x','keeplimits')
                xlim([min(mpr_mattime) max(mpr_mattime)])
                
                ylabel('Radiation for each \lambda (W m^{-2})')
                title('Dark Solar Standard')
                
                subplot(2,3,4,'replace'), hold on
                for j=20:10:130
                    plot(edd_data{end-1},edd_data{j},'.-','color',[0.4 0.4 0.4])
                end
                plot(edd_data{end-1},min(edd_wv_data,[],2),'k-','linewidth',2)
                plot(edd_data{end-1},max(edd_wv_data,[],2),'k-','linewidth',2)
                datetick('x','keeplimits')
                xlim([min(mpr_mattime) max(mpr_mattime)])
                
                ylabel('Radiation for each \lambda (W m^{-2})')
                title('Dark Downwelling')
                
                subplot(2,3,2,'replace'), hold on
                plot(esl_data{end},esl_PAR,'.-','color',[0 0.5 1])
                plot(edl_data{end},edl_PAR,'.-','color',[0.5 0.5 0.5])
                %plot(edl_data{end}(downcast_ind),edl_PAR(downcast_ind),'b.') %%%%%%%%%%%
                datetick('x','keeplimits')
                ylabel('PAR (\mumol/m^{2}/s)')
                title('Dark-corrected PAR')
                xlim([min(mpr_mattime) max(mpr_mattime)])
                
                subplot(2,3,5,'replace'), hold on
                plot(mpr_mattime,depth,'.-')
%                 if ~isempty(zz)
%                     plot(mpr_mattime(downcast_ind),tdepth(downcast_ind),'.-')
%                 end
                set(gca,'Ydir','reverse')
                ylabel('Depth (m)')
                xlim([min(mpr_mattime) max(mpr_mattime)])
                datetick('x','keeplimits')
                title('Cast Trajectory')
                
                %PAR with depth:
                subplot(2,3,3,'replace'), hold on
                plot(edl_PAR,depth(mpr2edl),'.')
                plot(esl_PAR,depth(mpr2esl),'g.')
%                 if ~isempty(zz)
%                     plot(edl_PAR(startcast:returncast),par_depths(startcast:returncast),'r.')
%                 end
                set(gca,'YDir','reverse')
                ylabel('Depth (m)')
                xlabel('PAR (\mumol/m^{2}/s)')
                title('PAR with Depth')
                legend('PAR depth','PAR solar std','location','SouthEast')
                
                %What would be used for an attentuation coefficient fit:
                subplot(2,3,6,'replace'), hold on
                plot(log(edl_PAR),depth(mpr2edl),'.'), set(gca,'Ydir','reverse')
                hold on
              % plot(log(edl_PAR(startcast:returncast)),par_depths(startcast:returncast),'.-')
              % plot(K(1)+K(2)*par_depths,par_depths,'-')
                
                xlabel('log(PAR) (\mumol/m^{2}/s)')
                ylabel('Depth (m)')
              %  title(['Attenuation coefficient check: ' num2str(K(2))])
              %  legend('Downwelling PAR','Downcast PAR','polynomial fit','location','northwest')
               
            end %plot_flag
            
            %keyboard
            
            data_comment = input('Good data? Enter a comment for future use :)\n');
            
            raw_data_log=[raw_data_log; [{datafolders{foldernum}} {sourcefiles(filenum).name} {'data present'} {data_comment}] ];
            
            %SAVE THAT DATA!
            %--------------------------------------------------------------------------------------------------
            eval(['save ' outfile ' *mpr* *edl* *esl* *edd* *esd* solarstd_flag *PAR* pressure_tare depth *time*'])
        else
            raw_data_log=[raw_data_log; [{datafolders{foldernum}} {sourcefiles(filenum).name} {'no data for at least one of the files'} {''}] ];
        end %empty_flag
        
        clearvars -except numfiles sourcefiles mastersourcepath masteroutpath calfiledir datafolders filenum
        %
    end;   % of filenum   
    
end;  %of foldernum
