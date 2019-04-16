% processPROII_MVCO
% Core pieces of this script that call the SatCon software were provided by Sam Laney at WHOI. 
% Other pieces, specific to MVCO processing, dark measurement subtraction,
% PAR calculation were contributed by KRHC 
% Script imports and processes data from 3 sensors on HyperPro, saves intermediate txt files


%addpath /Users/kristenhunter-cevera/Documents/MATLAB/mvco_tools/

clear all;
close all; %for any open files

plot_flag=1;
do_satcon=1; %if do not need to run satcon and want to only work wth preprocessed txt files, set to 0

%Load in initial screeing comments, if haven't made this file yet, script will ask you to screen now:
try
    load(fullfile('\\sosiknas1\lab_data\mvco\HyperPro_Radiometer\','initial_data_notes.mat'))
    comment_flag=0; %notes/comments have been made and loaded
catch
    fprintf(['Hmmm...looks like there are not any initial comments on the data...\n Initiating plots for screening and entering comments\n'])
    data_notes={};
    data_note_titles={'date' 'file name' 'data present' 'if file has data, comment on data'};
    comment_flag=1; %need to make comments on data
    plot_flag=1; %Automatically changed it
end

mastersourcepath=fullfile('\\sosiknas1\lab_data\mvco\HyperPro_Radiometer\raw_data\'); % path to folders with raw data...
%mastersourcepath=fullfile(pwd);%'/Volumes/Lab_data/SummerStudents/2013_Marco/MVCO_hyperpro'; %top level directory
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
% save_txtfiles=1; %if do or do not want to save the intermediate txt files
% that are generated from satcon, need at least a temporary one though to reload into matlab!

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
    
    
    % for each *.raw files in a directory and for each sensor that needs to
    % be imported:
    % 
    % convert raw into .txt file using appropriate cal file
    % read in that .txt file and convert to desired matlab variables
    % repeat for remaining sensors
    
    %below, we convert data from MPR (depth, time), 284 (downwelling irrad)
    % and 285 (solar reference)
    
    %From the SatCon Manual v1.5, pages 10-15: 
    
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
    
    
    for filenum = 1:length(sourcefiles)
        
        filename=sourcefiles(filenum).name;
        
        infile = fullfile(tempdatasource, filename);
        outfile = fullfile(matoutdir, [filename(1:end-3) 'mat']);
        
        fprintf('Processing %s into text into matlab files\n', sourcefiles(filenum).name);
        empty_flag=0;
        
        for sensortype = 1:5,
            
            switch sensortype
                
                case 1 %MPR sensor (depth, temperature, etc.)
                    calfilestr = [calfiledir '\MPR106a.cal'];
                    txtfile = fullfile(txtoutdir,[sourcefiles(filenum).name(1:end-4) '-MPR.txt']);
                    satconswitches = ' /x /c /tf /p10 /OverWrite=yes /EW=y /cl=nul';
                    sensor ='mpr';
                    
                case 2 %Solar standard, light measurments
                    calfilestr = [calfiledir '\Hse285c.cal'];
                    txtfile = fullfile(txtoutdir,[sourcefiles(filenum).name(1:end-4) '-ESTDL.txt']);
                    satconswitches = ' /x /c /tf /p10 /OverWrite=yes /EW=y /cl=nul';
                    sensor = 'esl';
                    
                case 3 %Solar standard, dark measurments
                    calfilestr = [calfiledir '\HED285c.cal'];
                    txtfile = fullfile(txtoutdir,[sourcefiles(filenum).name(1:end-4) '-ESTDD.txt']);
                    satconswitches = ' /x /c /tf /p10 /OverWrite=yes /EW=y /cl=nul';
                    sensor = 'esd';
                    
                case 4 %Downwelling radiation, light measurments
                    calfilestr = [calfiledir '\HPE284c.cal'];
                    txtfile = fullfile(txtoutdir,[sourcefiles(filenum).name(1:end-4) '-EDOWNL.txt']);
                    satconswitches = ' /x /ci /tf /p10 /OverWrite=yes /EW=y /cl=nul'; %want immersion corrected
                    sensor = 'edl';
                    
                case 5 %Downwelling radiation, dark measurments
                    calfilestr = [calfiledir '\PED284c.cal'];
                    txtfile = fullfile(txtoutdir,[sourcefiles(filenum).name(1:end-4) '-EDOWND.txt']);
                    satconswitches = ' /x /ci /tf /p10 /OverWrite=yes /EW=y /cl=nul'; %want immersion corrected
                    sensor = 'edd';
                    
                otherwise
                    fprintf('Error in switch\n');
                    fclose all;
                return;
                    
            end;
            
            
            % need to satcon these data into temp files -----------------------------------------------------------
            %fnamestr = ['"C:\Program Files (x86)\Satlantic\SatCon\SatCon.exe" ' infile ' ' txtfile ' ' calfilestr];
            if do_satcon==1
                cmdstr = ['"C:\Program Files (x86)\Satlantic\SatCon\SatCon.exe" ' infile ' ' txtfile ' ' calfilestr satconswitches];       
                [s,w] = system(cmdstr); %have matlab call the function
            end
            
            
            % ----------------REIMPORT DATA FROM TXT FILES ----------------------------------------------------
            
            %Once have made the txt files, reimport to matlab to save as .mat files:
            fprintf('reading in %s data\n', sensor);
            
            %could in theory find the formatSpec once for .txt import .... but doesn't look too taxing at the moment...
            
            %in the txtfiles, there are usually two header lines until you
            %get to the data columns, find this and then record the line
            %number and number of columns:
            
            %find out how many columns you have in the file:
            fid = fopen(txtfile);
            line = fgetl(fid); row=1;
            while isempty(strfind(line,'Index'))   %look for this key word - this denotes the line that has the column headers
                line = fgetl(fid);
                row=row+1;
            end
            %line = fgetl(fid);
            fclose(fid);
            row=row-1; %important - textscan doesn't seem to recognize an empty line with the below formatpsec
            
%           tmp = textscan(line, '%s'); %returns cell array
%           n = length(tmp{:}); %number of columns
            n = length(strsplit(line,'\t')); %number of columns
            
            %textscan magic:
            %first gather column header info
            fileID = fopen(txtfile,'r');
            formatSpec=repmat('%s',1,n);
            textscan(fileID, '%[^\n\r]', row-1, 'ReturnOnError', false); %reads in lines, 2 times, throw away until row 3
            column_headers = textscan(fileID, formatSpec, 2, 'Delimiter', '\t', 'ReturnOnError', false); %should be at column line, read in twice for column title and units
            
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
            
            eval([sensor '_data=dataArray;'])
            eval([sensor '_hdr=column_headers;'])
            
            clear dataArray column_headers txtfile
            
            %GATHER SATFILE HEADER INFORMATION: PRESSURE TARE, LAT/LON, OPERATOR
            
            if ftype==1 %just do this once :)
                
                delimiter = ' ';
                endRow = 19;
                formatSpec = '%q%q%q%q%q%q%q%[^\n\r]'; %'%q%q%q%q%q%f%q%[^\n\r]';
                
                fileID = fopen(infile,'r');
                dataArray = textscan(fileID, formatSpec, endRow, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'EmptyValue' ,NaN,'ReturnOnError', false);
                fclose(fileID);
                
                dataArray(6) = cellfun(@(x) num2cell(x), dataArray(6), 'UniformOutput', false);
                satfile_hdr = [dataArray{1:end-1}];
                
                % Clear temporary variables
                clearvars delimiter endRow formatSpec fileID dataArray ans
                
                tempstruc(filenum).file = filename;
                tempstruc(filenum).cruiseID = satfile_hdr{1,2};
                tempstruc(filenum).operator = [satfile_hdr{2,2} ' ' satfile_hdr{2,3}];
                tempstruc(filenum).lat = [satfile_hdr{3,2} ' ' satfile_hdr{3,3}];
                tempstruc(filenum).lon = [satfile_hdr{4,2} ' ' satfile_hdr{4,3}];
                tempstruc(filenum).timestamp=datenum([char(satfile_hdr{9,3}) '-' char(satfile_hdr{9,4}) '-' char(satfile_hdr{9,6}) ' ' char(satfile_hdr{9,5})]);
                
                if tempstruc(filenum).timestamp > datenum('1-0-2016')
                    disp('Ummm- somethings not right with the timestamp here...')
                    keyboard
                end
                
                %PRESSURE TARE - this is recorder in a header region of each .raw file:
                tempstruc(filenum).pressure_tare=str2num(satfile_hdr{11,2});
                
                if tempstruc(filenum).pressure_tare < 8
                    fprintf('Is this the right pressure tare? %2.3f',tempstruc(filenum).pressure_tare)
                    keyboard
                    %If indeed you find yourself here, what I did if there
                    %was no pressure tare, was to look at the pressures
                    %with record number and you should be able to see when
                    %it was bobbing at the surface. I just took an average
                    %of this as the pressure tare:
                    %figure, plot(mpr_data{:,7})
                    %set(gca,'YDir','reverse')
                    
                    %for 2008-255-183514.raw, this became:
                    %tempstruc(filenum).pressure_tare=nanmean(mpr_data{:,7}(2960:3060))
                    
                    %for 2007-346-191403.raw, 2008-267-163905.raw, 2009-116-214302.raw, 2009-116-214414.raw  these aren't a cast so just used:
                    %tempstruc(filenum).pressure_tare=min(mpr_data{:,7})
                end
                
                pressure_tare=tempstruc(filenum).pressure_tare; %for use later on in teh script
                
            end %if ftype==1
            
        end; % of ftype, and all available variables should be loaded in
        
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
        
        sensors={'mpr','edl','edd','esl','esd'};
        
        if empty_flag==0
            for q=2:5;
                eval(['hdr=' sensors{q} '_hdr;'])
                temp=regexp(hdr,'\d{3}\.\d{2}','match');
                ind=find(cellfun('isempty',temp)==0);
                temp=temp(ind);
                temp=[temp{:}]';
                
                eval([sensors{q} '_lambdas=str2num(char(temp));'])
                
                %extract just the wavelength measurments:
                eval([sensors{q} '_wv_data=cell2mat(' sensors{q} '_data(:,ind));'])
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
            
            depth=mpr_data{:,7}-pressure_tare;
            fprintf('Pressure tare: %2.2f         min deth: %2.2f max deth: %2.2f \n',pressure_tare,min(depth),max(depth));
            
            if plot_flag==1 || comment_flag==1
                
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
                xlim([min(mpr_jultime) max(mpr_jultime)])
                
                ylabel('Radiation for each \lambda (W m^{-2})')
                title('Dark Solar Standard')
                
                subplot(2,3,4,'replace'), hold on
                for j=20:10:130
                    plot(edd_data{end-1},edd_data{j},'.-','color',[0.4 0.4 0.4])
                end
                plot(edd_data{end-1},min(edd_wv_data,[],2),'k-','linewidth',2)
                plot(edd_data{end-1},max(edd_wv_data,[],2),'k-','linewidth',2)
                datetick('x','keeplimits')
                xlim([min(mpr_jultime) max(mpr_jultime)])
                
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
            
            %if have not already made this log:
            if comment_flag==1
                data_comment = input('Good data? Enter a comment for future use :)\n');
                data_notes=[data_notes; [{datafolders{foldernum}} {sourcefiles(filenum).name} {'data present'} {data_comment}] ];
                
                if strcmp(data_comment,'not a cast')
                    tempstruc(filenum).emptyflag = 2;
                else %good casts
                    tempstruc(filenum).emptyflag = 0;
                end
                
            else
                ii=find(strcmp(filename, raw_data_log(:,2))==1);
                
                if strcmp(raw_data_log{ii,4},'not a cast')
                    tempstruc(filenum).emptyflag = 2;
                else %good casts
                    tempstruc(filenum).emptyflag = 0;
                end
                
            end %end of comment_flag
            
            %SAVE THAT DATA!
            %--------------------------------------------------------------------------------------------------
            %eval(['save ' outfile ' *mpr* *edl* *esl* *edd* *esd* solarstd_flag *PAR* pressure_tare depth *time*'])
            eval(['save ' outfile ' mpr_data edl_data esl_data edd_data satfile_hdr'])    %save the rawer products as matrices
            
            %but the products in the structure?
            tempstruc(filenum).adj_esl=adj_esl;
            tempstruc(filenum).adj_edl=adj_edl;
            tempstruc(filenum).esl_PAR=esl_PAR;
            tempstruc(filenum).edl_PAR=edl_PAR;
            tempstruc(filenum).edl_ind=mpr2edl;
            tempstruc(filenum).esl_ind=mpr2esl;
            tempstruc(filenum).mprtime=mpr_mattime;
            tempstruc(filenum).solarflag=solarstd_flag;
            tempstruc(filenum).depth=depth;
            tempstruc(filenum).wavelen_solarst=lambdas_solarstd;
            tempstruc(filenum).wavelen_downwell=lambdas_downwell;
            
        else
            tempstruc(filenum).emptyflag = 1;
            if comment_flag==1
                data_notes=[data_notes; [{datafolders{foldernum}} {sourcefiles(filenum).name} {'no data for at least one of the files'} {''}] ];
            end
        end % empty_flag
        
        clearvars -except numfiles sourcefiles matoutdir txtoutdir numfiles mastersourcepath masteroutpath calfiledir datafolders filenum tempstruc foldernum tempdatasource plot_flag raw_data_log
        
    end;   % of filenum
    
    eval(['data_' datafolders{foldernum} '=tempstruc;'])
    eval(['save ' matoutdir '\data_' datafolders{foldernum} '.mat ' 'data_' datafolders{foldernum}])
    clear tempstruc
    
end;  %of foldernum

if comment_flag==1
    save initial_data_notesB data_note_titles data_notes
end