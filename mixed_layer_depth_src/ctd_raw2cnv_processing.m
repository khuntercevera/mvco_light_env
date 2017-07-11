
%Script to process CTD casts into raw .cnv files with Sea-Bird SBE Data
%Processing software:

CTDpath=fullfile('\\maddie\TaylorF\from_Samwise\data\MVCO\'); % path to folders with raw CTD data...
%mastersourcepath=fullfile('/Volumes/TaylorF/from_Samwise/data/MVCO/'); % path to folders with raw CTD data...

%%
%It appears that the CTD casts are stored in a few files:
%Those with the Tioga designation
%Survery Cruises
%todownload

d = dir(CTDpath);
isub = [d(:).isdir]; %# returns logical vector
foldernames = {d(isub).name}'; %names of folders
temp=regexp(foldernames,'Tioga_\d{3}'); %find all folders with this naming scheme
datafolders = foldernames(cellfun('isempty',temp)==0);

%Survey Cruises:
d = dir(fullfile([CTDpath '/SurveyCruises/']));
isub = [d(:).isdir]; %# returns logical vector
foldernames = {d(isub).name}'; %names of folders
temp=regexp(foldernames,'Tioga_\d{3}'); %find
ii=find(cellfun('isempty',temp)==0);
datafolders =[datafolders;  strcat(cellstr(repmat('/SurveyCruises/',length(ii),1)),foldernames(ii))];

%todownload folder:
d = dir(fullfile([CTDpath '/todownload/']));
isub = [d(:).isdir]; %# returns logical vector
foldernames = {d(isub).name}'; %names of folders
temp=regexp(foldernames,'(Tioga_\d{3})|(ti\d{3})'); %find
ii=find(cellfun('isempty',temp)==0);
datafolders =[datafolders;  strcat(cellstr(repmat('/todownload/',length(ii),1)),foldernames(ii))];

datafiles={};
%% Gather a complete list of data files by checking each of these folder locations for raw data files!

for j=1:length(datafolders)
    
    %check to see what files/folders are in each one of the data folders:
    %a bit more complicated because Matlab does not appear to have a
    %recursive search function....arghhh...
    
    filepath=fullfile(CTDpath,datafolders{j});
    filelist=dir(filepath);
    disp(['Searching folder: ' filepath])
    
    % keep on searching until have found a .hex or a .dat file...
    flag=1; %datafiles={};
    
    while flag
        
        rawind=find(cellfun('isempty',regexp({filelist(:).name}','(\w*\.hex)|(\w*_\w*\.dat)'))==0);
        
        if isempty(rawind)
            %disp('No data files here...looking in subdirectories...')
            
            subfolders=regexpi({filelist([filelist(:).isdir]).name}','(\w*ctd)|(original_files)','match');
            subind=find(cellfun('isempty',subfolders)==0);
            if length(subind)==1
                %disp(['found a ' char(subfolders{subind}) ' folder'])
                filepath=fullfile(filepath,subfolders{subind});
                filelist=dir(char(filepath));
            else
                'Uh-Oh! Two subfolders?'
                keyboard
                %if j=41:
                %subind=3;
                %               filepath=fullfile(filepath,subfolders{subind});
                %               filelist=dir(char(filepath));
                %               dbcont
            end
            
        else
            %disp(['I found ' num2str(length(rawind)) ' file(s)...adding to list'])
            datafiles=[datafiles; fullfile(filepath,{filelist(rawind).name}')];
            flag=0;
        end
    end
    
    
    
end

%%

for q=1:length(datafiles)
    
    disp(['Processing raw data file for: ' datafiles{q}])
    
    infile=datafiles{q};
    confile=regexprep(datafiles{q},'(\.hex)|(\.dat)','\.CON'); %use the corresponding .con file - this should be with the data file!
    
    outputdir=fullfile('\\sosiknas1\Lab_data\MVCO\','processed_CTD_casts');
    outtemp=strsplit(datafiles{q},'\');
    outfile=outtemp{end}(1:end-4);
    
    psafile='C:\Users\heidi\AppData\Local\Sea-Bird\SBEDataProcessing-Win32\DatCnv.psa';
    
    %/s = start processing now
    %/m = minimize window
    
    cmdstr = ['"C:\Program Files (x86)\Sea-Bird\SBEDataProcessing-Win32\datcnvw.exe" ' '/s /c' confile ' /i' infile ' /o' outputdir ' /f' outfile ' /p' psafile'];
    [s,w] = system(cmdstr);
    
end