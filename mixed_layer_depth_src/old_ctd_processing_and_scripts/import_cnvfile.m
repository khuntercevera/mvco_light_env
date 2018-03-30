function [lat,lon, UTC_time,upload_time,header,data]=import_cnvfile(filename)
% Script that uses MATLAB's textscan to import the .cnv file from the
% converted CTD files into:

%lat/long & time of cast
%name of data fields
%actual data

%Please note that this assumes a standard format to the CTD files, if
%software or file format changes, the script will have to be changed

fid= fopen(filename,'r'); %open the file with read access only

info=textscan(fid,'%s',40,'Delimiter','\b'); %import the first 40 lines, these should contain info on the data on lat/long/time and col headers
info=info{:};


%locate system upload time:
ind=strfind(info,'System UpLoad Time =');
ii= cellfun('isempty',ind)==0;
test=info{ii};
[a, b] = regexp(test,'System UpLoad Time = ');
upload_time=test(b+1:end);
upload_time=datenum(upload_time);

%locate UTC time:
ind=strfind(info,'UTC');
ii= find(cellfun('isempty',ind)==0);
if ~isempty(ii)
    test=info{ii};
    [a, b] = regexp(test,' = ');
    UTC_time=test(b+1:end);
    if ~isempty(regexp(UTC_time,'[A-Z][a-z]{2}')) %meaning that the value captured is an acutal time (and not a lat or lon - as sometimes happens)
        UTC_time=datenum(UTC_time);
    else
        UTC_time=[];
    end
else
    UTC_time=[];
end

%locate latitude:
ind=strfind(info,'Latitude');
ii=find(cellfun('isempty',ind)==0);
if ~isempty(ii)
    test=info{ii};
    lat=regexp(test,'(?<deg>\d{2,3})(\s)(?<min>\d{1,2})\.(?<sec>\d{2})','names');
else
    lat=[];
end

%locate longitude:
ind=strfind(info,'Longitude');
ii= find(cellfun('isempty',ind)==0);
if ~isempty(ii)
    test=info{ii};
    lon=regexp(test,'(?<deg>\d{2,3})(\s)(?<min>\d{1,2})\.(?<sec>\d{2})','names');
else
    lon=[];
end

%a double check for certain casts where the longitude is mixed up with the
%latutitude and lat mixed up with time:
if ~isempty(lat) & isempty(lon)
if ~isempty(strfind(lat.deg,'70'))
    
    lon=lat;
    jj=strfind(info,'NMEA');
    j2=find(cellfun('isempty',jj)==0);
    
     %now find latitude and time:  
    for w=1:length(j2)
        test=regexp(info{j2(w)},'(?<deg>\d{2,3})(\s)(?<min>\d{1,2})\.(?<sec>\d{2})\sN','names');
        test2=regexp(info{j2(w)},'[A-Z]{1}[a-z]{2}(\s)*\d{2}(\s*)\d{4}');

        if ~isempty(test)
            lat=test;
        end
        
        if ~isempty(test2)
            UTC_time=info{j2(w)}(test2:end);
            UTC_time=datenum(UTC_time);
        end
            
    end
end
end
%locate header names:
ind=strfind(info,'# name');
ii= cellfun('isempty',ind)==0;
header=info(ii);
for j=1:length(header)
    temp=header{j};
    temp2=regexprep(temp,'#\sname\s\d{1,2}\s=\s','');
    header{j}=temp2;
end

%stop at the position before the data:
t=textscan(fid,'%s',1,'Delimiter','\b','CommentStyle','#'); % find the '*END*' before the data

%gather as many data columns as there are headers:
FormatSpec=repmat('%f',1,length(header));
data=textscan(fid,FormatSpec,'Delimiter','\b','CommentStyle','#');

fclose(fid); %close the file!


%% if you know the exact format of your .cnv file:
%
% info=textscan(fid,'%s',19,'Delimiter','\b'); %scan first top for lat/lon/time:
% info=info{:};
% ind=strfind(info,'NMEA'); %this info will have a 'NMEA' out in front
% ii=find(cellfun('isempty',ind)==0);
%
% temp=info{ii(1)};
% lat=regexp(temp,'(?<=Latitude\s=\s)\w*\s*\w*.\w*','match');
%
% temp=info{ii(2)};
% lon=regexp(temp,'(?<=Longitude\s=\s)\w*\s*\w*.\w*','match');
%
% temp=info{ii(3)};
% time=regexp(temp,'(?<=\s=\s)\w{3}\s*\d{2}\s*\d{4}\s*\d{2}:\d{2}:\d{2}','match');
% time=datenum(time{:});
%
% %%
% datacols=textscan(fid,'# name%s',10,'Delimiter','\b');
% datacols=datacols{:};
%
% t=textscan(fid,'%s',1,'Delimiter','\b','CommentStyle','#'); % an '*END*' before the data
% data=textscan(fid,'%f%f%f%f%f%f%f%f%f%f','Delimiter','\b','CommentStyle','#');