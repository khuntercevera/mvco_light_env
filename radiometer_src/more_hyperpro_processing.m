%% to be incorporated in anotherscript...

%what we want to do in this script is for each "set" of data,
%see what the casts look like over the day
%calculate attentuation coefficient for each .raw file
%see what k looks like over the casts

%other fancy plots at the end :)

% first find all the folders with the processed matlab files...

sourcepath=fullfile('\\sosiknas1\lab_data\mvco\HyperPro_Radiometer\processed_radiomater_files\'); % path to folders with raw data...

%How many of these do we have?
d = dir(sourcepath);
isub = [d(:).isdir]; %# returns logical vector
foldernames = {d(isub).name}'; %names of folders
temp=regexp(foldernames,'\d{1,2}\w{3,4}20\d{2}'); %find all folders with this naming scheme
datafolders = foldernames(cellfun('isempty',temp)==0);


%can just manually enter in a foldernum for now .... put into a for loop later....
foldernum=2;

%find dat data!
matsource=fullfile(sourcepath,datafolders{foldernum});
sourcefiles=dir([tempdatasource '/mat_outfiles/*.mat']);
   
%where to save:
savepath=fullfile(sourcepath,datafolders{foldernum},'mat_outfiles'); %this directory should already exist!

%now go through each file:
%check the log to see if good casts are present
%if so, go ahead and calculate k
%save data in a structure format?

numfiles=length(sourcefiles);
for filenum=1:numfiles

    filename=sourcefiles(filenum).name;
    %check if good casts?
    ii=strcmp(filename(1:end-4),raw_data_log(:,2)); %find the file
    
    if isempty(raw_data_log(ii,4)) | strcmp(raw_data_log(ii,4),'not a cast') 
        disp(['Empty file or file with no valid casts - skipping: ' filename])
    else
        disp(['Processing: ' filename])
        
        
        
        
    

    end

end 

 % fill in info about cast....
            
            %FIND NUMBER OF CASTS IN FILE
            %--------------------------------------------------------------------------------------------------
            % split casts?
            
            %Some definitions:
            %full cast - one downward fall followed by one upward pull
            %half cast - either downward fall, no return or return but no fall...
            
            %find casts by seeing how many times there is a 'sign' switch from
            %an arbitrary depth (so as to avoid any bobbing at the surface or
            %depth:
            
            Q=quantile(depth,[0 .25 0.5 0.75 1]);
            dline=Q(2);
%             subplot(2,3,3,'replace')
%             plot([mpr_data{end-1}(1) mpr_data{end-1}(end)],[dline dline])
%             
            testd=depth-dline;
            testdB=zeros(length(depth),1);
            
            testdB(testd>1)=1;
            ii=find(diff(testdB~=0));
            
            if ~isempty(ii) & length(ii)>1
            zdir=testdB(ii) - testdB(ii+1);
            zz=find(zdir==-1); %profiler is falling
            
            
            fullcasts=0;
            halfcasts=zz(1)-1; %0 if this is the first fall, one if coming back up
            k=1;
            while k <= length(zz)
                fall=zz(k);
                if zz(k) == length(zdir) %meaning that no up pull was found
                    halfcasts=halfcasts+1;
                elseif zdir(fall+1) == 1 %profiler is being pulled up after being let down
                    fullcasts=fullcasts+1;
                end
                k=k+1;
            end
            
            fprintf('I think I found %1.0f full casts and %1.0f half-casts\n',fullcasts,halfcasts)
            castflag=0;
            else
                disp('Could not find any casts...')
                castflag=1;
            end
            %okay, so once we now what/how many casts we have - find the
            %downward casts only:
            
            %smooth depths
            % find all sign changes of difference in position
            % find sign changes that correspond to depth decrease
            
            if castflag==0
                smooth_depth=smooth2(tdepth,min(length(tdepth),15));
                sign_motion=sign(diff(smooth_depth));
                sign_motion(isnan(sign_motion))=0;
                zz=find(diff(sign_motion)~=0); %these should be all the changes vertical direction
            else
                zz=[]
            end
            %find the closest indexes to the max depth:
            %[~, bottom]=min(abs(smooth_depth(zz)-max(depth)));
            
            if ~isempty(zz)
                [~, is]=sort(abs(smooth_depth(zz)-max(depth)));
                
                bottoms=is(1:fullcasts);
                downcast_ind=[];
                for q=1:fullcasts
                    if bottoms(q)==1
                        downcast_ind=[downcast_ind 1:zz(bottoms(q))];
                    else
                        downcast_ind=[downcast_ind zz(bottoms(q)-1):zz(bottoms(q))];
                    end
                end
                
                downcast_ind=sort(downcast_ind); %because some bottoms may be found out of order???
                
                % find attenuation coefficient!
                % We are really only interested in the coefficient at 4m depth
                % (they can change with depth, and for each wavelength, and that's what the prosoft
                % does, but for correlations with other instruments, maybe just
                % one k value is the best??
                
                %And I guess we would really only be interested in the
                %downcast to calculate this...
                
                %             %Attenuation coefficient over the whole depth (assume one
                %             %value) with PAR:
                %             par_depths=depth(mpr2edl);
                %             mindepth=min(par_depths); jj=find(par_depths > mindepth+2);
                %             maxdepth=max(par_depths); ii=find(par_depths < maxdepth-2);
                %             qq=intersect(ii,jj);
                %
                %             %remove the first/last 2 meters, I think (this is a rough
                %             %guide, certain casts may not follow this:
                %             [B,~,~,~,STATS] = regress(log(edl_PAR(qq)),[ones(size(par_depths(qq))) par_depths(qq)]);
                %           %sync to light data time
                
                downcast_time_window = [mpr_data{end-1}(downcast_ind(1)) mpr_data{end-1}(downcast_ind(end))];
                [~, startcast]=min(abs(edl_data{end-1}-downcast_time_window(1)));
                [~, returncast]=min(abs(edl_data{end-1}-downcast_time_window(end)));
                par_depths=tdepth(mpr2edl);
                [K,~,~,~,STATS] = regress(log(edl_PAR(startcast:returncast)),[ones(size(par_depths(startcast:returncast))) par_depths(startcast:returncast)]);
                
                
                if STATS(1) < 0.9
                    disp('Hmmm...maybe not such a great regression fit to find k?')
                end
            end