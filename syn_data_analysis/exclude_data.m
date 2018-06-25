function to_use=exclude_data(cellresultsall,year2do)
%data to exclude for now (and ideally fix later)

%Issues include:
%beads classified as Syn
%merging issues due to baseline wobble
%just bad signatures
%strange acquisiiton times - not sure what exactly these are...
%syringe clogged - bad volume estimates

switch year2do
    case 2003
        
        ex_timepts={};
        
        ii=find(cellresultsall(:,1)>=datenum('5-10-03 21:00:00') & cellresultsall(:,1)<=datenum('5-10-03 23:00:00'));
        reason=repmat('Beads classified as Syn',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];
        
        i1=find(cellresultsall(:,1)>=datenum('5-17-03 03:45:00') & cellresultsall(:,1)<=datenum('5-17-03 04:15:00')); %really maybe only from 3:45-4:15?
        i2=find(cellresultsall(:,1)>=datenum('5-18-03 06:45:00') & cellresultsall(:,1)<=datenum('5-17-03 07:20:00')); %really maybe only from 3:45-4:15?
        i3=find(cellresultsall(:,1)>=datenum('5-18-03 11:40:00') & cellresultsall(:,1)<=datenum('5-17-03 12:00:00')); %really maybe only from 3:45-4:15?
        ii=[i1;i2;i3];
        reason=repmat('Too fast acquisition times?',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];
        
        ii=find(cellresultsall(:,1)>=datenum('06-08-03 22:00:00') & cellresultsall(:,1)<=datenum('06-08-03 23:00:00'));
        reason=repmat('Syn classified as picoeuks!',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];
        
        %         ii=cell2mat(ex_timepts(:,1));
        %         [~, ia, ib] =intersect(ii,cellresultsall(:,1));
        [bt_pts]=ismember(cellresultsall(:,1),cell2mat(ex_timepts(:,1))); %returns 0 if not in the time points to exclude!
        to_use=find(bt_pts==0);
    
    case 2004
        %Nothing! :)
        to_use=find(~isnan(cellresultsall(:,1)));
        
    case 2005
        
        ex_timepts={};
        
        ii=find(cellresultsall(:,1)>=datenum('4-20-05 01:00:00') & cellresultsall(:,1)<=datenum('4-20-05 03:00:00'));
        reason=repmat('Beads classified as Syn',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];
        
        ii=find(cellresultsall(:,1)>=datenum('10-24-05 23:00:00') & cellresultsall(:,1)<=datenum('10-25-05 00:00:00'));
        reason=repmat('Bad time file split?',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];
        
        ii=find(cellresultsall(:,1)>=datenum('10-30-05 23:00:00') & cellresultsall(:,1)<=datenum('10-31-05 00:00:00'));
        reason=repmat('Too long acquisition ties?',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];
        
        %         ii=cell2mat(ex_timepts(:,1));
        %         [~, ia, ib] =intersect(ii,cellresultsall(:,1));
        [bt_pts]=ismember(cellresultsall(:,1),cell2mat(ex_timepts(:,1))); %returns 0 if not in the time points to exclude!
        to_use=find(bt_pts==0);
        
    case  2006
        
        ex_timepts={};
        
        ii=find(cellresultsall(:,1)>=datenum('6-20-06 14:00:00') & cellresultsall(:,1)<=datenum('6-20-06 18:00:00'));
        reason=repmat('Too long acquisition times?',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];
        %
        %         ii=cell2mat(ex_timepts(:,1));
        %         [~, ia, ib] =intersect(ii,cellresultsall(:,1));
        
        [bt_pts]=ismember(cellresultsall(:,1),cell2mat(ex_timepts(:,1))); %returns 0 if not in the time points to exclude!
        to_use=find(bt_pts==0);
        
    case  2007
        ex_timepts={};
        
        ii=find(cellresultsall(:,1)>=datenum('5-19-07 11:00:00') & cellresultsall(:,1)<=datenum('5-19-07 14:00:00'));
        reason=repmat('Noise classified as Syn',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];
        
        ii=find(cellresultsall(:,1)>=datenum('12-13-07 17:00:00') & cellresultsall(:,1)<=datenum('12-13-07 18:00:00'));
        reason=repmat('Missed Syn patch',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];
      
        %         ii=cell2mat(ex_timepts(:,1));
        %         [~, ia, ib] =intersect(ii,cellresultsall(:,1));
        
        [bt_pts]=ismember(cellresultsall(:,1),cell2mat(ex_timepts(:,1))); %returns 0 if not in the time points to exclude!
        to_use=find(bt_pts==0);
        
    case 2008
        ex_timepts={};
        
        ii=find(cellresultsall(:,1)>=datenum('1-11-08 14:00:00') & cellresultsall(:,1)<=datenum('1-11-08 15:05:00'));
        reason=repmat('Strange acquisition times?',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];
        
        %         ii=cell2mat(ex_timepts(:,1));
        %         [~, ia, ib] =intersect(ii,cellresultsall(:,1));
        
        [bt_pts]=ismember(cellresultsall(:,1),cell2mat(ex_timepts(:,1))); %returns 0 if not in the time points to exclude!
        to_use=find(bt_pts==0);
        
    case 2009
        ex_timepts={};
        
        ii=find(cellresultsall(:,1)>=datenum('2-10-09 20:00:00') & cellresultsall(:,1)<=datenum('2-10-09 21:00:00'));
        reason=repmat('Strange acquisition times?',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];
        
        ii=find(cellresultsall(:,1)>=datenum('8-14-09 14:00:00') & cellresultsall(:,1)<=datenum('8-14-09 16:00:00'));
        reason=repmat('Noise classified as Syn',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];
        
        i1=find(cellresultsall(:,1)>=datenum('11-08-09 05:00:00') & cellresultsall(:,1)<=datenum('11-08-09 06:00:00'));
        i2=find(cellresultsall(:,1)>=datenum('11-08-09 20:00:00') & cellresultsall(:,1)<=datenum('11-09-09 08:00:00'));
        i3=find(cellresultsall(:,1)>=datenum('11-13-09 12:00:00') & cellresultsall(:,1)<=datenum('11-13-09 15:00:00')); %really 14:45
        i4=find(cellresultsall(:,1)>=datenum('11-14-09  03:00:00') & cellresultsall(:,1)<=datenum('11-14-09 06:00:00'));
        i5=find(cellresultsall(:,1)>=datenum('11-14-09  22:00:00') & cellresultsall(:,1)<=datenum('11-15-09 00:00:00'));
        i6=find(cellresultsall(:,1)>=datenum('11-15-09  02:00:00') & cellresultsall(:,1)<=datenum('11-15-09 03:00:00'));
        i7=find(cellresultsall(:,1)>=datenum('11-15-09  05:00:00') & cellresultsall(:,1)<=datenum('11-15-09 08:00:00'));
        i8=find(cellresultsall(:,1)>=datenum('11-15-09  19:00:00') & cellresultsall(:,1)<=datenum('11-15-09 21:00:00'));
        i9=find(cellresultsall(:,1)>=datenum('11-16-09  12:00:00') & cellresultsall(:,1)<=datenum('11-16-09 14:00:00'));
        i10=find(cellresultsall(:,1)>=datenum('11-17-09  22:00:00') & cellresultsall(:,1)<=datenum('11-18-09 00:00:00'));
        i11=find(cellresultsall(:,1)>=datenum('11-19-09  23:00:00') & cellresultsall(:,1)<=datenum('11-20-09 02:00:00'));
        ii=[i1;i2;i3;i4;i5;i6;i7;i8;9;i10;i11];
        reason=repmat('Slow acquisition times - inlet clog/sheath fill?',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];
        %
        %         ii=cell2mat(ex_timepts(:,1));
        %         [~, ia, ib] =intersect(ii,cellresultsall(:,1));
        [bt_pts]=ismember(cellresultsall(:,1),cell2mat(ex_timepts(:,1))); %returns 0 if not in the time points to exclude!
        to_use=find(bt_pts==0);
        
    case 2010
        
        %Possibly some dying syn around Dec 5-6, where distinct lower PE clusters
        %appear...
        ex_timepts={};
        i1=find(cellresultsall(:,1)>=datenum('8-21-10 4:45:00') & cellresultsall(:,1)<=datenum('8-21-10 06:00:00'));
        i2=find(cellresultsall(:,1)>=datenum('12-17-10 05:00:00') & cellresultsall(:,1)<=datenum('12-17-10 06:00:00'));
        %i3=find(cellresultsall(:,1)>=datenum('12-13-10 13:30:00') & cellresultsall(:,1)<=datenum('12-13-10 14:30:00')); %longer acq times - fix later!
        ii=[i1;i2];
        reason=repmat('Strange acquisition times?',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];
        
        ii=find(cellresultsall(:,1)>=datenum('8-22-10 6:00:00') & cellresultsall(:,1)<=datenum('8-22-10 07:00:00'));
        reason=repmat('Syn classified as picoeuks',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];
        
        i1=find(cellresultsall(:,1)>=datenum('2-13-10 17:00:00') & cellresultsall(:,1)<=datenum('2-13-10 19:00:00'));
        i2=find(cellresultsall(:,1)>=datenum('12-16-10 10:00:00') & cellresultsall(:,1)<=datenum('12-16-10 11:00:00'));
        ii=[i1;i2];
        reason=repmat('Strange data pattern - no detectable Syn?',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];   
     
        i1=find(cellresultsall(:,1)>=datenum('6-22-10 04:30:00') & cellresultsall(:,1)<=datenum('6-22-10 18:30:00'));
        i2=find(cellresultsall(:,1)>=datenum('6-22-10 20:30:00') & cellresultsall(:,1)<=datenum('6-23-10 02:00:00'));
        ii=[i1;i2];
        reason=repmat('Syring pump clogged and sucking sheath?',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];
        
        %         ii=cell2mat(ex_timepts(:,1));
        %         [~, ia, ib] =intersect(ii,cellresultsall(:,1));
        [bt_pts]=ismember(cellresultsall(:,1),cell2mat(ex_timepts(:,1))); %returns 0 if not in the time points to exclude!
        to_use=find(bt_pts==0);
        
    case 2011
        
        ex_timepts={};
        
        ii=find(floor(cellresultsall(:,1))==datenum('7-12-11 ')); ii=ii(3:end); %first two look okay
        reason=repmat('Merging issues',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];
        
%         ii=find(cellresultsall(:,1)>=datenum('7-31-11 04:30:00') & cellresultsall(:,1)<=datenum('7-31-11 6:10:00'));
%         reason=repmat('So few cells - clogged?',length(ii),1);
%         ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];
%        
        % maybe...
        % ii=find(floor(cellresultsall(:,1))==datenum('7-13-11'));
        % reason=repmat('Merging issues',length(ii),1);
        % ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];
        
        %maybe: August 3rd -August 6th - super smeared Syn on SSC - SSC noise being
        %captured
        
        i1=find(cellresultsall(:,1)>=datenum('8-06-11 21:00:00') & cellresultsall(:,1)<=datenum('8-07-11 4:00:00'));
        i2=find(cellresultsall(:,1)>=datenum('8-07-11 9:00:00') & cellresultsall(:,1)<=datenum('8-07-11 12:00:00'));
        i3=find(cellresultsall(:,1)>=datenum('8-10-11 3:00:00') & cellresultsall(:,1)<=datenum('8-10-11 4:00:00'));
        i4=find(cellresultsall(:,1)>=datenum('8-11-11 4:00:00') & cellresultsall(:,1)<=datenum('8-11-11 5:00:00'));
        ii=[i1;i2;i3;i4];
        reason=repmat('Noise classified as Syn',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];
        
        ii=find(cellresultsall(:,1)>=datenum('10-19-11 06:00:00') & cellresultsall(:,1)<=datenum('10-19-11 7:00:00'));
        reason=repmat('Syn patch classified as picoeuks',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];
        
        i1=find(cellresultsall(:,1)>=datenum('12-15-11 10:00:00') & cellresultsall(:,1)<=datenum('12-15-11 12:00:00'));
        i2=find(cellresultsall(:,1)>=datenum('12-15-11 12:00:00') & cellresultsall(:,1)<=datenum('12-15-11 14:00:00'));
        i3=find(cellresultsall(:,1)>=datenum('12-15-11 15:00:00') & cellresultsall(:,1)<=datenum('12-15-11 18:00:00'));
        i4=find(cellresultsall(:,1)>=datenum('12-15-11 20:00:00') & cellresultsall(:,1)<=datenum('12-15-11 22:00:00'));
        i5=find(cellresultsall(:,1)>=datenum('12-16-11 01:30:00') & cellresultsall(:,1)<=datenum('12-31-11 11:59:00'));
        ii=[i1;i2;i3;i4;i5];
        reason=repmat('Bad SSC data',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];
        
        [bt_pts]=ismember(cellresultsall(:,1),cell2mat(ex_timepts(:,1))); %returns 0 if not in the time points to exclude!
        to_use=find(bt_pts==0);
        
    case 2012
        
        ex_timepts={};
        ii=find(cellresultsall(:,1)>=datenum('Jan-09-12 20:00:00') & cellresultsall(:,1)<=datenum('Jan-09-12 22:00:00')); %start until 22:00 (Jan-9-2012)
        reason=repmat('Strang acquisition time pattern?',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];
        
        %         ii=cell2mat(ex_timepts(:,1));
        %         [~, ia, ib] =intersect(ii,cellresultsall(:,1));
        
        [bt_pts]=ismember(cellresultsall(:,1),cell2mat(ex_timepts(:,1))); %returns 0 if not in the time points to exclude!
        to_use=find(bt_pts==0);
        
    case 2013
        
        ex_timepts={};
        
        %note- there are more instances where classfier is grabbing noise, these
        %just looked the most severe:
        i1=find(cellresultsall(:,1)>=datenum('3-27-13 19:00:00') & cellresultsall(:,1)<=datenum('3-27-13 22:00:00'));
        i2=find(cellresultsall(:,1)>=datenum('4-3-13 20:00:00') & cellresultsall(:,1)<=datenum('4-3-13 22:00:00'));
        i3=find(cellresultsall(:,1)>=datenum('4-4-13 09:00:00') & cellresultsall(:,1)<=datenum('4-4-13 11:00:00'));
        i4=find(cellresultsall(:,1)>=datenum('4-5-13 09:00:00') & cellresultsall(:,1)<=datenum('4-5-13 11:00:00'));
        i5=find(cellresultsall(:,1)>=datenum('4-8-13 01:00:00') & cellresultsall(:,1)<=datenum('4-8-13 03:00:00'));
        i6=find(cellresultsall(:,1)>=datenum('4-13-13 07:00:00') & cellresultsall(:,1)<=datenum('4-13-13 08:00:00'));
        i7=find(cellresultsall(:,1)>=datenum('4-17-13 15:00:00') & cellresultsall(:,1)<=datenum('4-17-13 16:00:00'));
        i8=find(cellresultsall(:,1)>=datenum('4-27-13 16:00:00') & cellresultsall(:,1)<=datenum('4-27-13 17:00:00'));
        i9=find(cellresultsall(:,1)>=datenum('4-28-13 09:00:00') & cellresultsall(:,1)<=datenum('4-28-13 10:00:00'));
        i10=find(cellresultsall(:,1)>=datenum('5-18-13 17:00:00') & cellresultsall(:,1)<=datenum('5-18-13 19:00:00'));
        i11=find(cellresultsall(:,1)>=datenum('5-19-13 11:00:00') & cellresultsall(:,1)<=datenum('5-19-13 13:00:00'));
        ii=[i1;i2;i3;i4;i5;i6;i7;i8;i9;i10;i11];
        reason=repmat('Noise classified as Syn',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];
        
        i1=find(cellresultsall(:,1)>=datenum('5-19-13 04:00:00') & cellresultsall(:,1)<=datenum('5-19-13 5:00:00'));
        i2=find(cellresultsall(:,1)>=datenum('5-21-13 17:00:00') & cellresultsall(:,1)<=datenum('5-21-13 18:00:00'));
        i3=find(cellresultsall(:,1)>=datenum('11-21-13 18:00:00') & cellresultsall(:,1)<=datenum('11-21-13 19:00:00'));
        ii=[i1;i2;i3];
        reason=repmat('Strange data pattern - no detectable Syn?',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];
        
        ii=find(cellresultsall(:,1)>=datenum('5-27-13 00:00:00') & cellresultsall(:,1)<=datenum('5-27-13 01:00:00'));
        reason=repmat('Part of Syn patch classified as picoeuks',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];
        
        i1=find(cellresultsall(:,1)>=datenum('5-19-13 02:00:00') & cellresultsall(:,1)<=datenum('5-19-13 06:00:00'));
        i2=find(cellresultsall(:,1)>=datenum('6-21-13 12:00:00') & cellresultsall(:,1)<=datenum('6-21-13 20:00:00'));
        i3=find(cellresultsall(:,1)>=datenum('8-08-13 12:00:00') & cellresultsall(:,1)<=datenum('08-08-13 14:00:00'));
        i4=find(cellresultsall(:,1)>=datenum('11-23-13 12:00:00') & cellresultsall(:,1)<=datenum('11-23-13 13:00:00'));
        i5=find(cellresultsall(:,1)>=datenum('8-02-13 14:00:00') & cellresultsall(:,1)<=datenum('08-02-13 21:00:00'));
        ii=[i1;i2;i3;i4;i5];
        reason=repmat('Something not quite right with timing or volume?',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];
        
        %         ii=cell2mat(ex_timepts(:,1));
        %         [~, ia, ib] =intersect(ii,cellresultsall(:,1));
        %
        [bt_pts]=ismember(cellresultsall(:,1),cell2mat(ex_timepts(:,1))); %returns 0 if not in the time points to exclude!
        to_use=find(bt_pts==0);
        
    case 2014
        ex_timepts={};
        
        i1=find(cellresultsall(:,1)>=datenum('1-06-14 17:00:00') & cellresultsall(:,1)<=datenum('1-06-14 20:00:00'));
        i2=find(cellresultsall(:,1)>=datenum('1-10-14 06:00:00') & cellresultsall(:,1)<=datenum('1-10-14 07:00:00'));
        i3=find(cellresultsall(:,1)>=datenum('1-10-14 09:00:00') & cellresultsall(:,1)<=datenum('1-10-14 12:00:00'));
        i4=find(cellresultsall(:,1)>=datenum('5-09-14 10:00:00') & cellresultsall(:,1)<=datenum('5-09-14 11:00:00'));
        ii=[i1;i2;i3;i4];
        reason=repmat('Noise classified as Syn',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];
        
        i1=find(cellresultsall(:,1)>=datenum('1-08-14 01:00:00') & cellresultsall(:,1)<=datenum('1-08-14 06:00:00'));
        i2=find(cellresultsall(:,1)>=datenum('1-08-14 14:00:00') & cellresultsall(:,1)<=datenum('1-08-14 15:00:00'));
        i3=find(cellresultsall(:,1)>=datenum('1-14-14 13:00:00') & cellresultsall(:,1)<=datenum('1-14-14 15:00:00'));
        ii=[i1;i2;i3];
        reason=repmat('Bad signatures/data?',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];
        %
        ii=find(cellresultsall(:,1)>=datenum('Nov-11-14 09:00:00') & cellresultsall(:,1)<=datenum('Nov-11-14 10:00:00')); %after 9:30 - bad acq?
        reason=repmat('Bad acquisition times?',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];
        
        i1=find(cellresultsall(:,1)>=datenum('Apr-03-14 00:00:00') & cellresultsall(:,1)<=datenum('Apr-03-14 01:00:00'));
        i2=find(cellresultsall(:,1)>=datenum('Apr-04-14 23:00:00') & cellresultsall(:,1)<=datenum('Apr-05-14 00:00:00'));
        i3=find(cellresultsall(:,1)>=datenum('Apr-05-14 21:00:00') & cellresultsall(:,1)<=datenum('Apr-05-14 22:00:00'));
        i4=find(cellresultsall(:,1)>=datenum('Apr-07-14 22:00:00') & cellresultsall(:,1)<=datenum('Apr-07-14 23:00:00'));
        i5=find(cellresultsall(:,1)>=datenum('Apr-30-14 20:00:00') & cellresultsall(:,1)<=datenum('Apr-30-14 21:00:00'));
        ii=[i1;i2;i3;i4;i5];
        reason=repmat('Cryptos classified as Syn',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];
        
        ii=find(cellresultsall(:,1)>=datenum('Nov-10-14 19:00:00') & cellresultsall(:,1)<=datenum('Nov-10-14 20:00:00')); %from 19:00 - 20:30 - syr pump doesn't reach max
        reason=repmat('Clogged syringe?',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];
        
        %         ii=cell2mat(ex_timepts(:,1));
        %         [~, ia, ib] =intersect(ii,cellresultsall(:,1));
        %
        [bt_pts]=ismember(cellresultsall(:,1),cell2mat(ex_timepts(:,1))); %returns 0 if not in the time points to exclude!
        to_use=find(bt_pts==0);
        
    case 2015
        
        ex_timepts={};
        ii=find(cellresultsall(:,1)>=datenum('3-25-15 15:00:00') & cellresultsall(:,1)<=datenum('3-25-15 16:00:00'));
        reason=repmat('Beads classified as Syn',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];
        
        i1=find(cellresultsall(:,1)>=datenum('3-30-15 07:00:00') & cellresultsall(:,1)<=datenum('3-30-15 08:00:00'));
        i2=find(cellresultsall(:,1)>=datenum('6-08-15 07:00:00') & cellresultsall(:,1)<=datenum('6-08-15 08:00:00'));
        i3=find(cellresultsall(:,1)>=datenum('6-11-15 03:00:00') & cellresultsall(:,1)<=datenum('6-11-15 04:00:00'));
        i4=find(cellresultsall(:,1)>=datenum('6-11-15 09:00:00') & cellresultsall(:,1)<=datenum('6-11-15 11:00:00'));
        i5=find(cellresultsall(:,1)>=datenum('6-12-15 11:00:00') & cellresultsall(:,1)<=datenum('6-12-15 13:00:00'));
        i6=find(cellresultsall(:,1)>=datenum('6-12-15 14:00:00') & cellresultsall(:,1)<=datenum('6-12-15 15:00:00'));
        i7=find(cellresultsall(:,1)>=datenum('6-14-15 04:00:00') & cellresultsall(:,1)<=datenum('6-14-15 05:00:00'));
        ii=[i1;i2;i3;i4;i5;i6;i7];
        reason=repmat('Noise classified as Syn',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];
        
        ii=find(cellresultsall(:,1)>=datenum('5-17-15 01:00:00') & cellresultsall(:,1)<=datenum('5-17-15 02:00:00'));
        reason=repmat('Strange data patterns',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];
        
        i1=find(cellresultsall(:,1)>=datenum('11-18-15 19:00:00') & cellresultsall(:,1)<=datenum('11-18-15 20:00:00'));
        i2=find(cellresultsall(:,1)>=datenum('11-28-15 19:00:00') & cellresultsall(:,1)<=datenum('11-28-15 20:00:00'));
        ii=[i1;i2];
        reason=repmat('Strange acquisition times?',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];
        
        %         ii=cell2mat(ex_timepts(:,1));
        %         [~, ia, ib] =intersect(ii,cellresultsall(:,1));
        %
        [bt_pts]=ismember(cellresultsall(:,1),cell2mat(ex_timepts(:,1))); %returns 0 if not in the time points to exclude!
        to_use=find(bt_pts==0);
        
    case 2016
        
        ex_timepts={};
        
        %Jan 14-15, most bad data removed -some left in good files, remove:
        i1=find(cellresultsall(:,1)>=datenum('1-14-16 22:00:00') & cellresultsall(:,1)<=datenum('1-14-16 23:00:00'));
        i2=find(cellresultsall(:,1)>=datenum('1-15-16 16:00:00') & cellresultsall(:,1)<=datenum('1-15-16 17:50:00')); %really only up unitl '1-15-16 17:30:00'
        ii=[i1;i2];
        reason=repmat('Strange acquisition times?',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];
        
        ii=find(cellresultsall(:,1)>=datenum('Apr-5-16 11:00:00') & cellresultsall(:,1)<=datenum('Apr-5-16 12:00:00'));
        reason=repmat('Noise classified as Syn',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];
        
        
        i1=find(cellresultsall(:,1)>=datenum('May-9-16 04:00:00') & cellresultsall(:,1)<=datenum('May-9-16 06:00:00'));
        i2=find(cellresultsall(:,1)>=datenum('5-31-16 03:00:00') & cellresultsall(:,1)<=datenum('5-31-16 04:00:00')); %this one is probably tail end of bad data....
        i3=find(cellresultsall(:,1)>=datenum('6-15-16 12:00:00') & cellresultsall(:,1)<=datenum('6-19-16 05:00:00')); %mess of signatures - THERE MAYBE STILL GOOD FILES IN HERE THOUGH!

%         i3=find(cellresultsall(:,1)>=datenum('6-15-16 12:00:00') & cellresultsall(:,1)<=datenum('6-16-16 01:00:00')); %mess of signatures 
%         i4=find(cellresultsall(:,1)>=datenum('6-16-16 14:00:00') & cellresultsall(:,1)<=datenum('6-16-16 23:00:00')); %mess of signatures 
%         i5=find(cellresultsall(:,1)>=datenum('6-17-16 11:00:00') & cellresultsall(:,1)<=datenum('6-19-16 05:00:00')); %mess of signatures 

        ii=[i1;i2;i3];
        reason=repmat('Bad signatures',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];
        
        ii=find(cellresultsall(:,1)>=datenum('10-03-16 22:00:00') & cellresultsall(:,1)<=datenum('10-04-16 17:00:00')); 
        reason=repmat('SSC spread - bad signatures',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];
        
        ii=find(cellresultsall(:,1)>=datenum('10-11-16 05:00:00') & cellresultsall(:,1)<=datenum('10-11-16 15:00:00')); 
        reason=repmat('Bad signatures - spread on both PE and SSC',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];
       
        ii=find(cellresultsall(:,1)>=datenum('11-08-16 07:00:00') & cellresultsall(:,1)<=datenum('11-08-16 08:00:00')); 
        reason=repmat('Bottom end of Syn cloud classified as picoeuk',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];
 
        ii=find(cellresultsall(:,1)>=datenum('12-29-16 20:00:00') & cellresultsall(:,1)<=datenum('12-30-16 10:00:00')); 
        reason=repmat('Small SSC noise classified as Syn',length(ii),1);
        ex_timepts=[ex_timepts; num2cell(cellresultsall(ii,1)) cellstr(reason)];
 
        [bt_pts]=ismember(cellresultsall(:,1),cell2mat(ex_timepts(:,1))); %returns 0 if not in the time points to exclude!
        to_use=find(bt_pts==0);
        
end
