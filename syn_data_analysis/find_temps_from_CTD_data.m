% find corresponding temperature for days in a dataset:

%load in a datelist! or check MVCO samples!

daylist=datenum({'09-Aug-2016';		%MV366
	'19-Aug-2016';	 %MV367
	'24-Aug-2016'}); %MV368

load /Users/kristenhunter-cevera/MVCO_light_at_depth/mixed_layer_depth_src/QC_downcast.mat
mvco_ind=find(downcast_lat > 41.3 & downcast_lat < 41.35 & downcast_lon < -70.53 & downcast_lon > -70.60 & cellfun('isempty',{CTD_QC(:).flag})==1);
time_mvco=cell2mat({CTD_QC(mvco_ind).upload_time}');
time=cell2mat({CTD_QC.upload_time}');

%%

for j=1:lemght(daylist)
    jj=find(floor(time)==daylist(j))

end