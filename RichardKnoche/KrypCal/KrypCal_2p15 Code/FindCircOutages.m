function [ MFC1_time MFC1_date noflow_times ] = FindCircOutages()
%Finds circulation outages between Sep 1 2014, and Jan 1 2017
    
%Get MFC1 Settings
% disp('Starting MFC1 Data Collecting')
% start_time=(datenum(2014,09,01,0,0,0)-datenum(1970,1,1))*86400;
% end_time=(datenum(2014,10,01,0,0,0)-datenum(1970,1,1))*86400;
% querystr=sprintf('select * from control.sc_sens_ACRS_Flow where time > %d and time < %d order by time asc',start_time,end_time);
% query_result=MySQLQuery(querystr);
% MFC1_time=query_result.time;
% MFC1_value=query_result.value;
% sprintf('Finished Month of %d %d',09,2014)
% clear query_result
MFC1_time=[];
MFC1_value=[];
for i=9:11;
start_time=(datenum(2014,i,01,0,0,0)-datenum(1970,1,1))*86400;
end_time=(datenum(2014,i+1,01,0,0,0)-datenum(1970,1,1))*86400;
querystr=sprintf('select * from control.sc_sens_ACRS_Flow where time > %d and time < %d order by time asc',start_time,end_time);
query_result=MySQLQuery(querystr);
MFC1_time=[MFC1_time;query_result.time];
MFC1_value=[MFC1_value;query_result.value];
sprintf('Finished Month of %d %d',i,2014)
clear query_result
end
start_time=(datenum(2014,12,01,0,0,0)-datenum(1970,1,1))*86400;
end_time=(datenum(2015,1,01,0,0,0)-datenum(1970,1,1))*86400;
querystr=sprintf('select * from control.sc_sens_ACRS_Flow where time > %d and time < %d order by time asc',start_time,end_time);
query_result=MySQLQuery(querystr);
MFC1_time=[MFC1_time;query_result.time];
MFC1_value=[MFC1_value;query_result.value];
sprintf('Finished Month of %d %d',12, 2014)
clear query_result
for i=1:11;
start_time=(datenum(2015,i,01,0,0,0)-datenum(1970,1,1))*86400;
end_time=(datenum(2015,i+1,01,0,0,0)-datenum(1970,1,1))*86400;
querystr=sprintf('select * from control.sc_sens_ACRS_Flow where time > %d and time < %d order by time asc',start_time,end_time);
query_result=MySQLQuery(querystr);
MFC1_time=[MFC1_time;query_result.time];
MFC1_value=[MFC1_value;query_result.value];
sprintf('Finished Month of %d %d',i,2015)
clear query_result
end
start_time=(datenum(2015,12,01,0,0,0)-datenum(1970,1,1))*86400;
end_time=(datenum(2016,1,01,0,0,0)-datenum(1970,1,1))*86400;
querystr=sprintf('select * from control.sc_sens_ACRS_Flow where time > %d and time < %d order by time asc',start_time,end_time);
query_result=MySQLQuery(querystr);
MFC1_time=[MFC1_time;query_result.time];
MFC1_value=[MFC1_value;query_result.value];
sprintf('Finished Month of %d %d',12,2015)
clear query_result
for i=1:11;
start_time=(datenum(2016,i,01,0,0,0)-datenum(1970,1,1))*86400;
end_time=(datenum(2016,i+1,01,0,0,0)-datenum(1970,1,1))*86400;
querystr=sprintf('select * from control.sc_sens_ACRS_Flow where time > %d and time < %d order by time asc',start_time,end_time);
query_result=MySQLQuery(querystr);
MFC1_time=[MFC1_time;query_result.time];
MFC1_value=[MFC1_value;query_result.value];
sprintf('Finished Month of %d %d',i,2016)
clear query_result
end
start_time=(datenum(2016,12,01,0,0,0)-datenum(1970,1,1))*86400;
end_time=(datenum(2017,1,01,0,0,0)-datenum(1970,1,1))*86400;
querystr=sprintf('select * from control.sc_sens_ACRS_Flow where time > %d and time < %d order by time asc',start_time,end_time);
query_result=MySQLQuery(querystr);
MFC1_time=[MFC1_time;query_result.time];
MFC1_value=[MFC1_value;query_result.value];
sprintf('Finished Month of %d %d',12,2016)
clear query_result

 MFC1_date= (MFC1_time./86400 + datenum([1969 12 31 24 0 0])); %In eastern timezone
 noflow_times=(MFC1_value<15);

end

