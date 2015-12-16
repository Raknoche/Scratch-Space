%%%Examples of quering the LUG for IQs

%% Get electron lifetime and 1-sigma from lug

%get all the IQs that contain the term 'electron_lifetime' and store them as an XLM structure
out = MySQLQuery('select values_xml from control.lug_iqs where iq_type = "electron_lifetime" and strikeme = 0 and algorithm_version = 2.10;');

for j=1:size(out.values_xml,1); % 22(first entry) to last entry

%read in electron lifetime and 1 sigma, using the XMLParser function
value1 = XMLParser(out.values_xml{j});
lifetime(j)=value1.iq.correction.fit.electron_attenuation_us;
sigma_lifetime(j)=value1.iq.correction.fit.sigma_electron_attenuation_us;

%get the date from file name pref1ix (as Matlab time)
file_date=value1.iq.global.filename_prefix;
Date_MatTime(j)=datenum(file_date(7:19),'yyyymmddTHHMM');

end

%% Make a quick plot of the electron lifetime
errorbar(Date_MatTime,lifetime,sigma_lifetime,'.k')
datetick %formats the x axis to be a date label
myfigview(16); xlabel('Date'); ylabel('Electron Lifetime (\mus)');



%% Track SE Size from bottom array over time

clear all
%Query the LUG
out = MySQLQuery('select values_xml from control.lug_iqs where iq_type = "single_e_kr" and strikeme = 0 and algorithm_version = 2.10;');

%Loop through query return, using XML Parser to pull out IQs
for j=1:size(out.values_xml,1); % 
    
%read in SE Size and Error using the XMLParser function
value1 = XMLParser(out.values_xml{j});
SE_mean(j)=value1.iq.correction.fit.e_mean_bottom;
SE_error(j)=value1.iq.correction.fit.e_mean_err_bottom;

%get the date from file name prefix (as Matlab time)
file_date=value1.iq.global.filename_prefix;
Date_MatTime(j)=datenum(file_date(7:19),'yyyymmddTHHMM');

end

figure
errorbar(Date_MatTime,SE_mean,SE_error,'.k');
datetick;
myfigview(16); xlabel('Date'); ylabel('SE Size - Bottom Array (phe)');
