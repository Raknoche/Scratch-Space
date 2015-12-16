
%  path='C:\Users\Richard\Desktop\TestData\CH3T';
 path='C:\Users\Richard\Desktop\TestData\lux10_20150122T1337_cp18051';

rqs_to_load = {'pulse_area_phe','event_timestamp_samples'...
   ,'pulse_classification' ...
   ,'z_drift_samples' , 's1s2_pairing'...
   ,'top_bottom_ratio','x_cm','y_cm'...
   ,'full_evt_area_phe',...
   'event_number','chi2','prompt_fraction','aft_t1_samples','pulse_start_samples',...
   'pulse_end_samples','top_bottom_asymmetry','aft_t0_samples','aft_t2_samples',...
   'full_evt_area_phe','admin','Kr83fit_s1a_area_phe','Kr83fit_s1b_area_phe',...
   'hft_t50r_samples','hft_t50l_samples','Kr83fit_dt_samples'};

d = LUXLoadMultipleRQMs_framework(path,rqs_to_load);  

dp = d;
file_name='lux10_20150122T1337_cp18051';


%Get date of dataset
filename_prefix=file_name(1:19);
current_data_set_time=filename2epoch_framework(filename_prefix); %I think this is hard coded to some other time zone? +6 hours from MT?

%Load in Circ Outages - NOTE THIS ONLY WORKS FOR DATES BEFORE 11/1/2015!!
% load(OutagePath);
%noflow_times has timestamps of SC return

%If you want to run data after 11/1/2015, (and before 1/2017) refresh the noflow times list by running this (takes ~20 minutes)
%  [MFC1_time noflow_times_cut]=FindCircOutages; 
%  noflow_times=MFC1_time(noflow_times_cut);

%% Get all S2 corrections
out = MySQLQuery('select values_xml from control.lug_iqs where iq_type = "electron_lifetime" and strikeme = 0 and algorithm_version = 2.16;');


s2_norm_z_iqs = {};
s2_norm_iqs = {};
s2_bins_test = {};
s2_means_test = {};
iq_times = [];

for j=1:size(out.values_xml,1); % 22(first entry) to last entry
%read in electron lifetime and 1 sigma, using the XMLParser function
clear value1
value1 = XMLParser(out.values_xml{j});
           s2_norm_z_iqs = [s2_norm_z_iqs value1.iq.correction.fit.s2_both_norm_z_index];
           s2_norm_iqs = [s2_norm_iqs value1.iq.correction.fit.s2_both_norm_z];
           iq_times = [iq_times filename2epoch_framework(value1.iq.global.filename_prefix) ];
           
           s2_bins_test=[s2_bins_test value1.iq.correction.fit.bin_means_s2_both];
           s2_means_test=[s2_means_test value1.iq.correction.fit.means_s2_both];
end

 [s2_norm_z s2_norm]=InterpS2IQ_Run04(filename_prefix, iq_times, s2_norm_z_iqs, s2_norm_iqs);


%% Finding the S1 Z-dep values

z_dep_both_values = zeros(0,3);
z_dep_bottom_values = zeros(0,3);
iq_times = [];

clear out
out = MySQLQuery('select values_xml from control.lug_iqs where iq_type = "z_dep_s1_correction" and strikeme = 0 and algorithm_version = 2.16;');


for j=1:size(out.values_xml,1); % 22(first entry) to last entry
%read in electron lifetime and 1 sigma, using the XMLParser function
clear value1
value1 = XMLParser(out.values_xml{j});
           z_dep_both_values = vertcat(z_dep_both_values,value1.iq.correction.fit.s1_both_zdep_quad_fit);
           z_dep_bottom_values = vertcat(z_dep_bottom_values,value1.iq.correction.fit.s1_bottom_zdep_quad_fit);
           iq_times = [iq_times filename2epoch_framework(value1.iq.global.filename_prefix)];
end
 
   [index iq_time_zDep] = NearestIQ(filename_prefix,iq_times);
   
   z_dep_par_all = z_dep_both_values(index,:);
   z_dep_par_bot = z_dep_bottom_values(index,:);


%% Finding the S2 xy correction map values

s2_xy_index = [];
iq_times = [];


clear out
out = MySQLQuery('select values_xml from control.lug_iqs where iq_type = "s2_xy_correction" and strikeme = 0 and algorithm_version = 2.16;');


for j=1:size(out.values_xml,1); % 22(first entry) to last entry
%read in electron lifetime and 1 sigma, using the XMLParser function
clear value1
value1 = XMLParser(out.values_xml{j});
           s2_xy_index = [s2_xy_index j];
           iq_times = [iq_times filename2epoch_framework(value1.iq.global.filename_prefix)];
end
 
   [index iq_time_s2xyDep] =NearestIQ(filename_prefix,iq_times);
 
    clear value1
    value1 = XMLParser(out.values_xml{index});
   s2_x_bins = value1.iq.correction.fit.x_bin_center;
   s2_y_bins = value1.iq.correction.fit.y_bin_center;
   s2_map_all = value1.iq.correction.fit.norm_s2_both;
   s2_map_bottom = value1.iq.correction.fit.norm_s2_bottom;

%% Finding the S1 xy correction map values
   
s1_xy_index = [];
iq_times = [];


clear out
out = MySQLQuery('select values_xml from control.lug_iqs where iq_type = "s1_xy_correction" and strikeme = 0 and algorithm_version = 2.16;');


for j=1:size(out.values_xml,1); % 22(first entry) to last entry
%read in electron lifetime and 1 sigma, using the XMLParser function
clear value1
value1 = XMLParser(out.values_xml{j});
           s1_xy_index = [s1_xy_index j];
           iq_times = [iq_times filename2epoch_framework(value1.iq.global.filename_prefix)];
end

 
   [index iq_time]  = NearestIQ(filename_prefix,iq_times);
    clear value1
    value1 = XMLParser(out.values_xml{index});      
   s1_x_bins = value1.iq.correction.fit.x_bin_center;
   s1_y_bins = value1.iq.correction.fit.y_bin_center;
   s1_map_all = value1.iq.correction.fit.norm_s1_both;
   s1_map_bottom = value1.iq.correction.fit.norm_s1_bottom;
   
   %% Finding the detector center according to KrypCal
   
det_center_index = [];
iq_times = [];


clear out
out = MySQLQuery('select values_xml from control.lug_iqs where iq_type = "s1ab_xy" and strikeme = 0 and algorithm_version = 2.16;');


for j=1:size(out.values_xml,1); % 22(first entry) to last entry
%read in electron lifetime and 1 sigma, using the XMLParser function
clear value1
value1 = XMLParser(out.values_xml{j});
           det_center_index = [det_center_index j];
           iq_times = [iq_times filename2epoch_framework(value1.iq.global.filename_prefix)];
end
 
   [index iq_time_det_center] = NearestIQ(filename_prefix,iq_times);
    clear value1
   value1 = XMLParser(out.values_xml{index});      
     
   det_x_center = value1.iq.correction.fit.x_center;
   det_y_center = value1.iq.correction.fit.y_center;
   det_z_center = value1.iq.correction.fit.z_center;
   
   s1_ref_z_ns=det_z_center.*1000; %conver uSec to ns

   %%
dp.z_corrected_pulse_area_all_phe = dp.pulse_area_phe;
dp.xyz_corrected_pulse_area_all_phe = dp.pulse_area_phe;

dp.z_corrected_pulse_area_bot_phe = dp.pulse_area_phe./(1+dp.top_bottom_ratio);
dp.xyz_corrected_pulse_area_bot_phe = dp.pulse_area_phe./(1+dp.top_bottom_ratio);

N = length(dp.pulse_area_phe);

dp.correction_s2_z_dependence = zeros(10,N, 'single');
dp.correction_s1_z_dependence = zeros(10,N, 'single');
dp.correction_s1_xyz_dependence = zeros(10,N, 'single');
dp.correction_s2_xy_dependence = zeros(10,N, 'single');
% Initialize variable with negative value so that we're sure it goes through the code 
dp.s1_correction_is_3d= zeros(1,N,'int8');
dp.s1_correction_is_3d(:) = -99;

%% Finding the drift time, x and y associated with S1

% Finding the drift time and xy position from the largest S2 pulse in the
% pairing


[a b] = size(dp.z_drift_samples); % b is number of events in file ??

drift = dp.z_drift_samples;
drift(find(isnan(drift))) = 0.0;
s1_drift_ns = +(dp.pulse_classification==1);
s1_x_cm = +(dp.pulse_classification==1);
s1_y_cm = +(dp.pulse_classification==1);
s2_phe =  dp.pulse_area_phe.*(dp.pulse_classification==2);
s2_phe(find(isnan(s2_phe))) = 0.0;

%s1s = (sum(dp.pulse_classification==1)==1) & (sum(dp.pulse_classification==2)>0) & (sum(dp.s1s2_pairing==1)>1);

 s1s = (sum(dp.s1s2_pairing==1)>1); % tells you if its an event with at least one s1 and one s2.
  
 for i=1:length(s1s)  % for each event in the file
     if s1s(i)>0 % if the event had an S1
       if s1s(i)==1
           %... pair first S1 with first S2 -commented out by AD, July 21
           %2014
         %[r c v] = find(drift(:,i),1,'first');
         %[c1 r1 v1] = find(dp.pulse_classification(:,i)==1,1,'first');
         %[c2 r2 v2] = find(dp.pulse_classification(:,i)==2,1,'first');
         %s1_drift_ns(c1,i) = 10.*v;
         %s1_x_cm(c1,i) = dp.x_cm(c2,i);
         %s1_y_cm(c1,i) = dp.y_cm(c2,i);   
        
         [v r] = max(s2_phe(:,i));% value and index (r for row) of Max S2
         [c1 r1 v1] = find(dp.pulse_classification(:,i)==1,1,'first');
         %[c2 r2 v2] = find(dp.pulse_classification(:,i)==2,1,'first');
         s1_drift_ns(c1,i) = 10.*drift(r,i); %row c1,r and column i
         s1_x_cm(c1,i) = dp.x_cm(r,i);
         s1_y_cm(c1,i) = dp.y_cm(r,i);  
         
       else
          s1_drift_ns(:,i) = 0;
          s1_x_cm(:,i) = 0;
          s1_y_cm(:,i) = 0;
       end  
     end
 end

 
 
% fprintf('Done finding S1 depths :) \n');

s1_drift_ns(find((s1_drift_ns==0))) = nan;
s1_x_cm(find((s1_x_cm==0))) = nan;
s1_y_cm(find((s1_y_cm==0))) = nan;



%--------------------------------------------------------------------------
% Calculate Z corrections
%--------------------------------------------------------------------------

% Calculate S1 Z-correction

s1_z_correction = polyval(z_dep_par_all,s1_ref_z_ns./1000)./polyval(z_dep_par_all,s1_drift_ns./1000);
s1_z_correction(find(isnan(s1_z_correction))) = 1.0;

s1_z_correction_bot = polyval(z_dep_par_bot,s1_ref_z_ns./1000)./polyval(z_dep_par_bot,s1_drift_ns./1000);
s1_z_correction_bot(find(isnan(s1_z_correction_bot))) = 1.0;

% Calculate electron lifetime correction (S2 Z-correction)
% Reading the values of electron lifetime from the LUG and interpolating
% between the two either side

% electron_lifetime_correction = exp(double(dp.z_drift_samples)./(100.*electron_lifetime));
% electron_lifetime_correction(find(isnan(electron_lifetime_correction))) = 1.0;
s2_z_correction=interp1(s2_norm_z,s2_norm, double(dp.z_drift_samples)/100,'spline',1);
s2_z_correction(find(isnan(s2_z_correction))) = 1.0;

%--------------------------------------------------------------------------
% Calculate XY corrections
%--------------------------------------------------------------------------


% Calculate S1 XY-only corrections

s2x = dp.x_cm.*(+(dp.pulse_classification==2)); s2x(find(s2x==0)) = nan;
s2y = dp.y_cm.*(+(dp.pulse_classification==2)); s2y(find(s2y==0)) = nan;

s1xy_correction = interp2(s1_x_bins,s1_y_bins,s1_map_all,s1_x_cm,s1_y_cm,'spline',1);
s1xy_correction = s1xy_correction.*(+dp.pulse_classification==1);
s1xy_correction(find(s1xy_correction==0))=1.0;  s1xy_correction(find(isnan(s1xy_correction)))=1.0;

s1xy_correction_bot = interp2(s1_x_bins,s1_y_bins,s1_map_bottom,s1_x_cm,s1_y_cm,'spline',1);
s1xy_correction_bot = s1xy_correction_bot.*(+dp.pulse_classification==1);
s1xy_correction_bot(find(s1xy_correction_bot==0))=1.0;  s1xy_correction_bot(find(isnan(s1xy_correction_bot)))=1.0;

% Calculate S2 XY corrections

s2xy_correction = interp2(s2_x_bins,s2_y_bins,s2_map_all,dp.x_cm,dp.y_cm,'spline',1);
s2xy_correction = s2xy_correction.*(+dp.pulse_classification==2);
s2xy_correction(find(s2xy_correction==0))=1.0;  s2xy_correction(find(isnan(s2xy_correction)))=1.0;

s2xy_correction_bot = interp2(s2_x_bins,s2_y_bins,s2_map_bottom,dp.x_cm,dp.y_cm,'spline',1);
s2xy_correction_bot = s2xy_correction_bot.*(+dp.pulse_classification==2);
s2xy_correction_bot(find(s2xy_correction_bot==0))=1.0;  s2xy_correction_bot(find(isnan(s2xy_correction_bot)))=1.0;

