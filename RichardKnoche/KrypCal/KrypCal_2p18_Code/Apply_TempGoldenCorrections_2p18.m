function [ s2_z, s2_xyz, s1_z, s1_xyz ] = Apply_TempGoldenCorrections_2p18(s2,s1,s2x,s2y,drift_time,file_name,LoadLocation)
%Loads IQs from a .mat file provides corrected pulse quantities

%s2 and s1 should be uncorrected
%s2x and s2y should come from d.x_cm and d.y_cm
%drift time should be in uSec
%file_name should be in the format: 'lux10_20140903T1918_cp13848'
%LoadLocation is a string that specifies the IQ .mat file


%Apply like this:
%    IQ_location='C:\Program Files\MATLAB\R2012a\bin\LUXCode\Scratch\RichardKnoche\KrypCal\KrypCal_2p18_Code\IQ_Query_Results';
%     [ s2_phe_both_z, s2_phe_both_xyz, s1_phe_both_z, s1_phe_both_xyz ] = Apply_TempGoldenCorrections_2p18(s2_phe_both,s1_phe_both,s2x,s2y,drift_time,'lux10_20140903T1918_cp13848',IQ_location);
    


%Load in the list of IQs that is given by the query function
load(LoadLocation)
 

%Get date of dataset
filename_prefix=file_name(1:19);
current_data_set_time=filename2epoch_framework(filename_prefix); %I think this is hard coded to some other time zone? +6 hours from MT?

%S2 Z Corr
[s2_norm_z s2_norm]=InterpS2IQ_Run04(filename_prefix, s2z_iq_times, s2_norm_z_iqs, s2_norm_iqs);

%S1 Z Corr
[index iq_time_zDep] = NearestIQ(filename_prefix,s1z_iq_times);
   z_dep_par_all = z_dep_both_values(index,:);
   z_dep_par_bot = z_dep_bottom_values(index,:);

%S2 XY Corr
[index iq_time_s2xyDep] =NearestIQ(filename_prefix,s2xy_iq_times);
    clear value1
    value1 = XMLParser(s2xy_out.values_xml{index});
   s2_x_bins = value1.iq.correction.fit.x_bin_center;
   s2_y_bins = value1.iq.correction.fit.y_bin_center;
   s2_map_all = value1.iq.correction.fit.norm_s2_both;
   s2_map_bottom = value1.iq.correction.fit.norm_s2_bottom;
   
%S1 XY Corr
   [index iq_time]  = NearestIQ(filename_prefix,s1xy_iq_times);
    clear value1
    value1 = XMLParser(s1xy_out.values_xml{index});      
   s1_x_bins = value1.iq.correction.fit.x_bin_center;
   s1_y_bins = value1.iq.correction.fit.y_bin_center;
   s1_map_all = value1.iq.correction.fit.norm_s1_both;
   s1_map_bottom = value1.iq.correction.fit.norm_s1_bottom;

%S1 XYZ Corr

   [index iq_time_s1xyzDep] = NearestIQ(filename_prefix,s1xyz_iq_times);
      
   clear value1
   value1 = XMLParser(s1xyz_out.values_xml{index});      
   s1_xyz_x_bins = value1.iq.correction.fit.x_bin_center;
   s1_xyz_y_bins = value1.iq.correction.fit.y_bin_center;
   s1_xyz_z_bins = value1.iq.correction.fit.z_bin_center;
   
   xx = size(s1_xyz_x_bins); yy = size(s1_xyz_y_bins); zz = size(s1_xyz_z_bins);
   
   s1_xyz_map_all = value1.iq.correction.fit.norm_s1_both_xyz; s1_xyz_map_all = reshape(s1_xyz_map_all,xx(2),yy(2),zz(2));
   s1_xyz_map_bottom = value1.iq.correction.fit.norm_s1_bottom_xyz;   s1_xyz_map_bottom = reshape(s1_xyz_map_bottom,xx(2),yy(2),zz(2));

%Detector center
   [index iq_time_det_center] = NearestIQ(filename_prefix,detcenter_iq_times);
    clear value1
   value1 = XMLParser(detcenter_out.values_xml{index});      
     
   det_x_center = value1.iq.correction.fit.x_center;
   det_y_center = value1.iq.correction.fit.y_center;
   det_z_center = value1.iq.correction.fit.z_center;
   
   s1_ref_z_ns=det_z_center.*1000; %conver uSec to ns


%-----------------------------------------------------------------------------------------------------------------
%-----------------------------------------------------------------------------------------------------------------
%% Calculating corrections ----------------------------------------------------------------------------
%-----------------------------------------------------------------------------------------------------------------
%----------------------------------------------------------------------------------------------------------------- 

s1_drift_ns=drift_time*1000;


%--------------------------------------------------------------------------
% Calculate Z corrections
%--------------------------------------------------------------------------

% Calculate S1 Z-correction

s1_z_correction = polyval(z_dep_par_all,s1_ref_z_ns./1000)./polyval(z_dep_par_all,s1_drift_ns./1000);
s1_z_correction(find(isnan(s1_z_correction))) = 1.0;

s1_z_correction_bot = polyval(z_dep_par_bot,s1_ref_z_ns./1000)./polyval(z_dep_par_bot,s1_drift_ns./1000);
s1_z_correction_bot(find(isnan(s1_z_correction_bot))) = 1.0;

% electron_lifetime_correction = exp(double(dp.z_drift_samples)./(100.*electron_lifetime));
% electron_lifetime_correction(find(isnan(electron_lifetime_correction))) = 1.0;
s2_z_correction=RK_1DCubicInterp_LinearExtrap(s2_norm_z,s2_norm, drift_time);
s2_z_correction(find(isnan(s2_z_correction))) = 1.0;

%--------------------------------------------------------------------------
% Calculate XY corrections
%--------------------------------------------------------------------------

% Calculate S1 XY-only corrections
s1xy_correction = interp2(s1_x_bins,s1_y_bins,s1_map_all,s2x,s2y,'spline',1);
s1xy_correction(find(s1xy_correction==0))=1.0;  s1xy_correction(find(isnan(s1xy_correction)))=1.0;

s1xy_correction_bot = interp2(s1_x_bins,s1_y_bins,s1_map_bottom,s2x,s2y,'spline',1);
s1xy_correction_bot(find(s1xy_correction_bot==0))=1.0;  s1xy_correction_bot(find(isnan(s1xy_correction_bot)))=1.0;

% Calculate S2 XY corrections

s2xy_correction = interp2(s2_x_bins,s2_y_bins,s2_map_all,s2x,s2y,'spline',1);
s2xy_correction(find(s2xy_correction==0))=1.0;  s2xy_correction(find(isnan(s2xy_correction)))=1.0;

s2xy_correction_bot = interp2(s2_x_bins,s2_y_bins,s2_map_bottom,s2x,s2y,'spline',1);
s2xy_correction_bot(find(s2xy_correction_bot==0))=1.0;  s2xy_correction_bot(find(isnan(s2xy_correction_bot)))=1.0;

%--------------------------------------------------------------------------
% Calculate XYZ corrections
%--------------------------------------------------------------------------


% Use S1 XYZ correction if available, otherwise use S1 XY + Z correction
if abs(current_data_set_time - iq_time_s1xyzDep) <= abs(current_data_set_time - iq_time)    
    s1xyz_correction = interp3(s1_xyz_x_bins,s1_xyz_y_bins,s1_xyz_z_bins,s1_xyz_map_all,s2x,s2y,drift_time,'spline',1);
    s1xyz_correction(find(isnan(s1xyz_correction))) = 1.0;

    s1xyz_correction_bot = interp3(s1_xyz_x_bins,s1_xyz_y_bins,s1_xyz_z_bins,s1_xyz_map_bottom,s2x,s2y,drift_time,'spline',1);
    s1xyz_correction_bot(find(isnan(s1xyz_correction_bot))) = 1.0;

else
    s1xyz_correction=s1xy_correction.*s1_z_correction;
    s1xyz_correction_bot=s1xy_correction_bot.*s1_z_correction_bot;


end


%-----------------------------------------------------------------------------------------------------------------
%-----------------------------------------------------------------------------------------------------------------
%% Applying corrections ----------------------------------------------------------------------------
%-----------------------------------------------------------------------------------------------------------------
%----------------------------------------------------------------------------------------------------------------- 


s2_z=s2.*s2_z_correction;
s2_xyz=s2.*s2_z_correction.*s2xy_correction;

s1_z=s1.*s1_z_correction;
s1_xyz=s1.*s1xyz_correction;


end

