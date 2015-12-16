%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is used to apply S1 and S2 XYZ corrections 100 V/cm Run03 data
% Do NOT use the standard DP corrections for 100 V/cm.  The new field
% causes the S1 XYZ corrections to be wrong.
%
% Inputs: dp struct that contains z_drift_samples, pulse_classification,
% x_cm, y_cm, pulse_area_phe, s1s2_pairing, and top_bottom_ratio
%
% Outputs: the same dp struct with dp.z_corrected_pulse_area_all_phe,
% dp.z_corrected_pulse_area_bot_phe, dp.xyz_corrected_pulse_area_all_phe,
% dp.xyz_corrected_pulse_area_bot_phe added in
%
% Contact: Richard Knoche
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ dp ] = Apply_100Vcm_Corr( dp )

%Loading Corrections for 100 V/cm
load('KrypCal_100Vcm_Corr.mat')

%%
%Finding S1 Depths

[a b] = size(dp.z_drift_samples); % b is number of events in file ??

drift = dp.z_drift_samples; %in 10 nanoseconds
drift(find(isnan(drift))) = 0.0;
s1_drift_ns = +(dp.pulse_classification==1);
s1_x_cm = +(dp.pulse_classification==1);
s1_y_cm = +(dp.pulse_classification==1);
s2x = dp.x_cm.*(+(dp.pulse_classification==2)); s2x(find(s2x==0)) = nan;
s2y = dp.y_cm.*(+(dp.pulse_classification==2)); s2y(find(s2y==0)) = nan;
s2_phe =  dp.pulse_area_phe.*(dp.pulse_classification==2);
s2_phe(find(isnan(s2_phe))) = 0.0;

 s1s = (sum(dp.s1s2_pairing==1)>1); % tells you if its an event with at least one s1 and one s2.
  
 for i=1:length(s1s)  % for each event in the file
     if s1s(i)>0 % if the event had an S1
       if s1s(i)==1
        
         [v r] = max(s2_phe(:,i));% value and index (r for row) of Max S2
         [c1 r1 v1] = find(dp.pulse_classification(:,i)==1,1,'first');
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

s1_drift_ns(find((s1_drift_ns==0))) = nan;
s1_x_cm(find((s1_x_cm==0))) = nan;
s1_y_cm(find((s1_y_cm==0))) = nan;


%--------------------------------------------------------------------------
% Making Correction Matrices
%--------------------------------------------------------------------------

%S2 Z
electron_lifetime_correction = exp(double(dp.z_drift_samples)./(100.*electron_lifetime));
electron_lifetime_correction(find(isnan(electron_lifetime_correction))) = 1.0;

%S2 XY

s2xy_correction = interp2(s2x_bins,s2y_bins,norm_S2_both,dp.x_cm,dp.y_cm,'spline',1);
s2xy_correction = s2xy_correction.*(+dp.pulse_classification==2);
s2xy_correction(find(s2xy_correction==0))=1.0;  s2xy_correction(find(isnan(s2xy_correction)))=1.0;

s2xy_correction_bot = interp2(s2x_bins,s2y_bins,norm_S2_bot,dp.x_cm,dp.y_cm,'spline',1);
s2xy_correction_bot = s2xy_correction_bot.*(+dp.pulse_classification==2);
s2xy_correction_bot(find(s2xy_correction_bot==0))=1.0;  s2xy_correction_bot(find(isnan(s2xy_correction_bot)))=1.0;



%S1 Z
s1_z_correction = polyval(P_s1_both,(det_edge-4.32)/2)./polyval(P_s1_both,s1_drift_ns./1000);
s1_z_correction(find(isnan(s1_z_correction))) = 1.0;

s1_z_correction_bot = polyval(P_s1_bottom,(det_edge-4.32)/2)./polyval(P_s1_bottom,s1_drift_ns./1000);
s1_z_correction_bot(find(isnan(s1_z_correction_bot))) = 1.0;


%S1 XY
s1xy_correction = interp2(s1_x_bins,s1_y_bins,norm_S1_all,dp.x_cm,dp.y_cm,'spline',1);
s1xy_correction = s1xy_correction.*(+dp.pulse_classification==1);
s1xy_correction(find(s1xy_correction==0))=1.0;  s1xy_correction(find(isnan(s1xy_correction)))=1.0;

s1xy_correction_bot = interp2(s1_x_bins,s1_y_bins,norm_S1_bot,dp.x_cm,dp.y_cm,'spline',1);
s1xy_correction_bot = s1xy_correction_bot.*(+dp.pulse_classification==1);
s1xy_correction_bot(find(s1xy_correction_bot==0))=1.0;  s1xy_correction_bot(find(isnan(s1xy_correction_bot)))=1.0;


%S1 XYZ
s1xyz_correction = interp3(s1xbins,s1ybins,s1zbins,norm_s1_both_xyz,s1_x_cm,s1_y_cm,s1_drift_ns./1000,'spline');
s1xyz_correction(find(isnan(s1xyz_correction))) = 1.0;

s1xyz_correction_bot = interp3(s1xbins,s1ybins,s1zbins,norm_s1_bottom_xyz,s1_x_cm,s1_y_cm,s1_drift_ns./1000,'spline');
s1xyz_correction_bot(find(isnan(s1xyz_correction_bot))) = 1.0;

%--------------------------------------------------------------------------
% Apply Z corrections
%--------------------------------------------------------------------------

dp.z_corrected_pulse_area_all_phe = dp.pulse_area_phe;
dp.z_corrected_pulse_area_bot_phe = dp.pulse_area_phe./(1+dp.top_bottom_ratio);

dp.xyz_corrected_pulse_area_all_phe = dp.pulse_area_phe;
dp.xy_z_corrected_pulse_area_all_phe = dp.pulse_area_phe;
dp.xyz_corrected_pulse_area_bot_phe = dp.pulse_area_phe./(1+dp.top_bottom_ratio);


dp.z_corrected_pulse_area_all_phe = dp.z_corrected_pulse_area_all_phe.*electron_lifetime_correction.*s1_z_correction; 
dp.z_corrected_pulse_area_bot_phe = dp.z_corrected_pulse_area_bot_phe.*electron_lifetime_correction.*s1_z_correction_bot; 

%--------------------------------------------------------------------------
% Apply XYZ corrections
%--------------------------------------------------------------------------

dp.xyz_corrected_pulse_area_all_phe = dp.xyz_corrected_pulse_area_all_phe.*electron_lifetime_correction.*s2xy_correction.*s1xyz_correction;
dp.xyz_corrected_pulse_area_bot_phe = dp.xyz_corrected_pulse_area_bot_phe.*electron_lifetime_correction.*s2xy_correction_bot.*s1xyz_correction_bot;

% dp.xy_z_corrected_pulse_area_all_phe = dp.xy_z_corrected_pulse_area_all_phe.*electron_lifetime_correction.*s2xy_correction.*s1xy_correction.*s1_z_correction;

end

