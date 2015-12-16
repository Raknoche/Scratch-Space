function [] = Apply_TempGoldenCorrections_Query_2p16(SaveLocation)
%Can specify a location to save the query results as a .mat file, or leave the arguement blank to not save


%% Get all S2 corrections
out = MySQLQuery('select values_xml from control.lug_iqs where iq_type = "electron_lifetime" and strikeme = 0 and algorithm_version = 2.16;');


s2_norm_z_iqs = {};
s2_norm_iqs = {};
s2_bins_test = {};
s2_means_test = {};
s2z_iq_times = [];

for j=1:size(out.values_xml,1); % 22(first entry) to last entry
%read in electron lifetime and 1 sigma, using the XMLParser function
clear value1
value1 = XMLParser(out.values_xml{j});
           s2_norm_z_iqs = [s2_norm_z_iqs value1.iq.correction.fit.s2_both_norm_z_index];
           s2_norm_iqs = [s2_norm_iqs value1.iq.correction.fit.s2_both_norm_z];
           s2z_iq_times = [s2z_iq_times filename2epoch_framework(value1.iq.global.filename_prefix) ];
           
           s2_bins_test=[s2_bins_test value1.iq.correction.fit.bin_means_s2_both];
           s2_means_test=[s2_means_test value1.iq.correction.fit.means_s2_both];
end


%% Finding the S1 Z-dep values

z_dep_both_values = zeros(0,3);
z_dep_bottom_values = zeros(0,3);
s1z_iq_times = [];

clear out
out = MySQLQuery('select values_xml from control.lug_iqs where iq_type = "z_dep_s1_correction" and strikeme = 0 and algorithm_version = 2.16;');


for j=1:size(out.values_xml,1); % 22(first entry) to last entry
%read in electron lifetime and 1 sigma, using the XMLParser function
clear value1
value1 = XMLParser(out.values_xml{j});
           z_dep_both_values = vertcat(z_dep_both_values,value1.iq.correction.fit.s1_both_zdep_quad_fit);
           z_dep_bottom_values = vertcat(z_dep_bottom_values,value1.iq.correction.fit.s1_bottom_zdep_quad_fit);
           s1z_iq_times = [s1z_iq_times filename2epoch_framework(value1.iq.global.filename_prefix)];
end
   

%% Finding the S2 xy correction map values

s2_xy_index = [];
s2xy_iq_times = [];


clear out
out = MySQLQuery('select values_xml from control.lug_iqs where iq_type = "s2_xy_correction" and strikeme = 0 and algorithm_version = 2.16;');


for j=1:size(out.values_xml,1); % 22(first entry) to last entry
%read in electron lifetime and 1 sigma, using the XMLParser function
clear value1
value1 = XMLParser(out.values_xml{j});
           s2_xy_index = [s2_xy_index j];
           s2xy_iq_times = [s2xy_iq_times filename2epoch_framework(value1.iq.global.filename_prefix)];
end
 
s2xy_out=out;
  

%% Finding the S1 xy correction map values
   
s1_xy_index = [];
s1xy_iq_times = [];


clear out
out = MySQLQuery('select values_xml from control.lug_iqs where iq_type = "s1_xy_correction" and strikeme = 0 and algorithm_version = 2.16;');


for j=1:size(out.values_xml,1); % 22(first entry) to last entry
%read in electron lifetime and 1 sigma, using the XMLParser function
clear value1
value1 = XMLParser(out.values_xml{j});
           s1_xy_index = [s1_xy_index j];
           s1xy_iq_times = [s1xy_iq_times filename2epoch_framework(value1.iq.global.filename_prefix)];
end

s1xy_out=out;
 

%% Finding the S1 xyz correction map values
   
s1_xyz_index = [];
s1xyz_iq_times = [];

clear out
out = MySQLQuery('select values_xml from control.lug_iqs where iq_type = "s1_xyz_correction" and strikeme = 0 and algorithm_version = 2.16;');


for j=1:size(out.values_xml,1); % 22(first entry) to last entry
%read in electron lifetime and 1 sigma, using the XMLParser function
clear value1
value1 = XMLParser(out.values_xml{j});
           s1_xyz_index = [s1_xyz_index j];
           s1xyz_iq_times = [s1xyz_iq_times filename2epoch_framework(value1.iq.global.filename_prefix)];
end
 
s1xyz_out=out;
    
%% Finding the detector center according to KrypCal
   
det_center_index = [];
detcenter_iq_times = [];


clear out
out = MySQLQuery('select values_xml from control.lug_iqs where iq_type = "s1ab_xy" and strikeme = 0 and algorithm_version = 2.16;');


for j=1:size(out.values_xml,1); % 22(first entry) to last entry
%read in electron lifetime and 1 sigma, using the XMLParser function
clear value1
value1 = XMLParser(out.values_xml{j});
           det_center_index = [det_center_index j];
           detcenter_iq_times = [detcenter_iq_times filename2epoch_framework(value1.iq.global.filename_prefix)];
end
 
detcenter_out=out;

%saving .mat file
if nargin==1
save(SaveLocation,'s2z_iq_times', 's2_norm_z_iqs', 's2_norm_iqs',...
    's1z_iq_times','z_dep_both_values','z_dep_bottom_values',...
    's2xy_iq_times','s2xy_out',...
    's1xy_iq_times','s1xy_out',...
    's1xyz_iq_times','s1xyz_out',...
    'detcenter_iq_times','detcenter_out');
end

end

