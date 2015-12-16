%Find IQs that that are Kr-83 data

% % For TAGCB=1,7,1,10,2
% folders_to_load={'C:\Users\Richard\Desktop\LUX Analysis\Run04_CH3T_TAGCB0107011002\lux10_20140807T0906_cp10194','C:\Users\Richard\Desktop\LUX Analysis\Run04_CH3T_TAGCB0107011002\lux10_20140807T1345_cp10200','C:\Users\Richard\Desktop\LUX Analysis\Run04_CH3T_TAGCB0107011002\lux10_20140807T1839_cp10202','C:\Users\Richard\Desktop\LUX Analysis\Run04_CH3T_TAGCB0107011002\lux10_20140808T1400_cp10205','C:\Users\Richard\Desktop\LUX Analysis\Run04_CH3T_TAGCB0107011002\lux10_20140808T1849_cp10208',...
%     'C:\Users\Richard\Desktop\LUX Analysis\Run04_CH3T_TAGCB0107011002\lux10_20140809T0958_cp10219'}; %'C:\Users\Richard\Desktop\LUX Analysis\Run04_CH3T_TAGCB0107011002\lux10_21040808T0931_cp10204'
% %Background for TAGCB=1,7,1,10,2
% % folders_to_load={'C:\Users\Richard\Desktop\LUX Analysis\Run04_CH3T_TAGCB0107011002\lux10_20140802T2312_cp10117'};
% 


% For TAGCB=1,7,2,10,2
% folders_to_load={'C:\Users\Richard\Desktop\LUX Analysis\Run04_CH3T_TAGCB0107021002\lux10_20140812T1434_cp10482','C:\Users\Richard\Desktop\LUX Analysis\Run04_CH3T_TAGCB0107021002\lux10_20140812T2033_cp10358','C:\Users\Richard\Desktop\LUX Analysis\Run04_CH3T_TAGCB0107021002\lux10_20140813T0248_cp10356',...
%     'C:\Users\Richard\Desktop\LUX Analysis\Run04_CH3T_TAGCB0107021002\lux10_20140813T1444_cp10381','C:\Users\Richard\Desktop\LUX Analysis\Run04_CH3T_TAGCB0107021002\lux10_20140813T2055_cp10483','C:\Users\Richard\Desktop\LUX Analysis\Run04_CH3T_TAGCB0107021002\lux10_20140814T0310_cp10488',...
%     'C:\Users\Richard\Desktop\LUX Analysis\Run04_CH3T_TAGCB0107021002\lux10_20140814T0926_cp10503'};


%For TAGCB=1,7,1,8.5,2
folders_to_load={'C:\Users\Richard\Desktop\LUX Analysis\Run04_CH3T_TAGCB01070108502\lux10_20140905T1531_cp11010','C:\Users\Richard\Desktop\LUX Analysis\Run04_CH3T_TAGCB01070108502\lux10_20140905T2036_cp11012','C:\Users\Richard\Desktop\LUX Analysis\Run04_CH3T_TAGCB01070108502\lux10_20140906T0144_cp11014',...
    'C:\Users\Richard\Desktop\LUX Analysis\Run04_CH3T_TAGCB01070108502\lux10_20140906T0650_cp11016','C:\Users\Richard\Desktop\LUX Analysis\Run04_CH3T_TAGCB01070108502\lux10_20140906T1131_cp11005','C:\Users\Richard\Desktop\LUX Analysis\Run04_CH3T_TAGCB01070108502\lux10_20140906T1615_cp11018'};


rqs_to_load = {'event_number','pulse_area_phe','event_timestamp_samples','x_cm','y_cm'...
   ,'top_bottom_ratio','pulse_classification' ...
   ,'full_evt_area_phe'...
   ,'z_drift_samples' , 's1s2_pairing','golden'...
   ,'z_corrected_pulse_area_all_phe','xyz_corrected_pulse_area_all_phe','z_corrected_pulse_area_bot_phe'...
   ,'xyz_corrected_pulse_area_bot_phe'...
   ,'x_corrected','y_corrected'};
    
ii=1;

    path=folders_to_load{ii};
    
    d = LUXLoadMultipleRQMs_framework(path,rqs_to_load);
     
           %Correcting x and y

 myname = 'Corrections_PositionCorrection';           
 position_correction_path = 'C:\Program Files\MATLAB\R2012a\bin\LuxCode\Trunk\DataProcessing\MatlabModules\Corrections_PositionCorrection\Corrections_PositionCorrection.m';
 IQs = dir([position_correction_path(1:(end-numel(myname)-2)) '*Mercury.mat']);
 table_corrections = load(IQs(end).name);

 d = Corrections_PositionCorrection_FunctionUseThis(d,table_corrections);
           % defining cuts
        zcut_min = 30;%us
        zcut_max = 300;%us
        rcut_min = 0;
        rcut_max = 25;%cm
        s1area_bound_min = 0;%2 PMT hits
        s1area_bound_max = 10^6;%include Kr 83 data and gammas

        s2area_bound_min = 100;%100 ~ 0.5 keVee
        s2area_bound_max = 10^6;% %both PMT Arrays and Kr83

         
    d.s1s2_pairing(isnan(d.s1s2_pairing))=0; %by default it returns nan when there is no pair
    d.z_drift_samples(isnan(d.z_drift_samples)) = 0.0; % get rid of NaN   
    
    s1s2_pairing = d.s1s2_pairing;   
    num_paired_pulses = sum(s1s2_pairing);
    
    events=sort(size(d.pulse_area_phe)); %The first element is the number of pulses. sometimes it's 5, sometimes it's 10
          
    s1_area_cut= inrange(d.pulse_area_phe,[s1area_bound_min,s1area_bound_max]);
    s2_area_cut= inrange(d.pulse_area_phe,[s2area_bound_min,s2area_bound_max]);
    
    s1_class=d.pulse_classification==1 & d.golden==1 & s1_area_cut ;% & squeeze(sum((d.peak_area_phe>0.3),2)>1); %at least 3 PMTs fire
    s2_class=d.pulse_classification==2 & d.golden==1 & s2_area_cut ;      
    
     s1_before_s2_cut = repmat(num_paired_pulses ==2 ,events(1),1) ...
                         & repmat(sum(s1_class,1)>=1,events(1),1) & repmat(sum(s2_class,1)==1,events(1),1)...
                         & cumsum(s1_class+s2_class)<=2 ; % ignore all S1 that happen after the s2 !
        
    s1_single_cut = s1_class & s1_before_s2_cut; %
    s2_single_cut = s2_class & s1_before_s2_cut;
    
    drift_time = d.z_drift_samples(s2_single_cut)/100; %units of us
        
        
    d.phe_bottom=d.pulse_area_phe./(1+d.top_bottom_ratio); %bottom PMT pulse area
    
    s1_phe_both = d.pulse_area_phe(s1_single_cut);
    s1_phe_bottom = d.phe_bottom(s1_single_cut);
    s1_phe_both_z=d.z_corrected_pulse_area_all_phe(s1_single_cut);
    s1_phe_both_xyz=d.xyz_corrected_pulse_area_all_phe(s1_single_cut);

    s2_phe_both = d.pulse_area_phe(s2_single_cut);
    s2_phe_bottom = d.phe_bottom(s2_single_cut);
    s2_phe_bottom_z=d.z_corrected_pulse_area_bot_phe(s2_single_cut); 
    s2_phe_bottom_xyz=d.xyz_corrected_pulse_area_bot_phe(s2_single_cut); 

    s2x = d.x_cm(s2_single_cut);
    s2y = d.y_cm(s2_single_cut);
    s2radius = (s2x.^2+s2y.^2).^(0.5);
    s2x_c = d.x_corrected(s2_single_cut);
    s2y_c = d.y_corrected(s2_single_cut);
    s2radius_c = (s2x_c.^2+s2y_c.^2).^(0.5);
    %s2x_taxi = d.x_tmplt_corrected(s2_single_cut);
    %s2y_taxi = d.y_tmplt_corrected(s2_single_cut);
    
 
        timestamp_vec_1 = sort([d.livetime_latch_samples d.livetime_end_samples]);
    livetime_sec = sum(timestamp_vec_1(2:2:end) - timestamp_vec_1(1:2:end)) / 1e8;
    
    evt_cut=logical(sum(s2_single_cut));%Cut for all the events passing the single S1 & S2 cut
    event_number=d.event_number(evt_cut)';
    event_timestamp_samples=d.event_timestamp_samples(evt_cut)';
   
    
        livetime_samples_file(1) = sum(timestamp_vec_1(2:2:end) - timestamp_vec_1(1:2:end));
    livetime_frac_file(1)=livetime_samples_file(1)/((event_timestamp_samples(end)-event_timestamp_samples(1)));         
     event_time_min=flipdim(event_timestamp_samples,1)*10^(-8)/60;
     event_time_hour=event_time_min/60;
     bin_size=1;
     start_time=0;
     end_time=max(event_time_hour);
     time_x=bin_size/2+start_time:bin_size:end_time; %hours. HIST bin center
     count_evt_time=hist(event_time_hour,time_x);
     hist_evt_time=count_evt_time/bin_size/3600/livetime_frac_file(1); %number per second
     
     sigma_hist_evt_time=sqrt(count_evt_time)/bin_size/3600/livetime_frac_file(1);
     total_sigma_hist_evt_time=sigma_hist_evt_time;
     total_hist_evt_time=hist_evt_time;
     total_time_x=time_x;
            %time_since_1e5_phe = d.time_since_1e5_phe(evt_cut);
            %time_since_3e5_phe = d.time_since_3e5_phe(evt_cut);
            %time_since_5e5_phe = d.time_since_5e5_phe(evt_cut);
            %time_since_8e5_phe = d.time_since_8e5_phe(evt_cut);
            
    
    clean_cut=-(s1_phe_both+s2_phe_both)' + d.full_evt_area_phe(logical(squeeze(sum(s2_single_cut,1)))) < 100 ...
        & drift_time' > 0;%more than 1/2 the area is S1+S2. And S1 before S2
    clean_cut=clean_cut';
   
    %s1_double_hit_cut= sum(d.peak_area_phe(permute(repmat(s1_single_cut,[1,1,122]),[1 3 2])) >0.3, 2) > 2; 
    
    
%     s1_peak_area_phe = squeeze(sum(permute(repmat(s1_single_cut,[1 1 122]),[1 3 2]).*d.peak_area_phe,1));
%     s1_multiplicity_corr = sum(s1_peak_area_phe>0.25);
%     for pp=1:2:121
% 	s1_multiplicity_corr(sum(s1_peak_area_phe(pp:(pp+1),:)>0.25)==2 & s1_multiplicity_corr==2) = 1;
%     end
%     s1_multiplicity_corr=s1_multiplicity_corr(logical(sum(s1_single_cut)))'; % get back to s1_cut size
%     s1_multiplicity_corr_cut=s1_multiplicity_corr~=1;
    
    s1_min=0;
    s1_max=50;
    
    liquid_cut= inrange(drift_time,[5 320]) & inrange(s2radius_c,[0,24.5]) & inrange(s1_phe_both_xyz,[s1_min s1_max]) ...
    & inrange(s2_phe_bottom_xyz,[150 4000])  & clean_cut & s2_phe_both>100 ; % s1_multiplicity_corr_cut~=0
         
         s1_phe_both_2=s1_phe_both;
         s1_phe_both_z_2=s1_phe_both_z ;
         s1_phe_both_xyz_2=s1_phe_both_xyz ;        
         s2_phe_bottom_2=s2_phe_bottom ;
         s2_phe_bottom_xyz_2= s2_phe_bottom_xyz ;
         s2_phe_bottom_z_2= s2_phe_bottom_z ;
         s2_phe_both_2=s2_phe_both ;
         drift_time_2=drift_time ;
         event_timestamp_samples_2=event_timestamp_samples ;
         event_number_2=event_number;
         s2x_2=s2x ;
         s2y_2=s2y ;
         s2x_c_2 = s2x_c ;
         s2y_c_2 = s2y_c ;
%          s1_multiplicity_corr_cut_2=s1_multiplicity_corr_cut;
         clean_cut_2=clean_cut;
         liquid_cut_2=liquid_cut;
         %kr_cut_2=kr_cut;
        
 time_diff_samples_1=0; %no time difference from first file
    file_start_samples(ii)=time_diff_samples_1;
    file_end_samples(ii)=time_diff_samples_1+d.event_timestamp_samples(end);
    
for  ii=2:length(folders_to_load) % 
      tic; 
        clear d1;
        d1=d;
        clear d;
        
               path=folders_to_load{ii};
        d = LUXLoadMultipleRQMs_framework(path,rqs_to_load);

                   %Correcting x and y

 myname = 'Corrections_PositionCorrection';           
 position_correction_path = 'C:\Program Files\MATLAB\R2012a\bin\LuxCode\Trunk\DataProcessing\MatlabModules\Corrections_PositionCorrection\Corrections_PositionCorrection.m';
 IQs = dir([position_correction_path(1:(end-numel(myname)-2)) '*Mercury.mat']);
 table_corrections = load(IQs(end).name);

 d = Corrections_PositionCorrection_FunctionUseThis(d,table_corrections);
        
        time_diff_samples_1= time_diff_samples_1 + (datenum(d.admin.daq_settings.global.filename_prefix(7:end),'yyyymmddTHHMM')-datenum(d1.admin.daq_settings.global.filename_prefix(7:end),'yyyymmddTHHMM'))*24*60*60*1e8;
            file_start_samples(ii)=time_diff_samples_1;
            file_end_samples(ii)=time_diff_samples_1+d.event_timestamp_samples(end);
         

           % defining cuts
        zcut_min = 30;%us
        zcut_max = 300;%us
        rcut_min = 0;
        rcut_max = 25;%cm
        s1area_bound_min = 0;%2 PMT hits
        s1area_bound_max = 10^6;%include Kr 83 data and gammas

        s2area_bound_min = 100;%100 ~ 0.5 keVee
        s2area_bound_max = 10^6;% %both PMT Arrays and Kr83

         
    d.s1s2_pairing(isnan(d.s1s2_pairing))=0; %by default it returns nan when there is no pair
    d.z_drift_samples(isnan(d.z_drift_samples)) = 0.0; % get rid of NaN   
    
    s1s2_pairing = d.s1s2_pairing;   
    num_paired_pulses = sum(s1s2_pairing);
    
    events=sort(size(d.pulse_area_phe)); %The first element is the number of pulses. sometimes it's 5, sometimes it's 10
          
    s1_area_cut= inrange(d.pulse_area_phe,[s1area_bound_min,s1area_bound_max]);
    s2_area_cut= inrange(d.pulse_area_phe,[s2area_bound_min,s2area_bound_max]);
    
    s1_class=d.pulse_classification==1 & d.golden==1 & s1_area_cut ;% & squeeze(sum((d.peak_area_phe>0.3),2)>1); %at least 3 PMTs fire
    s2_class=d.pulse_classification==2 & d.golden==1 & s2_area_cut ;      
    
     s1_before_s2_cut = repmat(num_paired_pulses ==2 ,events(1),1) ...
                         & repmat(sum(s1_class,1)>=1,events(1),1) & repmat(sum(s2_class,1)==1,events(1),1)...
                         & cumsum(s1_class+s2_class)<=2 ; % ignore all S1 that happen after the s2 !
        
    s1_single_cut = s1_class & s1_before_s2_cut; %
    s2_single_cut = s2_class & s1_before_s2_cut;
    
    drift_time = d.z_drift_samples(s2_single_cut)/100; %units of us
        
        
    d.phe_bottom=d.pulse_area_phe./(1+d.top_bottom_ratio); %bottom PMT pulse area
    
    s1_phe_both = d.pulse_area_phe(s1_single_cut);
    s1_phe_bottom = d.phe_bottom(s1_single_cut);
    s1_phe_both_z=d.z_corrected_pulse_area_all_phe(s1_single_cut);
    s1_phe_both_xyz=d.xyz_corrected_pulse_area_all_phe(s1_single_cut);

    s2_phe_both = d.pulse_area_phe(s2_single_cut);
    s2_phe_bottom = d.phe_bottom(s2_single_cut);
    s2_phe_bottom_z=d.z_corrected_pulse_area_bot_phe(s2_single_cut); 
    s2_phe_bottom_xyz=d.xyz_corrected_pulse_area_bot_phe(s2_single_cut); 

    s2x = d.x_cm(s2_single_cut);
    s2y = d.y_cm(s2_single_cut);
    s2radius = (s2x.^2+s2y.^2).^(0.5);
    s2x_c = d.x_corrected(s2_single_cut);
    s2y_c = d.y_corrected(s2_single_cut);
    s2radius_c = (s2x_c.^2+s2y_c.^2).^(0.5);
    %s2x_taxi = d.x_tmplt_corrected(s2_single_cut);
    %s2y_taxi = d.y_tmplt_corrected(s2_single_cut);
    
 
        timestamp_vec_1 = sort([d.livetime_latch_samples d.livetime_end_samples]);
    livetime_sec = livetime_sec + sum(timestamp_vec_1(2:2:end) - timestamp_vec_1(1:2:end)) / 1e8;
    
    evt_cut=logical(sum(s2_single_cut));%Cut for all the events passing the single S1 & S2 cut
    event_number=d.event_number(evt_cut)';
    event_timestamp_samples=d.event_timestamp_samples(evt_cut)';
   
            %time_since_1e5_phe = d.time_since_1e5_phe(evt_cut);
            %time_since_3e5_phe = d.time_since_3e5_phe(evt_cut);
            %time_since_5e5_phe = d.time_since_5e5_phe(evt_cut);
            %time_since_8e5_phe = d.time_since_8e5_phe(evt_cut);
            livetime_samples_file(ii) = sum(timestamp_vec_1(2:2:end) - timestamp_vec_1(1:2:end));
    livetime_frac_file(ii)=livetime_samples_file(ii)/((event_timestamp_samples(end)-event_timestamp_samples(1)));         
     event_time_min=(flipdim(event_timestamp_samples+time_diff_samples_1,1))*10^(-8)/60;
     event_time_hour=event_time_min/60;
     bin_size=0.5;
     start_time=file_start_samples(ii)*10^(-8)/(60*60);
     end_time=file_end_samples(ii)*10^(-8)/(60*60);
     time_x=start_time:bin_size:end_time; %hours. HIST bin center
     count_evt_time=hist(event_time_hour,time_x);
     hist_evt_time=count_evt_time/bin_size/3600/livetime_frac_file(ii); %number per second
     sigma_hist_evt_time=sqrt(count_evt_time)/bin_size/3600/livetime_frac_file(ii);   
     total_sigma_hist_evt_time=[total_sigma_hist_evt_time, sigma_hist_evt_time];
     total_hist_evt_time=[total_hist_evt_time, hist_evt_time];
     total_time_x=[total_time_x, time_x];         
    
    clean_cut=-(s1_phe_both+s2_phe_both)' + d.full_evt_area_phe(logical(squeeze(sum(s2_single_cut,1)))) < 100 ...
        & drift_time' > 0;%more than 1/2 the area is S1+S2. And S1 before S2
    clean_cut=clean_cut';
   
    %s1_double_hit_cut= sum(d.peak_area_phe(permute(repmat(s1_single_cut,[1,1,122]),[1 3 2])) >0.3, 2) > 2; 
    
%     s1_peak_area_phe = squeeze(sum(permute(repmat(s1_single_cut,[1 1 122]),[1 3 2]).*d.peak_area_phe,1));
%     s1_multiplicity_corr = sum(s1_peak_area_phe>0.25);
%     for pp=1:2:121
% 	s1_multiplicity_corr(sum(s1_peak_area_phe(pp:(pp+1),:)>0.25)==2 & s1_multiplicity_corr==2) = 1;
%     end
%     s1_multiplicity_corr=s1_multiplicity_corr(logical(sum(s1_single_cut)))'; % get back to s1_cut size
%     s1_multiplicity_corr_cut=s1_multiplicity_corr~=1;
    
    s1_min=0;
    s1_max=50;
    
    liquid_cut= inrange(drift_time,[5 320]) & inrange(s2radius_c,[0,24.5]) & inrange(s1_phe_both_xyz,[s1_min s1_max]) ...
    & inrange(s2_phe_bottom_xyz,[150 4000]) & clean_cut & s2_phe_both>100 ; %& s1_multiplicity_corr_cut~=0

 
            %file_end_samples(ii)=time_diff_samples_1+d.event_timestamp_samples(end);

            
         s1_phe_both=[s1_phe_both ; s1_phe_both_2];
         s1_phe_both_z=[s1_phe_both_z ; s1_phe_both_z_2];
         s1_phe_both_xyz=[s1_phe_both_xyz ; s1_phe_both_xyz_2] ;        
         s2_phe_bottom=[s2_phe_bottom ; s2_phe_bottom_2] ;
         s2_phe_bottom_xyz= [s2_phe_bottom_xyz ; s2_phe_bottom_xyz_2];
         s2_phe_bottom_z= [s2_phe_bottom_z ; s2_phe_bottom_z_2];
         s2_phe_both=[s2_phe_both ; s2_phe_both_2];
         drift_time=[drift_time ; drift_time_2];
         event_timestamp_samples=[event_timestamp_samples_2 ; event_timestamp_samples+time_diff_samples_1;];
         event_number=[event_number; event_number_2];
         s2x=[s2x ; s2x_2];
         s2y=[s2y ; s2y_2];
         s2radius = (s2x.^2+s2y.^2).^(0.5);
         s2x_c = [s2x_c ; s2x_c_2];
         s2y_c = [s2y_c ; s2y_c_2];
         s2radius_c = (s2x_c.^2+s2y_c.^2).^(0.5);
%          s1_multiplicity_corr_cut=[s1_multiplicity_corr_cut ; s1_multiplicity_corr_cut_2];
         clean_cut=[clean_cut ; clean_cut_2 ];
         liquid_cut=[liquid_cut ; liquid_cut_2];
         %kr_cut=[kr_cut ; kr_cut_2 ];
         
         s1_phe_both_2=s1_phe_both;
         s1_phe_both_z_2=s1_phe_both_z ;
         s1_phe_both_xyz_2=s1_phe_both_xyz ;        
         s2_phe_bottom_2=s2_phe_bottom ;
         s2_phe_bottom_xyz_2= s2_phe_bottom_xyz ;
         s2_phe_bottom_z_2= s2_phe_bottom_z ;
         s2_phe_both_2=s2_phe_both ;
         drift_time_2=drift_time ;
         event_timestamp_samples_2=event_timestamp_samples ;
         event_number_2=event_number;
         s2x_2=s2x ;
         s2y_2=s2y ;
         s2x_c_2 = s2x_c ;
         s2y_c_2 = s2y_c ;
%          s1_multiplicity_corr_cut_2=s1_multiplicity_corr_cut;
         clean_cut_2=clean_cut;
         %kr_cut_2=kr_cut;
        
  save('CH3T_Sep2014_TAGCB_1_7_1_085_2_poscorrected','drift_time','s2radius','s1_phe_both_xyz','s2_phe_bottom_xyz','event_timestamp_samples'...
      ,'s2radius_c','s1_phe_both','s2_phe_bottom','event_number','s2x','s2y','s2x_c','s2y_c','livetime_sec','clean_cut'...
      ,'liquid_cut','s2_phe_both','total_hist_evt_time','total_time_x','livetime_frac_file','total_sigma_hist_evt_time');              

  fprintf(strcat(num2str(ii),' out of ', num2str(length(folders_to_load)),' data sets loaded') );
  
  toc;
  
end
    

    

fprintf('Finished');