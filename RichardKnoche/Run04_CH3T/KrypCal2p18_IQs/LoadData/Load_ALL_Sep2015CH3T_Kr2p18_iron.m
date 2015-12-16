%Find IQs that that are Kr-83 data

dir_path='/scratch3/LUXData/Run04_Sep2015_CH3T_2p18/';

potential_folders_to_load={' '};
folders_to_load={' '};

% Check that we are only loading folders with lux10_
% and that they contain at least 10 MB of data
list=dir(dir_path);

    k=1;
for ii=1:length(list)
    
    if length(list(ii).name) > 6
        
        if  strcmp(list(ii).name(1:6),'lux10_') 
            
            folder_contents=dir(strcat(dir_path,'/',list(ii).name));
            dir_size=sum([folder_contents.bytes]);
            
            if dir_size > 10^7 % 10 MB
             
                potential_folders_to_load{k}=list(ii).name;
                k=k+1;
            
            end
                
                
        end
    end
    
end


% Once we determined good folders to load, check their multiplicity
% Keep higher CP unless duplicated folder with a higher CP number contains no data
folder_names_char=char(potential_folders_to_load);
folder_name_no_cp=cellstr([folder_names_char(:,1:19)]); %convert the char back to a cell array of chars

k=1;
good_index=1;
while good_index<=length(potential_folders_to_load)
        
        
       folder_multiplicity=sum(strcmp(potential_folders_to_load{good_index}(1:19),folder_name_no_cp ));
       good_index=good_index+folder_multiplicity-1;
            
       folders_to_load{k}=potential_folders_to_load{good_index};
       
       k=k+1;     
       good_index=good_index+1;
    
end


rqs_to_load = {'pulse_area_phe','event_timestamp_samples'...
   ,'pulse_classification' ...
   ,'z_drift_samples' , 's1s2_pairing'...
   ,'z_corrected_pulse_area_all_phe','xyz_corrected_pulse_area_all_phe'...
   ,'xyz_corrected_pulse_area_bot_phe'...
   ,'top_bottom_ratio','x_cm','y_cm'...
   ,'z_corrected_pulse_area_bot_phe','xyz_corrected_pulse_area_bot_phe','full_evt_area_phe'...
   ,'x_cm_del','y_cm_del'};
    
ii=1;

    path=sprintf('%s%s',dir_path,folders_to_load{ii});
    
    d = LUXLoadMultipleRQMs_framework(path,rqs_to_load);

     
           % defining cuts
        zcut_min = 30;%us
        zcut_max = 300;%us
        rcut_min = 0;
        rcut_max = 25;%cm
        s1area_bound_min = 0;%2 PMT hits
        s1area_bound_max = 10^7;%include Kr 83 data and gammas

        s2area_bound_min = 0;%100 ~ 0.5 keVee
        s2area_bound_max = 10^7;% %both PMT Arrays and Kr83

         
    d.z_drift_samples(isnan(d.z_drift_samples)) = 0.0; % get rid of NaN      

    events=sort(size(d.pulse_area_phe)); %The first element is the number of pulses. sometimes it's 5, sometimes it's 10
          
  
cut_pulse_s1 = d.pulse_classification == 1;
cut_pulse_s2 = d.pulse_classification == 2;
cut_s2_with_threshold = d.pulse_area_phe.*cut_pulse_s2 > 100; % subset of cut_pulse_s2
cut_legit_s2_in_legit_event = d.s1s2_pairing.*cut_s2_with_threshold; % this should be used as s2 classification cuts
cut_golden_event = sum(cut_legit_s2_in_legit_event) == 1;
cut_s2_in_golden_events = logical(repmat(cut_golden_event,[10,1]).*cut_legit_s2_in_legit_event);
cut_s1_in_golden_events = logical(repmat(cut_golden_event,[10,1]).*cut_pulse_s1.*d.s1s2_pairing);
% select Kr83 events with cut on S2
cut_s2_area = inrange(d.pulse_area_phe, [s2area_bound_min, s2area_bound_max]);
cut_s1_area = inrange(d.pulse_area_phe, [s1area_bound_min, s1area_bound_max]);
cut_s2_for = cut_s2_in_golden_events.*cut_s2_area;
cut_s1_for = cut_s1_in_golden_events.*cut_s1_area;
cut_selected_events = sum(cut_s2_for) == 1 & sum(cut_s1_for) == 1;
s1_single_cut = logical(repmat(cut_selected_events,[10,1]).*cut_s1_in_golden_events);
s2_single_cut = logical(repmat(cut_selected_events,[10,1]).*cut_s2_in_golden_events);
    
    drift_time = d.z_drift_samples(s2_single_cut)/100; %units of us
        
        
    d.phe_bottom=d.pulse_area_phe./(1+d.top_bottom_ratio); %bottom PMT pulse area
    
    s1_phe_both = d.pulse_area_phe(s1_single_cut);
    s1_phe_bottom = d.phe_bottom(s1_single_cut);
    s1_phe_both_z=d.z_corrected_pulse_area_all_phe(s1_single_cut);
    s1_phe_both_xyz=d.xyz_corrected_pulse_area_all_phe(s1_single_cut);
    s2_phe_both = d.pulse_area_phe(s2_single_cut);
    s2_phe_both_z=d.z_corrected_pulse_area_all_phe(s2_single_cut); 
    s2_phe_both_xyz=d.xyz_corrected_pulse_area_all_phe(s2_single_cut); 
    s2_phe_bottom = d.phe_bottom(s2_single_cut);
    s2_phe_bottom_z=d.z_corrected_pulse_area_bot_phe(s2_single_cut); 
    s2_phe_bottom_xyz=d.xyz_corrected_pulse_area_bot_phe(s2_single_cut); 

    s2x = d.x_cm(s2_single_cut);
    s2y = d.y_cm(s2_single_cut);
    s2radius = (s2x.^2+s2y.^2).^(0.5);
%     s2x_c = d.x_corrected(s2_single_cut);
%     s2y_c = d.y_corrected(s2_single_cut);
%     s2radius_c = (s2x_c.^2+s2y_c.^2).^(0.5);
    %s2x_taxi = d.x_tmplt_corrected(s2_single_cut);
    %s2y_taxi = d.y_tmplt_corrected(s2_single_cut);
    s2x_del=d.x_cm_del(s2_single_cut);
    s2y_del=d.y_cm_del(s2_single_cut);
    s2radius_del=(s2x_del.^2+s2y_del.^2).^(0.5);    
 
    
     %Based on Rachel's work here: http://teacher.pas.rochester.edu:8080/wiki/pub/Lux/LuxDB00000361/BadAreaCut_summary.pdf
    good_area= (s1_phe_both+s2_phe_both)';
    bad_area=(d.full_evt_area_phe(logical(squeeze(sum(s2_single_cut,1)))) - good_area);
    clean_cut= ( bad_area < 80 & good_area < 253 ) | ( bad_area < (80.*( (10^-2.4).*good_area ).^0.4668) & good_area >= 253);
    clean_cut=clean_cut'; 
   
      
    s1_min=0;
    s1_max=110;
    
    max_r=25;
    liquid_cut= clean_cut & inrange(drift_time,[30 330]) & inrange(s2radius,[0,max_r]) & inrange(s1_phe_both_xyz,[s1_min s1_max]) ...
    & (log10(s1_phe_both_xyz)+0.5*log10(s2_phe_both_xyz)<4.2) & s2_phe_both>100 ; %used to be drift time 100 to 250

        timestamp_vec_1 = sort([d.livetime_latch_samples d.livetime_end_samples]);
    livetime_sec = sum(timestamp_vec_1(2:2:end) - timestamp_vec_1(1:2:end)) / 1e8;
    
    evt_cut=logical(sum(s2_single_cut));%Cut for all the events passing the single S1 & S2 cut
    event_number=d.event_number(evt_cut)';
    event_number=event_number(liquid_cut);
    
    event_timestamp_samples=d.event_timestamp_samples(evt_cut)';
        event_timestamp_samples=event_timestamp_samples(liquid_cut)';

    
        livetime_samples_file(1) = sum(timestamp_vec_1(2:2:end) - timestamp_vec_1(1:2:end));
    livetime_frac_file(1)=livetime_samples_file(1)/((event_timestamp_samples(end)-event_timestamp_samples(1)));         
     event_time_min=flipdim(event_timestamp_samples,1)*10^(-8)/60;
     event_time_hour=event_time_min/60;
     bin_size=0.5;
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
            
    
         
         s1_phe_both_2=s1_phe_both;
         s1_phe_both_z_2=s1_phe_both_z ;
         s1_phe_both_xyz_2=s1_phe_both_xyz ;        
         s2_phe_bottom_2=s2_phe_bottom ;
         s2_phe_bottom_xyz_2= s2_phe_bottom_xyz ;
         s2_phe_bottom_z_2= s2_phe_bottom_z ;
         s2_phe_both_2=s2_phe_both ;
         s2_phe_both_z_2=s2_phe_both_z;
         s2_phe_both_xyz_2=s2_phe_both_xyz;
         s1_phe_bottom_2 =  s1_phe_bottom;
         s2x_del_2 = s2x_del;
         s2y_del_2 = s2y_del;
         s2radius_del_2=s2radius_del;   
         
         drift_time_2=drift_time ;
         event_timestamp_samples_2=event_timestamp_samples ;
         event_number_2=event_number;
         s2x_2=s2x ;
         s2y_2=s2y ;
%          s2x_c_2 = s2x_c ;
%          s2y_c_2 = s2y_c ;
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
        
    path=sprintf('%s%s',dir_path,folders_to_load{ii});
        d = LUXLoadMultipleRQMs_framework(path,rqs_to_load);

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
        s2area_bound_max = 10^7;% %both PMT Arrays and Kr83

         
      d.z_drift_samples(isnan(d.z_drift_samples)) = 0.0; % get rid of NaN      

    events=sort(size(d.pulse_area_phe)); %The first element is the number of pulses. sometimes it's 5, sometimes it's 10
          
  
cut_pulse_s1 = d.pulse_classification == 1;
cut_pulse_s2 = d.pulse_classification == 2;
cut_s2_with_threshold = d.pulse_area_phe.*cut_pulse_s2 > 100; % subset of cut_pulse_s2
cut_legit_s2_in_legit_event = d.s1s2_pairing.*cut_s2_with_threshold; % this should be used as s2 classification cuts
cut_golden_event = sum(cut_legit_s2_in_legit_event) == 1;
cut_s2_in_golden_events = logical(repmat(cut_golden_event,[10,1]).*cut_legit_s2_in_legit_event);
cut_s1_in_golden_events = logical(repmat(cut_golden_event,[10,1]).*cut_pulse_s1.*d.s1s2_pairing);
% select Kr83 events with cut on S2
cut_s2_area = inrange(d.pulse_area_phe, [s2area_bound_min, s2area_bound_max]);
cut_s1_area = inrange(d.pulse_area_phe, [s1area_bound_min, s1area_bound_max]);
cut_s2_for = cut_s2_in_golden_events.*cut_s2_area;
cut_s1_for = cut_s1_in_golden_events.*cut_s1_area;
cut_selected_events = sum(cut_s2_for) == 1 & sum(cut_s1_for) == 1;
s1_single_cut = logical(repmat(cut_selected_events,[10,1]).*cut_s1_in_golden_events);
s2_single_cut = logical(repmat(cut_selected_events,[10,1]).*cut_s2_in_golden_events);
    
    drift_time = d.z_drift_samples(s2_single_cut)/100; %units of us
        
        
    d.phe_bottom=d.pulse_area_phe./(1+d.top_bottom_ratio); %bottom PMT pulse area
    
    s1_phe_both = d.pulse_area_phe(s1_single_cut);
    s1_phe_bottom = d.phe_bottom(s1_single_cut);
    s1_phe_both_z=d.z_corrected_pulse_area_all_phe(s1_single_cut);
    s1_phe_both_xyz=d.xyz_corrected_pulse_area_all_phe(s1_single_cut);

    s2_phe_both = d.pulse_area_phe(s2_single_cut);
    s2_phe_both_z=d.z_corrected_pulse_area_all_phe(s2_single_cut); 
    s2_phe_both_xyz=d.xyz_corrected_pulse_area_all_phe(s2_single_cut); 
    s2_phe_bottom = d.phe_bottom(s2_single_cut);
    s2_phe_bottom_z=d.z_corrected_pulse_area_bot_phe(s2_single_cut); 
    s2_phe_bottom_xyz=d.xyz_corrected_pulse_area_bot_phe(s2_single_cut); 

    s2x = d.x_cm(s2_single_cut);
    s2y = d.y_cm(s2_single_cut);
    s2radius = (s2x.^2+s2y.^2).^(0.5);
%     s2x_c = d.x_corrected(s2_single_cut);
%     s2y_c = d.y_corrected(s2_single_cut);
%     s2radius_c = (s2x_c.^2+s2y_c.^2).^(0.5);
    s2x_del=d.x_cm_del(s2_single_cut);
    s2y_del=d.y_cm_del(s2_single_cut);
    s2radius_del=(s2x_del.^2+s2y_del.^2).^(0.5);       
    %s2x_taxi = d.x_tmplt_corrected(s2_single_cut);
    %s2y_taxi = d.y_tmplt_corrected(s2_single_cut);
    
     
     %Based on Rachel's work here: http://teacher.pas.rochester.edu:8080/wiki/pub/Lux/LuxDB00000361/BadAreaCut_summary.pdf
    good_area= (s1_phe_both+s2_phe_both)';
    bad_area=(d.full_evt_area_phe(logical(squeeze(sum(s2_single_cut,1)))) - good_area);
    clean_cut= ( bad_area < 80 & good_area < 253 ) | ( bad_area < (80.*( (10^-2.4).*good_area ).^0.4668) & good_area >= 253);
    clean_cut=clean_cut'; 
   
      
    s1_min=0;
    s1_max=110;
    
    max_r=25;
    liquid_cut= clean_cut & inrange(drift_time,[30 330]) & inrange(s2radius,[0,max_r]) & inrange(s1_phe_both_xyz,[s1_min s1_max]) ...
    & (log10(s1_phe_both_xyz)+0.5*log10(s2_phe_both_xyz)<4.2) & s2_phe_both>100 ; %used to be drift time 100 to 250

        timestamp_vec_1 = sort([d.livetime_latch_samples d.livetime_end_samples]);
    livetime_sec = sum(timestamp_vec_1(2:2:end) - timestamp_vec_1(1:2:end)) / 1e8;
    
    evt_cut=logical(sum(s2_single_cut));%Cut for all the events passing the single S1 & S2 cut
    event_number=d.event_number(evt_cut)';
    event_number=event_number(liquid_cut);
    
    event_timestamp_samples=d.event_timestamp_samples(evt_cut)';
        event_timestamp_samples=event_timestamp_samples(liquid_cut)';

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
  
    
    %s1_double_hit_cut= sum(d.peak_area_phe(permute(repmat(s1_single_cut,[1,1,122]),[1 3 2])) >0.3, 2) > 2; 
    
%     s1_peak_area_phe = squeeze(sum(permute(repmat(s1_single_cut,[1 1 122]),[1 3 2]).*d.peak_area_phe,1));
%     s1_multiplicity_corr = sum(s1_peak_area_phe>0.25);
%     for pp=1:2:121
% 	s1_multiplicity_corr(sum(s1_peak_area_phe(pp:(pp+1),:)>0.25)==2 & s1_multiplicity_corr==2) = 1;
%     end
%     s1_multiplicity_corr=s1_multiplicity_corr(logical(sum(s1_single_cut)))'; % get back to s1_cut size
%     s1_multiplicity_corr_cut=s1_multiplicity_corr~=1;
    
   
 
            %file_end_samples(ii)=time_diff_samples_1+d.event_timestamp_samples(end);

            
         s1_phe_both=[s1_phe_both ; s1_phe_both_2];
         s1_phe_both_z=[s1_phe_both_z ; s1_phe_both_z_2];
         s1_phe_both_xyz=[s1_phe_both_xyz ; s1_phe_both_xyz_2] ;        
         s2_phe_bottom=[s2_phe_bottom ; s2_phe_bottom_2] ;
         s2_phe_bottom_xyz= [s2_phe_bottom_xyz ; s2_phe_bottom_xyz_2];
         s2_phe_bottom_z= [s2_phe_bottom_z ; s2_phe_bottom_z_2];
         s2_phe_both=[s2_phe_both ; s2_phe_both_2];
         s2_phe_both_z=[s2_phe_both_z ; s2_phe_both_z_2];
         s2_phe_both_xyz=[s2_phe_both_xyz ; s2_phe_both_xyz_2]; 
         drift_time=[drift_time ; drift_time_2];
         event_timestamp_samples=[event_timestamp_samples_2 , event_timestamp_samples+time_diff_samples_1;];
         event_number=[event_number; event_number_2];
         s2x=[s2x ; s2x_2];
         s2y=[s2y ; s2y_2];
         s2radius = (s2x.^2+s2y.^2).^(0.5);
%          s2x_c = [s2x_c ; s2x_c_2];
%          s2y_c = [s2y_c ; s2y_c_2];
%          s2radius_c = (s2x_c.^2+s2y_c.^2).^(0.5);
%          s1_multiplicity_corr_cut=[s1_multiplicity_corr_cut ; s1_multiplicity_corr_cut_2];
         clean_cut=[clean_cut ; clean_cut_2 ];
         liquid_cut=[liquid_cut ; liquid_cut_2];
         s1_phe_bottom=[s1_phe_bottom ; s1_phe_bottom_2];
         s2x_del=[s2x_del ; s2x_del_2];
         s2y_del=[s2y_del ; s2y_del_2];
         s2radius_del=[s2radius_del ; s2radius_del_2];
         %kr_cut=[kr_cut ; kr_cut_2 ];
         
         s1_phe_both_2=s1_phe_both;
         s1_phe_both_z_2=s1_phe_both_z ;
         s1_phe_both_xyz_2=s1_phe_both_xyz ;        
         s2_phe_bottom_2=s2_phe_bottom ;
         s2_phe_bottom_xyz_2= s2_phe_bottom_xyz ;
         s2_phe_bottom_z_2= s2_phe_bottom_z ;
         s2_phe_both_2=s2_phe_both ;
         s2_phe_both_z_2=s2_phe_both_z;
         s2_phe_both_xyz_2=s2_phe_both_xyz;
         s1_phe_bottom_2 =  s1_phe_bottom;
         s2x_del_2=s2x_del;
         s2y_del_2=s2y_del;
         s2radius_del_2=s2radius_del;
         
         drift_time_2=drift_time ;
         event_timestamp_samples_2=event_timestamp_samples ;
         event_number_2=event_number;
         s2x_2=s2x ;
         s2y_2=s2y ;
%          s2x_c_2 = s2x_c ;
%          s2y_c_2 = s2y_c ;
%          s1_multiplicity_corr_cut_2=s1_multiplicity_corr_cut;
         clean_cut_2=clean_cut;
         liquid_cut_2=liquid_cut;
         %kr_cut_2=kr_cut;
        
  save('CH3T_Sep2015_Kr2p18','drift_time','s2radius','s1_phe_both_xyz','s2_phe_bottom_xyz','event_timestamp_samples'...
      ,'s1_phe_both','s1_phe_bottom','s2_phe_bottom','event_number','s2x','s2y','livetime_sec','clean_cut'...
      ,'liquid_cut','s2_phe_both','total_hist_evt_time','total_time_x','livetime_frac_file','total_sigma_hist_evt_time','s2_phe_both_z','s2_phe_both_xyz'...
      ,'s2x_del','s2y_del','s2radius_del'); 
  
  fprintf(strcat(num2str(ii),' out of ', num2str(length(folders_to_load)),' data sets loaded') );
  
  toc;
  
end
    

    

fprintf('Finished');