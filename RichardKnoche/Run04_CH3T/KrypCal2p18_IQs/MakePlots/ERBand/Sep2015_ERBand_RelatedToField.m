%First measure the field in corrected space

path='C:\Users\Richard\Desktop\LUX Analysis\lux10_20150929T1905_cp18197\';

rqs_to_load = {'pulse_area_phe','event_timestamp_samples'...
   ,'pulse_classification' ...
   ,'z_drift_samples' , 's1s2_pairing'...
   ,'top_bottom_ratio','x_cm','y_cm','x_cm_del','y_cm_del'...
   ,'full_evt_area_phe',...
   'event_number','chi2','prompt_fraction','aft_t1_samples','pulse_start_samples',...
   'pulse_end_samples','top_bottom_asymmetry','aft_t0_samples','aft_t2_samples',...
   'full_evt_area_phe','admin','Kr83fit_s1a_area_phe','Kr83fit_s1b_area_phe',...
   'hft_t50r_samples','hft_t50l_samples','Kr83fit_dt_samples'};

d = LUXLoadMultipleRQMs_framework(path,rqs_to_load);  

file_id='lux10_20150929T1905_cp18197'; file_id_cp='lux10_20150929T1905_cp18197'; delay_min=0;

%%  Defining cuts and variables
       
grid_size=3; %cm XY plane
rcut_min = 0;
rcut_max = 25;%cm
s1area_bound_min = 50;%100
s1area_bound_max = 650;%600
s2area_bound_min = 1000;%200
s2area_bound_max = 32000;%30000
min_evts_per_bin = 300;
max_num_bins = 65;

S1xybinsize = grid_size;%cm
S1_xbin_min = -25;
S1_xbin_max = 25;
S1_ybin_min = -25;
S1_ybin_max = 25;


SE_xbin_min = -25;
SE_xbin_max = 25;
SE_ybin_min = -25;
SE_ybin_max = 25;
% SExybinsize = grid_size; Defined later, based on number of SE events

S2xybinsize = grid_size;%cm
S2_xbin_min = -25;
S2_xbin_max = 25;
S2_ybin_min = -25;
S2_ybin_max = 25;



%% Event selection with golden and Pulse Classification
  
    d.z_drift_samples(isnan(d.z_drift_samples)) = 0.0; % get rid of NaN        
    events=sort(size(d.pulse_area_phe)); %The first element is the number of pulses. sometimes it's 5, sometimes it's 10
          
    s1_area_cut= inrange(d.pulse_area_phe,[s1area_bound_min,s1area_bound_max]);
    s2_area_cut= inrange(d.pulse_area_phe,[s2area_bound_min,s2area_bound_max]);
    
    s1_class=(d.pulse_classification==1 )& s1_area_cut ; %area cut for Kr events
    s2_class=(d.pulse_classification==2) & s2_area_cut ;   
    s4_class=(d.pulse_classification==4) ;
  
events=sort(size(d.pulse_area_phe)); %The first element is the number of pulses. sometimes it's 5, sometimes it's 10
cut_pulse_s1 = d.pulse_classification == 1 | d.pulse_classification == 9;
cut_pulse_s2 = d.pulse_classification == 2;
cut_s2_with_threshold = d.pulse_area_phe.*cut_pulse_s2 > 100; % subset of cut_pulse_s2
cut_legit_s2_in_legit_event = d.s1s2_pairing.*cut_s2_with_threshold; % this should be used as s2 classification cuts
cut_golden_event = sum(cut_legit_s2_in_legit_event) == 1; %defines golden events to be events which have one and only one paired S2 above the threshold of 100 phe - there can be multiple S1s still
cut_s2_in_golden_events = logical(repmat(cut_golden_event,[10,1]).*cut_legit_s2_in_legit_event); %Selects S2 that is in a golden event
cut_s1_in_golden_events = logical(repmat(cut_golden_event,[10,1]).*cut_pulse_s1.*d.s1s2_pairing); %Selects first S1 that is in a golden event
% select Kr83 events with cut on S2
cut_s2_area = inrange(d.pulse_area_phe, [s2area_bound_min, s2area_bound_max]);
cut_s1_area = inrange(d.pulse_area_phe, [s1area_bound_min, s1area_bound_max]);
cut_s2_for = cut_s2_in_golden_events.*cut_s2_area; %Selects S2 that is in a golden event and in Kr area bounds
cut_s1_for = cut_s1_in_golden_events.*cut_s1_area; %Selects first S1 that is in a golden event and in Kr area bounds
cut_selected_events = sum(cut_s2_for) == 1 & sum(cut_s1_for) == 1 & sum(d.pulse_classification==1)==1; %Requires that "good" golden events have only one S1, that the S1 be within area bounds, and the S2 be within area bounds
%Note sum(cut_s1_for) == 1 above only requires that the first of the S1 in an event be within area bounds, since the S1S2pairing part of cut_s1_in_golden_events is 0 for all subsequent S1s in the events
s1_single_cut = logical(repmat(cut_selected_events,[10,1]).*cut_s1_in_golden_events);
s2_single_cut = logical(repmat(cut_selected_events,[10,1]).*cut_s2_in_golden_events);

 
%% Set up cuts and variables
    s1ab_cut=logical(sum(s1_single_cut(:,:),1));
    s1a_phe_both=d.Kr83fit_s1a_area_phe(s1ab_cut);
    s1b_phe_both=d.Kr83fit_s1b_area_phe(s1ab_cut);
    s1ab_timing=d.Kr83fit_dt_samples(s1ab_cut);
    s1ab_x=d.x_cm_del(s2_single_cut);
    s1ab_y=d.y_cm_del(s2_single_cut);
    s1ab_z=d.z_drift_samples(s2_single_cut)/100;
    s1ab_radius=sqrt(s1ab_x.^2+s1ab_y.^2);
    s1ab_timing_cut=s1ab_timing>13 & s1b_phe_both>30 & s1ab_radius.'<25;

    
    s1ab_x=s1ab_x(s1ab_timing_cut);
    s1ab_y=s1ab_y(s1ab_timing_cut);
    s1ab_z=s1ab_z(s1ab_timing_cut);
    s1a_phe_both=s1a_phe_both(s1ab_timing_cut);
    s1b_phe_both=s1b_phe_both(s1ab_timing_cut);
    s1ab_radius=s1ab_radius(s1ab_timing_cut);
       

%% S1a/b R^2 v Z map 
s1ab_r2bin_min=0;
s1ab_r2bin_max=25*25;
s1ab_zbin_min=10;
s1ab_zbin_max=330;

s1ab_r2z_binnum=(length(s1ab_z))/3000; %10k based on CH3T stats for ER band
s1ab_r2_binnum=floor(sqrt(s1ab_r2z_binnum));
s1ab_z_binnum=floor(sqrt(s1ab_r2z_binnum));

s1ab_r2_binsize=(s1ab_r2bin_max-s1ab_r2bin_min)/s1ab_r2_binnum;
s1ab_z_binsize=(s1ab_zbin_max-s1ab_zbin_min)/s1ab_z_binnum;

s1ab_r2bins=s1ab_r2bin_min+s1ab_r2_binsize/2:s1ab_r2_binsize:s1ab_r2bin_max;
s1ab_zbins=s1ab_zbin_min+s1ab_z_binsize/2:s1ab_z_binsize:s1ab_zbin_max; 

s1a_r2z_mean=zeros(length(s1ab_zbins),length(s1ab_r2bins));
s1b_r2z_mean=zeros(length(s1ab_zbins),length(s1ab_r2bins));
s1a_r2z_mean_err=zeros(length(s1ab_zbins),length(s1ab_r2bins));
s1b_r2z_mean_err=zeros(length(s1ab_zbins),length(s1ab_r2bins));
field_map_2014=zeros(length(s1ab_zbins),length(s1ab_r2bins));
field_map_2015=zeros(length(s1ab_zbins),length(s1ab_r2bins));
field_map_error=zeros(length(s1ab_zbins),length(s1ab_r2bins));

r2_counter=1;
 for r2_bin_max=s1ab_r2bin_min+s1ab_r2_binsize:s1ab_r2_binsize:s1ab_r2bin_max;
     z_counter=1;
        for z_bin_max=s1ab_zbin_min+s1ab_z_binsize:s1ab_z_binsize:s1ab_zbin_max;          
           
          bin_cut = s1a_phe_both>0 & s1b_phe_both>0 & inrange(s1ab_radius.^2,[r2_bin_max-s1ab_r2_binsize,r2_bin_max]).' & inrange(s1ab_z,[z_bin_max-s1ab_z_binsize,z_bin_max]).'; 
            if length(s1a_phe_both(bin_cut))>100;
                clear s1a_confint  s1b_confint
                s1a_r2z_fit=fit([0:2:400].',hist(s1a_phe_both(bin_cut & inrange(s1a_phe_both,[0 400])),[0:2:400]).','gauss1');
                s1a_r2z_mean(z_counter,r2_counter)=s1a_r2z_fit.b1;
                s1a_confint=confint(s1a_r2z_fit,0.68);
%               s1a_r2z_mean_err(z_counter,r2_counter)=s1a_r2z_fit.c1/sqrt(2)/length(s1a_phe_both_z(bin_cut.' & inrange(s1a_phe_both_z,[0 400])));
                s1a_r2z_mean_err(z_counter,r2_counter)=abs(s1a_confint(1,2)-s1a_r2z_fit.b1);
                s1b_r2z_fit=fit([0:2:300].',hist(s1b_phe_both(bin_cut & inrange(s1b_phe_both,[0 300])),[0:2:300]).','gauss1');
                s1b_r2z_mean(z_counter,r2_counter)=s1b_r2z_fit.b1;
                s1b_confint=confint(s1b_r2z_fit,0.68);
%                 s1b_xy_mean_err(y_count,x_count)=s1b_z_fit.c1/sqrt(2)/length(s1b_phe_both_z(bin_cut.' & inrange(s1b_phe_both_z,[0 300])));
                s1b_r2z_mean_err(z_counter,r2_counter)=abs(s1b_confint(1,2)-s1b_r2z_fit.b1);
            else
                s1a_r2z_mean(z_counter,r2_counter)=0;
                s1a_r2z_mean_err(z_counter,r2_counter)=0;
                s1b_r2z_mean(z_counter,r2_counter)=0;
                s1b_r2z_mean_err(z_counter,r2_counter)=0; 
                
            end
           z_counter=z_counter+1;     
        end
        r2_counter=r2_counter+1;
 end
     s1ab_r2z_mean=s1a_r2z_mean./s1b_r2z_mean;
     s1ab_r2z_mean_err=sqrt( (s1b_r2z_mean_err.*s1a_r2z_mean./(s1b_r2z_mean.^2)).^2 + (s1a_r2z_mean_err./s1b_r2z_mean).^2);
   
    %Plot s1a R2Z means
    color_range_max=max(max(s1a_r2z_mean));
    color_range_min=min(min(s1a_r2z_mean(s1a_r2z_mean>0)));
    vc_step=(color_range_max-color_range_min)/50;
    vc=color_range_min:vc_step:color_range_max;
    

    s1a_r2z_fig = figure;
    contourf(s1ab_r2bins,s1ab_zbins,s1a_r2z_mean,vc,'LineColor','none');
    xlabel('R^2 (cm^2)','fontsize', 18);
    ylabel('Z (uSec)','fontsize', 18);
    title(strcat(file_id_cp, '. S1a Mean vs. R^2 v Z.'),'fontsize',16,'Interpreter','none');
    caxis([color_range_min color_range_max])
    g=colorbar;
    ylabel(g,'S1a Mean (phe)','FontSize',16)
    set(g,'FontSize',18)
    hold on
    axis([0 625 0 330])
    myfigview(16);
    set(gca,'Ydir','reverse')

   
    %Plot s1a R2Z means   
    color_range_max=max(max(s1b_r2z_mean));
    color_range_min=min(min(s1b_r2z_mean(s1b_r2z_mean>0)));
    vc_step=(color_range_max-color_range_min)/50;
    vc=color_range_min:vc_step:color_range_max;
    
    s1b_r2z_fig = figure;
    contourf(s1ab_r2bins,s1ab_zbins,s1b_r2z_mean,vc,'LineColor','none');
    xlabel('R^2 (cm^2)','fontsize', 18);
    ylabel('Z (uSec)','fontsize', 18);
    title(strcat(file_id_cp, '. S1b Mean vs. R^2 v Z.'),'fontsize',16,'Interpreter','none');
    caxis([color_range_min color_range_max])
    g=colorbar;
    ylabel(g,'S1b Mean (phe)','FontSize',16)
    set(g,'FontSize',18)
    hold on
    axis([0 625 0 330])
    myfigview(16);
    set(gca,'Ydir','reverse')
 
    
    %s1a/b R2Z Ratio
    color_range_max=max(max(s1ab_r2z_mean));
    color_range_min=min(min(s1ab_r2z_mean(s1ab_r2z_mean>0)));
    vc_step=(color_range_max-color_range_min)/50;
    vc=color_range_min:vc_step:color_range_max;
    
    s1ab_r2z_fig = figure;
    contourf(s1ab_r2bins,s1ab_zbins,s1ab_r2z_mean,vc,'LineColor','none');
    xlabel('R^2 (cm^2)','fontsize', 18);
    ylabel('Z (uSec)','fontsize', 18);
    title(strcat(file_id_cp, '. S1a/S1b Mean vs. R^2 v Z.'),'fontsize',16,'Interpreter','none');
    caxis([color_range_min color_range_max])
    g=colorbar;
    ylabel(g,'S1a/b Ratio','FontSize',16)
    set(g,'FontSize',18)
    hold on
    axis([0 625 0 330])
    myfigview(16);
    set(gca,'Ydir','reverse')
 
 

%% Convert S1a/S1b R^2 v Z to field measurement based on Scott's Chi by Eye

load('C:\Program Files\MATLAB\R2012a\bin\LUXCode\Scratch\RichardKnoche\KrypCal\KrypCal_2p19_Code\2015_Field_to_S1aS1b.mat');
[field_map_2015, field_map_2015_staterr]=polyval(S1aS1b_fit_2015,s1ab_r2z_mean,S1aS1b_fit_2015_error);

load('C:\Program Files\MATLAB\R2012a\bin\LUXCode\Scratch\RichardKnoche\KrypCal\KrypCal_2p19_Code\2014_Field_to_S1aS1b.mat');
[field_map_2014, field_map_2014_staterr]=polyval(S1aS1b_fit_2014,s1ab_r2z_mean,S1aS1b_fit_2014_error);

field_map_syserr=((field_map_2015-(field_map_2014+field_map_2015)/2).^2 + (field_map_2014-(field_map_2014+field_map_2015)/2).^2).^(1/2);
mean_field_map=(field_map_2014+field_map_2015)/2;
mean_field_map_totalerr=sqrt(field_map_syserr.^2+((field_map_2014_staterr+field_map_2015_staterr)/2).^2);

    %Mean R^2vZ Field Map
    color_range_max=max(max(mean_field_map));
    color_range_min=min(min(mean_field_map(mean_field_map>0)));
    vc_step=(color_range_max-color_range_min)/50;
    vc=color_range_min:vc_step:color_range_max;
    
    field_r2z_fig = figure;
    contourf(s1ab_r2bins,s1ab_zbins,mean_field_map,vc,'LineColor','none');
    xlabel('R^2 (cm^2)','fontsize', 18);
    ylabel('Z (uSec)','fontsize', 18);
    title(strcat(file_id_cp, '. Electric Field vs. R^2 v Z.'),'fontsize',16,'Interpreter','none');
    caxis([color_range_min color_range_max])
    g=colorbar;
    ylabel(g,'Electric Field (V/cm)','FontSize',16)
    set(g,'FontSize',18)
    hold on
    axis([0 625 0 330])
    myfigview(16);
    set(gca,'Ydir','reverse')

    
    %Mean R^2vZ Field Map Error
    color_range_max=max(max(mean_field_map_totalerr));
    color_range_min=min(min(mean_field_map_totalerr(mean_field_map_totalerr>0)));
    vc_step=(color_range_max-color_range_min)/50;
    vc=color_range_min:vc_step:color_range_max;
    
    field_error_r2z_fig = figure;
    contourf(s1ab_r2bins,s1ab_zbins,mean_field_map_totalerr,vc,'LineColor','none');
    xlabel('R^2 (cm^2)','fontsize', 18);
    ylabel('Z (uSec)','fontsize', 18);
    title(strcat(file_id_cp, '. Electric Field Error vs. R^2 v Z.'),'fontsize',16,'Interpreter','none');
    caxis([color_range_min color_range_max])
    g=colorbar;
    ylabel(g,'Electric Field Error (V/cm)','FontSize',16)
    set(g,'FontSize',18)
    hold on
    axis([0 625 0 330])
    myfigview(16);
    set(gca,'Ydir','reverse')

  
    %Mean R^2vZ Field Map Fractional Error
    color_range_max=max(max(mean_field_map_totalerr./mean_field_map));
    color_range_min=min(min(mean_field_map_totalerr(mean_field_map_totalerr./mean_field_map>0)./mean_field_map(mean_field_map_totalerr./mean_field_map>0)));
    vc_step=(color_range_max-color_range_min)/50;
    vc=color_range_min:vc_step:color_range_max;
    
    field_frac_error_r2z_fig = figure;
    contourf(s1ab_r2bins,s1ab_zbins,mean_field_map_totalerr./mean_field_map,vc,'LineColor','none');
    xlabel('R^2 (cm^2)','fontsize', 18);
    ylabel('Z (uSec)','fontsize', 18);
    title(strcat(file_id_cp, '. Electric Field Fractional Error vs. R^2 v Z.'),'fontsize',16,'Interpreter','none');
    caxis([color_range_min color_range_max])
    g=colorbar;
    ylabel(g,'Electric Field Fractional Error','FontSize',16)
    set(g,'FontSize',18)
    hold on
    axis([0 625 0 330])
    myfigview(16);
    set(gca,'Ydir','reverse')


%% Next measure the ER band parameters

load('C:\Program Files\MATLAB\R2012a\bin\LUXCode\Scratch\RichardKnoche\Run04_CH3T\KrypCal2p18_IQs\MakePlots\CH3T_Sep2015_Kr2p18.mat');


s1_min=0;
s1_max=50;

%% First figure out where the outliers are in space

pulse_area_cut= (log10(s1_phe_both_xyz)+0.7*log10(s2_phe_both_xyz)<4.8) & (log10(s1_phe_both_xyz)+0.7*log10(s2_phe_both_xyz)>2.3) & ...
    (log10(s1_phe_both_xyz)-2*log10(s2_phe_both_xyz)<-5) & (log10(s1_phe_both_xyz)-2*log10(s2_phe_both_xyz)>-6.5);


max_r=25;

%% REDO THIS WITH SAME BINNING AS FIELD MAP

voxel_counter=1;

  ER_lower_power_fit_vox={};
  ER_upper_power_fit_vox={};
  ER_mean_power_fit_vox={};

 for r2_max=s1ab_r2bin_min+s1ab_r2_binsize:s1ab_r2_binsize:s1ab_r2bin_max;
        for z_max=s1ab_zbin_min+s1ab_z_binsize:s1ab_z_binsize:s1ab_zbin_max;   
            %note: voxels index in z, then y, then x
            
fiducial_cut=inrange(drift_time,[40 300]) & inrange(s2radius_del.^2,[0,r2_max]) & inrange(drift_time,[z_max-s1ab_z_binsize z_max]);  
tritium_cut= clean_cut & pulse_area_cut & fiducial_cut;

%% Make Band Plot
   band_sigma=1.28; %90% confidence bounds
    
             
            s2_s1_bin= 1:.05:5;
           ER_bin_size=2; %3 works
           index_i=1;
           ER_bins=ER_bin_size/2+s1_min:ER_bin_size:s1_max;
               lower_bound=ones(size(ER_bins));
               mean_ER_band=ones(size(ER_bins));
               upper_bound=ones(size(ER_bins));
                counts_ER_hist=ones(size(ER_bins));
               ER_sigma=ones(size(ER_bins));
               ER_sigma_sigma=ones(size(ER_bins));
        
%             temporary, for calculating correct E-lifetime   
        for ER_bin = (ER_bin_size/2)+s1_min:ER_bin_size:(s1_max); %HIST uses bin centers
   
        
            ER_fit_cut=tritium_cut & inrange(s1_phe_both_xyz,[ER_bin-ER_bin_size/2,ER_bin+ER_bin_size/2]);
             
            ER_band_hist=hist(log10(s2_phe_both_xyz(ER_fit_cut)./s1_phe_both_xyz(ER_fit_cut)),s2_s1_bin);
            counts_ER_hist(index_i)=sum(ER_band_hist);
           
            Fit_ER=fit(s2_s1_bin',ER_band_hist','gauss1');
             
            rms_fit=fit_ML_normal(log10(s2_phe_both_xyz(ER_fit_cut)./s1_phe_both_xyz(ER_fit_cut)),s2_s1_bin);
            
            %band=confint(Fit_ER,.1);
            %ER_sigma(index_i)= Fit_ER.c1/sqrt(2);
            ER_sigma(index_i)= sqrt(rms_fit.sig2);
            ER_sigma_sigma(index_i)=ER_sigma(index_i)/sqrt(2*counts_ER_hist(index_i));
            
            lower_bound(index_i)=Fit_ER.b1-ER_sigma(index_i)*band_sigma;
            mean_ER_band(index_i)=Fit_ER.b1;
            upper_bound(index_i)=Fit_ER.b1+ER_sigma(index_i)*band_sigma;
            
            %mean_ER_band_sim(index_i)=Fit_ER_sim.b1;
            
            index_i=index_i+1;
         end
        
         
         % Fit power law to the mean
         g=fittype('a*x^b');
         ER_mean_power_fit=fit(ER_bins(ER_bins>2)',mean_ER_band(ER_bins>2)',g,'startpoint',[ 2.5 -.1]);
         ER_mean_fit=ER_mean_power_fit.a.*ER_bins.^(ER_mean_power_fit.b);
            b=confint(ER_mean_power_fit,.68);
            sigma_a_ER_mean=(b(2,1)-b(1,1))/2;
            sigma_b_ER_mean=(b(2,2)-b(1,2))/2;
            
         ER_lower_power_fit=fit(ER_bins(ER_bins>2)',lower_bound(ER_bins>2)',g,'startpoint',[ 2.5 -.1]);
         ER_lower_fit=ER_lower_power_fit.a.*ER_bins.^(ER_lower_power_fit.b);
            b=confint(ER_lower_power_fit,.68);
            sigma_a_ER_lower=(b(2,1)-b(1,1))/2;
            sigma_b_ER_lower=(b(2,2)-b(1,2))/2;
            
         ER_upper_power_fit=fit(ER_bins(ER_bins>2)',upper_bound(ER_bins>2)',g,'startpoint',[ 2.5 -.1]);
         ER_upper_fit=ER_upper_power_fit.a.*ER_bins.^(ER_upper_power_fit.b);
            b=confint(ER_upper_power_fit,.68);
            sigma_a_ER_upper=(b(2,1)-b(1,1))/2;
            sigma_b_ER_upper=(b(2,2)-b(1,2))/2;

            
  ER_lower_power_fit_vox{voxel_counter}=ER_lower_power_fit;
  ER_upper_power_fit_vox{voxel_counter}=ER_upper_power_fit;
  ER_mean_power_fit_vox{voxel_counter}=ER_mean_power_fit;
  
  ER_lower_power_fit_vox_aerr{voxel_counter}=sigma_a_ER_lower;
  ER_upper_power_fit_vox_aerr{voxel_counter}=sigma_a_ER_upper;
  ER_mean_power_fit_vox_aerr{voxel_counter}=sigma_a_ER_mean;

  ER_lower_power_fit_vox_berr{voxel_counter}=sigma_b_ER_lower;
  ER_upper_power_fit_vox_berr{voxel_counter}=sigma_b_ER_upper;
  ER_mean_power_fit_vox_berr{voxel_counter}=sigma_b_ER_mean;
  
  r2_centers(voxel_counter)=r2_max-s1ab_r2_binsize/2;
  z_centers(voxel_counter)=z_max-s1ab_z_binsize/2;
  
  voxel_counter=voxel_counter+1;
        end
    end



%% Now relate the field to the ER band parameters

%Relating parameter a to field
temp_mean_field_map=mean_field_map(:);
temp_mean_field_map_totalerr=mean_field_map_totalerr(:);

for i=1:length(ER_mean_power_fit_vox)
temp_ER_mean_a(i)=ER_mean_power_fit_vox{i}.a;
temp_ER_mean_a_err(i)=ER_mean_power_fit_vox_aerr{i};

temp_ER_mean_b(i)=ER_mean_power_fit_vox{i}.b;
temp_ER_mean_b_err(i)=ER_mean_power_fit_vox_berr{i};

temp_ER_lower_a(i)=ER_lower_power_fit_vox{i}.a;
temp_ER_lower_a_err(i)=ER_lower_power_fit_vox_aerr{i};

temp_ER_lower_b(i)=ER_lower_power_fit_vox{i}.b;
temp_ER_lower_b_err(i)=ER_lower_power_fit_vox_berr{i};

temp_ER_upper_a(i)=ER_upper_power_fit_vox{i}.a;
temp_ER_upper_a_err(i)=ER_upper_power_fit_vox_aerr{i};

temp_ER_upper_b(i)=ER_upper_power_fit_vox{i}.b;
temp_ER_upper_b_err(i)=ER_upper_power_fit_vox_berr{i};
end


[field_to_mean_a_fit field_to_mean_a_fit_err]=polyfit(temp_mean_field_map,temp_ER_mean_a.',3);
[field_to_mean_b_fit field_to_mean_b_fit_err]=polyfit(temp_mean_field_map,temp_ER_mean_b.',3);
[field_to_lower_a_fit field_to_lower_a_fit_err]=polyfit(temp_mean_field_map,temp_ER_lower_a.',3);
[field_to_lower_b_fit field_to_lower_b_fit_err]=polyfit(temp_mean_field_map,temp_ER_lower_b.',3);
[field_to_upper_a_fit field_to_upper_a_fit_err]=polyfit(temp_mean_field_map,temp_ER_upper_a.',3);
[field_to_upper_b_fit field_to_upper_b_fit_err]=polyfit(temp_mean_field_map,temp_ER_upper_b.',3);

%% Plot coefficient relationship
figure
hold on;
xlabel('Electric Field (V/cm)'); ylabel('Power Law Coefficient');
myfigview(16);

[tempy tempyerr]=polyval(field_to_mean_a_fit,[0:1:600],field_to_mean_a_fit_err);
count=1;
for i=601:-1:1;
    tempfill(count)=tempy(i)-tempyerr(i);
    count=count+1;
end
fill([0:1:600,600:-1:0],[tempy+tempyerr,tempfill],[0.6 0.6 1],'EdgeColor',[0.6 0.6 1])
plot([0:1:600],tempy,'-b','LineWidth',2);
rkploterr(temp_mean_field_map,temp_ER_mean_a,temp_mean_field_map_totalerr,temp_ER_mean_a_err,[0 0 0],[],[],2)

%% Plot exponent relationship
figure
hold on;
xlabel('Electric Field (V/cm)'); ylabel('Power Law Exponent');
myfigview(16);

clear tempy tempyerr
[tempy tempyerr]=polyval(field_to_mean_b_fit,[0:1:600],field_to_mean_b_fit_err);
count=1;
for i=601:-1:1;
    tempfill(count)=tempy(i)-tempyerr(i);
    count=count+1;
end
fill([0:1:600,600:-1:0],[tempy+tempyerr,tempfill],[0.6 0.6 1],'EdgeColor',[0.6 0.6 1])
plot([0:1:600],tempy,'-b','LineWidth',2);
rkploterr(temp_mean_field_map,temp_ER_mean_b,temp_mean_field_map_totalerr,temp_ER_mean_b_err,[0 0 0],[],[],2)


save('Field_To_ERBand','field_to_mean_a_fit','field_to_mean_b_fit','field_to_lower_a_fit','field_to_lower_b_fit','field_to_upper_a_fit','field_to_upper_b_fit',...
    'field_to_mean_a_fit_err','field_to_mean_b_fit_err','field_to_lower_a_fit_err','field_to_lower_b_fit_err','field_to_upper_a_fit_err','field_to_upper_b_fit_err');