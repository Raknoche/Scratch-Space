%First measure the field in corrected space

path='C:\Users\Richard\Desktop\LUX Analysis\lux10_20141104T1748_cp??\';

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

file_id='lux10_20141104T1748_cp??'; file_id_cp='lux10_20141104T1748_cp??'; delay_min=0;

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
       

%% Find average S1a/S1b 

%% Convert S1a/S1b to field measurement - UPDATE THIS FOR AVERAGE

load('C:\Program Files\MATLAB\R2012a\bin\LUXCode\Scratch\RichardKnoche\KrypCal\KrypCal_2p19_Code\2015_Field_to_S1aS1b.mat');
[field_map_2015, field_map_2015_staterr]=polyval(S1aS1b_fit_2015,s1ab_r2z_mean,S1aS1b_fit_2015_error);

load('C:\Program Files\MATLAB\R2012a\bin\LUXCode\Scratch\RichardKnoche\KrypCal\KrypCal_2p19_Code\2014_Field_to_S1aS1b.mat');
[field_map_2014, field_map_2014_staterr]=polyval(S1aS1b_fit_2014,s1ab_r2z_mean,S1aS1b_fit_2014_error);

field_map_syserr=((field_map_2015-(field_map_2014+field_map_2015)/2).^2 + (field_map_2014-(field_map_2014+field_map_2015)/2).^2).^(1/2);
mean_field_map=(field_map_2014+field_map_2015)/2;
mean_field_map_totalerr=sqrt(field_map_syserr.^2+((field_map_2014_staterr+field_map_2015_staterr)/2).^2);

%% Make the ER band estimate from the Kr data set measurements

load('C:\Program Files\MATLAB\R2012a\bin\LUXCode\Scratch\RichardKnoche\Run04_CH3T\KrypCal2p18_IQs\MakePlots\ERBand\Field_To_ERBand.mat');

[powerlaw_a powerlaw_a_err]=polyval(field_to_mean_a_fit,average_field,field_to_mean_a_fit_err);
[powerlaw_b powerlaw_b_err]=polyval(field_to_mean_b_fit,average_field,field_to_mean_b_fit_err);

[powerlaw_a_upper powerlaw_a_upper_err]=polyval(field_to_upper_a_fit,average_field,field_to_upper_a_fit_err);
[powerlaw_b_upper powerlaw_b_upper_err]=polyval(field_to_upper_b_fit,average_field,field_to_upper_b_fit_err);

[powerlaw_a_lower powerlaw_a_lower_err]=polyval(field_to_lower_a_fit,average_field,field_to_lower_a_fit_err);
[powerlaw_b_lower powerlaw_b_lower_err]=polyval(field_to_lower_b_fit,average_field,field_to_lower_b_fit_err);

%% Make ER band from CH3T data and compare
load('C:\Program Files\MATLAB\R2012a\bin\LUXCode\Scratch\RichardKnoche\Run04_CH3T\KrypCal2p18_IQs\MakePlots\CH3T_Nov2014_Kr2p18.mat');

s1_min=0;
s1_max=50;

%% First figure out where the outliers are in space

pulse_area_cut= (log10(s1_phe_both_xyz)+0.7*log10(s2_phe_both_xyz)<4.8) & (log10(s1_phe_both_xyz)+0.7*log10(s2_phe_both_xyz)>2.3) & ...
    (log10(s1_phe_both_xyz)-2*log10(s2_phe_both_xyz)<-5) & (log10(s1_phe_both_xyz)-2*log10(s2_phe_both_xyz)>-6.5);


max_r=21;
fiducial_cut=inrange(drift_time,[40 300]) & inrange(s2radius_del,[0,21]);
   
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

            
  %Make Plot    

 no_fid_new=figure;
 hold on
 %plot(s1_phe_both_xyz(tritium_cut),log10(s2_phe_both_xyz(tritium_cut)./s1_phe_both_xyz(tritium_cut)),'ok','markersize',2)
 density(s1_phe_both_xyz(tritium_cut),log10(s2_phe_both_xyz(tritium_cut)./s1_phe_both_xyz(tritium_cut)),[0:.4:100],[0:.04:4]);
%  densityplot(s1_phe_both_xyz(tritium_cut),log10(s2_phe_both_xyz(tritium_cut)./s1_phe_both_xyz(tritium_cut)),'rectangle',[.5 .05],'lin');
 cc=colorbar;
    ylabel(cc,'Count','FontSize',16);
 colormap;
 ylabel('Corrected log10(S2/S1)'); xlabel('S1 corrected (Pulse Area Phe)'); 
 title(sprintf('Run 04 CH3T, R<%d cm',max_r),'fontsize',16,'Interpreter','none');
 myfigview(16);
 ylim([0 3.5]);xlim([0 s1_max]);
%  plot(ER_bins,lower_bound,'ok','markersize',6,'linewidth',2)
 plot(ER_bins,ER_lower_fit,'--k','markersize',14,'linewidth',2)
%  plot(ER_bins,mean_ER_band,'ok','markersize',6,'linewidth',2)%plot mean of ER band
 plot(ER_bins,ER_mean_fit,'-k','markersize',20,'linewidth',2)
%  plot(ER_bins,upper_bound,'ok','markersize',6,'linewidth',2)
 plot(ER_bins,ER_upper_fit,'--k','markersize',14,'linewidth',2)

 
     bin_step=2;

            total_count=length(log10(s2_phe_both_xyz(tritium_cut)));
 
    
 box;
  text(10,2.8,strcat('Total Count= ',num2str(total_count)),'fontsize',16);
  %text(5,.7,strcat('NEST\_SE_{bot} Fit=',num2str(ER_mean_NEST_fit.SE_ER_fit*73/64*27/max_r,3),'\pm ', num2str(sigma_SE_ER_fit*73/64*27/max_r,2),'[Phe]'),'fontsize',16,'color','m','interpreter','tex');
  %text(5,.8,'NEST Platinum V4-b', 'fontsize',16,'color','m','interpreter','tex' );
  text(5,.55,strcat('+',num2str(band_sigma),'\sigma = (',num2str(ER_upper_power_fit.a,3),'\pm', num2str(sigma_a_ER_upper,1) ,')\times S1_c ^{', num2str(ER_upper_power_fit.b,3),'\pm',num2str(sigma_b_ER_upper,1),'}'),'fontsize',16,'interpreter','tex');
  text(5,.3,strcat('Mean=(',num2str(ER_mean_power_fit.a,3),'\pm', num2str(sigma_a_ER_mean,1),')\times S1_c ^{', num2str(ER_mean_power_fit.b,3),'\pm', num2str(sigma_b_ER_mean,1),'}'),'fontsize',16,'interpreter','tex');
  text(5,.05,strcat('-',num2str(band_sigma),'\sigma = (',num2str(ER_lower_power_fit.a,3),'\pm', num2str(sigma_a_ER_lower,1),')\times S1_c ^{', num2str(ER_lower_power_fit.b,3),'\pm', num2str(sigma_b_ER_lower,1),'}'),'fontsize',16,'interpreter','tex');
  title('Nov2014 CH3T')
  
  s1_phe_both_xyz_all=s1_phe_both_xyz(tritium_cut & inrange(drift_time,[0,110]));
  log10_s2s1_xyz_all=log10(s2_phe_both_xyz(tritium_cut & inrange(drift_time,[0,110]))./s1_phe_both_xyz(tritium_cut & inrange(drift_time,[0,110])));
  ER_lower_fit_all=ER_lower_fit;
  ER_upper_fit_all=ER_upper_fit;
  ER_mean_fit_all=ER_mean_fit;
    
  
 
%% Overlay the estimate ER band
         ER_lower_fit=ER_lower_power_fit.a.*ER_bins.^(ER_lower_power_fit.b);

[powerlaw_a powerlaw_a_err]=polyval(field_to_mean_a_fit,average_field,field_to_mean_a_fit_err);
[powerlaw_b powerlaw_b_err]=polyval(field_to_mean_b_fit,average_field,field_to_mean_b_fit_err);

[powerlaw_a_upper powerlaw_a_upper_err]=polyval(field_to_upper_a_fit,average_field,field_to_upper_a_fit_err);
[powerlaw_b_upper powerlaw_b_upper_err]=polyval(field_to_upper_b_fit,average_field,field_to_upper_b_fit_err);

[powerlaw_a_lower powerlaw_a_lower_err]=polyval(field_to_lower_a_fit,average_field,field_to_lower_a_fit_err);
[powerlaw_b_lower powerlaw_b_lower_err]=polyval(field_to_lower_b_fit,average_field,field_to_lower_b_fit_err);

 plot(ER_bins,powerlaw_a.*ER_bins.^powerlaw_b,'-r','markersize',20,'linewidth',2)
 plot(ER_bins,powerlaw_a_upper.*ER_bins.^powerlaw_b_upper,'--r','markersize',14,'linewidth',2)
 plot(ER_bins,powerlaw_a_lower.*ER_bins.^powerlaw_b_lower,'--r','markersize',14,'linewidth',2)