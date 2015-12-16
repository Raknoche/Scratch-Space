path='C:\Users\Richard\Desktop\LUX Analysis\lux10_20140903T1918_cp13848';
% path='C:\Users\Richard\Desktop\lux10_20140903T1918_cp13943';
   
rqs_to_load = {'pulse_area_phe','event_timestamp_samples'...
   ,'pulse_classification' ...
   ,'z_drift_samples' , 's1s2_pairing','golden'...
   ,'selected_s1_s2','top_bottom_ratio','x_cm','y_cm'...
   ,'full_evt_area_phe',...
   'event_number','chi2','prompt_fraction','aft_t1_samples','pulse_start_samples',...
   'pulse_end_samples','top_bottom_asymmetry','aft_t0_samples','aft_t2_samples',...
   'full_evt_area_phe','admin','Kr83fit_s1a_area_phe','Kr83fit_s1b_area_phe',...
   'hft_t50r_samples','hft_t50l_samples','Kr83fit_dt_samples'};

d = LUXLoadMultipleRQMs_framework(path,rqs_to_load);  
d = Apply_TempCorrections_2p15(d,'lux10_20140903T1918_cp13848');

%%


%Hardcoding inputs for ease of use during testing
user_name='RichardKnoche';
dp_version='2.15_test';
algorithm_name='test_KrypCal';
submit=0;
delay_min=0;
inpaint_on=1;

        file=char(d.files_loaded(1));
        file_id=file(1:19);
        cp_start=strfind(file,'cp')+2;%find location of 'cp' in string
        cp=file(cp_start:cp_start+4);%cp number has four digits
        
                
                gs=0;
   
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

 
    d.z_drift_samples(isnan(d.z_drift_samples)) = 0.0; % get rid of NaN        
    events=sort(size(d.pulse_area_phe)); %The first element is the number of pulses. sometimes it's 5, sometimes it's 10
          
    s1_area_cut= inrange(d.pulse_area_phe,[s1area_bound_min,s1area_bound_max]);
    s2_area_cut= inrange(d.pulse_area_phe,[s2area_bound_min,s2area_bound_max]);
    
    s1_class=(d.pulse_classification==1) & s1_area_cut ; %area cut for Kr events
    s2_class=(d.pulse_classification==2) & s2_area_cut ;   
    s4_class=(d.pulse_classification==4) ;
   

%     s1_single_cut =logical( (s1_class & d.golden & d.selected_s1_s2).*repmat(sum(s2_class & d.golden & d.selected_s1_s2)==1,events(1),1) ); % was s1_before_s2_cut before using golden
%     s2_single_cut =logical( (s2_class & d.golden & d.selected_s1_s2).*repmat(sum(s1_class & d.golden & d.selected_s1_s2)==1,events(1),1) );


 
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

    drift_time = d.z_drift_samples(s2_single_cut)/100;  % us
    s2x=d.x_cm(s2_single_cut);
    s2y=d.y_cm(s1_single_cut);
    s2radius=sqrt(s2x.^2+s2y.^2);
    
    s1_phe_both_xyz=d.xyz_corrected_pulse_area_all_phe(s1_single_cut);
    s2_phe_both_xyz=d.xyz_corrected_pulse_area_all_phe(s2_single_cut);
    s1_phe_both=d.pulse_area_phe(s1_single_cut);
    s2_phe_both=d.pulse_area_phe(s2_single_cut);
    
    max_r=18;
    liquid_cut= inrange(drift_time,[30 300]) & inrange(s2radius,[0,max_r]) & inrange(s1_phe_both_xyz,[100 1000]) ...
    & (log10(s1_phe_both_xyz)+0.5*log10(s2_phe_both_xyz)>4.2) & s2_phe_both>100 ; %used to be drift time 100 to 250

    
%%
g1=0.092; g2=18.4;
energy_c=(1/73).*(s2_phe_both_xyz(liquid_cut)./(g2)+s1_phe_both_xyz(liquid_cut)./(g1));
energy_uc=(1/73).*(s2_phe_both(liquid_cut)./(g2)+s1_phe_both(liquid_cut)./(g1));


% energy_c=(1/73).*(d.xyz_corrected_pulse_area_all_phe(s2_single_cut)./(g2)+d.xyz_corrected_pulse_area_all_phe(s1_single_cut)./(g1));
% energy_uc=(1/73).*(d.pulse_area_phe(s2_single_cut)./(g2)+d.pulse_area_phe(s1_single_cut)./(g1));

figure
step(energy_c,[0:0.1:60],'-b')
hold on;
step(energy_uc,[0:0.1:60],'-r')
line([41.55 41.55],[0 30000],'Color','k','LineWidth',2)
xlabel('Pulse Area'); ylabel('Count');myfigview(16);
legend('Corrected','Uncorrected');

figure
step(s2_phe_both(liquid_cut),[0:200:30000],'-r')
hold on;
step(s2_phe_both_xyz(liquid_cut),[0:200:30000],'-b')
xlabel('S2');ylabel('Count');myfigview(16);
legend('DP, Uncorrected','DP, Corrected');

figure
step(s1_phe_both(liquid_cut),[0:2:600],'-r')
hold on;
step(s1_phe_both_xyz(liquid_cut),[0:2:600],'-b')
xlabel('S1');ylabel('Count');myfigview(16);
legend('DP, Uncorrected','DP, Corrected');