path='C:\Users\Richard\Desktop\TestData\CalibrationData\lux10_20141027T2102_cp14178'
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

inpaint_on=0;
delay_min=0;
file_id_cp='test';
file_id='test';

%%
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
        
    d.phe_bottom=d.pulse_area_phe./(1+d.top_bottom_ratio); %bottom PMT pulse area

    s1_phe_both = d.pulse_area_phe(s1_single_cut);
    s1_phe_bottom = d.phe_bottom(s1_single_cut);
    
    s2_phe_both = d.pulse_area_phe(s2_single_cut);
    s2_phe_bottom = d.phe_bottom(s2_single_cut);
    
    s2x = d.x_cm(s2_single_cut);
    s2y = d.y_cm(s2_single_cut);   

    d.livetime_sec=sum(d.livetime_end_samples-d.livetime_latch_samples)/1e8;
    evt_cut=logical(sum(s2_single_cut));%Cut for all the events passing the single S1 & S2 cut
    event_timestamp_samples=d.event_timestamp_samples(evt_cut);
    
    time_wait_cut=event_timestamp_samples/1e8/60 > delay_min; %allow x min for Kr mixing
    
    s2_width=(d.aft_t2_samples(s2_single_cut)-d.aft_t0_samples(s2_single_cut)); %cut above 800 samples
    s1_width=d.pulse_end_samples(s1_single_cut)-d.pulse_start_samples(s1_single_cut);
    
%     
%     clean_cut=-(s1_phe_both+s2_phe_both)' + d.full_evt_area_phe(logical(squeeze(sum(s2_single_cut,1)))) < 500 ...
%         & drift_time' > 0;%The WIMP clean cut at <100 Phe would remove 80% of the Kr83 events. 500 only removes 0.3%
%     clean_cut=clean_cut' & time_wait_cut';
%    
%          s1_phe_bottom=s1_phe_bottom(clean_cut);
%          s1_phe_both=s1_phe_both(clean_cut);
%          s2_phe_bottom=s2_phe_bottom(clean_cut);
%          s2_phe_both=s2_phe_both(clean_cut);
%          drift_time=drift_time(clean_cut);
%          s2x=s2x(clean_cut);
%          s2y=s2y(clean_cut);
         s2radius = (s2x.^2+s2y.^2).^(0.5);
%          event_timestamp_samples=event_timestamp_samples(clean_cut);
%          s2_width=s2_width(clean_cut);
%          s1_width=s1_width(clean_cut);


kr_energy_cut = inrange(s1_phe_both,[s1area_bound_min s1area_bound_max]) & inrange(s2_phe_both,[s2area_bound_min s2area_bound_max]);%for counting Kr events
Kr_events=length(drift_time(kr_energy_cut));% Are there enough events for the XY or XYZ map?

%%%% Detector Center Detection %%%%%%%
x_center=mean(s2x(drift_time>4)); %exclude extraction field region at uSec<4
y_center=mean(s2y(drift_time>4));
z_center=mean(drift_time(drift_time>4));
det_edge=330;
%%
    zcut_min = 10; %field changes at 4 uSec
    zcut_max = 0.95*det_edge;


    if Kr_events < 30000 %then create a 10x10 grid. Otherwise 25x25. Want about 30 evts per bin
           
            S1xybinsize = sqrt((50*50)/(Kr_events/150));
            S1_xbin_min = -25;
            S1_xbin_max = 25;
            S1_ybin_min = -25;
            S1_ybin_max = 25;
            SE_ybin_min = -25;
            SE_ybin_max = 25;

            S2xybinsize = sqrt((50*50)/(Kr_events/150));%cm
%             SExybinsize = 5; DEFINED LATER, based on number of SE events
            S2_xbin_min = -25;
            S2_xbin_max = 25;
            S2_ybin_min = -25;
            S2_ybin_max = 25;
            SE_xbin_min = -25;
            SE_xbin_max = 25;
            
                    
    end
    
    s1_xbins = (S1_xbin_max-S1_xbin_min)./S1xybinsize;
    s1_ybins = (S1_ybin_max-S1_ybin_min)./S1xybinsize;

    s2_xbins = (S2_xbin_max-S2_xbin_min)./S2xybinsize;
    s2_ybins = (S2_ybin_max-S2_ybin_min)./S2xybinsize;
        

 if Kr_events < 4000; % skip krpCal if less than 4000 events
     calibrations=zeros(6,1);
 end
 
  %% Kr S1a/S2b map

%Set up cuts and variables
    s1ab_cut=logical(sum(s1_single_cut(:,:),1));
    s1a_phe_both=d.Kr83fit_s1a_area_phe(s1ab_cut);
    s1b_phe_both=d.Kr83fit_s1b_area_phe(s1ab_cut);
    s1ab_timing=d.Kr83fit_dt_samples(s1ab_cut);
    s1ab_x=d.x_cm(s2_single_cut);
    s1ab_y=d.y_cm(s2_single_cut);
    s1ab_z=d.z_drift_samples(s2_single_cut)/100;
   
    s1ab_timing_cut=s1ab_timing>13 & s1b_phe_both>30;
       %plotting x=s1a v y=s1a/s1b cuts: y<0.035x-1, y>-0.035.*x+5; y>0.01.*x; %y<4.2; x<255
%     s1ab_areacut= (s1a_phe_both./s1b_phe_both-0.035.*s1a_phe_both<-1) & (s1a_phe_both./s1b_phe_both+0.035.*s1a_phe_both>5)...
%         & (s1a_phe_both./s1b_phe_both > 0.01.* s1a_phe_both) & (s1a_phe_both./s1b_phe_both<4.2) & (s1a_phe_both<255);
    
    s1ab_x=s1ab_x(s1ab_timing_cut);
    s1ab_y=s1ab_y(s1ab_timing_cut);
    s1ab_z=s1ab_z(s1ab_timing_cut);
    s1a_phe_both=s1a_phe_both(s1ab_timing_cut);
    s1b_phe_both=s1b_phe_both(s1ab_timing_cut);
    
    dT_step=det_edge/(length(s1ab_z)/300); %2000 events per z bin
    
    
    clear s1a_z_means s1a_z_means_err s1b_z_means s1b_z_means_err s1ab_bincenters s1ab_z_means s1ab_z_means_err
%Z Dependence of S1a/S1b
i=1;
    for z_max=10+dT_step:dT_step:det_edge;
        s1a_fit=fit([0:2:400].',hist(s1a_phe_both(s1a_phe_both>0 & s1a_phe_both<400 & inrange(s1ab_z,[z_max-dT_step,z_max]).'),[0:2:400]).','gauss1');
        s1b_fit=fit([0:2:300].',hist(s1b_phe_both(s1b_phe_both>0 & s1b_phe_both<300 & inrange(s1ab_z,[z_max-dT_step,z_max]).'),[0:2:300]).','gauss1');
        s1a_z_means(i)=s1a_fit.b1;
        s1a_z_means_err(i)=s1a_fit.c1/sqrt(2)/sqrt(length(s1a_phe_both(s1a_phe_both>0 & s1a_phe_both<400 & inrange(s1ab_z,[z_max-dT_step,z_max]).')));
        s1b_z_means(i)=s1b_fit.b1;
        s1b_z_means_err(i)=s1b_fit.c1/sqrt(2)/sqrt(length(s1b_phe_both(s1b_phe_both>0 & s1b_phe_both<300 & inrange(s1ab_z,[z_max-dT_step,z_max]).')));
        s1ab_bincenters(i)=z_max-dT_step/2;
        i=i+1;
    end

s1ab_z_means=s1a_z_means./s1b_z_means;
s1ab_z_means_err=sqrt( (s1b_z_means_err.*s1a_z_means./(s1b_z_means.^2)).^2 + (s1a_z_means_err./s1b_z_means).^2);
    
%fit polynomial to s1az and s1bz, to remove z dependence later in xy fit
[s1a_P, s1a_S]=polyfit(s1ab_bincenters,s1a_z_means,3);
[s1b_P, s1b_S]=polyfit(s1ab_bincenters,s1b_z_means,3);
[s1ab_P, s1ab_S]=polyfit(s1ab_bincenters,s1ab_z_means,3);

%S1a z dependence plot    
s1az_plot=figure;
errorbar(s1ab_bincenters,s1a_z_means,s1a_z_means_err,'.k')   
xlabel('Drift Time (uSec)');ylabel('S1a Mean (phe)'); title('S1a Z Dependence'); myfigview(16);
hold on;
plot([10:1:320],polyval(s1a_P,[10:1:320]),'-r','LineWidth',2)


%S1b z dependence plot    
s1bz_plot=figure;
errorbar(s1ab_bincenters,s1b_z_means,s1b_z_means_err,'.k')   
xlabel('Drift Time (uSec)');ylabel('S1b Mean (phe)'); title('S1b Z Dependence'); myfigview(16);
hold on;
plot([10:1:320],polyval(s1b_P,[10:1:320]),'-r','LineWidth',2)
 
%S1a/b z dependence plot
s1a_over_bz_plot=figure;
errorbar(s1ab_bincenters,s1ab_z_means,s1ab_z_means_err,'.k')   
xlabel('Drift Time (uSec)');ylabel('S1a/b'); title('S1 a/b Z Dependence'); myfigview(16);
hold on;
plot([10:1:320],polyval(s1ab_P,[10:1:320]),'-r','LineWidth',2)

%Correcting z dependence to get XY dependence in same manner that Kr/CH3T
%map does
    s1a_phe_both_z=s1a_phe_both.'.*polyval(s1a_P,z_center)./polyval(s1a_P,s1ab_z);
    s1b_phe_both_z=s1b_phe_both.'.*polyval(s1b_P,z_center)./polyval(s1b_P,s1ab_z);

%% S1a/b XY map (after z correction)    
s1ab_xbin_min=-25;
s1ab_ybin_min=-25;
s1ab_xbin_max=25;
s1ab_ybin_max=25;
s1ab_xybinsize=sqrt((50*50)/(length(s1ab_z)/300));
clear s1a_xy_mean s1a_xy_mean_err s1b_xy_mean s1b_xy_mean_err

s1ab_xbins=s1ab_xbin_min+s1ab_xybinsize/2:s1ab_xybinsize:s1ab_xbin_max;
s1ab_ybins=s1ab_ybin_min+s1ab_xybinsize/2:s1ab_xybinsize:s1ab_ybin_max; 

s1a_xy_mean=zeros(length(s1ab_xbins),length(s1ab_ybins));
s1b_xy_mean=zeros(length(s1ab_xbins),length(s1ab_ybins));

 for x_bin=s1ab_xbin_min:s1ab_xybinsize:(s1ab_xbin_max-s1ab_xybinsize/2);
        for y_bin=s1ab_ybin_min:s1ab_xybinsize:(s1ab_ybin_max-s1ab_xybinsize/2);          
          x_min = x_bin; x_max = x_bin+s1ab_xybinsize;
          y_min = y_bin; y_max = y_bin+s1ab_xybinsize;      

            x_count=int32(1+x_bin/s1ab_xybinsize+s1ab_xbin_max/s1ab_xybinsize);
            y_count=int32(1+y_bin/s1ab_xybinsize+s1ab_ybin_max/s1ab_xybinsize); 

          bin_cut = s1a_phe_both>0 & s1b_phe_both>0 & inrange(s1ab_x,[x_min,x_max]).' & inrange(s1ab_y,[y_min,y_max]).'; 
            if length(s1a_phe_both_z(bin_cut))>100;
                s1a_z_fit=fit([0:2:400].',hist(s1a_phe_both_z(bin_cut.' & inrange(s1a_phe_both_z,[0 400])),[0:2:400]).','gauss1');
                s1a_xy_mean(y_count,x_count)=s1a_z_fit.b1;
                s1a_xy_mean_err(y_count,x_count)=s1a_z_fit.c1/sqrt(2)/length(s1a_phe_both_z(bin_cut.' & inrange(s1a_phe_both_z,[0 400])));

                
                s1b_z_fit=fit([0:2:300].',hist(s1b_phe_both_z(bin_cut.' & inrange(s1b_phe_both_z,[0 300])),[0:2:300]).','gauss1');
                s1b_xy_mean(y_count,x_count)=s1b_z_fit.b1;
                s1b_xy_mean_err(y_count,x_count)=s1b_z_fit.c1/sqrt(2)/length(s1b_phe_both_z(bin_cut.' & inrange(s1b_phe_both_z,[0 300])));
            else
                s1a_xy_mean(y_count,x_count)=0;
                s1a_xy_mean_err(y_count,x_count)=0;
                s1b_xy_mean(y_count,x_count)=0;
                s1b_xy_mean_err(y_count,x_count)=0;   
            end
                
        end
 end
     s1ab_xy_mean=s1a_xy_mean./s1b_xy_mean;
     s1ab_xy_mean_err=sqrt( (s1b_xy_mean_err.*s1a_xy_mean./(s1b_xy_mean.^2)).^2 + (s1a_xy_mean_err./s1b_xy_mean).^2);
   

    %Plot s1a XY means
    color_range_max=max(max(s1a_xy_mean));
    color_range_min=min(min(s1a_xy_mean(s1a_xy_mean>0)));
    vc_step=(color_range_max-color_range_min)/50;
    vc=color_range_min:vc_step:color_range_max;
    
    s1a_xy_fig = figure;
    contourf(s1ab_xbins,s1ab_ybins,s1a_xy_mean,vc,'LineColor','none');
    xlabel('x (cm)','fontsize', 18);
    ylabel('y (cm)','fontsize', 18);
    title(strcat(file_id_cp, '. S1a Mean vs. XY.'),'fontsize',16,'Interpreter','none');
    caxis([color_range_min color_range_max])
    g=colorbar;
    ylabel(g,'S1a Mean (phe)','FontSize',16)
    set(g,'FontSize',18)
    hold on
    LUXPlotTopPMTs
    axis([-25 25 -25 25])
    myfigview(16);

    %Plot s1a XY means
    color_range_max=max(max(s1b_xy_mean));
    color_range_min=min(min(s1b_xy_mean(s1b_xy_mean>0)));
    vc_step=(color_range_max-color_range_min)/50;
    vc=color_range_min:vc_step:color_range_max;
    
    s1b_xy_fig = figure;
    contourf(s1ab_xbins,s1ab_ybins,s1b_xy_mean,vc,'LineColor','none');
    xlabel('x (cm)','fontsize', 18);
    ylabel('y (cm)','fontsize', 18);
    title(strcat(file_id_cp, '. S1b Mean vs. XY.'),'fontsize',16,'Interpreter','none');
    caxis([color_range_min color_range_max])
    g=colorbar;
    ylabel(g,'S1b Mean (phe)','FontSize',16)
    set(g,'FontSize',18)
    hold on
    LUXPlotTopPMTs
    axis([-25 25 -25 25])
    myfigview(16);

    %s1a/b XY Ratio
    color_range_max=max(max(s1ab_xy_mean));
    color_range_min=min(min(s1ab_xy_mean(s1ab_xy_mean>0)));
    vc_step=(color_range_max-color_range_min)/50;
    vc=color_range_min:vc_step:color_range_max;
    
    s1ab_xy_fig = figure;
    contourf(s1ab_xbins,s1ab_ybins,s1ab_xy_mean,vc,'LineColor','none');
    xlabel('x (cm)','fontsize', 18);
    ylabel('y (cm)','fontsize', 18);
    title(strcat(file_id_cp, '. S1a/b Ratio vs. XY.'),'fontsize',16,'Interpreter','none');
    caxis([color_range_min color_range_max])
    g=colorbar;
    ylabel(g,'S1a/b Ratio','FontSize',16)
    set(g,'FontSize',18)
    hold on
    LUXPlotTopPMTs
    axis([-25 25 -25 25])
    myfigview(16);

    %% 3D Kr s1a/s1b map

if length(s1ab_z)>100000; %require 100,000 events to do 3D binning
s1ab_xyz_numbins=floor((length(s1ab_z)/200)^(1/3)); %number of bins in one direction
s1ab_xyz_zstep=det_edge/s1ab_xyz_numbins;
s1ab_xyz_xstep=50/s1ab_xyz_numbins;
s1ab_xyz_ystep=50/s1ab_xyz_numbins;
s1ab_xyz_xmax=25;
r_max=25;

s1ab_xyz_zbins=10+s1ab_xyz_zstep/2:s1ab_xyz_zstep:10+s1ab_xyz_zstep/2+s1ab_xyz_zstep*s1ab_xyz_numbins; %20 us
s1ab_xyz_xbins=(-r_max+s1ab_xyz_xstep/2):s1ab_xyz_xstep:r_max;
s1ab_xyz_ybins=(-r_max+s1ab_xyz_ystep/2):s1ab_xyz_ystep:r_max;

s1a_xyz_mean=zeros(length(s1ab_xyz_xbins),length(s1ab_xyz_ybins),length(s1ab_xyz_zbins));
s1a_xyz_mean_err=zeros(length(s1ab_xyz_xbins),length(s1ab_xyz_ybins),length(s1ab_xyz_zbins));
s1b_xyz_mean=zeros(length(s1ab_xyz_xbins),length(s1ab_xyz_ybins),length(s1ab_xyz_zbins));
s1b_xyz_mean_err=zeros(length(s1ab_xyz_xbins),length(s1ab_xyz_ybins),length(s1ab_xyz_zbins));

for k = s1ab_xyz_zbins; %s1zbins are the center of the bin    
    for i = (-r_max+s1ab_xyz_xstep):s1ab_xyz_xstep:r_max %map to rows going down (y_cm) -- s1x bins are the max edge of the bin
        for j = (-r_max+s1ab_xyz_ystep):s1ab_xyz_ystep:r_max %map columns across (x_cm) -- s1y bins are the max edge of the bin
                       
            l=int8(i/s1ab_xyz_xstep+s1ab_xyz_xmax/s1ab_xyz_xstep);
            m=int8(j/s1ab_xyz_ystep+s1ab_xyz_xmax/s1ab_xyz_ystep);
            n=int8((k-(10+s1ab_xyz_zstep/2))/s1ab_xyz_zstep + 1);
            
           %sort, and make the cut. using the variable q
           q = s1ab_x<j & s1ab_x>(j-s1ab_xyz_xstep) & s1ab_y<i & s1ab_y>(i-s1ab_xyz_ystep) & inrange(s1ab_z,[(k-s1ab_xyz_zstep/2) , (k+s1ab_xyz_zstep/2)]); %no 0th element!
                    %swap x,y due to matrix x-row, y-column definition

            %Count the number of events per bin
            Count_S1ab_3D(l,m,n)=length(s1a_phe_both(q));
            
            if (Count_S1ab_3D(l,m,n) >= 100) % at least 100 counts before fitting. 
                s1a_z_fit=fit([0:2:400].',hist(s1a_phe_both(q.' & inrange(s1a_phe_both,[0 400])),[0:2:400]).','gauss1');
                s1a_xyz_mean(l,m,n)=s1a_z_fit.b1;
                s1a_xyz_mean_err(l,m,n)=s1a_z_fit.c1/sqrt(2)/length(s1a_phe_both(q.' & inrange(s1a_phe_both,[0 400])));

                
                s1b_z_fit=fit([0:2:300].',hist(s1b_phe_both(q.' & inrange(s1b_phe_both,[0 300]) ),[0:2:300]).','gauss1');
                s1b_xyz_mean(l,m,n)=s1b_z_fit.b1;
                s1b_xyz_mean_err(l,m,n)=s1b_z_fit.c1/sqrt(2)/length(s1b_phe_both(q.' & inrange(s1b_phe_both,[0 300])));
             else %not enough stats to do the fit
                s1a_xyz_mean(l,m,n)=0;
                s1a_xyz_mean_err(l,m,n)=0;
                s1b_xyz_mean(l,m,n)=0;
                s1b_xyz_mean_err(l,m,n)=0;
            end
                                      
        end
    end
k
end

     s1ab_xyz_mean=s1a_xyz_mean./s1b_xyz_mean;
     s1ab_xyz_mean_err=sqrt( (s1b_xyz_mean_err.*s1a_xyz_mean./(s1b_xyz_mean.^2)).^2 + (s1a_xyz_mean_err./s1b_xyz_mean).^2);

end
 

%% REMOVING FIELD


 %Hardcoded numbers from Energy Chi2 Fit
           
s1ab_s2_xyz_fit_map.p1=-0.8722;
s1ab_s2_xyz_fit_map.p2=3.4276;
c1=1-s1ab_s2_xyz_fit_map.p1.*2.7347.^2-2.7347.*s1ab_s2_xyz_fit_map.p2;
s1ab_s2_xyz_fit_map.p3=c1;

s1ab_s1_xyz_fit_map.p1=0.46;
b2=1-2.7347.*s1ab_s1_xyz_fit_map.p1;
s1ab_s1_xyz_fit_map.p2=b2;
 
 %If we can use 3D map, do so.  If there aren't enough events do 1D then 2D
 if length(s1ab_z)>100000; %require 100,000 events to do 3D binning
 
%Removing zeros from map with nearest neighbor method
s1ab_xyz_mean_old=s1ab_xyz_mean;
temp_map=s1ab_xyz_mean;
q= ~isnan(temp_map);
r=nearestpoint(find(~q),find(q));
s1ab_xyz_mean=temp_map; temp_mapQ=temp_map(q); s1ab_xyz_mean(~q)=temp_mapQ(r);s1ab_xyz_mean=reshape(s1ab_xyz_mean,size(s1ab_xyz_mean_old,1),size(s1ab_xyz_mean_old,2),size(s1ab_xyz_mean_old,3));

%interpolating with cubic so that we get NaN when extrapolating
 s1as1b_values=interp3(s1ab_xyz_xbins,s1ab_xyz_ybins,s1ab_xyz_zbins,s1ab_xyz_mean,s2x,s2y,drift_time,'cubic');
 s1as1b_at_center=interp3(s1ab_xyz_xbins,s1ab_xyz_ybins,s1ab_xyz_zbins,s1ab_xyz_mean,x_center,y_center,z_center,'cubic');
 s2_3Dfieldremoval= (s1ab_s2_xyz_fit_map.p1.*s1as1b_values.^2+s1ab_s2_xyz_fit_map.p2.*s1as1b_values+s1ab_s2_xyz_fit_map.p3)./(s1ab_s2_xyz_fit_map.p1.*s1as1b_at_center.^2+s1ab_s2_xyz_fit_map.p2.*s1as1b_at_center+s1ab_s2_xyz_fit_map.p3); 
 s1_3Dfieldremoval= (s1ab_s1_xyz_fit_map.p1.*s1as1b_values+s1ab_s1_xyz_fit_map.p2)./(s1ab_s1_xyz_fit_map.p1.*s1as1b_at_center+s1ab_s1_xyz_fit_map.p2); 
 
%Filling in NaN values with the s1a/s1b Z map polynomial 
  s1as1b_zvalues= polyval(s1ab_P,drift_time);
  s1as1b_z_at_center=polyval(s1ab_P,z_center);   
  s2_zfieldremoval= (s1ab_s2_xyz_fit_map.p1.*s1as1b_zvalues.^2+s1ab_s2_xyz_fit_map.p2.*s1as1b_zvalues+s1ab_s2_xyz_fit_map.p3)./(s1ab_s2_xyz_fit_map.p1.*s1as1b_z_at_center.^2+s1ab_s2_xyz_fit_map.p2.*s1as1b_z_at_center+s1ab_s2_xyz_fit_map.p3); 
  s1_zfieldremoval= (s1ab_s1_xyz_fit_map.p1.*s1as1b_zvalues+s1ab_s1_xyz_fit_map.p2)./(s1ab_s1_xyz_fit_map.p1.*s1as1b_z_at_center+s1ab_s1_xyz_fit_map.p2); 
  s2_3Dfieldremoval(isnan(s2_3Dfieldremoval))=s2_zfieldremoval(isnan(s2_3Dfieldremoval));
  s1_3Dfieldremoval(isnan(s1_3Dfieldremoval))=s1_zfieldremoval(isnan(s1_3Dfieldremoval));
 
s2_phe_bottom=s2_phe_bottom./s2_3Dfieldremoval;
s2_phe_both=s2_phe_both./s2_3Dfieldremoval;
s1_phe_bottom=s1_phe_bottom./s1_3Dfieldremoval;
s1_phe_both=s1_phe_both./s1_3Dfieldremoval;

 else 
%Removing zeros from map with nearest neighbor method
s1ab_xy_mean_old=s1ab_xy_mean;
temp_map=s1ab_xy_mean;
q= ~isnan(temp_map);
r=nearestpoint(find(~q),find(q));
s1ab_xy_mean=temp_map; temp_mapQ=temp_map(q); s1ab_xy_mean(~q)=temp_mapQ(r);s1ab_xy_mean=reshape(s1ab_xy_mean,size(s1ab_xy_mean_old,1),size(s1ab_xy_mean_old,2));

%interpolating with cubic so that we get NaN when extrapolating
 s1as1b_zvalues= polyval(s1ab_P,drift_time);
 s1as1b_zcorr_xyvalues=interp2(s1ab_xbins,s1ab_ybins,s1ab_xy_mean,s2x,s2y,'cubic');
 s1as1b_z_at_center=polyval(s1ab_P,z_center); 
 s1as1b_xy_at_center=interp2(s1ab_xbins,s1ab_ybins,s1ab_xy_mean,x_center,y_center,'cubic');
 
 s2_zfieldremoval= (s1ab_s2_xyz_fit_map.p1.*s1as1b_zvalues.^2+s1ab_s2_xyz_fit_map.p2.*s1as1b_zvalues+s1ab_s2_xyz_fit_map.p3)./(s1ab_s2_xyz_fit_map.p1.*s1as1b_z_at_center.^2+s1ab_s2_xyz_fit_map.p2.*s1as1b_z_at_center+s1ab_s2_xyz_fit_map.p3);
 s2_xyfieldremoval= (s1ab_s2_xyz_fit_map.p1.*s1as1b_zcorr_xyvalues.^2+s1ab_s2_xyz_fit_map.p2.*s1as1b_zcorr_xyvalues+s1ab_s2_xyz_fit_map.p3)./(s1ab_s2_xyz_fit_map.p1.*s1as1b_xy_at_center.^2+s1ab_s2_xyz_fit_map.p2.*s1as1b_xy_at_center+s1ab_s2_xyz_fit_map.p3); 
 s1_zfieldremoval= (s1ab_s1_xyz_fit_map.p1.*s1as1b_zvalues+s1ab_s1_xyz_fit_map.p2)./(s1ab_s1_xyz_fit_map.p1.*s1as1b_z_at_center+s1ab_s1_xyz_fit_map.p2); 
 s1_xyfieldremoval= (s1ab_s1_xyz_fit_map.p1.*s1as1b_zcorr_xyvalues+s1ab_s1_xyz_fit_map.p2)./(s1ab_s1_xyz_fit_map.p1.*s1as1b_xy_at_center+s1ab_s1_xyz_fit_map.p2); 
 
 %Filling in NaN values with no XY correction.  (Z dependence should never be NaN)
 s2_xyfieldremoval(isnan(s2_xyfieldremoval))=1;
 s2_zfieldremoval(isnan(s2_zfieldremoval))=1;
 s1_xyfieldremoval(isnan(s1_xyfieldremoval))=1;
 s1_zfieldremoval(isnan(s1_zfieldremoval))=1;

 
s2_phe_bottom=s2_phe_bottom./(s2_zfieldremoval.*s2_xyfieldremoval);
s2_phe_both=s2_phe_both./(s2_zfieldremoval.*s2_xyfieldremoval);
s1_phe_bottom=s1_phe_bottom./(s1_zfieldremoval.*s1_xyfieldremoval);
s1_phe_both=s1_phe_both./(s1_zfieldremoval.*s1_xyfieldremoval);

 end
 
%   figure
% hold on;
% plot(drift_time,s2_zfieldremoval.*s2_xyfieldremoval,'.r')
% plot(drift_time,s2_zfieldremoval,'.k')
% 
% figure
% plot(drift_time,s2_3Dfieldremoval,'.k')
% hold on;
% plot(drift_time,s2_zfieldremoval,'.g')


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Calculate Single Electron Size 
%   Only using events between good Kr S1 and S2
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %   Calculate Single Electron Size 
    s4_class=(d.pulse_classification==4) ; %Use or class=2 to fix 33 phe cut off... have to edit se_wrt_first_s1_func tho
    d = se_wrt_first_s1_func(d);
    cut_se_s1_z = d.t_btw_se_first_s1 > 100*10;

    SE_good_cut =logical( cut_se_s1_z & s4_class & ( cumsum(s1_single_cut)-cumsum(s2_single_cut) ) ); % only events between good S1 and S2s    

    SE_phe_both=d.pulse_area_phe(SE_good_cut);
    SE_phe_bot=d.phe_bottom(SE_good_cut);
    SE_phe_top=SE_phe_both-SE_phe_bot;
    
    SE_x=d.x_cm(SE_good_cut);
    SE_y=d.y_cm(SE_good_cut);
    SE_radius=(SE_x.^2+SE_y.^2).^(1/2);
    
    %Both PMT first fit
    SE_fit=fit([0:1:60]',hist(SE_phe_both(SE_radius<17 & inrange(SE_phe_both,[0 60])),[0:1:60])','gauss1');
    SE_Size=SE_fit.b1;
    SE_sig=SE_fit.c1/sqrt(2);
    SE_both_conf=confint(SE_fit,0.683);
    SE_sig_mu_both=abs(SE_Size-SE_both_conf(1,2)); % one sigma
    SE_sig_sig_both=abs(SE_sig-SE_both_conf(1,3)/sqrt(2)); % one sigma
    
    %Both PMT Second fit
    n_sig_fit=2;
    
    fit_start=SE_Size-n_sig_fit*SE_sig;
    fit_end=SE_Size+n_sig_fit*SE_sig;
    xfit=fit_start:0.1:fit_end;
    
    SE_cut=SE_phe_both(SE_radius<17);
    [xxo, yyo, xo, no] = step(SE_cut(inrange(SE_cut,[fit_start fit_end])),[0:0.2:60]);
    s = fitoptions('Method','NonlinearLeastSquares','Robust','On','Algorithm','Trust-Region','Lower',[-5,0,0,0],'Upper',[5,Inf,Inf,Inf],...
        'Startpoint',[1,max(no),10,5]);
    skewGauss1 = fittype('a1/c1*exp(-((x-b1)/c1)^2)*(1+erf(a*(x-b1)/c1))','options',s);

    [f_both,g_both] = fit(xo',no',skewGauss1);
    delta_both = f_both.a/sqrt(1+f_both.a^2);
    SE_mu_both = f_both.b1 + f_both.c1/sqrt(2)*delta_both*sqrt(2/pi);
    SE_sig_both = sqrt((f_both.c1/sqrt(2))^2 * (1 - 2 * (delta_both)^2/pi));
    skew_both = (4-pi)/2 * (delta_both * sqrt(2/pi))^3 / (1 - 2 * delta_both^2/pi)^(3/2);
    kurt_both = 2*(pi-3) * (delta_both * sqrt(2/pi))^3 / (1 - 2 * delta_both^2/pi)^2;
    ci_both = confint(f_both,0.683);
    y_fit_both = f_both(xfit);   
   
    both_SE_fig=figure;
    hold on;
    step(SE_cut,[0:0.2:60],'k');  
    h_fit=plot(xfit,y_fit_both,'-r','linewidth',2);
    box on;
    myfigview(18);
    xlabel('Pulse Area Phe (Phe) '); ylabel('Count/Phe ');
    legend([h_fit],strcat('\mu = ',num2str(SE_mu_both,'%2.2f'),'\newline \sigma= ',num2str(SE_sig_both,'%2.2f'),'\newline skew=', num2str(skew_both,'%2.2f') ), 'location','northeast');


   %% SE Both XY Map
    SExybinsize=sqrt((50*50)/(length(SE_phe_both)/100));
    num_SE_events=length(SE_phe_both);
   se_xbins=SE_xbin_min+SExybinsize/2:SExybinsize:SE_xbin_max;
   se_ybins=SE_ybin_min+SExybinsize/2:SExybinsize:SE_ybin_max;
   
   mean_SE_both_xy = zeros(length(se_xbins),length(se_ybins));
   mean_err_SE_both_xy = zeros(length(se_xbins),length(se_ybins));
   skew_SE_both_xy = zeros(length(se_xbins),length(se_ybins));
   sigma_SE_both_xy = zeros(length(se_xbins),length(se_ybins));
   kurt_SE_both_xy = zeros(length(se_xbins),length(se_ybins));
   sigma_err_SE_both_xy = zeros(length(se_xbins),length(se_ybins));
   
  for x_bin=SE_xbin_min:SExybinsize:(SE_xbin_max-SExybinsize);
    for y_bin=SE_ybin_min:SExybinsize:(SE_ybin_max-SExybinsize);          
      x_min = x_bin; x_max = x_bin+SExybinsize;
      y_min = y_bin; y_max = y_bin+SExybinsize;      
      
        x_count=int32(1+x_bin/SExybinsize+SE_xbin_max/SExybinsize);
        y_count=int32(1+y_bin/SExybinsize+SE_ybin_max/SExybinsize); 
      
      
      bin_cut = (SE_phe_both>0) & inrange(SE_x,[x_min,x_max]) & inrange(SE_y,[y_min,y_max]);            

    if length(SE_phe_both(bin_cut)) > 60; 
        %Both PMT first fit
        clear temp_SE
        temp_SE=SE_phe_both(bin_cut);
        SE_fit_xy=fit([0:1:60]',hist(temp_SE(inrange(temp_SE,[0 60])),[0:1:60])','gauss1');
        SE_Size_xy=SE_fit_xy.b1;
        SE_sig_xy=SE_fit_xy.c1/sqrt(2);
        SE_both_conf_xy=confint(SE_fit_xy,0.683);
        SE_sig_mu_both_xy=abs(SE_Size_xy-SE_both_conf_xy(1,2)); % one sigma
        SE_sig_sig_both_xy=abs(SE_sig_xy-SE_both_conf_xy(1,3)/sqrt(2)); % one sigma
    
        %Both PMT Second fit
        n_sig_fit=2;
    
        fit_start=SE_Size_xy-n_sig_fit*SE_sig_xy;
        fit_end=SE_Size_xy+n_sig_fit*SE_sig_xy;
        xfit=fit_start:0.1:fit_end;
    
        [xxo, yyo, xo, no] = step(temp_SE(inrange(temp_SE,[fit_start fit_end])),[0:0.2:60]);
        s = fitoptions('Method','NonlinearLeastSquares','Robust','On','Algorithm','Trust-Region','Lower',[-5,0,0,0],'Upper',[5,Inf,Inf,Inf],...
        'Startpoint',[1,max(no),10,5]);
        skewGauss1 = fittype('a1/c1*exp(-((x-b1)/c1)^2)*(1+erf(a*(x-b1)/c1))','options',s);

    [f_both_xy,g_both_xy] = fit(xo',no',skewGauss1);
    delta_both_xy = f_both_xy.a/sqrt(1+f_both_xy.a^2);
    SE_mu_both_xy = f_both_xy.b1 + f_both_xy.c1/sqrt(2)*delta_both_xy*sqrt(2/pi);
    SE_sig_both_xy = sqrt((f_both_xy.c1/sqrt(2))^2 * (1 - 2 * (delta_both_xy)^2/pi));
    skew_both_xy = (4-pi)/2 * (delta_both_xy * sqrt(2/pi))^3 / (1 - 2 * delta_both_xy^2/pi)^(3/2);
    kurt_both_xy = 2*(pi-3) * (delta_both_xy * sqrt(2/pi))^3 / (1 - 2 * delta_both_xy^2/pi)^2;
    ci_both_xy = confint(f_both_xy,0.683);
    y_fit_both_xy = f_both_xy(xfit);   
      
    mean_SE_both_xy(y_count,x_count)=SE_mu_both_xy;  
    skew_SE_both_xy(y_count,x_count)=skew_both_xy; 
    sigma_SE_both_xy(y_count,x_count)=SE_sig_both_xy;
    kurt_SE_both_xy(y_count,x_count)=kurt_both_xy;
    mean_err_SE_both_xy(y_count,x_count)=SE_sig_mu_both_xy;
    sigma_err_SE_both_xy(y_count,x_count)=SE_sig_sig_both_xy;    
    

      else
          
        mean_SE_both_xy(y_count,x_count)=0;  
        skew_SE_both_xy(y_count,x_count)=0; 
        sigma_SE_both_xy(y_count,x_count)=0;
        kurt_SE_both_xy(y_count,x_count)=0;
        mean_err_SE_both_xy(y_count,x_count)=0;
        sigma_err_SE_both_xy(y_count,x_count)=0;    
    
      end           
    end    
  end

%%Plot the mean SE_XY both %%%%%%%%%%%%%%%%
color_range_max=max(max(mean_SE_both_xy));
color_range_min=min(min(mean_SE_both_xy(mean_SE_both_xy>0)));
vc_step=(color_range_max-color_range_min)/50;
vc=color_range_min:vc_step:color_range_max;
 
SE_xy_both_fig = figure;
contourf(se_xbins,se_ybins,mean_SE_both_xy,vc,'LineColor','none');
xlabel('x (cm)','fontsize', 18);
ylabel('y (cm)','fontsize', 18);
title(strcat(file_id_cp, '. SE Mean (Both PMT) vs. XY.'),'fontsize',16,'Interpreter','none');
caxis([color_range_min color_range_max])
g=colorbar;
ylabel(g,'phe','FontSize',16)
set(g,'FontSize',18)
hold on
LUXPlotTopPMTs
axis([-25 25 -25 25])
myfigview(16);

%%Plot the 1 sigma of the mean

color_range_max=max(max(mean_err_SE_both_xy));
color_range_min=min(min(mean_err_SE_both_xy(mean_err_SE_both_xy>0)));
vc_step=(color_range_max-color_range_min)/50;
vc=color_range_min:vc_step:color_range_max;
 
sigma_SE_xy_both_fig = figure;
contourf(se_xbins,se_ybins,mean_err_SE_both_xy,vc,'LineColor','none');
xlabel('x (cm)','fontsize', 18);
ylabel('y (cm)','fontsize', 18);
title(strcat(file_id_cp, '. 1-Sigma SE Mean(Both PMT) vs. XY.'),'fontsize',16,'Interpreter','none');
caxis([color_range_min color_range_max])
g=colorbar;
ylabel(g,'phe','FontSize',16)
set(g,'FontSize',18)
hold on
LUXPlotTopPMTs
axis([-25 25 -25 25])
myfigview(16);


%Calculate corrections matrix and 1 sigma corrections matrix for SE_both
% if inpaint_on==1;  mean_SE_both_xy(mean_SE_both_xy==0)=nan; mean_SE_both_xy=inpaint_nans(mean_SE_both_xy,3); end;
SE_both_mu_center=interp2(se_xbins,se_ybins,mean_SE_both_xy,x_center,y_center,'spline');%Normalize to the center (x=y=0)
norm_SE_both_mu=SE_both_mu_center./mean_SE_both_xy; 
norm_SE_both_mu(isinf(norm_SE_both_mu))=1;%no infinity! Do not correct events outside r=25cm
norm_SE_both_mu(isnan(norm_SE_both_mu))=1;

%Calculate error on SE center
% if inpaint_on==1;  mean_err_SE_both_xy(mean_err_SE_both_xy==0)=nan; mean_err_SE_both_xy=inpaint_nans(mean_err_SE_both_xy,3); end;
error_SE_both_mu_center=interp2(se_xbins,se_ybins,mean_err_SE_both_xy,x_center,y_center,'spline');%Normalize to the center (x=y=0)
error_SE_both_mu_norm=sqrt((error_SE_both_mu_center./mean_SE_both_xy).^2+(mean_err_SE_both_xy.*SE_both_mu_center./mean_SE_both_xy.^2).^2);
error_SE_both_mu_norm(isinf(error_SE_both_mu_norm))=1;%no infinity. Set Sigma=1, which is 100% uncertainty
error_SE_both_mu_norm(isnan(error_SE_both_mu_norm))=1;%no nan. Set Sigma=1, which is 100% uncertainty

%Calculate kurt,sigma, and skew at center
% if inpaint_on==1;  kurt_SE_both_xy(kurt_SE_both_xy==0)=nan; kurt_SE_both_xy=inpaint_nans(kurt_SE_both_xy,3); end;
% if inpaint_on==1;  skew_SE_both_xy(skew_SE_both_xy==0)=nan; skew_SE_both_xy=inpaint_nans(skew_SE_both_xy,3); end;
% if inpaint_on==1;  sigma_SE_both_xy(sigma_SE_both_xy==0)=nan; sigma_SE_both_xy=inpaint_nans(sigma_SE_both_xy,3); end;
SE_both_kurt_center=interp2(se_xbins,se_ybins,kurt_SE_both_xy,x_center,y_center,'spline');%Normalize to the center (x=y=0)
SE_both_skew_center=interp2(se_xbins,se_ybins,skew_SE_both_xy,x_center,y_center,'spline');%Normalize to the center (x=y=0)
SE_both_sigma_center=interp2(se_xbins,se_ybins,sigma_SE_both_xy,x_center,y_center,'spline');%Normalize to the center (x=y=0)

    
    %%
    
    %Bottom SE Size
    SE_fit=fit([0:0.2:30]',hist(SE_phe_bot(SE_radius<17 & inrange(SE_phe_bot,[0 30])),[0:0.2:30])','gauss1');
    SE_Size_bot=SE_fit.b1;
    SE_sig_bot=SE_fit.c1/sqrt(2);
    SE_bot_conf=confint(SE_fit,0.683);
    SE_sig_mu_bot=abs(SE_Size_bot-SE_bot_conf(1,2)); % one sigma
    SE_sig_sig_bot=abs(SE_sig_bot-SE_bot_conf(1,3)/sqrt(2)); % one sigma
    
        %Bottom PMT Second fit
    n_sig_fit=2;
    
    fit_start=SE_Size_bot-n_sig_fit*SE_sig;
    fit_end=SE_Size_bot+n_sig_fit*SE_sig;
    xfit=fit_start:0.1:fit_end;
        
    SE_cut_bot=SE_phe_bot(SE_radius<17);

    [xxo, yyo, xo, no] = step(SE_cut_bot(inrange(SE_cut_bot,[fit_start fit_end])),[0:0.2:30]);
    s = fitoptions('Method','NonlinearLeastSquares','Robust','On','Algorithm','Trust-Region','Lower',[-5,0,0,0],'Upper',[5,Inf,Inf,Inf],...
        'Startpoint',[1,max(no),10,5]);
    skewGauss1 = fittype('a1/c1*exp(-((x-b1)/c1)^2)*(1+erf(a*(x-b1)/c1))','options',s);

    [f_bot,g_bot] = fit(xo',no',skewGauss1);
    delta_bot = f_bot.a/sqrt(1+f_bot.a^2);
    SE_mu_bot = f_bot.b1 + f_bot.c1/sqrt(2)*delta_bot*sqrt(2/pi);
    SE_sig_bot = sqrt((f_bot.c1/sqrt(2))^2 * (1 - 2 * (delta_bot)^2/pi));
    skew_bot = (4-pi)/2 * (delta_bot * sqrt(2/pi))^3 / (1 - 2 * delta_bot^2/pi)^(3/2);
    kurt_bot = 2*(pi-3) * (delta_bot * sqrt(2/pi))^3 / (1 - 2 * delta_bot^2/pi)^2;
    ci_bot = confint(f_bot);
    y_fit_bot = f_bot(xfit);  

    
    bottom_SE_fig=figure;
    hold on;
    step(SE_cut_bot,[0:0.2:30],'k');  
    h_fit=plot(xfit,y_fit_bot,'-r','linewidth',2);
    box on;
    myfigview(18);
    xlabel('Pulse Area Phe (Phe) '); ylabel('Count/Phe ');
    legend([h_fit],strcat('\mu = ',num2str(SE_mu_bot,'%2.2f'),'\newline \sigma= ',num2str(SE_sig_bot,'%2.2f'),'\newline skew=', num2str(skew_bot,'%2.2f') ), 'location','northeast'); 

    
    %% SE Bot XY Map

   mean_SE_bot_xy = zeros(length(se_xbins),length(se_ybins));
   mean_err_SE_bot_xy = zeros(length(se_xbins),length(se_ybins));
   skew_SE_bot_xy = zeros(length(se_xbins),length(se_ybins));
   sigma_SE_bot_xy = zeros(length(se_xbins),length(se_ybins));
   kurt_SE_bot_xy = zeros(length(se_xbins),length(se_ybins));
   sigma_err_SE_bot_xy = zeros(length(se_xbins),length(se_ybins));
   
  for x_bin=SE_xbin_min:SExybinsize:(SE_xbin_max-SExybinsize);
    for y_bin=SE_ybin_min:SExybinsize:(SE_ybin_max-SExybinsize);          
      x_min = x_bin; x_max = x_bin+SExybinsize;
      y_min = y_bin; y_max = y_bin+SExybinsize;      
      
        x_count=int32(1+x_bin/SExybinsize+SE_xbin_max/SExybinsize);
        y_count=int32(1+y_bin/SExybinsize+SE_ybin_max/SExybinsize); 
      
      
      bin_cut = (SE_phe_bot>0) & inrange(SE_x,[x_min,x_max]) & inrange(SE_y,[y_min,y_max]);            

    if length(SE_phe_bot(bin_cut)) > 60; 
        %Both PMT first fit
        clear temp_SE
        temp_SE=SE_phe_bot(bin_cut);
        SE_fit_xy=fit([0:1:60]',hist(temp_SE(inrange(temp_SE,[0 60])),[0:1:60])','gauss1');
        SE_Size_xy=SE_fit_xy.b1;
        SE_sig_xy=SE_fit_xy.c1/sqrt(2);
        SE_conf_xy=confint(SE_fit_xy,0.683);
        SE_sig_mu_xy=abs(SE_Size_xy-SE_conf_xy(1,2)); % one sigma
        SE_sig_sig_xy=abs(SE_sig_xy-SE_conf_xy(1,3)/sqrt(2)); % one sigma
    
        %Both PMT Second fit
        n_sig_fit=2;
    
        fit_start=SE_Size_xy-n_sig_fit*SE_sig_xy;
        fit_end=SE_Size_xy+n_sig_fit*SE_sig_xy;
        xfit=fit_start:0.1:fit_end;
    
        [xxo, yyo, xo, no] = step(temp_SE(inrange(temp_SE,[fit_start fit_end])),[0:0.2:60]);
        s = fitoptions('Method','NonlinearLeastSquares','Robust','On','Algorithm','Trust-Region','Lower',[-5,0,0,0],'Upper',[5,Inf,Inf,Inf],...
        'Startpoint',[1,max(no),10,5]);
        skewGauss1 = fittype('a1/c1*exp(-((x-b1)/c1)^2)*(1+erf(a*(x-b1)/c1))','options',s);

    [f_xy,g_xy] = fit(xo',no',skewGauss1);
    delta_xy = f_xy.a/sqrt(1+f_xy.a^2);
    SE_mu_xy = f_xy.b1 + f_xy.c1/sqrt(2)*delta_xy*sqrt(2/pi);
    SE_sig_xy = sqrt((f_xy.c1/sqrt(2))^2 * (1 - 2 * (delta_xy)^2/pi));
    skew_xy = (4-pi)/2 * (delta_xy * sqrt(2/pi))^3 / (1 - 2 * delta_xy^2/pi)^(3/2);
    kurt_xy = 2*(pi-3) * (delta_xy * sqrt(2/pi))^3 / (1 - 2 * delta_xy^2/pi)^2;
    ci_xy = confint(f_xy,0.683);
    y_fit_xy = f_xy(xfit);   
      
    mean_SE_bot_xy(y_count,x_count)=SE_mu_xy;  
    skew_SE_bot_xy(y_count,x_count)=skew_xy; 
    sigma_SE_bot_xy(y_count,x_count)=SE_sig_xy;
    kurt_SE_bot_xy(y_count,x_count)=kurt_xy;
    mean_err_SE_bot_xy(y_count,x_count)=SE_sig_mu_xy;
    sigma_err_SE_bot_xy(y_count,x_count)=SE_sig_sig_xy;    
    

      else
          
        mean_SE_bot_xy(y_count,x_count)=0;  
        skew_SE_bot_xy(y_count,x_count)=0; 
        sigma_SE_bot_xy(y_count,x_count)=0;
        kurt_SE_bot_xy(y_count,x_count)=0;
        mean_err_SE_bot_xy(y_count,x_count)=0;
        sigma_err_SE_bot_xy(y_count,x_count)=0;    
    
      end           
    end    
  end

%%Plot the mean SE_XY bot %%%%%%%%%%%%%%%%
color_range_max=max(max(mean_SE_bot_xy));
color_range_min=min(min(mean_SE_bot_xy(mean_SE_bot_xy>0)));
vc_step=(color_range_max-color_range_min)/50;
vc=color_range_min:vc_step:color_range_max;
 
SE_xy_bot_fig = figure;
contourf(se_xbins,se_ybins,mean_SE_bot_xy,vc,'LineColor','none');
xlabel('x (cm)','fontsize', 18);
ylabel('y (cm)','fontsize', 18);
title(strcat(file_id_cp, '. SE Mean (Bot PMT) vs. XY.'),'fontsize',16,'Interpreter','none');
caxis([color_range_min color_range_max])
g=colorbar;
ylabel(g,'phe','FontSize',16)
set(g,'FontSize',18)
hold on
LUXPlotTopPMTs
axis([-25 25 -25 25])
myfigview(16);


%%Plot the 1 sigma of the mean

color_range_max=max(max(mean_err_SE_bot_xy));
color_range_min=min(min(mean_err_SE_bot_xy(mean_err_SE_bot_xy>0)));
vc_step=(color_range_max-color_range_min)/50;
vc=color_range_min:vc_step:color_range_max;
 
sigma_SE_xy_bot_fig = figure;
contourf(se_xbins,se_ybins,mean_err_SE_bot_xy,vc,'LineColor','none');
xlabel('x (cm)','fontsize', 18);
ylabel('y (cm)','fontsize', 18);
title(strcat(file_id_cp, '. 1-Sigma SE Mean(Bot PMT) vs. XY.'),'fontsize',16,'Interpreter','none');
caxis([color_range_min color_range_max])
g=colorbar;
ylabel(g,'phe','FontSize',16)
set(g,'FontSize',18)
hold on
LUXPlotTopPMTs
axis([-25 25 -25 25])
myfigview(16);
%Calculate corrections matrix and 1 sigma corrections matrix for SE_both
% if inpaint_on==1;  mean_SE_bot_xy(mean_SE_bot_xy==0)=nan; mean_SE_bot_xy=inpaint_nans(mean_SE_bot_xy,3); end;
SE_bot_mu_center=interp2(se_xbins,se_ybins,mean_SE_bot_xy,x_center,y_center,'spline');%Normalize to the center (x=y=0)
norm_SE_bot_mu=SE_bot_mu_center./mean_SE_bot_xy; 
norm_SE_bot_mu(isinf(norm_SE_bot_mu))=1;%no infinity! Do not correct events outside r=25cm
norm_SE_bot_mu(isnan(norm_SE_bot_mu))=1;

%Calculate error on SE center
% if inpaint_on==1;  mean_err_SE_bot_xy(mean_err_SE_bot_xy==0)=nan; mean_err_SE_bot_xy=inpaint_nans(mean_err_SE_bot_xy,3); end;
error_SE_bot_mu_center=interp2(se_xbins,se_ybins,mean_err_SE_bot_xy,x_center,y_center,'spline');%Normalize to the center (x=y=0)
error_SE_bot_mu_norm=sqrt((error_SE_bot_mu_center./mean_SE_bot_xy).^2+(mean_err_SE_bot_xy.*SE_bot_mu_center./mean_SE_bot_xy.^2).^2);
error_SE_bot_mu_norm(isinf(error_SE_bot_mu_norm))=1;%no infinity. Set Sigma=1, which is 100% uncertainty
error_SE_bot_mu_norm(isnan(error_SE_bot_mu_norm))=1;%no nan. Set Sigma=1, which is 100% uncertainty

%Calculate kurt,sigma, and skew at center
% if inpaint_on==1;  kurt_SE_bot_xy(kurt_SE_bot_xy==0)=nan; kurt_SE_bot_xy=inpaint_nans(kurt_SE_bot_xy,3); end;
% if inpaint_on==1;  skew_SE_bot_xy(skew_SE_bot_xy==0)=nan; skew_SE_bot_xy=inpaint_nans(skew_SE_bot_xy,3); end;
% if inpaint_on==1;  sigma_SE_bot_xy(sigma_SE_bot_xy==0)=nan; sigma_SE_bot_xy=inpaint_nans(sigma_SE_bot_xy,3); end;
SE_bot_kurt_center=interp2(se_xbins,se_ybins,kurt_SE_bot_xy,x_center,y_center,'spline');%Normalize to the center (x=y=0)
SE_bot_skew_center=interp2(se_xbins,se_ybins,skew_SE_bot_xy,x_center,y_center,'spline');%Normalize to the center (x=y=0)
SE_bot_sigma_center=interp2(se_xbins,se_ybins,sigma_SE_bot_xy,x_center,y_center,'spline');%Normalize to the center (x=y=0)

   
    %%
    %Top SE Size
    SE_fit=fit([0:0.2:30]',hist(SE_phe_top(SE_radius<17 & inrange(SE_phe_top,[0 30])),[0:0.2:30])','gauss1');
    SE_Size_top=SE_fit.b1;
    SE_sig_top=SE_fit.c1/sqrt(2);
    SE_top_conf=confint(SE_fit,0.683);
    SE_sig_mu_top=abs(SE_Size_top-SE_top_conf(1,2)); % one sigma
    SE_sig_sig_top=abs(SE_sig_top-SE_top_conf(1,3)/sqrt(2)); % one sigma
    
        %Top PMT Second fit
    n_sig_fit=2;
    
    fit_start=SE_Size_top-n_sig_fit*SE_sig;
    fit_end=SE_Size_top+n_sig_fit*SE_sig;
    xfit=fit_start:0.1:fit_end;
        
    SE_cut_top=SE_phe_top(SE_radius<17);

    [xxo, yyo, xo, no] = step(SE_cut_top(inrange(SE_cut_top,[fit_start fit_end])),[0:0.2:30]);
    s = fitoptions('Method','NonlinearLeastSquares','Robust','On','Algorithm','Trust-Region','Lower',[-5,0,0,0],'Upper',[5,Inf,Inf,Inf],...
        'Startpoint',[1,max(no),10,5]);
    skewGauss1 = fittype('a1/c1*exp(-((x-b1)/c1)^2)*(1+erf(a*(x-b1)/c1))','options',s);

    [f_top,g_top] = fit(xo',no',skewGauss1);
    delta_top = f_top.a/sqrt(1+f_top.a^2);
    SE_mu_top = f_top.b1 + f_top.c1/sqrt(2)*delta_top*sqrt(2/pi);
    SE_sig_top = sqrt((f_top.c1/sqrt(2))^2 * (1 - 2 * (delta_top)^2/pi));
    skew_top = (4-pi)/2 * (delta_top * sqrt(2/pi))^3 / (1 - 2 * delta_top^2/pi)^(3/2);
    kurt_top = 2*(pi-3) * (delta_top * sqrt(2/pi))^3 / (1 - 2 * delta_top^2/pi)^2;
    ci_top = confint(f_top);
    y_fit_top = f_top(xfit);  

    
    
    top_SE_fig=figure;
    hold on;
    step(SE_cut_top,[0:0.2:30],'k');
    h_fit=plot(xfit,y_fit_top,'-r','linewidth',2);
    box on;
    myfigview(18);
    xlabel('Pulse Area Phe (Phe) '); ylabel('Count/Phe ');
    legend([h_fit],strcat('\mu = ',num2str(SE_mu_top,'%2.2f'),'\newline \sigma= ',num2str(SE_sig_top,'%2.2f'),'\newline skew=', num2str(skew_top,'%2.2f') ), 'location','northeast');

    %% SE Top XY Map

   mean_SE_top_xy = zeros(length(se_xbins),length(se_ybins));
   mean_err_SE_top_xy = zeros(length(se_xbins),length(se_ybins));
   skew_SE_top_xy = zeros(length(se_xbins),length(se_ybins));
   sigma_SE_top_xy = zeros(length(se_xbins),length(se_ybins));
   kurt_SE_top_xy = zeros(length(se_xbins),length(se_ybins));
   sigma_err_SE_top_xy = zeros(length(se_xbins),length(se_ybins));
   
  for x_bin=SE_xbin_min:SExybinsize:(SE_xbin_max-SExybinsize);
    for y_bin=SE_ybin_min:SExybinsize:(SE_ybin_max-SExybinsize);          
      x_min = x_bin; x_max = x_bin+SExybinsize;
      y_min = y_bin; y_max = y_bin+SExybinsize;      
      
        x_count=int32(1+x_bin/SExybinsize+SE_xbin_max/SExybinsize);
        y_count=int32(1+y_bin/SExybinsize+SE_ybin_max/SExybinsize); 
      
      
      bin_cut = (SE_phe_top>0) & inrange(SE_x,[x_min,x_max]) & inrange(SE_y,[y_min,y_max]);            

    if length(SE_phe_top(bin_cut)) > 60; 
        %top PMT first fit
        clear temp_SE
        temp_SE=SE_phe_top(bin_cut);
        SE_fit_xy=fit([0:1:60]',hist(temp_SE(inrange(temp_SE,[0 60])),[0:1:60])','gauss1');
        SE_Size_xy=SE_fit_xy.b1;
        SE_sig_xy=SE_fit_xy.c1/sqrt(2);
        SE_conf_xy=confint(SE_fit_xy,0.683);
        SE_sig_mu_xy=abs(SE_Size_xy-SE_conf_xy(1,2)); % one sigma
        SE_sig_sig_xy=abs(SE_sig_xy-SE_conf_xy(1,3)/sqrt(2)); % one sigma
    
        %top PMT Second fit
        n_sig_fit=2;
    
        fit_start=SE_Size_xy-n_sig_fit*SE_sig_xy;
        fit_end=SE_Size_xy+n_sig_fit*SE_sig_xy;
        xfit=fit_start:0.1:fit_end;
    
        [xxo, yyo, xo, no] = step(temp_SE(inrange(temp_SE,[fit_start fit_end])),[0:0.2:60]);
        s = fitoptions('Method','NonlinearLeastSquares','Robust','On','Algorithm','Trust-Region','Lower',[-5,0,0,0],'Upper',[5,Inf,Inf,Inf],...
        'Startpoint',[1,max(no),10,5]);
        skewGauss1 = fittype('a1/c1*exp(-((x-b1)/c1)^2)*(1+erf(a*(x-b1)/c1))','options',s);

    [f_xy,g_xy] = fit(xo',no',skewGauss1);
    delta_xy = f_xy.a/sqrt(1+f_xy.a^2);
    SE_mu_xy = f_xy.b1 + f_xy.c1/sqrt(2)*delta_xy*sqrt(2/pi);
    SE_sig_xy = sqrt((f_xy.c1/sqrt(2))^2 * (1 - 2 * (delta_xy)^2/pi));
    skew_xy = (4-pi)/2 * (delta_xy * sqrt(2/pi))^3 / (1 - 2 * delta_xy^2/pi)^(3/2);
    kurt_xy = 2*(pi-3) * (delta_xy * sqrt(2/pi))^3 / (1 - 2 * delta_xy^2/pi)^2;
    ci_xy = confint(f_xy,0.683);
    y_fit_xy = f_xy(xfit);   
      
    mean_SE_top_xy(y_count,x_count)=SE_mu_xy;  
    skew_SE_top_xy(y_count,x_count)=skew_xy; 
    sigma_SE_top_xy(y_count,x_count)=SE_sig_xy;
    kurt_SE_top_xy(y_count,x_count)=kurt_xy;
    mean_err_SE_top_xy(y_count,x_count)=SE_sig_mu_xy;
    sigma_err_SE_top_xy(y_count,x_count)=SE_sig_sig_xy;    
    

      else
          
        mean_SE_top_xy(y_count,x_count)=0;  
        skew_SE_top_xy(y_count,x_count)=0; 
        sigma_SE_top_xy(y_count,x_count)=0;
        kurt_SE_top_xy(y_count,x_count)=0;
        mean_err_SE_top_xy(y_count,x_count)=0;
        sigma_err_SE_top_xy(y_count,x_count)=0;    
    
      end           
    end    
  end

%%Plot the mean SE_XY top %%%%%%%%%%%%%%%%
color_range_max=max(max(mean_SE_top_xy));
color_range_min=min(min(mean_SE_top_xy(mean_SE_top_xy>0)));
vc_step=(color_range_max-color_range_min)/50;
vc=color_range_min:vc_step:color_range_max;
 
SE_xy_top_fig = figure;
contourf(se_xbins,se_ybins,mean_SE_top_xy,vc,'LineColor','none');
xlabel('x (cm)','fontsize', 18);
ylabel('y (cm)','fontsize', 18);
title(strcat(file_id_cp, '. SE Mean (top PMT) vs. XY.'),'fontsize',16,'Interpreter','none');
caxis([color_range_min color_range_max])
g=colorbar;
ylabel(g,'phe','FontSize',16)
set(g,'FontSize',18)
hold on
LUXPlotTopPMTs
axis([-25 25 -25 25])
myfigview(16);


%%Plot the 1 sigma of the mean

color_range_max=max(max(mean_err_SE_top_xy));
color_range_min=min(min(mean_err_SE_top_xy(mean_err_SE_top_xy>0)));
vc_step=(color_range_max-color_range_min)/50;
vc=color_range_min:vc_step:color_range_max;
 
sigma_SE_xy_top_fig = figure;
contourf(se_xbins,se_ybins,mean_err_SE_top_xy,vc,'LineColor','none');
xlabel('x (cm)','fontsize', 18);
ylabel('y (cm)','fontsize', 18);
title(strcat(file_id_cp, '. 1-Sigma SE Mean(top PMT) vs. XY.'),'fontsize',16,'Interpreter','none');
caxis([color_range_min color_range_max])
g=colorbar;
ylabel(g,'phe','FontSize',16)
set(g,'FontSize',18)
hold on
LUXPlotTopPMTs
axis([-25 25 -25 25])
myfigview(16);

%Calculate corrections matrix and 1 sigma corrections matrix for SE_top
% if inpaint_on==1;  mean_SE_top_xy(mean_SE_top_xy==0)=nan; mean_SE_top_xy=inpaint_nans(mean_SE_top_xy,3); end;
SE_top_mu_center=interp2(se_xbins,se_ybins,mean_SE_top_xy,x_center,y_center,'spline');%Normalize to the center (x=y=0)
norm_SE_top_mu=SE_top_mu_center./mean_SE_top_xy; 
norm_SE_top_mu(isinf(norm_SE_top_mu))=1;%no infinity! Do not correct events outside r=25cm
norm_SE_top_mu(isnan(norm_SE_top_mu))=1;

%Calculate error on SE center
% if inpaint_on==1;  mean_err_SE_top_xy(mean_err_SE_top_xy==0)=nan; mean_err_SE_top_xy=inpaint_nans(mean_err_SE_top_xy,3); end;
error_SE_top_mu_center=interp2(se_xbins,se_ybins,mean_err_SE_top_xy,x_center,y_center,'spline');%Normalize to the center (x=y=0)
error_SE_top_mu_norm=sqrt((error_SE_top_mu_center./mean_SE_top_xy).^2+(mean_err_SE_top_xy.*SE_top_mu_center./mean_SE_top_xy.^2).^2);
error_SE_top_mu_norm(isinf(error_SE_top_mu_norm))=1;%no infinity. Set Sigma=1, which is 100% uncertainty
error_SE_top_mu_norm(isnan(error_SE_top_mu_norm))=1;%no nan. Set Sigma=1, which is 100% uncertainty

%Calculate kurt,sigma, and skew at center
% if inpaint_on==1;  kurt_SE_top_xy(kurt_SE_top_xy==0)=nan; kurt_SE_top_xy=inpaint_nans(kurt_SE_top_xy,3); end;
% if inpaint_on==1;  skew_SE_top_xy(skew_SE_top_xy==0)=nan; skew_SE_top_xy=inpaint_nans(skew_SE_top_xy,3); end;
% if inpaint_on==1;  sigma_SE_top_xy(sigma_SE_top_xy==0)=nan; sigma_SE_top_xy=inpaint_nans(sigma_SE_top_xy,3); end;
SE_top_kurt_center=interp2(se_xbins,se_ybins,kurt_SE_top_xy,x_center,y_center,'spline');%Normalize to the center (x=y=0)
SE_top_skew_center=interp2(se_xbins,se_ybins,skew_SE_top_xy,x_center,y_center,'spline');%Normalize to the center (x=y=0)
SE_top_sigma_center=interp2(se_xbins,se_ybins,sigma_SE_top_xy,x_center,y_center,'spline');%Normalize to the center (x=y=0)



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% Calculating the S1 Z-dependence (Both PMT arrays)
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    
    clear y_s1_both y_s1_bottom Fit_s1_both Fit_s1_bottom mean_s1_both mean_s1_both
    clear sigma_s1_both sigma_s1_both mean_index means mean_index means_error x temp_hist
    clear hist_s1_both x xfit S
    
    %set up energy bins for the gaussian fit
    bin_s1=8;%5
    s1_bin_min=100;
    s1_bin_max=500;
    x=s1_bin_min:bin_s1:s1_bin_max; %150:bin_s1:500;
    x=x';
        
    cut = inrange(s2radius,[rcut_min,rcut_max]) & inrange(s1_phe_both,[s1_bin_min,s1_bin_max]) & inrange(s2_width,[100,1000]) & inrange(s1_width,[30,400]);   
    
    A = sort(size(s1_phe_both(cut))); %%sort, incase 1xn or nx1 vector is defined.
    bins = ceil(A(2)./min_evts_per_bin); %get the number of points, adjust binning if stats are low.
    if bins>max_num_bins
       bins = max_num_bins; % usually we will have 65 bin. 0:5:325
    end
%     bin_size = (488/bins);%span the entire detector. ~49.0 cm. Bin centers every 4\mus
     bin_size = (floor(det_edge/bins*4)/4);%span the entire detector. 325us. typically = 5us. Round to nearest 0.25.

   
    hist_s1_both = zeros(length(x),bins);
    Fit_s1_both = cell(1,bins);
    means = zeros(bins,1);
    means_error = zeros(bins,1);
    mean_index = zeros(bins,1);
    
      
    for bin=1:bins

     binStart = ((bin-1).*bin_size);  binEnd = (bin.*bin_size);   

%      time_cut = inrange(drift_distance_\mus,[binStart,binEnd]) & cut==1;
       time_cut = inrange(drift_time,[binStart,binEnd]) & cut==1;
     
        hist_s1_both(:,bin)= hist(s1_phe_both(time_cut),x)'/bin_s1;
        temp_hist=hist_s1_both(:,bin);   
            amp_start=max(temp_hist);
            mean_start=sum(temp_hist.*x)/sum(temp_hist);
            sigma_start=std(temp_hist.*x);
        Fit_s1_both{bin}=fit(x(2:end-1),temp_hist(2:end-1),'gauss1');%remove bin edges

     means(bin)=Fit_s1_both{bin}.b1;
     means_error(bin) = Fit_s1_both{bin}.c1/sqrt(2)/sqrt(sum(hist_s1_both(:,bin))*bin_s1); % == sigma/sqrt(N);
     mean_index(bin) = (binStart+binEnd)/2;     

    end
    
 % //////////////////////////////Make Plots/////////////////////////////////
 
    dT_range= inrange( mean_index,[zcut_min, zcut_max] );
    dT_fit=0:1:det_edge;    
    [P_s1_both S] = polyfit(mean_index(dT_range),means(dT_range),2);
    sigma_P_s1_both = sqrt(diag(inv(S.R)*inv(S.R')).*S.normr.^2./S.df)';    
    
    light_correction_fig = figure;
    errorbar(mean_index,means,means_error,'.k');
    hold on     
    [yplot_s1_both] = polyval(P_s1_both,dT_fit);
    [center_s1_both, center_s1_both_err] = polyval(P_s1_both,z_center, S); %160[us]*1.51[\mus/us]
    plot(dT_fit,yplot_s1_both,'-r','LineWidth',1);    
    xlabel('Depth (\mus)', 'FontSize',16); ylabel('Mean (Pulse Area Phe)', 'FontSize', 16);
    title(strcat(file_id_cp,'.Kr83m. S1 Mean vs. Z'), 'FontSize', 16,'Interpreter','none');
    line([0 det_edge],[center_s1_both center_s1_both],'linestyle','--');
    line([z_center; z_center;],[min(means) max(means)],'linestyle','--');
    legend('S1\_both Mean', strcat('y=', num2str(P_s1_both(1),'%10.2e'), '*Z^2 + ', num2str(P_s1_both(2),'%10.2e')...
        , '*Z + ', num2str(P_s1_both(3),4)),strcat('Det center = ',num2str(center_s1_both,4), ' [Phe]'),'location','northwest')
    
    box on;
    myfigview(16);
    xlim([0 det_edge]);
    ylim([0.8*min(means) 1.15*max(means)]);
  
    %Plot 7 of the histograms     
    hist_fig = figure;
    hold on
    plot_step = floor ((bins-1)/7);    
    j=1;
    t = {num2str(zeros(7,1))};
    for i=2:plot_step:bins  %first bin contains anode... skip it. start at i=2.
        if j<= 7
        xfit=100:2:500;
        yfit=Fit_s1_both{i}.a1.*exp(-((xfit-Fit_s1_both{i}.b1)./Fit_s1_both{i}.c1).^2);
            
            if j==1
                mark = '+'; color = 'k'; t{j}=num2str(mean_index(i),4);
            elseif j==2
                mark = '.'; color = 'b'; t{j}=num2str(mean_index(i),4);    
            elseif j==3
                mark = '*'; color = 'r'; t{j}=num2str(mean_index(i),4);
            elseif j==4
                mark = 'x'; color = 'g'; t{j}=num2str(mean_index(i),4);
            elseif j==5
                mark = '+'; color = 'm'; t{j}=num2str(mean_index(i),4);
            elseif j==6
                mark = '.'; color = 'c'; t{j}=num2str(mean_index(i),4);                
            elseif j==7
                mark = '*'; color = 'y'; t{j}=num2str(mean_index(i),4);                
            end                  
                
        plot(x,hist_s1_both(:,i),strcat(mark,color));
        plot(xfit,yfit,color);
        end
        
        j=j+1;
    end
        title(strcat(file_id_cp,'.Kr83m. S1_both Hist Data'), 'FontSize', 16,'Interpreter','none');
        legend(strcat(t{1},' \mus'),strcat(t{1},' \mus Fit'),strcat(t{2},' \mus'),strcat(t{2},' \mus Fit')...
            ,strcat(t{3},' \mus'),strcat(t{3},' \mus Fit'),strcat(t{4},' \mus'),strcat(t{4},' \mus Fit')...
            ,strcat(t{5},' \mus'),strcat(t{5},' \mus Fit'),strcat(t{6},' \mus'),strcat(t{6},' \mus Fit')...
            ,strcat(t{7},' \mus'),strcat(t{7},' \mus Fit') ); 
        xlabel('Pulse Area Phe','FontSize',16); ylabel('Counts/Pulse-Area-Phe','FontSize',16);
    myfigview(16);
    xlim([s1_bin_min-100 s1_bin_max+100]);
    box on;
    
    s1_both_means=means;
    s1_both_means_sigma=means_error;
    s1_both_means_bin=mean_index;
    
    
%%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% Calculating the S1 Z-dependence (Bottom PMT array)
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    clear y_s1_bottom y_s1_bottom Fit_s1_bottom Fit_s1_bottom mean_s1_bottom mean_s1_bottom
    clear sigma_s1_both sigma_s1_bottom mean_index means mean_index means_error x temp_hist
    clear hist_s1_both dT_range dT_fit S x xfit
    
    %set up energy bins for the gaussian fit
    bin_s1=8;
    s1_bin_min=50;
    s1_bin_max=400;
    x=s1_bin_min:bin_s1:s1_bin_max;
    x=x';
        
    
    cut = inrange(s2radius,[rcut_min,rcut_max]) & inrange(s1_phe_bottom,[s1_bin_min,s1_bin_max]) & inrange(s2_width,[100,1000]) & inrange(s1_width,[30,400]);
    A = sort(size(s1_phe_bottom(cut))); %%sort, incase 1xn or nx1 vector is defined.
    bins = ceil(A(2)./min_evts_per_bin); %get the number of points, adjust binning if stats are low.
    if bins>max_num_bins
       bins = max_num_bins; % usually we will have 65 bin. 0:5:325
    end
%     bin_size = (488/bins);%span the entire detector. 325us. typically = 5us
    bin_size = (floor(det_edge/bins*4)/4);%span the entire detector. 325us. typically = 5us. Round to nearest 0.25.
    
        hist_s1_bottom = zeros(length(x),bins);
        Fit_s1_bottom = cell(1,bins);
        means = zeros(bins,1);
        means_error = zeros(bins,1);
        mean_index = zeros(bins,1);

    for bin=1:bins

     binStart = ((bin-1).*bin_size);  binEnd = (bin.*bin_size);   

%      time_cut = inrange(drift_distance_\mus,[binStart,binEnd]) & cut==1;
       time_cut = inrange(drift_time,[binStart,binEnd]) & cut==1;
       
     %[A Y chisq dof J hb hy] = IterFHGauss(s1area(time_cut)); 
     
        hist_s1_bottom(:,bin)= hist(s1_phe_bottom(time_cut),x)'/bin_s1;
        temp_hist=hist_s1_bottom(:,bin);
            amp_start=max(temp_hist);
            mean_start=sum(temp_hist.*x)/sum(temp_hist);
            sigma_start=std(temp_hist.*x);
        Fit_s1_bottom{bin}=fit(x(2:end-1),temp_hist(2:end-1),'gauss1');%remove bin edges
     
      
     means(bin)=Fit_s1_bottom{bin}.b1;
     means_error(bin) = Fit_s1_bottom{bin}.c1/sqrt(2)/sqrt(sum(hist_s1_bottom(:,bin))*bin_s1); % == sigma/sqrt(N);
     mean_index(bin) = (binStart+binEnd)/2;     

    end
    
 % //////////////////////////////Make Plots/////////////////////////////////
 
    dT_range= inrange( mean_index,[zcut_min, zcut_max] );
    dT_fit=0:1:det_edge;    
    [P_s1_bottom S] = polyfit(mean_index(dT_range),means(dT_range),2);
    sigma_P_s1_bottom = sqrt(diag(inv(S.R)*inv(S.R')).*S.normr.^2./S.df)';
 
    light_correction_fig_bottom = figure;
    errorbar(mean_index,means,means_error,'.k');
    hold on     
    [yplot_s1_bottom] = polyval(P_s1_bottom,dT_fit);
    [center_s1_bottom, center_s1_bottom_err] = polyval(P_s1_bottom,z_center,S); %160[us]*1.51[\mus/us]
    plot(dT_fit,yplot_s1_bottom,'-r','LineWidth',1);    
    xlabel('Depth (\mus)', 'FontSize',16); ylabel('Mean (Pulse Area Phe)', 'FontSize', 16);
    title(strcat(file_id_cp,'.Kr83m. S1_bottom Mean vs. Z'), 'FontSize', 16,'Interpreter','none');
    line([0 det_edge],[center_s1_bottom center_s1_bottom],'linestyle','--');
    line([z_center; z_center;],[min(means) max(means)],'linestyle','--');
    legend('S1\_bottom Mean', strcat('y=', num2str(P_s1_bottom(1),'%10.2e'), '*Z^2 + ', num2str(P_s1_bottom(2),'%10.2e')...
        , '*Z + ', num2str(P_s1_bottom(3),4)),strcat('Det center = ',num2str(center_s1_bottom,4), ' [Phe]'),'location','northwest')
    box on;
    myfigview(16);
    xlim([0 det_edge]);
    ylim([0.8*min(means) 1.15*max(means)]);

    %Plot 7 of the histograms     
    hist_s1_bottom_fig = figure;
    hold on
    plot_step = floor((bins-1)/7);    
    j=1;
    t = {num2str(zeros(7,1))};
    for i=2:plot_step:bins  %first bin contains anode... skip it. start at i=2.
        if j<=7
        xfit=50:2:500;
        yfit=Fit_s1_bottom{i}.a1.*exp(-((xfit-Fit_s1_bottom{i}.b1)./Fit_s1_bottom{i}.c1).^2);
            
            if j==1
                mark = '+'; color = 'k'; t{j}=num2str(mean_index(i),4);
            elseif j==2
                mark = '.'; color = 'b'; t{j}=num2str(mean_index(i),4);    
            elseif j==3
                mark = '*'; color = 'r'; t{j}=num2str(mean_index(i),4);
            elseif j==4
                mark = 'x'; color = 'g'; t{j}=num2str(mean_index(i),4);
            elseif j==5
                mark = '+'; color = 'm'; t{j}=num2str(mean_index(i),4);
            elseif j==6
                mark = '.'; color = 'c'; t{j}=num2str(mean_index(i),4);                
            elseif j==7
                mark = '*'; color = 'y'; t{j}=num2str(mean_index(i),4);                
            end                  
                
        plot(x,hist_s1_bottom(:,i),strcat(mark,color));
        plot(xfit,yfit,color);
        end
        j=j+1;
    end
        title(strcat(file_id_cp,'.Kr83m. S1_bottom Hist Data'), 'FontSize', 16,'Interpreter','none');
        legend(strcat(t{1},' \mus'),strcat(t{1},' \mus Fit'),strcat(t{2},' \mus'),strcat(t{2},' \mus Fit')...
            ,strcat(t{3},' \mus'),strcat(t{3},' \mus Fit'),strcat(t{4},' \mus'),strcat(t{4},' \mus Fit')...
            ,strcat(t{5},' \mus'),strcat(t{5},' \mus Fit'),strcat(t{6},' \mus'),strcat(t{6},' \mus Fit')...
            ,strcat(t{7},' \mus'),strcat(t{7},' \mus Fit') ); 
        xlabel('Pulse Area Phe','FontSize',16); ylabel('Counts/Pulse-Area-Phe','FontSize',16);
    myfigview(16);
    xlim([s1_bin_min-50 s1_bin_max+50]);
    box on;
      
    s1_bottom_means=means;
    s1_bottom_means_sigma=means_error;
    s1_bottom_means_bin=mean_index;
     



%%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% Calculating the S1 Z-dependence (Top PMT array)
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dT_range= inrange( mean_index,[zcut_min, zcut_max] );
    
    [P_s1_top S] = polyfit(s1_both_means_bin(dT_range),s1_both_means(dT_range)-s1_bottom_means(dT_range),2);
    sigma_P_s1_top = sqrt(diag(inv(S.R)*inv(S.R')).*S.normr.^2./S.df)';
    [center_s1_top, center_s1_top_err] = polyval(P_s1_top,z_center,S); %160[us]*1.51[\mus/us]



%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%   Calculating the S2 Z-dependence (i.e. electron lifetime)
%   Using Bottom PMT array
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    
    
%     s1_z_correction=polyval(P_s1_both,241.6)./polyval(P_s1_both,drift_distance_\mus);%normalize s1 to 160 ns.
    s1_z_correction=polyval(P_s1_both,z_center)./polyval(P_s1_both,drift_time);
    
    %s1_both_z_correction=interp1(s1_both_means_bin,s1_both_means,160)./interp1(s1_both_means_bin,s1_both_means,drift_time);
    s1_phe_both_z=s1_phe_both.*s1_z_correction;

    clear y_s2_bottom y_s2_bottom Fit_s2_bottom Fit_EL xEL yfitEL xfit
    clear sigma_s2_bottom sigma_s2_bottom mean_index means mean_index means_error x temp_hist
    clear hist_s2_bottom dT_range dT_fit S x2
    
    %set up energy bins for the gaussian fit. Using bottom S2 only
   bin_s2=300;
   s2_bin_min=50;
   s2_bin_max=10000;
   x2=s2_bin_min:bin_s2:s2_bin_max;
   x2=x2';
   
   %x3=2:1:30; %%for S2/S1
   %x3=x3';
        
    
    cut = inrange(s2radius,[rcut_min,rcut_max]) & inrange(s2_phe_bottom,[s2_bin_min,s2_bin_max]) & inrange(s2_width,[100,1000]) & inrange(s1_width,[30,400]) ;
    A = sort(size(s2_phe_bottom(cut))); %%sort, incase 1xn or nx1 vector is defined.
    bins = ceil(A(2)./min_evts_per_bin); %get the number of points, adjust binning if stats are low.
    if bins>max_num_bins
       bins = max_num_bins; % usually we will have 65 bin. 0:5:325
    end
%     bin_size = (488/bins);%span the entire detector. 325us. typically = 5us. 488 \mus
   bin_size = (floor(det_edge/bins*4)/4);%span the entire detector. 325us. typically = 5us. Round to nearest 0.25.

    
        hist_s2_bottom = zeros(length(x2),bins);
        Fit_s2_bottom = cell(1,bins);
        means=zeros(bins,1);
        means_error = zeros(bins,1);
        mean_index= zeros(bins,1);                   
        
    for bin=1:bins

     binStart = ((bin-1).*bin_size);  binEnd = (bin.*bin_size);   

%      time_cut = inrange(drift_distance_\mus,[binStart,binEnd]) & cut==1;
     time_cut = inrange(drift_time,[binStart,binEnd]) & cut==1;
     
        hist_s2_bottom(:,bin)= hist(s2_phe_bottom(time_cut),x2)'/bin_s2;
        temp_hist=hist_s2_bottom(:,bin);
            amp_start=max(temp_hist);
            mean_start=sum(temp_hist.*x2)/sum(temp_hist);
            sigma_start=std(temp_hist.*x2);
        Fit_s2_bottom{bin}=fit(x2(2:end-1),temp_hist(2:end-1),'gauss1');%remove bin edges
          

         means(bin)=Fit_s2_bottom{bin}.b1;
         means_error(bin) = Fit_s2_bottom{bin}.c1/sqrt(2)/sqrt(sum(hist_s2_bottom(:,bin))*bin_s2); % == sigma/sqrt(N);
         mean_index(bin) = (binStart+binEnd)/2;     
     
     %%calculate S2/S1 also, just to double check the lifetime from S2 only
     % hist_s2_bottom_s1(:,bin)= hist(s2_phe_bottom(time_cut)./s1_phe_both_z(time_cut),x3)';
     % temp_hist=hist_s2_bottom_s1(:,bin);   
      %Fit_s2_bottom_s1{bin}=fit(x3(2:end-1),temp_hist(2:end-1),'gauss1');%remove bin edges        
       % means_s2_s1(bin)=Fit_s2_bottom_s1{bin}.b1;
       % means_error_s2_s1(bin) = Fit_s2_bottom_s1{bin}.c1/2/sqrt(sum(hist_s2_bottom_s1(:,bin))); % == sigma/sqrt(N);
       % mean_index_s2_s1(bin) = (binStart+binEnd)/2;     

    end
    
 % //////////////////////////////Make Plots/////////////////////////////////
 
     dT_range= inrange( mean_index,[5,500] );%from 50 to 300us for the fit

 
    s2_bottom_norm_z_index=mean_index(dT_range);
    s2_bottom_norm_z_means=means(dT_range);
    s2_at_top_bottom=interp1(s2_bottom_norm_z_index,s2_bottom_norm_z_means,4,'spline'); %normalize to right below the gate
    s2_bottom_norm_z=interp1(s2_bottom_norm_z_index,s2_bottom_norm_z_means,4,'spline')./s2_bottom_norm_z_means;
  
    %%Tracker for pseudo-lifetime
    s2_at_bottom_bottom=interp1(s2_bottom_norm_z_index,s2_bottom_norm_z_means,320,'spline');
    pseudo_lifetime_bottom=-320/log(s2_at_bottom_bottom/s2_at_top_bottom);
 
    dT_fit=0:1:det_edge;    
    
    Fit_EL= fit(mean_index(dT_range),means(dT_range),'exp1');
    yfitEL = Fit_EL.a.*exp(Fit_EL.b.*dT_fit);
    
    e_lifetime=-1/Fit_EL.b;
    s2_z0_bottom=Fit_EL.a;
       
    b=confint(Fit_EL, 0.683);%1 sigma matrix
    sigma_e_lifetime=(1/b(1,2)-1/b(2,2))/2;
    sigma_s2_z0_bottom=(b(2,1)-b(1,1))/2;
 
    lifetime_fig_bottom = figure;
    errorbar(mean_index(2:end),means(2:end),means_error(2:end),'.k');  %first bin is anode, skip it
    hold on     
    plot(dT_fit,yfitEL,'-r','LineWidth',1);   
    xlabel('Depth (\mus)', 'FontSize',16); ylabel('Mean S2 (Pulse Area Phe)', 'FontSize', 16);
    title(strcat(file_id_cp,'.Kr83m. S2_bottom Electron Lifetime'), 'FontSize', 16,'Interpreter','none');
    legend('S2\_bottom Mean', strcat( '\lambda= ', num2str(e_lifetime,4), ' \pm ', num2str(sigma_e_lifetime,2), ' \mus' )...
        ,'location','northeast');
    myfigview(16);
    xlim([0 det_edge]);
        
    %Plot 7 of the histograms     
    hist_s2_bottom_fig = figure;
    hold on
    plot_step = floor((bins-1)/7);    
    j=1;
    t = {num2str(zeros(7,1))};
    for i=2:plot_step:bins  %first bin contains anode... skip it. start at i=2.
        if j<=7
        xfit=100:50:10000;
        yfit=Fit_s2_bottom{i}.a1.*exp(-((xfit-Fit_s2_bottom{i}.b1)./Fit_s2_bottom{i}.c1).^2);
            
            if j==1
                mark = '+'; color = 'k'; t{j}=num2str(mean_index(i),4);
            elseif j==2
                mark = '.'; color = 'b'; t{j}=num2str(mean_index(i),4);    
            elseif j==3
                mark = '*'; color = 'r'; t{j}=num2str(mean_index(i),4);
            elseif j==4
                mark = 'x'; color = 'g'; t{j}=num2str(mean_index(i),4);
            elseif j==5
                mark = '+'; color = 'm'; t{j}=num2str(mean_index(i),4);
            elseif j==6
                mark = '.'; color = 'c'; t{j}=num2str(mean_index(i),4);                
            elseif j==7
                mark = '*'; color = 'y'; t{j}=num2str(mean_index(i),4);                
            end                  
                
        plot(x2,hist_s2_bottom(:,i)./d.livetime_sec,strcat(mark,color));
        plot(xfit,yfit/d.livetime_sec,color);
        end
        j=j+1;
    end
        title(strcat(file_id_cp,'.Kr83m.S2_Bottom'), 'FontSize', 16,'Interpreter','none');
        legend(strcat(t{1},' \mus'),strcat(t{1},' \mus Fit'),strcat(t{2},' \mus'),strcat(t{2},' \mus Fit')...
            ,strcat(t{3},' \mus'),strcat(t{3},' \mus Fit'),strcat(t{4},' \mus'),strcat(t{4},' \mus Fit')...
            ,strcat(t{5},' \mus'),strcat(t{5},' \mus Fit'),strcat(t{6},' \mus'),strcat(t{6},' \mus Fit')...
            ,strcat(t{7},' \mus'),strcat(t{7},' \mus Fit')); 
        xlabel('Pulse Area Phe','FontSize',16); ylabel('Counts/Phe/sec','FontSize',16);
    myfigview(16);
    xlim([s2_bin_min-50 s2_bin_max+500]);
   
    s2_bottom_means=means;
    s2_bottom_means_sigma=means_error;
    s2_bottom_means_bin=mean_index;      
    

%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%   Calculating the S2 Z-dependence (i.e. electron lifetime)
%   Using Both PMT arrays
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    clear y_s1_both y_s1_both Fit_s1_both Fit_EL xEL yfitEL xfit
    clear sigma_s1_both sigma_s1_both mean_index means mean_index means_error x temp_hist
    clear hist_s1_both dT_range dT_fit S x2
    
    %set up energy bins for the gaussian fit. Using both S2 only
   bin_s2=500;
   s2_bin_min=1000;
   s2_bin_max=40000;
   x2=s2_bin_min:bin_s2:s2_bin_max;
   x2=x2';
        
    cut = inrange(s2radius,[rcut_min,rcut_max]) & inrange(s2_phe_both,[s2_bin_min,s2_bin_max]) & inrange(s2_width,[100,1000]) & inrange(s1_width,[30,400]) ;
    A = sort(size(s2_phe_both(cut))); %%sort, incase 1xn or nx1 vector is defined.
    bins = ceil(A(2)./min_evts_per_bin); %get the number of points, adjust binning if stats are low.
    if bins>max_num_bins
       bins = max_num_bins; % usually we will have 65 bin. 0:5:325
    end
%     bin_size = (488/bins);%span the entire detector. 325us. typically = 5us
    bin_size = (floor(det_edge/bins*4)/4);%span the entire detector. 325us. typically = 5us. Round to nearest 0.25.
   
        hist_s2_both = zeros(length(x2),bins);
        Fit_s2_both = cell(1,bins);
        means=zeros(bins,1);
        means_error = zeros(bins,1);
        mean_index= zeros(bins,1);
        
    for bin=1:bins

     binStart = ((bin-1).*bin_size);  binEnd = (bin.*bin_size);   

%      time_cut = inrange(drift_distance_\mus,[binStart,binEnd]) & cut==1;
       time_cut = inrange(drift_time,[binStart,binEnd]) & cut==1; 
       
        hist_s2_both(:,bin)= hist(s2_phe_both(time_cut),x2)'/bin_s2;
        temp_hist=hist_s2_both(:,bin);
            amp_start=max(temp_hist);
            mean_start=sum(temp_hist.*x2)/sum(temp_hist);
            sigma_start=std(temp_hist.*x2);
        Fit_s2_both{bin}=fit(x2(2:end-1),temp_hist(2:end-1),'gauss1');%remove bin edges
          

         means(bin)=Fit_s2_both{bin}.b1;
         means_error(bin) = Fit_s2_both{bin}.c1/sqrt(2)/sqrt(sum(hist_s2_both(:,bin))*bin_s2); % == sigma/sqrt(N);
         mean_index(bin) = (binStart+binEnd)/2;     
     
     
     %%calculate S2/S1 also, just to double check the lifetime from S2 only
     % hist_s2_both_s1(:,bin)= hist(s2_phe_both(time_cut)./s1_phe_both_z(time_cut),x3)';
     % temp_hist=hist_s2_both_s1(:,bin);   
      %Fit_s2_both_s1{bin}=fit(x3(2:end-1),temp_hist(2:end-1),'gauss1');%remove bin edges        
       % means_s2_s1(bin)=Fit_s2_both_s1{bin}.b1;
       % means_error_s2_s1(bin) = Fit_s2_both_s1{bin}.c1/2/sqrt(sum(hist_s2_both_s1(:,bin))); % == sigma/sqrt(N);
       % mean_index_s2_s1(bin) = (binStart+binEnd)/2;     

    end
    
 % //////////////////////////////Make Plots/////////////////////////////////
 
  
     dT_range= inrange( mean_index,[5,500] );%from 50 to 300us for the fit


    s2_both_norm_z_index=mean_index(dT_range);
    s2_both_norm_z_means=means(dT_range);
    s2_at_top_both=interp1(s2_both_norm_z_index,s2_both_norm_z_means,4,'spline'); %normalize to right below the gate
    s2_both_norm_z=interp1(s2_both_norm_z_index,s2_both_norm_z_means,4,'spline')./s2_both_norm_z_means;
  
    %%Tracker for pseudo-lifetime
    s2_at_bottom_both=interp1(s2_both_norm_z_index,s2_both_norm_z_means,320,'spline');
    pseudo_lifetime_both=-320/log(s2_at_bottom_both/s2_at_top_both);
 
    dT_fit=0:1:det_edge;    
    
    Fit_EL= fit(mean_index(dT_range),means(dT_range),'exp1');
    yfitEL = Fit_EL.a.*exp(Fit_EL.b.*dT_fit);
    
    e_lifetime_both=-1/Fit_EL.b;
    s2_z0_both=Fit_EL.a;
       
    b=confint(Fit_EL, 0.683);%1 sigma matrix
    sigma_e_lifetime_both=(1/b(1,2)-1/b(2,2))/2;
    sigma_s2_z0_both=(b(2,1)-b(1,1))/2;
 
    lifetime_fig_both = figure;
    errorbar(mean_index(2:end),means(2:end),means_error(2:end),'.k');  %first bin is anode, skip it
    hold on     
    plot(dT_fit,yfitEL,'-r','LineWidth',1);   
    xlabel('Time (\mus)', 'FontSize',16); ylabel('Mean S2 (Pulse Area Phe)', 'FontSize', 16);
    title(strcat(file_id_cp,'.Kr83m. S2_both Electron Lifetime'), 'FontSize', 16,'Interpreter','none');
    legend('S2\_both Mean', strcat( '\lambda= ', num2str(e_lifetime_both,4), ' \pm ', num2str(sigma_e_lifetime_both,2), ' \mus' )...
        ,'location','northeast');
    myfigview(16);
    xlim([0 det_edge]);
   
    %Plot 7 of the histograms     
    hist_s2_both_fig = figure;
    hold on
    plot_step = floor((bins-1)/7);    
    j=1;
    t = {num2str(zeros(7,1))};
    for i=2:plot_step:bins  %first bin contains anode... skip it. start at i=2.
        if j<=7
        xfit=s2_bin_min:bin_s2:s2_bin_max;
        yfit=Fit_s2_both{i}.a1.*exp(-((xfit-Fit_s2_both{i}.b1)./Fit_s2_both{i}.c1).^2);
            
            if j==1
                mark = '+'; color = 'k'; t{j}=num2str(mean_index(i),4);
            elseif j==2
                mark = '.'; color = 'b'; t{j}=num2str(mean_index(i),4);    
            elseif j==3
                mark = '*'; color = 'r'; t{j}=num2str(mean_index(i),4);
            elseif j==4
                mark = 'x'; color = 'g'; t{j}=num2str(mean_index(i),4);
            elseif j==5
                mark = '+'; color = 'm'; t{j}=num2str(mean_index(i),4);
            elseif j==6
                mark = '.'; color = 'c'; t{j}=num2str(mean_index(i),4);                
            elseif j==7
                mark = '*'; color = 'y'; t{j}=num2str(mean_index(i),4);                
            end                  
                
        plot(x2,hist_s2_both(:,i)./d.livetime_sec,strcat(mark,color));
        plot(xfit,yfit./d.livetime_sec,color);
        end
        j=j+1;
    end
        title(strcat(file_id_cp,'.Kr83m. S2_both'), 'FontSize', 16,'Interpreter','none');
        legend(strcat(t{1},' \mus'),strcat(t{1},' \mus Fit'),strcat(t{2},' \mus'),strcat(t{2},' \mus Fit')...
            ,strcat(t{3},' \mus'),strcat(t{3},' \mus Fit'),strcat(t{4},' \mus'),strcat(t{4},' \mus Fit')...
            ,strcat(t{5},' \mus'),strcat(t{5},' \mus Fit'),strcat(t{6},' \mus'),strcat(t{6},' \mus Fit')...
            ,strcat(t{7},' \mus'),strcat(t{7},' \mus Fit') ); 
        xlabel('Pulse Area Phe','FontSize',16); ylabel('Counts/Phe/sec','FontSize',16);
    myfigview(16);
    xlim([s2_bin_min-50 s2_bin_max+500]);
   
    s2_both_means=means;
    s2_both_means_sigma=means_error;
    s2_both_means_bin=mean_index;      
    
%% Getting S2 both (z_corrected) mean and error
% s2_phe_both_z=s2_phe_both.*exp(drift_time/e_lifetime_both);
s2_phe_both_z=s2_phe_both.*interp1(s2_both_norm_z_index,s2_both_norm_z, drift_time,'spline',1);
        Fit_s2_both_z=fit(x2,hist(s2_phe_both_z(inrange(s2_phe_both,[s2_bin_min,s2_bin_max])),x2)'/bin_s2,'gauss1');%remove bin edges
        s2_both_center_z=Fit_s2_both_z.b1;
        temp_conf=confint(Fit_s2_both_z);
        s2_both_center_error_z=abs(s2_both_center_z-temp_conf(1,2));
    
    
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%  Mapping S2 in XY (after correcting for Z-dep). Both PMT
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    clear x2 hist_bin

    %set up energy bins for the gaussian fit. Using both S2 only
   bin_s2=500;
   s2_bin_min=1000;
   s2_bin_max=40000;
   x2=s2_bin_min:bin_s2:s2_bin_max;
   x2=x2';

%set up matricies to be filled
mean_S2_both_xy = zeros(floor(s2_ybins),floor(s2_xbins));
Sigma_S2_both = zeros(floor(s2_ybins),floor(s2_xbins));
Count_S2_both = zeros(floor(s2_ybins),floor(s2_xbins));

time_energy_cut = inrange(drift_time,[zcut_min,zcut_max]) & inrange(s1_phe_both_z,[100,s1area_bound_max]) & inrange(s2_phe_both_z,[s2_bin_min,s2_bin_max]) & inrange(s2radius,[rcut_min,rcut_max]) ...
    & inrange(s2_width,[100,1000]) & inrange(s1_width,[30,400]);

%%%%variables for uber turbo boost!
s2_phe_both_z_2=s2_phe_both_z(time_energy_cut);
s2x_2=s2x(time_energy_cut);
s2y_2=s2y(time_energy_cut);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   for x_bin=S2_xbin_min:S2xybinsize:(S2_xbin_max-S2xybinsize);
      %fprintf('X = %d \n',x_bin)
    for y_bin=S2_ybin_min:S2xybinsize:(S2_ybin_max-S2xybinsize);          
      x_min = x_bin; x_max = x_bin+S2xybinsize;
      y_min = y_bin; y_max = y_bin+S2xybinsize;      
      
        x_count=int32(1+x_bin/S2xybinsize+S2_xbin_max/S2xybinsize);
        y_count=int32(1+y_bin/S2xybinsize+S2_ybin_max/S2xybinsize); 
      
      
        bin_cut = inrange(s2x_2,[x_min,x_max]) & inrange(s2y_2,[y_min,y_max]);            
      
       hist_bin = hist(s2_phe_both_z_2(bin_cut) ,x2)'/bin_s2; 
       Count_S2_both(y_count,x_count)=sum(hist_bin)*bin_s2; % X and Y flip ... because MATLAB
       
       %%%%%%%%Truncate variables after binning is complete
        s2_phe_both_z_2=s2_phe_both_z_2(~bin_cut);
        s2x_2=s2x_2(~bin_cut);
        s2y_2=s2y_2(~bin_cut);
        clear bin_cut;
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
      if Count_S2_both(y_count,x_count)>30;                         
                
                   amp_start=max(hist_bin);
                   mean_start=sum(hist_bin.*x2)/sum(hist_bin);
                   sigma_start=std(hist_bin.*x2);    
                Fit_S2_both_xy= fit(x2(2:end-1),hist_bin(2:end-1),'gauss1'); %remove edges
                
                %Peak location (mean)
                mean_S2_both_xy(y_count,x_count)=Fit_S2_both_xy.b1; % X and Y flip ... because MATLAB
                
                %1-sigma of peak position
                Sig_b=Fit_S2_both_xy.c1/sqrt(2)/sqrt(Count_S2_both(y_count,x_count)); %sigma/sqrt(N)
   

                    %Fit error checking
                        if(isnan(Sig_b))
                         mean_S2_both_xy(y_count,x_count)=0;
                         Sig_b=0;
                        end

                    %uncertainty in mean
                    Sigma_S2_both(y_count,x_count)=Sig_b;% X and Y flip ... because MATLAB
       
      else
          
        mean_S2_both_xy(y_count,x_count)=0;  %don't so gaussin fit if less than 30 events
       
      end             
    end     
   end
   
   mean_S2_both_xy(isnan(mean_S2_both_xy)) = 0;
       
   
%%Plot the mean S2_XY both %%%%%%%%%%%%%%%%

s2xbins=(S2_xbin_min+S2xybinsize/2):S2xybinsize:(S2_xbin_max-S2xybinsize/2);
s2ybins=(S2_ybin_min+S2xybinsize/2):S2xybinsize:(S2_ybin_max-S2xybinsize/2);

color_range_max=max(max(mean_S2_both_xy));
color_range_min=min(min(mean_S2_both_xy(mean_S2_both_xy>0)));
vc_step=(color_range_max-color_range_min)/50;

vc=color_range_min:vc_step:color_range_max;
 
s2_xy_both_fig = figure;
contourf(s2xbins,s2ybins,mean_S2_both_xy,vc,'LineColor','none');
xlabel('x (cm)','fontsize', 18);
ylabel('y (cm)','fontsize', 18);
title(strcat(file_id_cp, '. S2 Mean (Both PMT) vs. XY.'),'fontsize',16,'Interpreter','none');
caxis([color_range_min color_range_max])
g=colorbar;
ylabel(g,'phe','FontSize',16)
set(g,'FontSize',18)
hold on
LUXPlotTopPMTs
axis([-25 25 -25 25])
myfigview(16);

%%Plot the 1 sigma of the mean

color_range_max=max(max(Sigma_S2_both));
color_range_min=min(min(Sigma_S2_both(Sigma_S2_both>0)));
vc_step=(color_range_max-color_range_min)/50;

vc=color_range_min:vc_step:color_range_max;
 
sigma_s2_xy_both_fig = figure;
contourf(s2xbins,s2ybins,Sigma_S2_both,vc,'LineColor','none');
xlabel('x (cm)','fontsize', 18);
ylabel('y (cm)','fontsize', 18);
title(strcat(file_id_cp, '. 1-Sigma S2 Mean(Both PMT) vs. XY.'),'fontsize',16,'Interpreter','none');
caxis([color_range_min color_range_max])
g=colorbar;
ylabel(g,'phe','FontSize',16)
set(g,'FontSize',18)
hold on
LUXPlotTopPMTs
axis([-25 25 -25 25])
myfigview(16);

%Calculate corrections matrix and 1 sigma corrections matrix
center=interp2(s2xbins,s2ybins,mean_S2_both_xy,x_center,y_center,'spline');%Normalize to the center (x=y=0)
norm_S2_all=center./mean_S2_both_xy; 
norm_S2_all(isinf(norm_S2_all))=1;%no infinity! Do not correct events outside r=25cm
norm_S2_all(isnan(norm_S2_all))=1;

sigma_center=interp2(s2xbins,s2ybins,Sigma_S2_both,x_center,y_center,'spline');%Normalize to the center (x=y=0)
sigma_norm_S2_all=sqrt((sigma_center./mean_S2_both_xy).^2+(Sigma_S2_both.*center./mean_S2_both_xy.^2).^2);
sigma_norm_S2_all(isinf(sigma_norm_S2_all))=1;%no infinity. Set Sigma=1, which is 100% uncertainty
sigma_norm_S2_all(isnan(sigma_norm_S2_all))=1;%no nan. Set Sigma=1, which is 100% uncertainty


s2_both_center_xyz=center;
s2_both_center_xyz_error=sigma_center;

s2xy_correction = interp2(s2xbins,s2ybins,norm_S2_all,s2x,s2y,'spline',1);
