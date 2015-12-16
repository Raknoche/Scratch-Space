%Result for DP2.1 : g2_bot= 6.47, g2_both=15.0

dir_path='C:\Users\Richard\Desktop\LUX Analysis\';
user_name='RichardKnoche';
dp_version='2.0';
algorithm_name='LUXkrypCal';
use_ccv_flag=1;

folders_to_load{1}='lux10_20140903T1918_cp12334';

 rqs_to_load = {'event_number','pulse_area_phe','event_timestamp_samples','x_cm','y_cm','chi2'...
   ,'prompt_fraction','aft_t1_samples','pulse_start_samples','pulse_end_samples'...
   , 'top_bottom_asymmetry','top_bottom_ratio','pulse_classification' ...
   ,'aft_t0_samples','aft_t2_samples','full_evt_area_phe'...
  ,'admin','z_drift_samples' , 's1s2_pairing','golden','x_corrected','y_corrected'};

%


file_id_cp='lux10_20140903T1918_cp11020';
path=strcat(dir_path,'/',folders_to_load{1});         
d = LUXLoadMultipleRQMs_framework(path,rqs_to_load);                
%           
% %Correcting x and y
%  myname = 'Corrections_PositionCorrection';           
%  position_correction_path = 'C:\Program Files\MATLAB\R2012a\bin\LuxCode\Trunk\DataProcessing\MatlabModules\Corrections_PositionCorrection\Corrections_PositionCorrection.m';
%  IQs = dir([position_correction_path(1:(end-numel(myname)-2)) '*Mercury.mat']);
%  table_corrections = load(IQs(end).name);
% 
%  d = Corrections_PositionCorrection_FunctionUseThis(d,table_corrections);
%%  Defining cuts and variables
clearvars -except d dir_path user_name dp_version algorithm_name use_ccv_flag folders_to_load file_id_cp rqs_to_load path d 

submit=0;
delay_min=0;
grid_size=2; %cm XY plane

% defining cuts
% zcut_min = 50;%\mus NOW DEFINED AFTER EDGE IS FOUND
% zcut_max = 460;%\mus

rcut_min = 0;
rcut_max = 25;%cm
s1area_bound_min = 20;%100
s1area_bound_max = 500;%600

s2area_bound_min =    5600;%set to this to remove CH3T, used to be 200
s2area_bound_max = 30000;%30000

min_evts_per_bin = 200;
max_num_bins = 65;

S1xybinsize = grid_size;%cm
S1_xbin_min = -25;
S1_xbin_max = 25;
S1_ybin_min = -25;
S1_ybin_max = 25;

S2xybinsize = grid_size;%cm
S2_xbin_min = -25;
S2_xbin_max = 25;
S2_ybin_min = -25;
S2_ybin_max = 25;


% Pulse Classification
 Pulse_Class=1;
  
   
    d.s1s2_pairing(isnan(d.s1s2_pairing))=0; %by default it returns nan when there is no pair
    d.z_drift_samples(isnan(d.z_drift_samples)) = 0.0; % get rid of NaN   
    
    s1s2_pairing = d.s1s2_pairing;
    num_paired_pulses = sum(s1s2_pairing); % When num_paired_pulses == 2 there is an S1 paired with an S2.
%     num_S1_before_S2 = sum(d.s1_before_s2); % We only want 1 S1 before the S2.
    
    
    events=sort(size(d.pulse_area_phe)); %The first element is the number of pulses. sometimes it's 5, sometimes it's 10
          
    s1_area_cut= inrange(d.pulse_area_phe,[s1area_bound_min,s1area_bound_max]);
    s2_area_cut= inrange(d.pulse_area_phe,[s2area_bound_min,s2area_bound_max]);
    
    s1_class=(d.pulse_classification==1 )& s1_area_cut ; %area cut for Kr events
    s2_class=(d.pulse_classification==2) & s2_area_cut ;   
    s4_class=(d.pulse_classification==4) ;

                           
                         s1_before_s2_cut=   repmat(sum(s1_class,1)==1,events(1),1) & repmat(sum(s2_class,1)==1,events(1),1)...
                            & cumsum(s1_class+s2_class)<=2 ; % ignore all S1 that happen after the s2 !

    s1_single_cut = s1_class & s1_before_s2_cut ; %
    s2_single_cut = s2_class & s1_before_s2_cut;
    SE_good_cut = s4_class & cumsum(s1_before_s2_cut& s1_single_cut)-cumsum(s1_before_s2_cut & s2_single_cut) ;
    
    drift_time = d.z_drift_samples(s2_single_cut)/100;  % us
        
    d.phe_bottom=d.pulse_area_phe./(1+d.top_bottom_ratio); %bottom PMT pulse area
    
    s1_phe_both = d.pulse_area_phe(s1_single_cut);
    s1_phe_bottom = d.phe_bottom(s1_single_cut);
    
    s2_phe_both = d.pulse_area_phe(s2_single_cut);
    s2_phe_bottom = d.phe_bottom(s2_single_cut);
    
    s2x = d.x_corrected(s2_single_cut);
    s2y = d.y_corrected(s2_single_cut);   

    d.livetime_sec=sum(d.livetime_end_samples-d.livetime_latch_samples)/1e8;
    evt_cut=logical(sum(s2_single_cut));%Cut for all the events passing the single S1 & S2 cut
    event_number=d.event_number(evt_cut);
    event_timestamp_samples=d.event_timestamp_samples(evt_cut);
    
    time_wait_cut=event_timestamp_samples/1e8/60 > delay_min; %allow x min for Kr mixing
    
    s2_width=(d.aft_t2_samples(s2_single_cut)-d.aft_t0_samples(s2_single_cut)); %cut above 800 samples
        %s1_width=d.aft_t2_samples(s1_single_cut)-d.aft_t0_samples(s1_single_cut); %length of S1, cut above 150 samples
    s1_width=d.pulse_end_samples(s1_single_cut)-d.pulse_start_samples(s1_single_cut);
    
    
    clean_cut=-(s1_phe_both+s2_phe_both)' + d.full_evt_area_phe(logical(squeeze(sum(s2_single_cut,1)))) < 500 ...
        & drift_time' > 0;%The WIMP clean cut at <100 Phe would remove 80% of the Kr83 events. 500 only removes 0.3%
    clean_cut=clean_cut' & time_wait_cut';
   
         s1_phe_bottom=s1_phe_bottom(clean_cut);
         s1_phe_both=s1_phe_both(clean_cut);
         s2_phe_bottom=s2_phe_bottom(clean_cut);
         s2_phe_both=s2_phe_both(clean_cut);
         drift_time=drift_time(clean_cut);
         s2x=s2x(clean_cut);
         s2y=s2y(clean_cut);
         s2radius = (s2x.^2+s2y.^2).^(0.5);
         event_timestamp_samples=event_timestamp_samples(clean_cut);
         s2_width=s2_width(clean_cut);
         s1_width=s1_width(clean_cut);

   %%%%%%%%%%%%%
   
kr_energy_cut = inrange(s1_phe_both,[s1area_bound_min s1area_bound_max]) & inrange(s2_phe_both,[s2area_bound_min s2area_bound_max]);%for counting Kr events
Kr_events=length(drift_time(kr_energy_cut));% Are there enough events for the XY or XYZ map?


 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

kr_energy_cut = inrange(s1_phe_both,[s1area_bound_min s1area_bound_max]) & inrange(s2_phe_both,[s2area_bound_min s2area_bound_max]);%for counting Kr events
Kr_events=length(drift_time(kr_energy_cut));% Are there enough events for the XY or XYZ map?

    
            %%%% Detector Edge Detection %%%%%%%

            %Finding R_cut
            edge_r_max=sqrt(4000*(25^2)/length(drift_time));

            if edge_r_max <2
                edge_r_max=2;
            end

            edge_r_cut=inrange(s2radius,[0,edge_r_max]);
            bin_size=2;

            dT_bins=0:bin_size:1000;
            dT_hist=hist(drift_time(edge_r_cut),dT_bins);

            %Finds first 10 consecutive zeros... replace this with first 10 less than
            %10% of average at start
            zero_count=zeros(length(dT_hist));
            mean_dT=mean(dT_hist(30:80));

            for i=2:length(dT_hist)
                if dT_hist(i-1)< 0.50*mean_dT;
                    zero_count(i)=1+zero_count(i-1);
                    if zero_count(i)==10;
                        det_edge=(i-11)*bin_size;
                    end
                end
            end


            %%%%%%%%%%%%%%%%%%%

    zcut_min = 0.1*det_edge;
    zcut_max = 0.95*det_edge;


    if Kr_events < 30000 %then create a 10x10 grid. Otherwise 25x25. Want about 30 evts per bin
           
            S1xybinsize = 5;%cm
            S1_xbin_min = -25;
            S1_xbin_max = 25;
            S1_ybin_min = -25;
            S1_ybin_max = 25;

            S2xybinsize = 5;%cm
            S2_xbin_min = -25;
            S2_xbin_max = 25;
            S2_ybin_min = -25;
            S2_ybin_max = 25;
                    
    end
    
    s1_xbins = (S1_xbin_max-S1_xbin_min)./S1xybinsize;
    s1_ybins = (S1_ybin_max-S1_ybin_min)./S1xybinsize;

    s2_xbins = (S2_xbin_max-S2_xbin_min)./S2xybinsize;
    s2_ybins = (S2_ybin_max-S2_ybin_min)./S2xybinsize;
        

 if Kr_events < 4000; % skip krpCal if less than 4000 events
     calibrations=zeros(6,1);
 end
     
 if Kr_events > 300000; % Do a 3D map if we have 400k events
     calibrations(5)=1;
 else
     calibrations(5)=0; % else skip it
 end 
    

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Calculate Single Electron Size 
%   Only using events between good Kr S1 and S2
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    n_sig=2; % will refit to -2sigma to +2 sigma around the first mean.
    
    SE_phe_both=d.pulse_area_phe(SE_good_cut);
    SE_phe_bottom=d.phe_bottom(SE_good_cut);
    SE_phe_top=SE_phe_both-SE_phe_bottom;
    
    [a]=mle(SE_phe_both); % retuns mu, sigma and 95% CF bounds
    SE_mu_both=a(1);
    SE_sig_both=a(2);
    % refit to -2 sigma and +2 sigma around the mean
    [a b]=mle(SE_phe_both(inrange(SE_phe_both,[SE_mu_both-n_sig*SE_sig_both SE_mu_both+n_sig*SE_sig_both])) );
    SE_mu_both=a(1);
    SE_sig_both=a(2);
    SE_sig_mu_both=(b(2,1)-b(1,1) )/4 ; % one sigma
    SE_sig_sig_both=(b(2,2)-b(1,2) )/4 ; % one sigma
    
    [a]=mle(SE_phe_bottom); % retuns mu, sigma and 95% CF bounds
    SE_mu_bottom=a(1);
    SE_sig_bottom=a(2);
    % refit to -2 sigma and +2 sigma around the mean
    [a b]=mle(SE_phe_bottom(inrange(SE_phe_bottom,[SE_mu_bottom-n_sig*SE_sig_bottom SE_mu_bottom+n_sig*SE_sig_bottom])) );
    SE_mu_bottom=a(1);
    SE_sig_bottom=a(2);
    SE_sig_mu_bottom=(b(2,1)-b(1,1) )/4 ; % one sigma
    SE_sig_sig_bottom=(b(2,2)-b(1,2) )/4 ; % one sigma
      
    [a]=mle(SE_phe_top); % retuns mu, sigma and 95% CF bounds
    SE_mu_top=a(1);
    SE_sig_top=a(2);
    % refit to -2 sigma and +2 sigma around the mean
    [a b]=mle(SE_phe_top(inrange(SE_phe_top,[SE_mu_top-n_sig*SE_sig_top SE_mu_top+n_sig*SE_sig_top])) );
    SE_mu_top=a(1);
    SE_sig_top=a(2);
    SE_sig_mu_top=(b(2,1)-b(1,1) )/4 ; % one sigma
    SE_sig_sig_top=(b(2,2)-b(1,2) )/4 ; % one sigma
    
    
 % Calculate fits and make the plots for SE
    
  %Both PMT
    both_SE_fig=figure;
    hold on;
    n_sig_fit=3;
    
    fit_start=SE_mu_both-n_sig_fit*SE_sig_both;
    fit_end=SE_mu_both+n_sig_fit*SE_sig_both;
    xfit=fit_start:0.1:fit_end;
    
    [xxo, yyo, xo, no] = step(SE_phe_both(inrange(SE_phe_both,[fit_start fit_end])));
    s = fitoptions('Method','NonlinearLeastSquares','Robust','On','Algorithm','Trust-Region','Lower',[-5,0,0,0],'Upper',[5,Inf,Inf,Inf],...
        'Startpoint',[1,max(no),10,5]);
    skewGauss1 = fittype('a1/c1*exp(-((x-b1)/c1)^2)*(1+erf(a*(x-b1)/c1))','options',s);

    [f_bot,g_bot] = fit(xo',no',skewGauss1);
    delta_bot = f_bot.a/sqrt(1+f_bot.a^2);
    SE_mu_both = f_bot.b1 + f_bot.c1/sqrt(2)*delta_bot*sqrt(2/pi);
    SE_sig_both = sqrt((f_bot.c1/sqrt(2))^2 * (1 - 2 * (delta_bot)^2/pi));
    skew_both = (4-pi)/2 * (delta_bot * sqrt(2/pi))^3 / (1 - 2 * delta_bot^2/pi)^(3/2);
    kurt_both = 2*(pi-3) * (delta_bot * sqrt(2/pi))^3 / (1 - 2 * delta_bot^2/pi)^2;
    ci_both = confint(f_bot);
    y_fit_both = f_bot(xfit);
    
    
    step(SE_phe_both(inrange(SE_phe_both,[fit_start fit_end])),100,'k');
   
    h_fit=plot(xfit,y_fit_both,'-r','linewidth',2);
    box on;
    myfigview(18);
    xlabel('Pulse Area Phe (Phe) '); ylabel('Count/Phe ');
    legend([h_fit],strcat('\mu = ',num2str(SE_mu_both,'%2.2f'),'\newline \sigma= ',num2str(SE_sig_both,'%2.2f'),'\newline skew=', num2str(skew_both,'%2.2f') ), 'location','northeast');
%     
%     saveas(both_SE_fig,strcat('LUX_corrections/',file_id_cp,'/both_SE_',file_id_cp),'jpg');
%     saveas(both_SE_fig,strcat('LUX_corrections/',file_id_cp,'/both_SE_',file_id_cp),'fig'); 
    
 % For Bottom PMT use skew Gauss.
    bottom_SE_fig=figure;
    hold on;
    n_sig_fit=3;
    
    fit_start=SE_mu_bottom-n_sig_fit*SE_sig_bottom;
    fit_end=SE_mu_bottom+n_sig_fit*SE_sig_bottom;
    xfit=fit_start:0.1:fit_end;
    
    [xxo, yyo, xo, no] = step(SE_phe_bottom(inrange(SE_phe_bottom,[fit_start fit_end])));
    s = fitoptions('Method','NonlinearLeastSquares','Robust','On','Algorithm','Trust-Region','Lower',[-5,0,0,0],'Upper',[5,Inf,Inf,Inf],...
        'Startpoint',[1,max(no),10,5]);
    skewGauss1 = fittype('a1/c1*exp(-((x-b1)/c1)^2)*(1+erf(a*(x-b1)/c1))','options',s);

    [f_bot,g_bot] = fit(xo',no',skewGauss1);
    delta_bot = f_bot.a/sqrt(1+f_bot.a^2);
    SE_mu_bottom = f_bot.b1 + f_bot.c1/sqrt(2)*delta_bot*sqrt(2/pi);
    SE_sig_bottom = sqrt((f_bot.c1/sqrt(2))^2 * (1 - 2 * (delta_bot)^2/pi));
    skew_bottom = (4-pi)/2 * (delta_bot * sqrt(2/pi))^3 / (1 - 2 * delta_bot^2/pi)^(3/2);
    kurt_bottom = 2*(pi-3) * (delta_bot * sqrt(2/pi))^3 / (1 - 2 * delta_bot^2/pi)^2;
    ci_bottom = confint(f_bot);
    y_fit_bottom = f_bot(xfit);
    
    
    step(SE_phe_bottom(inrange(SE_phe_bottom,[fit_start fit_end])),100,'k');
   
    h_fit=plot(xfit,y_fit_bottom,'-r','linewidth',2);
    box on;
    myfigview(18);
    xlabel('Pulse Area Phe (Phe) '); ylabel('Count/Phe ');
    legend([h_fit],strcat('\mu = ',num2str(SE_mu_bottom,'%2.2f'),'\newline \sigma= ',num2str(SE_sig_bottom,'%2.2f'),'\newline skew=', num2str(skew_bottom,'%2.2f') ), 'location','northeast');
    
%     saveas(bottom_SE_fig,strcat('LUX_corrections/',file_id_cp,'/bottom_SE_',file_id_cp),'jpg');
%     saveas(bottom_SE_fig,strcat('LUX_corrections/',file_id_cp,'/bottom_SE_',file_id_cp),'fig');   
 
 % For Top PMT use skew Gauss.
    top_SE_fig=figure;
    hold on;
    n_sig_fit=3;
    
    fit_start=SE_mu_top-n_sig_fit*SE_sig_top;
    fit_end=SE_mu_top+n_sig_fit*SE_sig_top;
    xfit=fit_start:0.1:fit_end;
    
    [xxo, yyo, xo, no] = step(SE_phe_top(inrange(SE_phe_top,[fit_start fit_end])));
    s = fitoptions('Method','NonlinearLeastSquares','Robust','On','Algorithm','Trust-Region','Lower',[-5,0,0,0],'Upper',[5,Inf,Inf,Inf],...
        'Startpoint',[1,max(no),10,5]);
    skewGauss1 = fittype('a1/c1*exp(-((x-b1)/c1)^2)*(1+erf(a*(x-b1)/c1))','options',s);

    [f_bot,g_bot] = fit(xo',no',skewGauss1);
    delta_bot = f_bot.a/sqrt(1+f_bot.a^2);
    SE_mu_top = f_bot.b1 + f_bot.c1/sqrt(2)*delta_bot*sqrt(2/pi);
    SE_sig_top = sqrt((f_bot.c1/sqrt(2))^2 * (1 - 2 * (delta_bot)^2/pi));
    skew_top = (4-pi)/2 * (delta_bot * sqrt(2/pi))^3 / (1 - 2 * delta_bot^2/pi)^(3/2);
    kurt_top = 2*(pi-3) * (delta_bot * sqrt(2/pi))^3 / (1 - 2 * delta_bot^2/pi)^2;
    ci_top = confint(f_bot);
    y_fit_top = f_bot(xfit);
    
    
    step(SE_phe_top(inrange(SE_phe_top,[fit_start fit_end])),100,'k');
   
    h_fit=plot(xfit,y_fit_top,'-r','linewidth',2);
    box on;
    myfigview(18);
    xlabel('Pulse Area Phe (Phe) '); ylabel('Count/Phe ');
    legend([h_fit],strcat('\mu = ',num2str(SE_mu_top,'%2.2f'),'\newline \sigma= ',num2str(SE_sig_top,'%2.2f'),'\newline skew=', num2str(skew_top,'%2.2f') ), 'location','northeast');
%     
%     saveas(top_SE_fig,strcat('LUX_corrections/',file_id_cp,'/top_SE_',file_id_cp),'jpg');
%     saveas(top_SE_fig,strcat('LUX_corrections/',file_id_cp,'/top_SE_',file_id_cp),'fig');

    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% Correcting S1 XYZ
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%////////////////Set Up Variables

s2_cut=inrange(s2_phe_both,[200 30000]);%cut on appropriate S2 size

s1=s1_phe_both(s2_cut); %No Corrections. Use Bottom, Top or Both PMT.
s1_bottom=s1_phe_bottom(s2_cut);
xc=s2x(s2_cut);
yc=s2y(s2_cut);
dT=drift_time(s2_cut);

Kr_events=length(s1);


%////////////////

%s1=s1_xyz;

bin=5;
x=100:bin:500; %Set up S1 energy bins
x=x';

cut_fit=x>110 & x<490; %remove bin edges

xy_step=4;% 2 cm bins
r_max=25; % max radius

z_step=ceil( (0.95*det_edge-0.05*det_edge)/15 );  %xy_step size about = 30 mm, to create 16 z bins
s1zbins=floor(0.05*det_edge):z_step:floor(0.05*det_edge)+z_step*15; %20 us


clear 'q';
clear 'hist_q';
clear 'hist_q_bottom';
clear 'yqfit';
clear 'yqfit_bottom';
clear 'Mean_S1_3D_both';
clear 'Mean_S1_3D_bottom';
clear 'Sig';
clear 'Sig_b';
clear 'Sigma_S1_3D_both';
clear 'Sigma_S1_3D_bottom';


Count_S1_3D=zeros(12,12,16);
Count_S1_3D_bottom=zeros(12,12,16);
Mean_S1_3D_both=zeros(12,12,16);
Mean_S1_3D_bottom=zeros(12,12,16);
Sigma_S1_3D_both=zeros(12,12,16);
Sigma_S1_3D_bottom=zeros(12,12,16);

for k = s1zbins; % out to 485 mm % start at about 20 mm bin center. xy_steps of 30 mm
    
    tic;
    
    for i = (-r_max+xy_step):xy_step:r_max %map to rows going down (y_cm)
        for j = (-r_max+xy_step):xy_step:r_max %map columns across (x_cm)
                       
            l=int8(i/xy_step+(r_max)/xy_step); %real y from 1:25
            m=int8(j/xy_step+(r_max)/xy_step); % real x from 1:25
            n=(k-floor(0.05*det_edge))/z_step + 1; %z from 1:16
            
           %sort, and make the cut. using the variable q
           q = xc<j & xc>(j-xy_step) & yc<i & yc>(i-xy_step) & inrange(dT,[(k-z_step/2) , (k+z_step/2)]); %no 0th element!
                    %swap x,y due to matrix x-row, y-column definition 
                    
           % Hist 
            hist_q = hist(s1(q) ,x)'/bin;
            hist_q_bottom = hist(s1_bottom(q) ,x)'/bin;
            
            %memory clean up. Improve speed by cutting down the size of S1      
            s1=s1(~q); %remove used portion of s1
            s1_bottom=s1_bottom(~q);
            xc=xc(~q);
            yc=yc(~q);
            dT=dT(~q);
            clear 'q'; %free up memory
            %///////////////////////////////////////////////////////
                        
            %Count the number of events per bin
            Count_S1_3D(l,m,n)=sum(hist_q)*bin;
            Count_S1_3D_bottom(l,m,n)=sum(hist_q_bottom)*bin;
            
            if (Count_S1_3D(l,m,n) >= 30) % at least 30 counts before fitting. 
                
                yqfit=hist_q;
                    amp_start=max(yqfit(cut_fit));
                    mean_start=sum(yqfit(cut_fit).*x(cut_fit))/sum(yqfit(cut_fit));
                    sigma_start=std(yqfit(cut_fit).*x(cut_fit));    
                Fit_q= fit(x(cut_fit),yqfit(cut_fit),'gauss1','start',[amp_start mean_start sigma_start]);
                
                %Peak location (mean)
                Mean_S1_3D_both(l,m,n)=Fit_q.b1;
                
                %1-sigma of peak position
                Sig_b=Fit_q.c1/sqrt(2)/sqrt(Count_S1_3D(l,m,n));
   

                %error checking
                   if(strcmp(num2str(Sig_b),'NaN'))
                     Mean_S1_3D_both(l,m,n)=0;
                     Sig_b=0;
                   end
                  
                %uncertainty in mean
                Sigma_S1_3D_both(l,m,n)=Sig_b;
                
             %end IF
            
             else %not enough stats to do the fit
                
                 Fit_q= 0;
                %Find Peak location
                Mean_S1_3D_both(l,m,n)=0;
                %1-sigma
                 Sigma_S1_3D_both(l,m,n)=0;
            end
             
  
                                  
              
        end
    end
toc;
k

end

s1xbins=-r_max+xy_step/2:xy_step:r_max-xy_step/2;
s1ybins=-r_max+xy_step/2:xy_step:r_max-xy_step/2;


%% S1 3D correction. Both PMT (normalized to dT=160)

center_phe_both_xyz=interp3(s1xbins,s1ybins,s1zbins,Mean_S1_3D_both,0,0,(det_edge-4)/2,'spline');%Normalize to the center (x=y=0. dT=160us)
norm_s1_both_xyz=center_phe_both_xyz./Mean_S1_3D_both; %Normalize to center (x=y=0. dT=160)
norm_s1_both_xyz(isinf(norm_s1_both_xyz))=1;%remove infinity. no correction outside 25cm
norm_s1_both_xyz(isnan(norm_s1_both_xyz))=1;%remove nan. no correction outside 25cm

%Get 1 sigma of the corrections matrix.
sigma_center_phe_both_xyz=interp3(s1xbins,s1ybins,s1zbins,Sigma_S1_3D_both,0,0,(det_edge-4)/2,'spline');%Normalize to the center (x=y=0. dT=160us)
sigma_norm_s1_both_xyz=sqrt((sigma_center_phe_both_xyz./Mean_S1_3D_both).^2+(Sigma_S1_3D_both.*center_phe_both_xyz./Mean_S1_3D_both.^2).^2);
sigma_norm_s1_both_xyz(isinf(sigma_norm_s1_both_xyz))=1;%no infinity. Set Sigma=1, which is 100% uncertainty
sigma_norm_s1_both_xyz(isnan(sigma_norm_s1_both_xyz))=1;%no nan. Set Sigma=1, which is 100% uncertainty


mean_center_3D_both=interp3(s1xbins,s1ybins,s1zbins,Mean_S1_3D_both,0,0,(det_edge-4)/2,'spline');
sigma_mean_center_3D_both=interp3(s1xbins,s1ybins,s1zbins,Sigma_S1_3D_both,0,0,(det_edge-4)/2,'spline');
%to Apply the correction use the following.
s1_phe_both_xyz=s1_phe_both.*interp3(s1xbins,s1ybins,s1zbins,norm_s1_both_xyz,s2x,s2y,drift_time,'spline');
    

    %%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%   Calculating the S2 Z-dependence (i.e. electron lifetime)
%   Using Bottom PMT array
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    


    clear y_s2_bottom y_s2_bottom Fit_s2_bottom Fit_EL xEL yfitEL xfit
    clear sigma_s2_bottom sigma_s2_bottom mean_index means mean_index means_error x temp_hist
    clear hist_s2_bottom dT_range dT_fit S x2
    
    %set up energy bins for the gaussian fit. Using bottom S2 only
   bin_s2=50;
   x2=50:bin_s2:10000;
   x2=x2';
   
   %x3=2:1:30; %%for S2/S1
   %x3=x3';
        
    
    cut = inrange(s2radius,[rcut_min,rcut_max]) & inrange(s2_phe_both,[s2area_bound_min,s2area_bound_max]) & inrange(s2_width,[100,800]) & inrange(s1_width,[30,400]) ;
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
        Fit_s2_bottom{bin}=fit(x2(2:end-1),temp_hist(2:end-1),'gauss1','start',[amp_start mean_start sigma_start]);%remove bin edges
          

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
 
    dT_range= inrange( mean_index,[zcut_min, zcut_max] );%from 30 to 300us for the fit
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
    
%     saveas(lifetime_fig_bottom,strcat('LUX_corrections/',file_id_cp,'/S2_bottom_lifetime_',file_id_cp),'jpg');
%     saveas(lifetime_fig_bottom,strcat('LUX_corrections/',file_id_cp,'/S2_bottom_lifetime_',file_id_cp),'fig');
%     
    %Plot 7 of the histograms     
    hist_s2_bottom_fig = figure;
    hold on
    plot_step = ceil((bins-1)/7);    
    j=1;
    t = {num2str(zeros(7,1))};
    for i=2:plot_step:bins-1  %first bin contains anode... skip it. start at i=2.
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
        
        j=j+1;
    end
        title(strcat(file_id_cp,'.Kr83m.S2_Bottom'), 'FontSize', 16,'Interpreter','none');
        legend(strcat(t{1},' \mus'),strcat(t{1},' \mus Fit'),strcat(t{2},' \mus'),strcat(t{2},' \mus Fit')...
            ,strcat(t{3},' \mus'),strcat(t{3},' \mus Fit'),strcat(t{4},' \mus'),strcat(t{4},' \mus Fit')...
            ,strcat(t{5},' \mus'),strcat(t{5},' \mus Fit'),strcat(t{6},' \mus'),strcat(t{6},' \mus Fit')...
            ,strcat(t{7},' \mus'),strcat(t{7},' \mus Fit')); 
        xlabel('Pulse Area Phe','FontSize',16); ylabel('Counts/Phe/sec','FontSize',16);
    myfigview(16);
    xlim([200 8500]);
    
%     saveas(hist_s2_bottom_fig,strcat('LUX_corrections/',file_id_cp,'/S2_bottom_hist_EL_',file_id_cp),'jpg');
%     saveas(hist_s2_bottom_fig,strcat('LUX_corrections/',file_id_cp,'/S2_bottom_hist_EL_',file_id_cp),'fig');
%     
    s2_bottom_means=means;
    s2_bottom_means_sigma=means_error;
    s2_bottom_means_bin=mean_index;      
    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%  Mapping S2 in XY (after correcting for Z-dep). Bottom PMT
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear x2 hist_bin

s2_bin=100;  
x2=1000:s2_bin:12000; % Set up energy bin size. Bottom PMT
x2=x2';

%set up matricies to be filled
mean_S2_bottom_xy = zeros(s2_ybins,s2_xbins);
Sigma_S2_bottom = zeros(s2_ybins,s2_xbins);
Count_S2_bottom = zeros(s2_ybins,s2_xbins);
hist_bin_S2_bottom = zeros(length(x2),s2_ybins,s2_xbins);

s1_phe_both_z=s1_phe_both_xyz;
% s2_phe_bottom_z=s2_phe_bottom.*exp(drift_distance_\mus/e_lifetime); %Correct S2 for lifeimte
s2_phe_bottom_z=s2_phe_bottom.*exp(drift_time/e_lifetime); %Correct S2 for lifeimte


% time_energy_cut = inrange(drift_distance_\mus,[zcut_min,zcut_max]) & inrange(s1_phe_both_z,[200,500]) & inrange(s2_phe_both_z,[1000,20000]) & inrange(s2radius,[rcut_min,rcut_max]) ...
%     & inrange(s2_width,[100,800]) & inrange(s1_width,[30,400]);

time_energy_cut = inrange(drift_time,[zcut_min,zcut_max]) & inrange(s1_phe_both_z,[200,s1area_bound_max]) & inrange(s2_phe_bottom_z,[1000,12000]) & inrange(s2radius,[rcut_min,rcut_max]) ...
    & inrange(s2_width,[100,800]) & inrange(s1_width,[30,400]);

%%%%variables for uber turbo boost!
s2_phe_bottom_z_2=s2_phe_bottom_z(time_energy_cut);
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
      
       hist_bin = hist(s2_phe_bottom_z_2(bin_cut) ,x2)'/s2_bin; 
       hist_bin_S2_bottom(:,y_count,x_count)= hist_bin;
       Count_S2_bottom(y_count,x_count)=sum(hist_bin)*s2_bin; % X and Y flip ... because MATLAB        
        
       %%%%%%%%Truncate variables after binning is complete
        s2_phe_bottom_z_2=s2_phe_bottom_z_2(~bin_cut);
        s2x_2=s2x_2(~bin_cut);
        s2y_2=s2y_2(~bin_cut);
        clear bin_cut;
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
      if Count_S2_bottom(y_count,x_count)>40; % At least 40 events for a gaussian fit. Else no good
          
                    amp_start=max(hist_bin);
                    mean_start=sum(hist_bin.*x2)/sum(hist_bin);
                    sigma_start=std(hist_bin.*x2);    
                Fit_S2_bottom_xy= fit(x2(2:end-1),hist_bin(2:end-1),'gauss1','start',[amp_start mean_start sigma_start]); %remove edges                                                 
                %Fit_S2_bottom{:,y_count,x_count}=Fit_S2_bottom_xy;
                %Peak location (mean)
                mean_S2_bottom_xy(y_count,x_count)=Fit_S2_bottom_xy.b1; % X and Y flip ... because MATLAB             
               %1-sigma of peak position
                Sig_b=Fit_S2_bottom_xy.c1/sqrt(2)/sqrt(Count_S2_bottom(y_count,x_count)); %sigma/sqrt(N)
   
                    %Fit error checking
                        if(isnan(Sig_b))
                         mean_S2_bottom_xy(y_count,x_count)=0;
                         Sig_b=0;
                        end
                        
                    %uncertainty in mean
                    Sigma_S2_bottom(y_count,x_count)=Sig_b;% X and Y flip ... because MATLAB      
     
      else %less than 10 events in a bin... no correction
          
        mean_S2_bottom_xy(y_count,x_count)=0;     
        
      end
    end
   end
   
   mean_S2_bottom_xy(isnan(mean_S2_bottom_xy)) = 0;
     
   
%%Plot the correction S2_XY both %%%%%%%%%%%%%%%

    s2xbins=(S2_xbin_min+S2xybinsize/2):S2xybinsize:(S2_xbin_max-S2xybinsize/2);
    s2ybins=(S2_ybin_min+S2xybinsize/2):S2xybinsize:(S2_ybin_max-S2xybinsize/2);

    
    color_range_max=max(max(mean_S2_bottom_xy));
    color_range_min=min(min(mean_S2_bottom_xy(mean_S2_bottom_xy>0)));
    vc_step=(color_range_max-color_range_min)/50;
    
    vc=color_range_min:vc_step:color_range_max;

    s2_xy_bottom_fig = figure;
    contourf(s2xbins,s2ybins,mean_S2_bottom_xy,vc,'LineColor','none');
    xlabel('x (cm)','fontsize', 18);
    ylabel('y (cm)','fontsize', 18);
    title(strcat(file_id_cp,'. S2 Mean (Bottom PMT) vs. XY.'),'fontsize',16,'Interpreter','none');
    caxis([color_range_min color_range_max])
    g=colorbar;
    ylabel(g,'Phe','FontSize',16)
    set(g,'FontSize',18)
    hold on
    LUXPlotTopPMTs
    axis([-25 25 -25 25])
    myfigview(16);
% 
% saveas(s2_xy_bottom_fig,strcat('LUX_corrections/',file_id_cp,'/S2_XY_bottom_',file_id_cp),'jpg');
% saveas(s2_xy_bottom_fig,strcat('LUX_corrections/',file_id_cp,'/S2_XY_bottom_',file_id_cp),'fig');

%%Plot the 1 sigma of the mean

    color_range_max=max(max(Sigma_S2_bottom));
    color_range_min=min(min(Sigma_S2_bottom(Sigma_S2_bottom>0)));
    vc_step=(color_range_max-color_range_min)/50;

    vc=color_range_min:vc_step:color_range_max;

    sigma_s2_xy_bottom_fig = figure;
    contourf(s2xbins,s2ybins,Sigma_S2_bottom,vc,'LineColor','none');
    xlabel('x (cm)','fontsize', 18);
    ylabel('y (cm)','fontsize', 18);
    title(strcat(file_id_cp, '. 1-Sigma S2 Mean (Bottom PMT) vs. XY.'),'fontsize',16,'Interpreter','none');
    caxis([color_range_min color_range_max])
    g=colorbar;
    ylabel(g,'Phe','FontSize',16)
    set(g,'FontSize',18)
    hold on
    LUXPlotTopPMTs
    axis([-25 25 -25 25])
    myfigview(16);
% 
% saveas(sigma_s2_xy_bottom_fig,strcat('LUX_corrections/',file_id_cp,'/Sigma_S2_XY_bottom_',file_id_cp),'jpg');
% saveas(sigma_s2_xy_bottom_fig,strcat('LUX_corrections/',file_id_cp,'/Sigma_S2_XY_bottom_',file_id_cp),'fig');

%Calculate corrections matrix and 1 sigma corrections matrix
center=interp2(s2xbins,s2ybins,mean_S2_bottom_xy,0,0,'spline',1);%Normalize to the center (x=y=0)
norm_S2_bot=center./mean_S2_bottom_xy; 
norm_S2_bot(isinf(norm_S2_bot))=1;%no infinity! Do not correct events outside r=25cm
norm_S2_bot(isnan(norm_S2_bot))=1;

sigma_center=interp2(s2xbins,s2ybins,Sigma_S2_bottom,0,0,'spline',1);%Normalize to the center (x=y=0)
sigma_norm_S2_bot=sqrt((sigma_center./mean_S2_bottom_xy).^2+(Sigma_S2_bottom.*center./mean_S2_bottom_xy.^2).^2);
sigma_norm_S2_bot(isinf(sigma_norm_S2_bot))=1;%no infinity. Set Sigma=1, which is 100% uncertainty
sigma_norm_S2_bot(isnan(sigma_norm_S2_bot))=1;%no nan. Set Sigma=1, which is 100% uncertainty

s2_phe_bottom_xyz=s2_phe_bottom_z.*interp2(s2xbins,s2ybins,norm_S2_bot,s2x,s2y,'spline');



%% S2 both
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
        
    cut = inrange(s2radius,[rcut_min,rcut_max]) & inrange(s2_phe_both,[s2_bin_min,s2_bin_max]) & inrange(s2_width,[100,800]) & inrange(s1_width,[30,400]) ;
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
 
    dT_range= inrange( mean_index,[zcut_min, zcut_max] );%from 30 to 300us for the fit
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
    
    %saveas(lifetime_fig_both,strcat('LUX_corrections/',file_id_cp,'/S2_both_lifetime_',file_id_cp),'jpg');
    %saveas(lifetime_fig_both,strcat('LUX_corrections/',file_id_cp,'/S2_both_lifetime_',file_id_cp),'fig');
    
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
    
    %saveas(hist_s2_both_fig,strcat('LUX_corrections/',file_id_cp,'/S2_both_hist_EL_',file_id_cp),'jpg');
    %saveas(hist_s2_both_fig,strcat('LUX_corrections/',file_id_cp,'/S2_both_hist_EL_',file_id_cp),'fig');
    
    s2_both_means=means;
    s2_both_means_sigma=means_error;
    s2_both_means_bin=mean_index;      




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
mean_S2_both_xy = zeros(s2_xbins,s2_ybins);
Sigma_S2_both = zeros(s2_xbins,s2_ybins);
Count_S2_both = zeros(s2_xbins,s2_ybins);

s2_phe_both_z=s2_phe_both.*exp(drift_time/e_lifetime_both);

time_energy_cut = inrange(drift_time,[zcut_min,zcut_max]) & inrange(s1_phe_both_z,[100,s1area_bound_max]) & inrange(s2_phe_both_z,[s2_bin_min,s2_bin_max]) & inrange(s2radius,[rcut_min,rcut_max]) ...
    & inrange(s2_width,[100,800]) & inrange(s1_width,[30,400]);

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

%saveas(s2_xy_both_fig,strcat('LUX_corrections/',file_id_cp,'/S2_XY_both_',file_id_cp),'jpg');
%saveas(s2_xy_both_fig,strcat('LUX_corrections/',file_id_cp,'/S2_XY_both_',file_id_cp),'fig');

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

%saveas(sigma_s2_xy_both_fig,strcat('LUX_corrections/',file_id_cp,'/Sigma_S2_XY_both_',file_id_cp),'jpg');
%saveas(sigma_s2_xy_both_fig,strcat('LUX_corrections/',file_id_cp,'/Sigma_S2_XY_both_',file_id_cp),'fig');

%Calculate corrections matrix and 1 sigma corrections matrix
center=interp2(s2xbins,s2ybins,mean_S2_both_xy,0,0,'cubic');%Normalize to the center (x=y=0)
norm_S2_all=center./mean_S2_both_xy; 
norm_S2_all(isinf(norm_S2_all))=1;%no infinity! Do not correct events outside r=25cm
norm_S2_all(isnan(norm_S2_all))=1;

sigma_center=interp2(s2xbins,s2ybins,Sigma_S2_both,0,0,'cubic');%Normalize to the center (x=y=0)
sigma_norm_S2_all=sqrt((sigma_center./mean_S2_both_xy).^2+(Sigma_S2_both.*center./mean_S2_both_xy.^2).^2);
sigma_norm_S2_all(isinf(sigma_norm_S2_all))=1;%no infinity. Set Sigma=1, which is 100% uncertainty
sigma_norm_S2_all(isnan(sigma_norm_S2_all))=1;%no nan. Set Sigma=1, which is 100% uncertainty


s2_phe_both_xyz=s2_phe_both_z.*interp2(s2xbins,s2ybins,norm_S2_all,s2x,s2y,'spline');


save('run04_G1G2')



%% Fixing Cs peak method

clearvars g1 g2
current_best=900;

%  g1=0.097;
g1=0.123;
i=1;
    for g2=5.4:0.01:8
        gs(i)=g2;
        energy=(1/73)*(s1_phe_both_xyz./g1 + s2_phe_bottom_xyz./g2);
        krfit=fit([0:.1:300].',hist(energy,[0:.1:300]).','gauss1');
        kr_peak=krfit.b1;
        gs_diff(i)=abs(kr_peak-41.6);
        gs_mean(i)=kr_peak;

        if abs(kr_peak-41.6) < abs(current_best-41.6)
            current_best=kr_peak;
            best_g2=g2;
            best_energy=energy;
        end
        i=i+1;
    end

 extrac_eff=best_g2/SE_mu_bottom;

 figure
 plot(gs,gs_diff,'.k')
myfigview(16);
xlabel('G2 Value');ylabel('Energy Peak - 41.6 keV');
title('Finding G2 by matching Kr Peak');


%% s2 both g2

clearvars g1 g2 gs gs_diff gs_mean
current_best=900;

%  g1=0.097;
g1=0.123;
i=1;
    for g2=13.5:0.01:19
        gs(i)=g2;
        energy=(1/73)*(s1_phe_both_xyz./g1 + s2_phe_both_xyz./g2);
        krfit=fit([0:.1:300].',hist(energy,[0:.1:300]).','gauss1');
        kr_peak=krfit.b1;
        gs_diff(i)=abs(kr_peak-41.6);
        gs_mean(i)=kr_peak;

        if abs(kr_peak-41.6) < abs(current_best-41.6)
            current_best=kr_peak;
            best_g2=g2;
            best_energy=energy;
        end
        i=i+1;
    end

 extrac_eff=best_g2/SE_mu_bottom;

 figure
 plot(gs,gs_diff,'.k')
myfigview(16);
xlabel('G2 Value');ylabel('Energy Peak - 41.6 keV');
title('Finding G2 by matching Kr Peak');
