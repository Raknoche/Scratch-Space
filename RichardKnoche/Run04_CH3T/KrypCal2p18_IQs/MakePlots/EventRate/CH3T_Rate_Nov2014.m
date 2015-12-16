%% run after load all CH3T

% d.livetime_sec=sum(d.livetime_end_samples-d.livetime_latch_samples)/1e8;
% live_fraction=livetime_sec/((event_timestamp_samples(end)-event_timestamp_samples(1))/1e8);
%        live_fraction=0.97;
%       event_time_min=flipdim(event_timestamp_samples,1)*10^(-8)/60;
%       event_time_hour=event_time_min/60;
%      
%     start_time = 0; %0 
%     offset_time=.5;%890;%hours
%     bin_size=1;%hours
%     end_time=60; %23 hours
%     time_x=bin_size/2+start_time:bin_size:end_time; %hours. HIST bin centers
%     fit_start_t = offset_time+2;
%     fit_end_t= offset_time+13; %was 23
% 
% 
%     count_evt_time=hist(event_time_hour,time_x);
%     hist_evt_time=count_evt_time/bin_size/3600/live_fraction; %number per second
%     sigma_hist_evt_time=sqrt(count_evt_time)/bin_size/3600/live_fraction;

load('C:\Program Files\MATLAB\R2012a\bin\LUXCode\Scratch\RichardKnoche\Run04_CH3T\KrypCal2p18_IQs\MakePlots\CH3T_Nov2014_Kr2p18.mat');


    %Take out ugly data for august 2014 data sets
     total_time_x([1,2,3,7,8,15,16,23,24,25,31,32,33,38,39,40,47,48,49,56,57,58,68,66,67,68,69,76,77,78,85,86,87,88])=[];
     total_hist_evt_time([1,2,3,7,8,15,16,23,24,25,31,32,33,38,39,40,47,48,49,56,57,58,68,66,67,68,69,76,77,78,85,86,87,88])=[];
     total_sigma_hist_evt_time([1,2,3,7,8,15,16,23,24,25,31,32,33,38,39,40,47,48,49,56,57,58,68,66,67,68,69,76,77,78,85,86,87,88])=[];

    
    bin_size=1;
offset_time=0;
    rate_plot=figure;
    hold on
    stairs(total_time_x-bin_size/2-offset_time,total_hist_evt_time,'linewidth',2)
    xlabel('Time (hours) ');ylabel('Rate (Hz)'); title('CH3T R~23.7 cm','interpreter','none');
    myfigview(16)

    fit_start_t=total_time_x(1);
    fit_end_t=total_time_x(30);
     t_cut=inrange(total_time_x,[fit_start_t fit_end_t]);

        Fit_CH3T= fit(total_time_x(t_cut)',total_hist_evt_time(t_cut)','exp1','WEIGHTS',total_sigma_hist_evt_time(t_cut));
        yfitCH3T = Fit_CH3T.a.*exp(Fit_CH3T.b.*total_time_x);

        time_const=(-1/Fit_CH3T.b);

        b=confint(Fit_CH3T, 0.683);%1 sigma matrix
        sigma_time_const=(1/b(1,2)-1/b(2,2))/2;

        chi_sqrd=sum((total_hist_evt_time(t_cut)-yfitCH3T(t_cut)).^2/total_sigma_hist_evt_time(t_cut).^2);

%nominal background for August 2014 (first week, during different trigger settings) is 0.4376 on y axis
     plot(total_time_x(t_cut)-offset_time,yfitCH3T(t_cut),'--r','linewidth',3);
%      plot([-10 end_time],[0.4376 0.4376],'--g')
     errorbar(total_time_x-offset_time,total_hist_evt_time,total_sigma_hist_evt_time,'.k','linewidth',1);  


        legend('CH3T Rate', strcat( '\tau= ', num2str(time_const,2), ' \pm ', num2str(sigma_time_const,1), ' Hours' ),...
            'location','northeast');
        box;

        logy;
        grid on

            