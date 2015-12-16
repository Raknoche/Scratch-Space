%% Load Data

%First load Sep2015 CH3T Calibration
load('C:\Program Files\MATLAB\R2012a\bin\LUXCode\Scratch\RichardKnoche\Run04_CH3T\KrypCal2p18_IQs\MakePlots\ERBand\ERBandFits\Sep2015_2p18_ERFit')

Sep2015_mean_ER_band=mean_ER_band;
Sep2015_ER_sigma=ER_sigma;
Sep2015_ER_bins=ER_bins;
Sep2015_ER_mean_power_fit=ER_mean_power_fit;
Sep2015_ER_upper_power_fit=ER_upper_power_fit;
Sep2015_ER_lower_power_fit=ER_lower_power_fit;
clearvars -except Sep2015_mean_ER_band Sep2015_ER_sigma Sep2015_ER_bins Sep2015_ER_lower_power_fit Sep2015_ER_upper_power_fit Sep2015_ER_mean_power_fit

%Next load Sep2014 CH3T Calibration
load('C:\Program Files\MATLAB\R2012a\bin\LUXCode\Scratch\RichardKnoche\Run04_CH3T\KrypCal2p18_IQs\MakePlots\ERBand\ERBandFits\Sep2014_2p18_ERFit')

Sep2014_mean_ER_band=mean_ER_band;
Sep2014_ER_sigma=ER_sigma;
Sep2014_ER_bins=ER_bins;

clearvars -except Sep2015_mean_ER_band Sep2015_ER_sigma Sep2015_ER_bins Sep2015_ER_lower_power_fit Sep2015_ER_upper_power_fit Sep2015_ER_mean_power_fit ...
    Sep2014_mean_ER_band Sep2014_ER_sigma Sep2014_ER_bins


%Next load Nov2014 CH3T Calibration
load('C:\Program Files\MATLAB\R2012a\bin\LUXCode\Scratch\RichardKnoche\Run04_CH3T\KrypCal2p18_IQs\MakePlots\ERBand\ERBandFits\Nov2014_2p18_ERFit')

Nov2014_mean_ER_band=mean_ER_band;
Nov2014_ER_sigma=ER_sigma;
Nov2014_ER_bins=ER_bins;

clearvars -except Sep2015_mean_ER_band Sep2015_ER_sigma Sep2015_ER_bins Sep2015_ER_lower_power_fit Sep2015_ER_upper_power_fit Sep2015_ER_mean_power_fit ...
    Sep2014_mean_ER_band Sep2014_ER_sigma Sep2014_ER_bins ...
    Nov2014_mean_ER_band Nov2014_ER_sigma Nov2014_ER_bins

%Next load Feb2015 CH3T Calibration
load('C:\Program Files\MATLAB\R2012a\bin\LUXCode\Scratch\RichardKnoche\Run04_CH3T\KrypCal2p18_IQs\MakePlots\ERBand\ERBandFits\Feb2015_2p18_ERFit')

Feb2015_mean_ER_band=mean_ER_band;
Feb2015_ER_sigma=ER_sigma;
Feb2015_ER_bins=ER_bins;

clearvars -except Sep2015_mean_ER_band Sep2015_ER_sigma Sep2015_ER_bins Sep2015_ER_lower_power_fit Sep2015_ER_upper_power_fit Sep2015_ER_mean_power_fit...
    Sep2014_mean_ER_band Sep2014_ER_sigma Sep2014_ER_bins ...
    Nov2014_mean_ER_band Nov2014_ER_sigma Nov2014_ER_bins ...
    Feb2015_mean_ER_band Feb2015_ER_sigma Feb2015_ER_bins

%% Draw from the gaussians to make the full band
% should weight this by livetime eventually, but for now they have equal weight
            s2_s1_bin= 1:.05:5;
band_sigma=1.28; %90 bands

figure
subplot(5,5,1)
myfigview(16)
hold on;
for i=1:length(Sep2015_ER_bins)
ER_montecarlo_Sep2015(i,:)=normrnd(Sep2015_mean_ER_band(i),Sep2015_ER_sigma(i),10000,1);
ER_montecarlo_Sep2014(i,:)=normrnd(Sep2014_mean_ER_band(i),Sep2014_ER_sigma(i),10000,1);
ER_montecarlo_Feb2015(i,:)=normrnd(Feb2015_mean_ER_band(i),Feb2015_ER_sigma(i),10000,1);
ER_montecarlo_Nov2014(i,:)=normrnd(Nov2014_mean_ER_band(i),Nov2014_ER_sigma(i),10000,1);
total_ER_montecarlo(i,:)=[ER_montecarlo_Sep2015(i,:),ER_montecarlo_Sep2014(i,:),ER_montecarlo_Feb2015(i,:),ER_montecarlo_Nov2014(i,:)];

%fit the total_ER_montecarlo
clear tempfit
tempfit=fit(s2_s1_bin.',hist(total_ER_montecarlo(i,:),s2_s1_bin).','gauss1');
mean_ER_band(i)=tempfit.b1;
lower_bound(i)=tempfit.b1-(tempfit.c1/sqrt(2))*band_sigma;
upper_bound(i)=tempfit.b1+(tempfit.c1/sqrt(2))*band_sigma;


%Make plot of gaussianity
subplot(5,5,i,'position',[2 2 2 2])
% pause(0.5)
step(total_ER_montecarlo(i,:)-mean_ER_band(i),[-0.6:0.02:0.6],'-k')
hold on;
x=[-0.6:0.02:0.6];
y=(0.2/0.5).*tempfit.a1.*exp(-((x)./tempfit.c1).^2);
plot(x,y,'-b','LineWidth',2)
logy
ylim([1 10^4.2])

    pos = get(gca, 'Position');
    if rem(i,5)==0
        xadd=0.19*5;
    else
        xadd=0;
    end
    pos(1) = 0.19*(i-(floor(i/5)*5))-0.16+xadd;
    pos(2) = 0.8-(ceil(i/5)-1)*0.18;
    pos(3) = 0.18;
    pos(4) = 0.17;
    set(gca,'yticklabel',{}) 
    set(gca,'xticklabel',{}) 
    set(gca,'FontSize',12)
    set(gca, 'Position', pos)

    if i>=21 
        format_ticks(gca,{'-0.5','0','0.5'},{},[-0.5:0.5:0.6])
    end
    
    if rem(i-1,5)==0
        format_ticks(gca,{},{'$10^0$','$10^1$','$10^2$','$10^3$','$10^4$'},[],[10^0 10^1 10^2 10^3 10^4])
    end
    
end

ER_bins=Sep2015_ER_bins;
   

%% plot total band  
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
max_r=21; %parameters from individual ER band code
s1_max=50;
s1_min=0;

 no_fid_new=figure;
 hold on
 ylabel('Corrected log10(S2/S1)'); xlabel('S1 corrected (Pulse Area Phe)'); 
 title(sprintf('Run 04 Combined CH3T, R<%d cm',max_r),'fontsize',16,'Interpreter','none');
 myfigview(16);
 ylim([1.8 3.5]);xlim([0 s1_max]);
 plot(ER_bins,ER_lower_fit,'--k','markersize',14,'linewidth',2)
 plot(ER_bins,ER_mean_fit,'-k','markersize',20,'linewidth',2)
 plot(ER_bins,ER_upper_fit,'--k','markersize',14,'linewidth',2)

 
  %Make Plot with Sep2015 overlay for comparison
max_r=21; %parameters from individual ER band code
s1_max=50;
s1_min=0;

        
         % Fit power law to the mean
         Sep2015_ER_mean_fit=Sep2015_ER_mean_power_fit.a.*ER_bins.^(Sep2015_ER_mean_power_fit.b);
         Sep2015_ER_lower_fit=Sep2015_ER_lower_power_fit.a.*ER_bins.^(Sep2015_ER_lower_power_fit.b);           
         Sep2015_ER_upper_fit=Sep2015_ER_upper_power_fit.a.*ER_bins.^(Sep2015_ER_upper_power_fit.b);



 no_fid_new=figure;
 hold on
 ylabel('Corrected log10(S2/S1)'); xlabel('S1 corrected (Pulse Area Phe)'); 
 title(sprintf('Run 04 Combined CH3T, R<%d cm',max_r),'fontsize',16,'Interpreter','none');
 myfigview(16);
 ylim([1.8 3.5]);xlim([0 s1_max]);
 plot(ER_bins,ER_lower_fit,'--k','markersize',14,'linewidth',2)
 plot(ER_bins,ER_mean_fit,'-k','markersize',20,'linewidth',2)
 plot(ER_bins,ER_upper_fit,'--k','markersize',14,'linewidth',2)
 plot(Sep2015_ER_bins,Sep2015_ER_lower_fit,'--r','markersize',14,'linewidth',2)
 plot(Sep2015_ER_bins,Sep2015_ER_mean_fit,'-r','markersize',20,'linewidth',2)
 plot(Sep2015_ER_bins,Sep2015_ER_upper_fit,'--r','markersize',14,'linewidth',2)
