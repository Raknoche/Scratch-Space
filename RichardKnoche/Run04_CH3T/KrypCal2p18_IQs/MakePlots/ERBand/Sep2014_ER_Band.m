load('C:\Program Files\MATLAB\R2012a\bin\LUXCode\Scratch\RichardKnoche\Run04_CH3T\KrypCal2p18_IQs\MakePlots\CH3T_Sep2014_Kr2p18.mat');


s1_min=0;
s1_max=50;

%% First make plot of CH3T selection

figure
density(log10(s2_phe_both_xyz),log10(s1_phe_both_xyz),[2.5:0.01:5],[0:0.01:3])
cc=colorbar;
ylabel(cc,'Count','FontSize',16);
caxis
hold on;
xlabel('log10(Corrected S2)'); ylabel ('log10(Corrected S1)'); myfigview(16)
x=[3.629:0.001:4.185];
y=-0.7.*x+4.8;
plot(x,y,'-r','LineWidth',2);
x=[2.703:0.001:3.26];
y=-0.7.*x+2.3;
plot(x,y,'-r','LineWidth',2);
x=[2.703:0.001:3.629];
y=2.*x-5;
plot(x,y,'-r','LineWidth',2);
x=[3.26:0.001:4.185];
y=2.*x-6.5;
plot(x,y,'-r','LineWidth',2);
title('Sep2014 CH3T');

%% First figure out where the outliers are in space

pulse_area_cut= (log10(s1_phe_both_xyz)+0.7*log10(s2_phe_both_xyz)<4.8) & (log10(s1_phe_both_xyz)+0.7*log10(s2_phe_both_xyz)>2.3) & ...
    (log10(s1_phe_both_xyz)-2*log10(s2_phe_both_xyz)<-5) & (log10(s1_phe_both_xyz)-2*log10(s2_phe_both_xyz)>-6.5);


max_r=21;
fiducial_cut=inrange(drift_time,[40 300]) & inrange(s2radius,[0,10]); %delenser not working for this
   
tritium_cut= clean_cut & pulse_area_cut & fiducial_cut;

%% Next, make plot of fiducial cut

figure
density(s2radius_del.^2,drift_time,[0:10:1000],[0:2:450])
cc=colorbar;
ylabel(cc,'Count','FontSize',16);
set(gca,'Ydir','reverse')
hold on;
line([441 441],[40 300],'Color','r','LineWidth',2)
line([0 441],[300 300],'Color','r','LineWidth',2)
line([0 441],[40 40],'Color','r','LineWidth',2)
xlabel('Delensed Radius^2 (cm^2)'); ylabel('Drift Time (uSec)'); myfigview(16);
% plot(s2radius_del(~fiducial_cut).^2,drift_time(~fiducial_cut),'.k','MarkerSize',2)
title('Sep2014 CH3T');




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
  title('Sep2014 CH3T');

  s1_phe_both_xyz_all=s1_phe_both_xyz(tritium_cut & inrange(drift_time,[0,110]));
  log10_s2s1_xyz_all=log10(s2_phe_both_xyz(tritium_cut & inrange(drift_time,[0,110]))./s1_phe_both_xyz(tritium_cut & inrange(drift_time,[0,110])));
  ER_lower_fit_all=ER_lower_fit;
  ER_upper_fit_all=ER_upper_fit;
  ER_mean_fit_all=ER_mean_fit;
    
  

            
  %Save ER Band fits
  save('Sep2014_2p18_ERFit','ER_bins','mean_ER_band','ER_sigma','ER_lower_power_fit','ER_upper_power_fit','ER_mean_power_fit','ER_lower_fit','ER_upper_fit','ER_mean_fit');
   