load('CH3T_Sep2014_TAGCB_1_7_1_085_2_poscorrected.mat');


s1_min=0;
s1_max=50;
    
    liquid_cut= inrange(drift_time,[30 300]) & inrange(s2radius_c,[0,18]) & inrange(s1_phe_both_xyz,[s1_min s1_max]) ...
    & (log10(s1_phe_both_xyz)+0.5*log10(s2_phe_bottom_xyz)<3.8) & s2_phe_both>100 ; %used to be drift time 100 to 250

   
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
   
        
            ER_fit_cut=liquid_cut & inrange(s1_phe_both_xyz,[ER_bin-ER_bin_size/2,ER_bin+ER_bin_size/2]);
             
            ER_band_hist=hist(log10(s2_phe_bottom_xyz(ER_fit_cut)./s1_phe_both_xyz(ER_fit_cut)),s2_s1_bin);
            counts_ER_hist(index_i)=sum(ER_band_hist);
           
            Fit_ER=fit(s2_s1_bin',ER_band_hist','gauss1');
             
            rms_fit=fit_ML_normal(log10(s2_phe_bottom_xyz(ER_fit_cut)./s1_phe_both_xyz(ER_fit_cut)),s2_s1_bin);
            
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
 %plot(s1_phe_both_xyz(liquid_cut),log10(s2_phe_bottom_xyz(liquid_cut)./s1_phe_both_xyz(liquid_cut)),'ok','markersize',2)
 %density(s1_phe_both_xyz(liquid_cut),log10(s2_phe_bottom_xyz(liquid_cut)./s1_phe_both_xyz(liquid_cut)),[0:.5:100],[0:.05:3]);
 densityplot(s1_phe_both_xyz(liquid_cut),log10(s2_phe_bottom_xyz(liquid_cut)./s1_phe_both_xyz(liquid_cut)),'rectangle',[.5 .05],'lin');
 cc=colorbar;
    ylabel(cc,'Count','FontSize',16);
 colormap;
 ylabel('Corrected log10(S2_{b}/S1)'); xlabel('S1 corrected (Pulse Area Phe)'); 
 title('CH3T, R<24.5 cm','fontsize',16,'Interpreter','none');
 myfigview(16);
 ylim([0 3]);xlim([0 s1_max]);
%  plot(ER_bins,lower_bound,'ok','markersize',6,'linewidth',2)
 plot(ER_bins,ER_lower_fit,'--k','markersize',14,'linewidth',2)
%  plot(ER_bins,mean_ER_band,'ok','markersize',6,'linewidth',2)%plot mean of ER band
 plot(ER_bins,ER_mean_fit,'-k','markersize',20,'linewidth',2)
%  plot(ER_bins,upper_bound,'ok','markersize',6,'linewidth',2)
 plot(ER_bins,ER_upper_fit,'--k','markersize',14,'linewidth',2)

 
     bin_step=2;

            total_count=length(log10(s2_phe_bottom_xyz(liquid_cut)));
 
    
 box;
  text(10,2.8,strcat('Total Count= ',num2str(total_count)),'fontsize',16);
  %text(5,.7,strcat('NEST\_SE_{bot} Fit=',num2str(ER_mean_NEST_fit.SE_ER_fit*73/64*27/24.5,3),'\pm ', num2str(sigma_SE_ER_fit*73/64*27/24.5,2),'[Phe]'),'fontsize',16,'color','m','interpreter','tex');
  %text(5,.8,'NEST Platinum V4-b', 'fontsize',16,'color','m','interpreter','tex' );
  text(5,.55,strcat('+',num2str(band_sigma),'\sigma = (',num2str(ER_upper_power_fit.a,3),'\pm', num2str(sigma_a_ER_upper,1) ,')\times S1_c ^{', num2str(ER_upper_power_fit.b,3),'\pm',num2str(sigma_b_ER_upper,1),'}'),'fontsize',16,'interpreter','tex');
  text(5,.3,strcat('Mean=(',num2str(ER_mean_power_fit.a,3),'\pm', num2str(sigma_a_ER_mean,1),')\times S1_c ^{', num2str(ER_mean_power_fit.b,3),'\pm', num2str(sigma_b_ER_mean,1),'}'),'fontsize',16,'interpreter','tex');
  text(5,.05,strcat('-',num2str(band_sigma),'\sigma = (',num2str(ER_lower_power_fit.a,3),'\pm', num2str(sigma_a_ER_lower,1),')\times S1_c ^{', num2str(ER_lower_power_fit.b,3),'\pm', num2str(sigma_b_ER_lower,1),'}'),'fontsize',16,'interpreter','tex');
  
  s1_phe_both_xyz_all=s1_phe_both_xyz(liquid_cut & inrange(drift_time,[0,110]));
  log10_s2s1_xyz_all=log10(s2_phe_bottom_xyz(liquid_cut & inrange(drift_time,[0,110]))./s1_phe_both_xyz(liquid_cut & inrange(drift_time,[0,110])));
  ER_lower_fit_all=ER_lower_fit;
  ER_upper_fit_all=ER_upper_fit;
  ER_mean_fit_all=ER_mean_fit;
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
 % Z SLICE ER BAND
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%
 
            ER_bin_size=2; %3 works
           index_i=1;
           ER_bins=ER_bin_size/2+s1_min:ER_bin_size:s1_max;
               lower_bound=ones(size(ER_bins));
               mean_ER_band=ones(size(ER_bins));
               upper_bound=ones(size(ER_bins));
                counts_ER_hist=ones(size(ER_bins));
               ER_sigma=ones(size(ER_bins));
               ER_sigma_sigma=ones(size(ER_bins));
        
               
        for ER_bin = (ER_bin_size/2)+s1_min:ER_bin_size:(s1_max); %HIST uses bin centers
   
        
            ER_fit_cut=liquid_cut & inrange(s1_phe_both_xyz,[ER_bin-ER_bin_size/2,ER_bin+ER_bin_size/2]) & inrange(drift_time,[0,110]);
             
            ER_band_hist=hist(log10(s2_phe_bottom_xyz(ER_fit_cut)./s1_phe_both_xyz(ER_fit_cut)),s2_s1_bin);
            counts_ER_hist(index_i)=sum(ER_band_hist);
           
            Fit_ER=fit(s2_s1_bin',ER_band_hist','gauss1');
             
            rms_fit=fit_ML_normal(log10(s2_phe_bottom_xyz(ER_fit_cut)./s1_phe_both_xyz(ER_fit_cut)),s2_s1_bin);
            
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
  
  s1_phe_both_xyz_110=s1_phe_both_xyz(liquid_cut & inrange(drift_time,[0,110]));
  s2_phe_bottom_xyz_110=s2_phe_bottom_xyz(liquid_cut & inrange(drift_time,[0,110]));
  log10_s2s1_xyz_110=log10(s2_phe_bottom_xyz(liquid_cut & inrange(drift_time,[0,110]))./s1_phe_both_xyz(liquid_cut & inrange(drift_time,[0,110])));
  ER_lower_fit_110=ER_lower_fit;
  ER_upper_fit_110=ER_upper_fit;
  ER_mean_fit_110=ER_mean_fit;

 no_fid_new=figure;
 hold on
 %plot(s1_phe_both_xyz(liquid_cut),log10(s2_phe_bottom_xyz(liquid_cut)./s1_phe_both_xyz(liquid_cut)),'ok','markersize',2)
 %density(s1_phe_both_xyz(liquid_cut),log10(s2_phe_bottom_xyz(liquid_cut)./s1_phe_both_xyz(liquid_cut)),[0:.5:100],[0:.05:3]);
 densityplot(s1_phe_both_xyz(liquid_cut & inrange(drift_time,[0,110])),log10(s2_phe_bottom_xyz(liquid_cut & inrange(drift_time,[0,110]))./s1_phe_both_xyz(liquid_cut & inrange(drift_time,[0,110]))),'rectangle',[.5 .05],'lin');
 cc=colorbar;
    ylabel(cc,'Count','FontSize',16);
 colormap;
 ylabel('Corrected log10(S2_{b}/S1)'); xlabel('S1 corrected (Pulse Area Phe)'); 
 title('0-110 uSec CH3T, R<24.5 cm','fontsize',16,'Interpreter','none');
 myfigview(16);
 ylim([0 3]);xlim([0 s1_max]);
%  plot(ER_bins,lower_bound,'ok','markersize',6,'linewidth',2)
 plot(ER_bins,ER_lower_fit,'--k','markersize',14,'linewidth',2)
%  plot(ER_bins,mean_ER_band,'ok','markersize',6,'linewidth',2)%plot mean of ER band
 plot(ER_bins,ER_mean_fit,'-k','markersize',20,'linewidth',2)
%  plot(ER_bins,upper_bound,'ok','markersize',6,'linewidth',2)
 plot(ER_bins,ER_upper_fit,'--k','markersize',14,'linewidth',2)

 
     bin_step=2;
            total_count=length(log10(s2_phe_bottom_xyz(liquid_cut & inrange(drift_time,[0,110]))));
    
 box;
  text(10,2.8,strcat('Total Count= ',num2str(total_count)),'fontsize',16);
  %text(5,.7,strcat('NEST\_SE_{bot} Fit=',num2str(ER_mean_NEST_fit.SE_ER_fit*73/64*27/24.5,3),'\pm ', num2str(sigma_SE_ER_fit*73/64*27/24.5,2),'[Phe]'),'fontsize',16,'color','m','interpreter','tex');
  %text(5,.8,'NEST Platinum V4-b', 'fontsize',16,'color','m','interpreter','tex' );
  text(5,.55,strcat('+',num2str(band_sigma),'\sigma = (',num2str(ER_upper_power_fit.a,3),'\pm', num2str(sigma_a_ER_upper,1) ,')\times S1_c ^{', num2str(ER_upper_power_fit.b,3),'\pm',num2str(sigma_b_ER_upper,1),'}'),'fontsize',16,'interpreter','tex');
  text(5,.3,strcat('Mean=(',num2str(ER_mean_power_fit.a,3),'\pm', num2str(sigma_a_ER_mean,1),')\times S1_c ^{', num2str(ER_mean_power_fit.b,3),'\pm', num2str(sigma_b_ER_mean,1),'}'),'fontsize',16,'interpreter','tex');
  text(5,.05,strcat('-',num2str(band_sigma),'\sigma = (',num2str(ER_lower_power_fit.a,3),'\pm', num2str(sigma_a_ER_lower,1),')\times S1_c ^{', num2str(ER_lower_power_fit.b,3),'\pm', num2str(sigma_b_ER_lower,1),'}'),'fontsize',16,'interpreter','tex');
  
  
  %%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%
 
            ER_bin_size=2; %3 works
           index_i=1;
           ER_bins=ER_bin_size/2+s1_min:ER_bin_size:s1_max;
               lower_bound=ones(size(ER_bins));
               mean_ER_band=ones(size(ER_bins));
               upper_bound=ones(size(ER_bins));
                counts_ER_hist=ones(size(ER_bins));
               ER_sigma=ones(size(ER_bins));
               ER_sigma_sigma=ones(size(ER_bins));
        
               
        for ER_bin = (ER_bin_size/2)+s1_min:ER_bin_size:(s1_max); %HIST uses bin centers
   
        
            ER_fit_cut=liquid_cut & inrange(s1_phe_both_xyz,[ER_bin-ER_bin_size/2,ER_bin+ER_bin_size/2]) & inrange(drift_time,[110,220]);
             
            ER_band_hist=hist(log10(s2_phe_bottom_xyz(ER_fit_cut)./s1_phe_both_xyz(ER_fit_cut)),s2_s1_bin);
            counts_ER_hist(index_i)=sum(ER_band_hist);
           
            Fit_ER=fit(s2_s1_bin',ER_band_hist','gauss1');
             
            rms_fit=fit_ML_normal(log10(s2_phe_bottom_xyz(ER_fit_cut)./s1_phe_both_xyz(ER_fit_cut)),s2_s1_bin);
            
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
  
  s1_phe_both_xyz_220=s1_phe_both_xyz(liquid_cut & inrange(drift_time,[110,220]));
  s2_phe_bottom_xyz_220=s2_phe_bottom_xyz(liquid_cut & inrange(drift_time,[110,220]));
  log10_s2s1_xyz_220=log10(s2_phe_bottom_xyz(liquid_cut & inrange(drift_time,[110,220]))./s1_phe_both_xyz(liquid_cut & inrange(drift_time,[110,220])));
  ER_lower_fit_220=ER_lower_fit;
  ER_upper_fit_220=ER_upper_fit;
  ER_mean_fit_220=ER_mean_fit;

 no_fid_new=figure;
 hold on
 %plot(s1_phe_both_xyz(liquid_cut),log10(s2_phe_bottom_xyz(liquid_cut)./s1_phe_both_xyz(liquid_cut)),'ok','markersize',2)
 %density(s1_phe_both_xyz(liquid_cut),log10(s2_phe_bottom_xyz(liquid_cut)./s1_phe_both_xyz(liquid_cut)),[0:.5:100],[0:.05:3]);
 densityplot(s1_phe_both_xyz(liquid_cut & inrange(drift_time,[110,220])),log10(s2_phe_bottom_xyz(liquid_cut & inrange(drift_time,[110,220]))./s1_phe_both_xyz(liquid_cut & inrange(drift_time,[110,220]))),'rectangle',[.5 .05],'lin');
 cc=colorbar;
    ylabel(cc,'Count','FontSize',16);
 colormap;
 ylabel('Corrected log10(S2_{b}/S1)'); xlabel('S1 corrected (Pulse Area Phe)'); 
 title('110-220 uSec CH3T, R<24.5 cm','fontsize',16,'Interpreter','none');
 myfigview(16);
 ylim([0 3]);xlim([0 s1_max]);
%  plot(ER_bins,lower_bound,'ok','markersize',6,'linewidth',2)
 plot(ER_bins,ER_lower_fit,'--k','markersize',14,'linewidth',2)
%  plot(ER_bins,mean_ER_band,'ok','markersize',6,'linewidth',2)%plot mean of ER band
 plot(ER_bins,ER_mean_fit,'-k','markersize',20,'linewidth',2)
%  plot(ER_bins,upper_bound,'ok','markersize',6,'linewidth',2)
 plot(ER_bins,ER_upper_fit,'--k','markersize',14,'linewidth',2)

 
     bin_step=2;
            total_count=length(log10(s2_phe_bottom_xyz(liquid_cut & inrange(drift_time,[110,220]))));
    
 box;
  text(10,2.8,strcat('Total Count= ',num2str(total_count)),'fontsize',16);
  %text(5,.7,strcat('NEST\_SE_{bot} Fit=',num2str(ER_mean_NEST_fit.SE_ER_fit*73/64*27/24.5,3),'\pm ', num2str(sigma_SE_ER_fit*73/64*27/24.5,2),'[Phe]'),'fontsize',16,'color','m','interpreter','tex');
  %text(5,.8,'NEST Platinum V4-b', 'fontsize',16,'color','m','interpreter','tex' );
  text(5,.55,strcat('+',num2str(band_sigma),'\sigma = (',num2str(ER_upper_power_fit.a,3),'\pm', num2str(sigma_a_ER_upper,1) ,')\times S1_c ^{', num2str(ER_upper_power_fit.b,3),'\pm',num2str(sigma_b_ER_upper,1),'}'),'fontsize',16,'interpreter','tex');
  text(5,.3,strcat('Mean=(',num2str(ER_mean_power_fit.a,3),'\pm', num2str(sigma_a_ER_mean,1),')\times S1_c ^{', num2str(ER_mean_power_fit.b,3),'\pm', num2str(sigma_b_ER_mean,1),'}'),'fontsize',16,'interpreter','tex');
  text(5,.05,strcat('-',num2str(band_sigma),'\sigma = (',num2str(ER_lower_power_fit.a,3),'\pm', num2str(sigma_a_ER_lower,1),')\times S1_c ^{', num2str(ER_lower_power_fit.b,3),'\pm', num2str(sigma_b_ER_lower,1),'}'),'fontsize',16,'interpreter','tex');
  
  
  %%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%
 
            ER_bin_size=2; %3 works
           index_i=1;
           ER_bins=ER_bin_size/2+s1_min:ER_bin_size:s1_max;
               lower_bound=ones(size(ER_bins));
               mean_ER_band=ones(size(ER_bins));
               upper_bound=ones(size(ER_bins));
                counts_ER_hist=ones(size(ER_bins));
               ER_sigma=ones(size(ER_bins));
               ER_sigma_sigma=ones(size(ER_bins));
        
               
        for ER_bin = (ER_bin_size/2)+s1_min:ER_bin_size:(s1_max); %HIST uses bin centers
   
        
            ER_fit_cut=liquid_cut & inrange(s1_phe_both_xyz,[ER_bin-ER_bin_size/2,ER_bin+ER_bin_size/2]) & inrange(drift_time,[220,330]);
             
            ER_band_hist=hist(log10(s2_phe_bottom_xyz(ER_fit_cut)./s1_phe_both_xyz(ER_fit_cut)),s2_s1_bin);
            counts_ER_hist(index_i)=sum(ER_band_hist);
           
            Fit_ER=fit(s2_s1_bin',ER_band_hist','gauss1');
             
            rms_fit=fit_ML_normal(log10(s2_phe_bottom_xyz(ER_fit_cut)./s1_phe_both_xyz(ER_fit_cut)),s2_s1_bin);
            
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
  
  s1_phe_both_xyz_330=s1_phe_both_xyz(liquid_cut & inrange(drift_time,[220,330]));
  s2_phe_bottom_xyz_330=s2_phe_bottom_xyz(liquid_cut & inrange(drift_time,[220,330]));
  log10_s2s1_xyz_330=log10(s2_phe_bottom_xyz(liquid_cut & inrange(drift_time,[220,330]))./s1_phe_both_xyz(liquid_cut & inrange(drift_time,[220,330])));
  ER_lower_fit_330=ER_lower_fit;
  ER_upper_fit_330=ER_upper_fit;
  ER_mean_fit_330=ER_mean_fit;

 no_fid_new=figure;
 hold on
 %plot(s1_phe_both_xyz(liquid_cut),log10(s2_phe_bottom_xyz(liquid_cut)./s1_phe_both_xyz(liquid_cut)),'ok','markersize',2)
 %density(s1_phe_both_xyz(liquid_cut),log10(s2_phe_bottom_xyz(liquid_cut)./s1_phe_both_xyz(liquid_cut)),[0:.5:100],[0:.05:3]);
 densityplot(s1_phe_both_xyz(liquid_cut & inrange(drift_time,[220,330])),log10(s2_phe_bottom_xyz(liquid_cut & inrange(drift_time,[220,330]))./s1_phe_both_xyz(liquid_cut & inrange(drift_time,[220,330]))),'rectangle',[.5 .05],'lin');
 cc=colorbar;
    ylabel(cc,'Count','FontSize',16);
 colormap;
 ylabel('Corrected log10(S2_{b}/S1)'); xlabel('S1 corrected (Pulse Area Phe)'); 
 title('220-330 uSec CH3T, R<24.5 cm','fontsize',16,'Interpreter','none');
 myfigview(16);
 ylim([0 3]);xlim([0 s1_max]);
%  plot(ER_bins,lower_bound,'ok','markersize',6,'linewidth',2)
 plot(ER_bins,ER_lower_fit,'--k','markersize',14,'linewidth',2)
%  plot(ER_bins,mean_ER_band,'ok','markersize',6,'linewidth',2)%plot mean of ER band
 plot(ER_bins,ER_mean_fit,'-k','markersize',20,'linewidth',2)
%  plot(ER_bins,upper_bound,'ok','markersize',6,'linewidth',2)
 plot(ER_bins,ER_upper_fit,'--k','markersize',14,'linewidth',2)

 
     bin_step=2;
            total_count=length(log10(s2_phe_bottom_xyz(liquid_cut & inrange(drift_time,[220,330]))));
    
 box;
  text(10,2.8,strcat('Total Count= ',num2str(total_count)),'fontsize',16);
  %text(5,.7,strcat('NEST\_SE_{bot} Fit=',num2str(ER_mean_NEST_fit.SE_ER_fit*73/64*27/24.5,3),'\pm ', num2str(sigma_SE_ER_fit*73/64*27/24.5,2),'[Phe]'),'fontsize',16,'color','m','interpreter','tex');
  %text(5,.8,'NEST Platinum V4-b', 'fontsize',16,'color','m','interpreter','tex' );
  text(5,.55,strcat('+',num2str(band_sigma),'\sigma = (',num2str(ER_upper_power_fit.a,3),'\pm', num2str(sigma_a_ER_upper,1) ,')\times S1_c ^{', num2str(ER_upper_power_fit.b,3),'\pm',num2str(sigma_b_ER_upper,1),'}'),'fontsize',16,'interpreter','tex');
  text(5,.3,strcat('Mean=(',num2str(ER_mean_power_fit.a,3),'\pm', num2str(sigma_a_ER_mean,1),')\times S1_c ^{', num2str(ER_mean_power_fit.b,3),'\pm', num2str(sigma_b_ER_mean,1),'}'),'fontsize',16,'interpreter','tex');
  text(5,.05,strcat('-',num2str(band_sigma),'\sigma = (',num2str(ER_lower_power_fit.a,3),'\pm', num2str(sigma_a_ER_lower,1),')\times S1_c ^{', num2str(ER_lower_power_fit.b,3),'\pm', num2str(sigma_b_ER_lower,1),'}'),'fontsize',16,'interpreter','tex');
  
%%%%%%%%%%%% 
% Overlay Z-Slice ER Bands
%%%%%%%%%%%
  
  zslice_ER_overlay=figure;
  plot(ER_bins,ER_mean_fit_110,'-r','markersize',20,'linewidth',2);
  hold on;
  plot(ER_bins,ER_mean_fit_220,'-b','markersize',20,'linewidth',2);
  plot(ER_bins,ER_mean_fit_330,'-g','markersize',20,'linewidth',2);
  plot(ER_bins,ER_mean_fit_all,'-k','markersize',20,'linewidth',2);
  legend('0-110 uSec','110-220 uSec','220-330 uSec');
  plot(ER_bins,ER_lower_fit_110,'--r','markersize',14,'linewidth',2)
  plot(ER_bins,ER_upper_fit_110,'--r','markersize',14,'linewidth',2)
  plot(ER_bins,ER_lower_fit_220,'--b','markersize',14,'linewidth',2)
  plot(ER_bins,ER_upper_fit_220,'--b','markersize',14,'linewidth',2)  
  plot(ER_bins,ER_lower_fit_330,'--g','markersize',14,'linewidth',2)
  plot(ER_bins,ER_upper_fit_330,'--g','markersize',14,'linewidth',2)
  plot(ER_bins,ER_lower_fit_all,'--k','markersize',14,'linewidth',2)
  plot(ER_bins,ER_upper_fit_all,'--k','markersize',14,'linewidth',2)
  title('Z-Slice CH3T ER Bands, R<24.5 cm','fontsize',16,'Interpreter','none');
   ylabel('Corrected log10(S2_{b}/S1)'); xlabel('S1 corrected (Pulse Area Phe)'); 
   myfigview(16);
   
   %%%%%%%%%%%
   % Overlay Z-slice s1 and s2 spectra
   %%%%%%%%%%
   
   s1_phe_both_xyz_110_2=s1_phe_both_xyz(inrange(drift_time,[0,110]));
   s1_phe_both_xyz_220_2=s1_phe_both_xyz(inrange(drift_time,[110,220]));
   s1_phe_both_xyz_330_2=s1_phe_both_xyz(inrange(drift_time,[220,330]));
   
   zslice_S1_overlay=figure;
   step(s1_phe_both_xyz_110_2,[0:1:200],'r');
   hold on;
   step(s1_phe_both_xyz_220_2,[0:1:200],'b');
   step(s1_phe_both_xyz_330_2,[0:1:200],'g');
   legend('0-110 uSec','110-220 uSec','220-330 uSec')
   title('Z-Slice S1 Spectra');xlabel('Corrected S1 (phe)'); ylabel('Count');
   myfigview(16);
   logy;
   
   
 

   
   s2_phe_bottom_xyz_110_2=s2_phe_bottom_xyz(inrange(drift_time,[0,110]));
   s2_phe_bottom_xyz_220_2=s2_phe_bottom_xyz(inrange(drift_time,[110,220]));
   s2_phe_bottom_xyz_330_2=s2_phe_bottom_xyz(inrange(drift_time,[220,330]));

   
   zslice_S2_overlay=figure;
   step(s2_phe_bottom_xyz_110_2,[0:100:10000],'r');
   hold on;
   step(s2_phe_bottom_xyz_220_2,[0:100:10000],'b');
   step(s2_phe_bottom_xyz_330_2,[0:100:10000],'g');
   legend('0-110 uSec','110-220 uSec','220-330 uSec')
   title(strcat('Z-Slice S2 Spectra'));xlabel('Corrected S2 (phe)'); ylabel('Count');
   myfigview(16);
   logy;
 


%%%%%%%%%%%%%%%%%%%%%
% Normalized Gaussian fit
%%%%%%%%%%%%%%%%%%%%%
  
  
  
    clear y_ER_norm Fit_ER Fit_ER_norm mean_ER_band_norm gauss_fit_bins normalized_s2_s1
    normalized_s2_s1=log10(s2_phe_bottom_xyz(liquid_cut)./s1_phe_both_xyz(liquid_cut))-interp1(ER_bins,ER_mean_fit,s1_phe_both_xyz(liquid_cut)) ;
            
           ER_bin_start=-.51;
           ER_bin_end=0.51;
           ER_bin_size=.01;
           gauss_fit_bins=ER_bin_start-ER_bin_size/2:ER_bin_size:ER_bin_end+ER_bin_size/2;
           mean_ER_band_norm=hist(normalized_s2_s1,gauss_fit_bins)/ER_bin_size;
           Fit_ER_norm=fit(gauss_fit_bins',mean_ER_band_norm','gauss1');
           y_ER_norm=Fit_ER_norm.a1.*exp(-((gauss_fit_bins-Fit_ER_norm.b1)./Fit_ER_norm.c1).^2);
           sigma_error=confint(Fit_ER_norm,0.683);
           sigma_error_upper=sigma_error(2,3)/sqrt(2);
           sigma_error_lower=sigma_error(1,3)/sqrt(2);
           
     ER_band_gauss_fig=figure;

     stairs(gauss_fit_bins-ER_bin_size/2,mean_ER_band_norm,'linewidth',2)
     hold on
     plot(gauss_fit_bins,y_ER_norm,'-r','linewidth',2);
     xlim([-.5 .5])
     set(gca, 'XTickMode', 'manual');
        set(gca, 'XTick', [-0.5:.1:.5]);
     xlabel('Corrected log10(S2_b/S1)-Expected ');ylabel('Count/0.01 Phe');
     title(strcat('Run 04 CH3T, R<24.5 cm'),'fontsize',16,'Interpreter','none');
     legend('CH3T data', strcat('\sigma =',num2str(Fit_ER_norm.c1/sqrt(2),4),'(',num2str(sigma_error_lower,4),',',num2str(sigma_error_upper,4),')')); %was divided by 2
     grid
     myfigview(16)
     
     %%%%%%%%%%%%%%%%%%
  
  
  ER_bins_new=ER_bins;
  ER_mean_fit_new=ER_mean_fit;
  ER_upper_fit_new=ER_upper_fit;
  ER_lower_fit_new=ER_lower_fit;
  s1_phe_both_xyz_new=s1_phe_both_xyz;
  liquid_cut_new=liquid_cut;
  s2_phe_bottom_xyz_new=s2_phe_bottom_xyz;
  total_count_new=total_count;
  ER_upper_power_fit_new=ER_upper_power_fit;
  ER_lower_power_fit_new=ER_lower_power_fit;
  ER_mean_power_fit_new=ER_mean_power_fit;
  sigma_a_ER_upper_new=sigma_a_ER_upper;
  sigma_b_ER_upper_new=sigma_b_ER_upper;
  sigma_a_ER_lower_new=sigma_a_ER_lower;
  sigma_b_ER_lower_new=sigma_b_ER_lower; 
  sigma_a_ER_mean_new=sigma_a_ER_mean;
  sigma_b_ER_mean_new=sigma_b_ER_mean;
  mean_ER_band_new=mean_ER_band;
  upper_bound_new=upper_bound;
  lower_bound_new=lower_bound;
  mean_ER_band_new=mean_ER_band;
  
  
  
  s1_phe_both_xyz_new=s1_phe_both_xyz;
  s2_phe_bottom_xyz_new=s2_phe_bottom_xyz;
  
%   
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   % OVERLAY WITH Old Krypcal %
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
%    
%     band_sigma=1.28; %90% confidence bounds
%      
%              
%             s2_s1_bin= 1:.05:3;
%            ER_bin_size=2;
%            index_i=1;
%            ER_bins=ER_bin_size/2+s1_min:ER_bin_size:s1_max;
%                lower_bound=ones(size(ER_bins));
%                mean_ER_band=ones(size(ER_bins));
%                upper_bound=ones(size(ER_bins));
%                 counts_ER_hist=ones(size(ER_bins));
%                ER_sigma=ones(size(ER_bins));
%                ER_sigma_sigma=ones(size(ER_bins));
%         
%                
%         for ER_bin = (ER_bin_size/2)+s1_min:ER_bin_size:(s1_max); %HIST uses bin centers
%    
%         
%             ER_fit_cut=liquid_cut & inrange(s1_phe_both_xyz_old,[ER_bin-ER_bin_size/2,ER_bin+ER_bin_size/2]);
%              
%             ER_band_hist=hist(log10(s2_phe_bottom_xyz_old(ER_fit_cut)./s1_phe_both_xyz_old(ER_fit_cut)),s2_s1_bin);
%             counts_ER_hist(index_i)=sum(ER_band_hist);
%            
%             Fit_ER=fit(s2_s1_bin',ER_band_hist','gauss1');
%              
%             rms_fit=fit_ML_normal(log10(s2_phe_bottom_xyz_old(ER_fit_cut)./s1_phe_both_xyz_old(ER_fit_cut)),s2_s1_bin);
%             
%             %band=confint(Fit_ER,.1);
%             %ER_sigma(index_i)= Fit_ER.c1/sqrt(2);
%             ER_sigma(index_i)= sqrt(rms_fit.sig2);
%             ER_sigma_sigma(index_i)=ER_sigma(index_i)/sqrt(2*counts_ER_hist(index_i));
%             
%             lower_bound(index_i)=Fit_ER.b1-ER_sigma(index_i)*band_sigma;
%             mean_ER_band(index_i)=Fit_ER.b1;
%             upper_bound(index_i)=Fit_ER.b1+ER_sigma(index_i)*band_sigma;
%             
%             %mean_ER_band_sim(index_i)=Fit_ER_sim.b1;
%             
%             index_i=index_i+1;
%          end
%         
%          
%          % Fit power law to the mean
%          g=fittype('a*x^b');
%          ER_mean_power_fit=fit(ER_bins(ER_bins>2)',mean_ER_band(ER_bins>2)',g,'startpoint',[ 2.5 -.1]);
%          ER_mean_fit=ER_mean_power_fit.a.*ER_bins.^(ER_mean_power_fit.b);
%             b=confint(ER_mean_power_fit,.68);
%             sigma_a_ER_mean=(b(2,1)-b(1,1))/2;
%             sigma_b_ER_mean=(b(2,2)-b(1,2))/2;
%             
%          ER_lower_power_fit=fit(ER_bins(ER_bins>2)',lower_bound(ER_bins>2)',g,'startpoint',[ 2.5 -.1]);
%          ER_lower_fit=ER_lower_power_fit.a.*ER_bins.^(ER_lower_power_fit.b);
%             b=confint(ER_lower_power_fit,.68);
%             sigma_a_ER_lower=(b(2,1)-b(1,1))/2;
%             sigma_b_ER_lower=(b(2,2)-b(1,2))/2;
%             
%          ER_upper_power_fit=fit(ER_bins(ER_bins>2)',upper_bound(ER_bins>2)',g,'startpoint',[ 2.5 -.1]);
%          ER_upper_fit=ER_upper_power_fit.a.*ER_bins.^(ER_upper_power_fit.b);
%             b=confint(ER_upper_power_fit,.68);
%             sigma_a_ER_upper=(b(2,1)-b(1,1))/2;
%             sigma_b_ER_upper=(b(2,2)-b(1,2))/2;
% 
%  no_fid=figure;
%  hold on
%  %plot(s1_phe_both_xyz(liquid_cut),log10(s2_phe_bottom_xyz(liquid_cut)./s1_phe_both_xyz(liquid_cut)),'ok','markersize',2)
%  %density(s1_phe_both_xyz(liquid_cut),log10(s2_phe_bottom_xyz(liquid_cut)./s1_phe_both_xyz(liquid_cut)),[0:.5:100],[0:.05:3]);
%  densityplot(s1_phe_both_xyz_old(liquid_cut),log10(s2_phe_bottom_xyz_old(liquid_cut)./s1_phe_both_xyz_old(liquid_cut)),'rectangle',[.5 .05],'lin');
%  cc=colorbar;
%     ylabel(cc,'Count','FontSize',16);
%  colormap;
%  ylabel('Corrected log10(S2_{b}/S1)'); xlabel('S1 corrected (Pulse Area Phe)'); 
%  title('CH3T, R<24.5 cm','fontsize',16,'Interpreter','none');
%  myfigview(16);
%  ylim([0 3]);xlim([0 s1_max]);
% %  plot(ER_bins,lower_bound,'ok','markersize',6,'linewidth',2)
%  plot(ER_bins,ER_lower_fit,'--k','markersize',14,'linewidth',2)
% %  plot(ER_bins,mean_ER_band,'ok','markersize',6,'linewidth',2)%plot mean of ER band
%  plot(ER_bins,ER_mean_fit,'-k','markersize',20,'linewidth',2)
% %  plot(ER_bins,upper_bound,'ok','markersize',6,'linewidth',2)
%  plot(ER_bins,ER_upper_fit,'--k','markersize',14,'linewidth',2)
% 
%  
%      bin_step=2;
%             total_count=length(s2_phe_bottom_xyz_old(liquid_cut));
%     
%  box;
%   text(10,2.8,strcat('Total Count= ',num2str(total_count)),'fontsize',16);
%   %text(5,.7,strcat('NEST\_SE_{bot} Fit=',num2str(ER_mean_NEST_fit.SE_ER_fit*73/64*27/24.5,3),'\pm ', num2str(sigma_SE_ER_fit*73/64*27/24.5,2),'[Phe]'),'fontsize',16,'color','m','interpreter','tex');
%   %text(5,.8,'NEST Platinum V4-b', 'fontsize',16,'color','m','interpreter','tex' );
%   text(5,.55,strcat('+',num2str(band_sigma),'\sigma = (',num2str(ER_upper_power_fit.a,3),'\pm', num2str(sigma_a_ER_upper,1) ,')\times S1_c ^{', num2str(ER_upper_power_fit.b,3),'\pm',num2str(sigma_b_ER_upper,1),'}'),'fontsize',16,'interpreter','tex');
%   text(5,.3,strcat('Mean=(',num2str(ER_mean_power_fit.a,3),'\pm', num2str(sigma_a_ER_mean,1),')\times S1_c ^{', num2str(ER_mean_power_fit.b,3),'\pm', num2str(sigma_b_ER_mean,1),'}'),'fontsize',16,'interpreter','tex');
%   text(5,.05,strcat('-',num2str(band_sigma),'\sigma = (',num2str(ER_lower_power_fit.a,3),'\pm', num2str(sigma_a_ER_lower,1),')\times S1_c ^{', num2str(ER_lower_power_fit.b,3),'\pm', num2str(sigma_b_ER_lower,1),'}'),'fontsize',16,'interpreter','tex');
%   
% 
%   
%  %%%%%%%
%  %Overlay Plot
%  %%%%%%%%%
%  
%  figure
%   hold on;
%  scatter(s1_phe_both_xyz_old(liquid_cut),log10(s2_phe_bottom_xyz_old(liquid_cut)./s1_phe_both_xyz_old(liquid_cut)),2,[0.6,0.3,0.3],'filled');
%  scatter(s1_phe_both_xyz_new(liquid_cut_new),log10(s2_phe_bottom_xyz_new(liquid_cut_new)./s1_phe_both_xyz_new(liquid_cut_new)),2,[0.3,0.3,0.6],'filled');
%   plot(ER_bins,ER_lower_fit,'--r','markersize',14,'linewidth',2)
% %  plot(ER_bins,mean_ER_band,'ok','markersize',6,'linewidth',2)%plot mean of ER band
%  plot(ER_bins,ER_mean_fit,'-r','markersize',20,'linewidth',2)
% %  plot(ER_bins,upper_bound,'ok','markersize',6,'linewidth',2)
%  plot(ER_bins,ER_upper_fit,'--r','markersize',14,'linewidth',2)
%  hold on;
%   plot(ER_bins_new,ER_lower_fit_new,'--b','markersize',14,'linewidth',2)
% %  plot(ER_bins,mean_ER_band,'ok','markersize',6,'linewidth',2)%plot mean of ER band
%  plot(ER_bins_new,ER_mean_fit_new,'-b','markersize',20,'linewidth',2)
% %  plot(ER_bins,upper_bound,'ok','markersize',6,'linewidth',2)
%  plot(ER_bins_new,ER_upper_fit_new,'--b','markersize',14,'linewidth',2)
%  ylabel('Corrected log10(S2_{b}/S1)'); xlabel('S1 corrected (Pulse Area Phe)'); legend('Run03 ER Band','Run04 ER Band');
%  title('CH3T ER Band Comparison, R<24.5 cm','fontsize',16,'Interpreter','none');
%  myfigview(16);
%  
%   box;
%   text(10,2.8,strcat('Run03 Total Count= ',num2str(total_count)),'fontsize',16,'color','r');
%   %text(5,.7,strcat('NEST\_SE_{bot} Fit=',num2str(ER_mean_NEST_fit.SE_ER_fit*73/64*27/24.5,3),'\pm ', num2str(sigma_SE_ER_fit*73/64*27/24.5,2),'[Phe]'),'fontsize',16,'color','m','interpreter','tex');
%   %text(5,.8,'NEST Platinum V4-b', 'fontsize',16,'color','m','interpreter','tex' );
%   text(5,.55,strcat('+',num2str(band_sigma),'\sigma = (',num2str(ER_upper_power_fit.a,3),'\pm', num2str(sigma_a_ER_upper,1) ,')\times S1_c ^{', num2str(ER_upper_power_fit.b,3),'\pm',num2str(sigma_b_ER_upper,1),'}'),'fontsize',16,'color','r','interpreter','tex');
%   text(5,.3,strcat('Mean=(',num2str(ER_mean_power_fit.a,3),'\pm', num2str(sigma_a_ER_mean,1),')\times S1_c ^{', num2str(ER_mean_power_fit.b,3),'\pm', num2str(sigma_b_ER_mean,1),'}'),'fontsize',16,'color','r','interpreter','tex');
%   text(5,.05,strcat('-',num2str(band_sigma),'\sigma = (',num2str(ER_lower_power_fit.a,3),'\pm', num2str(sigma_a_ER_lower,1),')\times S1_c ^{', num2str(ER_lower_power_fit.b,3),'\pm', num2str(sigma_b_ER_lower,1),'}'),'fontsize',16,'color','r','interpreter','tex');
%   
%    box;
%   text(10,2.8,strcat('Run04 Total Count= ',num2str(total_count_new)),'fontsize',16,'color','b');
%   %text(5,.7,strcat('NEST\_SE_{bot} Fit=',num2str(ER_mean_NEST_fit.SE_ER_fit*73/64*27/24.5,3),'\pm ', num2str(sigma_SE_ER_fit*73/64*27/24.5,2),'[Phe]'),'fontsize',16,'color','m','interpreter','tex');
%   %text(5,.8,'NEST Platinum V4-b', 'fontsize',16,'color','m','interpreter','tex' );
%   text(5,.55,strcat('+',num2str(band_sigma),'\sigma = (',num2str(ER_upper_power_fit_new.a,3),'\pm', num2str(sigma_a_ER_upper_new,1) ,')\times S1_c ^{', num2str(ER_upper_power_fit_new.b,3),'\pm',num2str(sigma_b_ER_upper_new,1),'}'),'fontsize',16,'color','b','interpreter','tex');
%   text(5,.3,strcat('Mean=(',num2str(ER_mean_power_fit_new.a,3),'\pm', num2str(sigma_a_ER_mean_new,1),')\times S1_c ^{', num2str(ER_mean_power_fit_new.b,3),'\pm', num2str(sigma_b_ER_mean_new,1),'}'),'fontsize',16,'color','b','interpreter','tex');
%   text(5,.05,strcat('-',num2str(band_sigma),'\sigma = (',num2str(ER_lower_power_fit_new.a,3),'\pm', num2str(sigma_a_ER_lower_new,1),')\times S1_c ^{', num2str(ER_lower_power_fit_new.b,3),'\pm', num2str(sigma_b_ER_lower_new,1),'}'),'fontsize',16,'color','b','interpreter','tex');
%   
%   %%%%%%%%%%
%   
%   
%     clear y_ER_norm Fit_ER Fit_ER_norm mean_ER_band_norm gauss_fit_bins normalized_s2_s1
%     normalized_s2_s1=log10(s2_phe_bottom_xyz_old(liquid_cut)./s1_phe_both_xyz_old(liquid_cut))-interp1(ER_bins,ER_mean_fit,s1_phe_both_xyz_old(liquid_cut)) ;
%             
%            ER_bin_start=-.51;
%            ER_bin_end=0.51;
%            ER_bin_size=.01;
%            gauss_fit_bins=ER_bin_start-ER_bin_size/2:ER_bin_size:ER_bin_end+ER_bin_size/2;
%            mean_ER_band_norm=hist(normalized_s2_s1,gauss_fit_bins)/ER_bin_size;
%            Fit_ER_norm=fit(gauss_fit_bins',mean_ER_band_norm','gauss1');
%            y_ER_norm=Fit_ER_norm.a1.*exp(-((gauss_fit_bins-Fit_ER_norm.b1)./Fit_ER_norm.c1).^2);
%            sigma_error=confint(Fit_ER_norm,0.683);
%            sigma_error_upper=sigma_error(2,3)/sqrt(2);
%            sigma_error_lower=sigma_error(1,3)/sqrt(2);
%            
%      ER_band_gauss_fig=figure;
% 
%      stairs(gauss_fit_bins-ER_bin_size/2,mean_ER_band_norm,'linewidth',2)
%      hold on
%      plot(gauss_fit_bins,y_ER_norm,'-r','linewidth',2);
%      xlim([-.5 .5])
%      set(gca, 'XTickMode', 'manual');
%         set(gca, 'XTick', [-0.5:.1:.5]);
%      xlabel('Corrected log10(S2_b/S1)-Expected ');ylabel('Count/0.01 Phe');
%      title(strcat('Run 03 CH3T, R<24.5 cm'),'fontsize',16,'Interpreter','none');
%      legend('CH3T data', strcat('\sigma =',num2str(Fit_ER_norm.c1/sqrt(2),4),'(',num2str(sigma_error_lower,4),',',num2str(sigma_error_upper,4),')')); %was divided by 2
%      grid
%      myfigview(16)
%   
%   
%   %%%% S1 v S2 scatter overlay
% %   
% %         
% %     figure
% %     scatter(s1_phe_both_xyz,s2_phe_bottom_xyz,'.r');
% %     hold on;
% %     scatter(s1_phe_both_xyz_new,s2_phe_bottom_xyz_new,'.g');
% %     xlabel('S1');ylabel('S2');legend('Run03','Run04');
%   eff_boost=(10.^mean_ER_band_new)./(10.^mean_ER_band);
%   
% figure
% plot(ER_bins,eff_boost)
% mean(eff_boost)
% xlabel('Corrected S1 (phe)')
% ylabel('Shift in Run04 Centroid from Run03')
% myfigview(16)
% title('ER Centroid Shift from Run03 to Run04')
% 
%   
%   relative_width_new= (10.^ER_upper_fit_new-10.^mean_ER_band_new)./(10.^mean_ER_band_new);
%   relative_width_old= (10.^ER_upper_fit-10.^mean_ER_band)./(10.^mean_ER_band);
