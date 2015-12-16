       %% Kr83m at 42.6 keV and 100 V/cm  -- NEEDS ITS OWN XYZ Corrections (LOAD DATA GIVES NO XYZ CORRECTION)
         clear radius s2x s2y drift_time s1_phe_both_xyz s2_phe_bottom_xyz S2_numelectrons_VUV
load('C:\Program Files\MATLAB\R2012a\bin\Run03Kr_100Vcm_BotAndBoth.mat'); %Load Patched Kr Data

%% Setting variables
s2radius=(s2x.^2+s2y.^2).^(1/2);
grid_size=2;
rcut_min = 0;
rcut_max = 25;%cm
s1area_bound_min = 20;%100
s1area_bound_max = 500;%600

s2area_bound_min = 100;%200
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

det_edge=351;

    zcut_min = 0.1*det_edge;
    zcut_max = 0.95*det_edge;


 
    s1_xbins = (S1_xbin_max-S1_xbin_min)./S1xybinsize;
    s1_ybins = (S1_ybin_max-S1_ybin_min)./S1xybinsize;

    s2_xbins = (S2_xbin_max-S2_xbin_min)./S2xybinsize;
    s2_ybins = (S2_ybin_max-S2_ybin_min)./S2xybinsize;
        
s1_x_bins=(S1_xbin_min+S1xybinsize/2):S1xybinsize:(S1_xbin_max-S1xybinsize/2);
s1_y_bins=(S1_ybin_min+S1xybinsize/2):S1xybinsize:(S1_ybin_max-S1xybinsize/2);

s1_phe_both=s1_phe_both_xyz_VUV;
s1_phe_bottom=s1_phe_bot_xyz_VUV;
s2_phe_both=s2_phe_both_xyz_VUV;
s2_phe_bottom=s2_phe_bot_xyz_VUV;


%%
%Getting s2 z correction
    
dT_step=15;
dT_max=30+dT_step;
i=1;
det_edge=351;

for dT_max=30+dT_step:dT_step:310-dT_step;
s2z_mean_fit=fit([0:100:15000]',hist(s2_phe_both_xyz_VUV(inrange(drift_time,[dT_max-dT_step dT_max])),[0:100:15000])','gauss1'); 
s2z_mean(i)=s2z_mean_fit.b1;
s2z_dT_loc(i)=dT_max-dT_step/2;
i=i+1;
end

s2z_mean_fit=fit(s2z_dT_loc.',s2z_mean.','exp1');
electron_lifetime=-1/s2z_mean_fit.b;
s2_phe_both_z=s2_phe_both_xyz_VUV.*exp(drift_time./electron_lifetime);
s2_phe_bot_z=s2_phe_bot_xyz_VUV.*exp(drift_time./electron_lifetime);


%%
n_sig=3;
     radius=(s2x.^2+s2y.^2).^(1/2);
  
    clear y_s1_both y_s1_bottom Fit_s1_both Fit_s1_bottom mean_s1_both mean_s1_both
    clear sigma_s1_both sigma_s1_both mean_index means mean_index means_error x temp_hist
    clear hist_s1_both x xfit S
    
    %set up energy bins for the gaussian fit
    bin_s1=5;%5
    x=50:bin_s1:500; %150:bin_s1:500;
    x=x';
        
    cut = inrange(s2radius,[rcut_min,rcut_max]) & inrange(s1_phe_both,[s1area_bound_min,s1area_bound_max]);% & inrange(s2_width,[100,800]) & inrange(s1_width,[30,400]);   
    
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
        Fit_s1_both{bin}=fit(x(2:end-1),temp_hist(2:end-1),'gauss1','start',[amp_start mean_start sigma_start]);%remove bin edges

     means(bin)=Fit_s1_both{bin}.b1;
     means_error(bin) = Fit_s1_both{bin}.c1/sqrt(2)/sqrt(sum(hist_s1_both(:,bin))*bin_s1); % == sigma/sqrt(N);
     mean_index(bin) = (binStart+binEnd)/2;     

    end
    
    s1_both_means=means;
    s1_both_means_sigma=means_error;
    s1_both_means_bin=mean_index;
    
    
    dT_range= inrange( mean_index,[zcut_min, zcut_max] );
    dT_fit=0:1:det_edge;    
    [P_s1_both S] = polyfit(mean_index(dT_range),means(dT_range),2);
    sigma_P_s1_both = sqrt(diag(inv(S.R)*inv(S.R')).*S.normr.^2./S.df)';    
    
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
    bin_s1=5;
    x=50:bin_s1:400;
    x=x';
        
    
    cut = inrange(s2radius,[rcut_min,rcut_max]) & inrange(s1_phe_both,[s1area_bound_min,s1area_bound_max]);% & inrange(s2_width,[100,800]) & inrange(s1_width,[30,400]);
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
        Fit_s1_bottom{bin}=fit(x(2:end-1),temp_hist(2:end-1),'gauss1','start',[amp_start mean_start sigma_start]);%remove bin edges
     
      
     means(bin)=Fit_s1_bottom{bin}.b1;
     means_error(bin) = Fit_s1_bottom{bin}.c1/sqrt(2)/sqrt(sum(hist_s1_bottom(:,bin))*bin_s1); % == sigma/sqrt(N);
     mean_index(bin) = (binStart+binEnd)/2;     

    end

    s1_bottom_means=means;
    s1_bottom_means_sigma=means_error;
    s1_bottom_means_bin=mean_index;
     
    
    dT_range= inrange( mean_index,[zcut_min, zcut_max] );
    dT_fit=0:1:det_edge;    
    [P_s1_bottom S] = polyfit(mean_index(dT_range),means(dT_range),2);
    sigma_P_s1_bottom = sqrt(diag(inv(S.R)*inv(S.R')).*S.normr.^2./S.df)';
 
    

%% S1 XY

clear hist_bin x

s1_bin=5;
x=100:s1_bin:500;
x=x';

%set up matricies to be filled
mean_S1_both_xy = zeros(s1_xbins,s1_ybins);
Sigma_S1_both = zeros(s1_xbins,s1_ybins);
Count_S1_both = zeros(s1_xbins,s1_ybins);

% s1_z_correction=polyval(P_s1_both,241.6)./polyval(P_s1_both,drift_distance_\mus);%normalize s1 to 160 ns.
s1_z_correction=polyval(P_s1_both,(det_edge-4.32)/2)./polyval(P_s1_both,drift_time);%normalize s1 to 160 ns.

s1_phe_both_z=s1_phe_both.*s1_z_correction;
s2_phe_both_z=s2_phe_both.*exp(drift_time/electron_lifetime); %Correct S2 for lifeimte  %s2_phe_both_z[5000,18000]



time_energy_cut = inrange(drift_time,[zcut_min,zcut_max]) & inrange(s1_phe_both_z,[100,s1area_bound_max]) & inrange(s2_phe_both_z,[1000,s2area_bound_max]) & inrange(s2radius,[rcut_min,rcut_max]);


%%%%variables for uber turbo boost!
s1_phe_both_z_2=s1_phe_both_z(time_energy_cut);
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
      
       hist_bin = hist(s1_phe_both_z_2(bin_cut) ,x)'/s1_bin; 
       Count_S1_both(y_count,x_count)=sum(hist_bin)*s1_bin; % X and Y flip ... because MATLAB
       
       %%%%%%%%Truncate variables after binning is complete
        s1_phe_both_z_2=s1_phe_both_z_2(~bin_cut);
        s2x_2=s2x_2(~bin_cut);
        s2y_2=s2y_2(~bin_cut);
        clear bin_cut;
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
      if Count_S1_both(y_count,x_count)>30; 
                        
                    amp_start=max(hist_bin);
                    mean_start=sum(hist_bin.*x)/sum(hist_bin);
                    sigma_start=std(hist_bin.*x);                     
                Fit_S1_both_xy= fit(x(2:end-1),hist_bin(2:end-1),'gauss1','start',[amp_start mean_start sigma_start]); %remove edges
                
                %Peak location (mean)
                mean_S1_both_xy(y_count,x_count)=Fit_S1_both_xy.b1; % X and Y flip ... because MATLAB
                
                %1-sigma of peak position
                Sig_b=Fit_S1_both_xy.c1/sqrt(2)/sqrt(Count_S1_both(y_count,x_count)); %sigma/sqrt(N)
   
                    %Fit error checking
                        if(isnan(Sig_b))
                          mean_S1_both_xy(y_count,x_count)=0;
                         Sig_b=0;
                        end
                                
                    %uncertainty in mean
                    Sigma_S1_both(y_count,x_count)=Sig_b;% X and Y flip ... because MATLAB   
     
      else
          
        mean_S1_both_xy(y_count,x_count)=0;  
       
      end            
    end    
   end
   
   mean_S1_both_xy(isnan(mean_S1_both_xy)) = 0;
   
  
%Calculate corrections matrix and 1 sigma corrections matrix
center=interp2(s1_x_bins,s1_y_bins,mean_S1_both_xy,0,0,'cubic');%Normalize to the center (x=y=0)
norm_S1_all=center./mean_S1_both_xy; 
norm_S1_all(isinf(norm_S1_all))=1;%no infinity! Do not correct events outside r=25cm
norm_S1_all(isnan(norm_S1_all))=1;

sigma_center=interp2(s1_x_bins,s1_y_bins,Sigma_S1_both,0,0,'cubic');%Normalize to the center (x=y=0)
sigma_norm_S1_all=sqrt((sigma_center./mean_S1_both_xy).^2+(Sigma_S1_both.*center./mean_S1_both_xy.^2).^2);
sigma_norm_S1_all(isinf(sigma_norm_S1_all))=1;%no infinity. Set Sigma=1, which is 100% uncertainty
sigma_norm_S1_all(isnan(sigma_norm_S1_all))=1;%no nan. Set Sigma=1, which is 100% uncertainty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%  Mapping S1 in XY (after correction for S1 Z-dependence). Bottom PMT
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear x hist_bin

s1_bin=5;
x=50:s1_bin:400;
x=x';

%set up matricies to be filled
mean_S1_bottom_xy = zeros(s1_xbins,s1_ybins);
Sigma_S1_bottom = zeros(s1_xbins,s1_ybins);
Count_S1_bottom = zeros(s1_xbins,s1_ybins);

% s1_z_correction_bottom=polyval(P_s1_bottom,241.6)./polyval(P_s1_bottom,drift_distance_\mus);%normalize s1 to 160 ns.
s1_z_correction_bottom=polyval(P_s1_bottom,(det_edge-4.32)/2)./polyval(P_s1_bottom,drift_time);%normalize s1 to 160 ns.

s1_phe_bottom_z=s1_phe_bottom.*s1_z_correction_bottom;
% s2_phe_both_z=s2_phe_both.*exp(drift_distance_\mus/electron_lifetime); %Correct S2 for lifeimte
s2_phe_both_z=s2_phe_both.*exp(drift_time/electron_lifetime); %Correct S2 for lifeimte


% time_energy_cut = inrange(drift_distance_\mus,[zcut_min,zcut_max]) & inrange(s1_phe_both_z,[200,500]) & inrange(s2_phe_both_z,[1000,20000]) & inrange(s2radius,[rcut_min,rcut_max]) ...
%     &inrange(s2_width,[100 800]) & inrange(s1_width,[30 400]);
time_energy_cut = inrange(drift_time,[zcut_min,zcut_max]) & inrange(s1_phe_both_z,[100,s1area_bound_max]) & inrange(s2_phe_both_z,[1000,s2area_bound_max]) & inrange(s2radius,[rcut_min,rcut_max]);


%%%%variables for uber turbo boost!
s1_phe_bottom_z_2=s1_phe_bottom_z(time_energy_cut);
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
      
       hist_bin = hist(s1_phe_bottom_z_2(bin_cut) ,x)'/s1_bin; 
       Count_S1_bottom(y_count,x_count)=sum(hist_bin)*s1_bin; % X and Y flip ... because MATLAB
       
       %%%%%%%%Truncate variables after binning is complete
        s1_phe_bottom_z_2=s1_phe_bottom_z_2(~bin_cut);
        s2x_2=s2x_2(~bin_cut);
        s2y_2=s2y_2(~bin_cut);
        clear bin_cut;
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
      if Count_S1_bottom(y_count,x_count)>30; 
                        
                    amp_start=max(hist_bin);
                    mean_start=sum(hist_bin.*x)/sum(hist_bin);
                    sigma_start=std(hist_bin.*x);    
                Fit_S1_bottom_xy= fit(x(2:end-1),hist_bin(2:end-1),'gauss1','start',[amp_start mean_start sigma_start]); %remove edges
                
                %Peak location (mean)
                mean_S1_bottom_xy(y_count,x_count)=Fit_S1_bottom_xy.b1; % X and Y flip ... because MATLAB
                
                %1-sigma of peak position
                Sig_b=Fit_S1_bottom_xy.c1/sqrt(2)/sqrt(Count_S1_bottom(y_count,x_count)); %sigma/sqrt(N)
   
                    %Fit error checking
                        if(isnan(Sig_b))
                         mean_S1_bottom_xy(y_count,x_count)=0;
                         Sig_b=0;
                        end                    
              
                    %uncertainty in mean
                    Sigma_S1_bottom(y_count,x_count)=Sig_b;% X and Y flip ... because MATLAB        
      else
          
        mean_S1_bottom_xy(y_count,x_count)=0;  
       
      end           
    end    
   end
   
   mean_S1_bottom_xy(isnan(mean_S1_bottom_xy)) = 0;
       

%Calculate corrections matrix and 1 sigma corrections matrix
center=interp2(s1_x_bins,s1_y_bins,mean_S1_bottom_xy,0,0,'cubic');%Normalize to the center (x=y=0)
norm_S1_bot=center./mean_S1_bottom_xy; 
norm_S1_bot(isinf(norm_S1_bot))=1;%no infinity! Do not correct events outside r=25cm
norm_S1_bot(isnan(norm_S1_bot))=1;%no NaN

sigma_center=interp2(s1_x_bins,s1_y_bins,Sigma_S1_bottom,0,0,'cubic');%Normalize to the center (x=y=0)
sigma_norm_S1_bot=sqrt((sigma_center./mean_S1_bottom_xy).^2+(Sigma_S1_bottom.*center./mean_S1_bottom_xy.^2).^2);
sigma_norm_S1_bot(isinf(sigma_norm_S1_bot))=1;%no infinity. Set Sigma=1, which is 100% uncertainty
sigma_norm_S1_bot(isnan(sigma_norm_S1_bot))=1;%no nan. Set Sigma=1, which is 100% uncertainty


%%
%Getting s2 xy correction

clear x2 hist_bin

xx_step=2;
yy_step=2;
s2_both_xymeans = zeros(floor(50/yy_step),floor(50/xx_step));
s2_bot_xymeans = zeros(floor(50/yy_step),floor(50/xx_step));

i=1;
   for x_bin=-25+xx_step/2:xx_step:25-xx_step/2
       j=1;
    for y_bin=-25+yy_step/2:yy_step:25-yy_step/2              
      if length(s2_phe_both_z(inrange(s2x,[x_bin-xx_step x_bin+xx_step]) & inrange(s2y,[y_bin-yy_step y_bin+yy_step]))) > 300;
            s2_both_xymeans_fit=fit([0:200:20000]',hist(s2_phe_both_z(inrange(s2x,[x_bin-xx_step x_bin+xx_step]) & inrange(s2y,[y_bin-yy_step y_bin+yy_step])),[0:200:20000])','gauss1'); 
            s2_both_xymeans(j,i)=s2_both_xymeans_fit.b1;
            
            s2_bot_xymeans_fit=fit([0:200:20000]',hist(s2_phe_bot_z(inrange(s2x,[x_bin-xx_step x_bin+xx_step]) & inrange(s2y,[y_bin-yy_step y_bin+yy_step])),[0:200:20000])','gauss1'); 
            s2_bot_xymeans(j,i)=s2_bot_xymeans_fit.b1;
           else 
          
        s2_both_xymeans(j,i)=0;     
        s2_bot_xymeans(j,i)=0;     
        
      end
      j=j+1;
    end
    i=i+1;
   end
   
s2_both_xymeans(isnan(s2_both_xymeans)) = 0;
s2_bot_xymeans(isnan(s2_bot_xymeans)) = 0;


s2x_bins=-25+xx_step/2:xx_step:25-xx_step/2;
s2y_bins=-25+yy_step/2:yy_step:25-yy_step/2;

    %Calculate corrections matrix and 1 sigma corrections matrix
center=interp2(s2x_bins,s2y_bins,s2_both_xymeans,0,0,'spline',1);%Normalize to the center (x=y=0)
norm_S2_both=center./s2_both_xymeans; 
norm_S2_both(isinf(norm_S2_both))=1;%no infinity! Do not correct events outside r=25cm
norm_S2_both(isnan(norm_S2_both))=1;

center_bot=interp2(s2x_bins,s2y_bins,s2_bot_xymeans,0,0,'spline',1);%Normalize to the center (x=y=0)
norm_S2_bot=center_bot./s2_bot_xymeans; 
norm_S2_bot(isinf(norm_S2_bot))=1;%no infinity! Do not correct events outside r=25cm
norm_S2_bot(isnan(norm_S2_bot))=1;

%% S1 XYZ correction
%bin the detector in 30mm for z

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
x=10:bin:500; %Set up S1 energy bins
x=x';

cut_fit=x>80 & x<490; %remove bin edges

xy_step=2;% 2 cm bins
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


Count_S1_3D=zeros(25,25,16);
Count_S1_3D_bottom=zeros(25,25,16);
Mean_S1_3D_both=zeros(25,25,16);
Mean_S1_3D_bottom=zeros(25,25,16);
Sigma_S1_3D_both=zeros(25,25,16);
Sigma_S1_3D_bottom=zeros(25,25,16);

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
             
   %///////Now do the Bottom PMT arraycalculation//////////////////////////////
           
            clear yqfit Sig_b
   
            if (Count_S1_3D_bottom(l,m,n) >= 30) % at least 10 counts before fitting. 
                
                yqfit_bottom=hist_q_bottom;
                    amp_start=max(yqfit_bottom(cut_fit));
                    mean_start=sum( yqfit_bottom(cut_fit).*x(cut_fit))/sum( yqfit_bottom(cut_fit));
                    sigma_start=std( yqfit_bottom(cut_fit).*x(cut_fit));    
                Fit_q_bottom= fit(x(cut_fit),yqfit_bottom(cut_fit),'gauss1','start',[amp_start mean_start sigma_start]);
                
                %Peak location (mean)
                Mean_S1_3D_bottom(l,m,n)=Fit_q_bottom.b1;
                
                %1-sigma of peak position
                Sig_b=Fit_q_bottom.c1/sqrt(2)/sqrt(Count_S1_3D_bottom(l,m,n));
   

                %error checking
                   if(strcmp(num2str(Sig_b),'NaN'))
                     Mean_S1_3D_bottom(l,m,n)=0;
                     Sig_b=0;
                   end
                  
                %uncertainty in mean
                Sigma_S1_3D_bottom(l,m,n)=Sig_b;
                
             %end IF
            
             else %not enough stats to do the fit
                
                Fit_q_bottom= 0;
                %Find Peak location
                Mean_S1_3D_bottom(l,m,n)=0;
                %1-sigma
                Sigma_S1_3D_bottom(l,m,n)=0;
            end          
                                  
              
        end
    end
    %get 3D raw-count matrix
    %Count_q(:,:,n)=squeeze(sum(squeeze(hist_q(:,:,:,n)),1))*bin;
toc;
k

end

s1xbins=-r_max+xy_step/2:xy_step:r_max-xy_step/2;
s1ybins=-r_max+xy_step/2:xy_step:r_max-xy_step/2;


%% S1 3D correction. Both PMT (normalized to dT=160)

center_phe_both_xyz=interp3(s1xbins,s1ybins,s1zbins,Mean_S1_3D_both,0,0,(det_edge-4.32)/2,'spline');%Normalize to the center (x=y=0. dT=160us)
norm_s1_both_xyz=center_phe_both_xyz./Mean_S1_3D_both; %Normalize to center (x=y=0. dT=160)
norm_s1_both_xyz(isinf(norm_s1_both_xyz))=1;%remove infinity. no correction outside 25cm
norm_s1_both_xyz(isnan(norm_s1_both_xyz))=1;%remove nan. no correction outside 25cm

%Get 1 sigma of the corrections matrix.
sigma_center_phe_both_xyz=interp3(s1xbins,s1ybins,s1zbins,Sigma_S1_3D_both,0,0,(det_edge-4.32)/2,'spline');%Normalize to the center (x=y=0. dT=160us)
sigma_norm_s1_both_xyz=sqrt((sigma_center_phe_both_xyz./Mean_S1_3D_both).^2+(Sigma_S1_3D_both.*center_phe_both_xyz./Mean_S1_3D_both.^2).^2);
sigma_norm_s1_both_xyz(isinf(sigma_norm_s1_both_xyz))=1;%no infinity. Set Sigma=1, which is 100% uncertainty
sigma_norm_s1_both_xyz(isnan(sigma_norm_s1_both_xyz))=1;%no nan. Set Sigma=1, which is 100% uncertainty


mean_center_3D_both=interp3(s1xbins,s1ybins,s1zbins,Mean_S1_3D_both,0,0,(det_edge-4.32)/2,'spline');
sigma_mean_center_3D_both=interp3(s1xbins,s1ybins,s1zbins,Sigma_S1_3D_both,0,0,(det_edge-4.32)/2,'spline');
%to Apply the correction use the following.
%s1_phe_both_xyz=s1_phe_both.*interp3(xx,yy,zz,norm_s1_both_xyz,s2x,s2y,drift_time,'spline');


%% S1 3D correction. Bottom PMT (normalized to dT=160)

center_phe_bottom_xyz=interp3(s1xbins,s1ybins,s1zbins,Mean_S1_3D_bottom,0,0,(det_edge-4.32)/2,'spline');%Normalize to the center (x=y=0. dT=160us)
norm_s1_bottom_xyz=center_phe_bottom_xyz./Mean_S1_3D_bottom; %Normalize to center (x=y=0. dT=160)
norm_s1_bottom_xyz(isinf(norm_s1_bottom_xyz))=1;%remove infinity. no correction outside 25cm
norm_s1_bottom_xyz(isnan(norm_s1_bottom_xyz))=1;%remove nan. no correction outside 25cm

%Get 1 sigma of the corrections matrix.
sigma_center_phe_bottom_xyz=interp3(s1xbins,s1ybins,s1zbins,Sigma_S1_3D_bottom,0,0,(det_edge-4.32)/2,'spline');%Normalize to the center (x=y=0. dT=160us)
sigma_norm_s1_bottom_xyz=sqrt((sigma_center_phe_bottom_xyz./Mean_S1_3D_bottom).^2+(Sigma_S1_3D_bottom.*center_phe_bottom_xyz./Mean_S1_3D_bottom.^2).^2);
sigma_norm_s1_bottom_xyz(isinf(sigma_norm_s1_bottom_xyz))=1;%no infinity. Set Sigma=1, which is 100% uncertainty
sigma_norm_s1_bottom_xyz(isnan(sigma_norm_s1_bottom_xyz))=1;%no nan. Set Sigma=1, which is 100% uncertainty

mean_center_3D_bottom=interp3(s1xbins,s1ybins,s1zbins,Mean_S1_3D_bottom,0,0,(det_edge-4.32)/2,'spline');
sigma_mean_center_3D_bottom=interp3(s1xbins,s1ybins,s1zbins,Sigma_S1_3D_bottom,0,0,(det_edge-4.32)/2,'spline');

%to Apply the correction use the following.
%s1_phe_bottom_xyz=s1_phe_bottom.*interp3(xx,yy,zz,norm_s1_bottom_xyz,s2x,s2y,drift_time,'spline');
% s1_phe_both_xyz_z=s1_phe_both_xyz_VUV.*polyval(mean_fit_s1z,(det_edge-4.32)/2)./polyval(mean_fit_s1z,drift_time);
% s1_phe_both_xyz_VUV=s1_phe_both_xyz_z.*interp2(s2x_bins,s2y_bins,norm_S1_both,s2x,s2y,'spline');
% s2_phe_both_z=s2_phe_both_xyz_VUV.*exp(drift_time./electron_lifetime);
% s2_phe_both_xyz_VUV=s2_phe_both_z.*interp2(s2x_bins,s2y_bins,norm_S2_bot,s2x,s2y,'spline');



save('100Vcm_Corr','P_s1_both','norm_s1_bottom_xyz','norm_s1_both_xyz','s1xbins','s1ybins','s1zbins','P_s1_bottom','det_edge','norm_S1_all','norm_S1_bot','s1_x_bins','s1_y_bins','s2x_bins','s2y_bins','electron_lifetime','norm_S2_both','norm_S2_bot');