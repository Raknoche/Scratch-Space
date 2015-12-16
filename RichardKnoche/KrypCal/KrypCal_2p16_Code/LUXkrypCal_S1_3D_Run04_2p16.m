function [lug_val] = LUXkrypCal_S1_3D_2p16(s1_phe_both,s1_phe_bottom,s2_phe_both,s2x,s2y,drift_time,x_center,y_center,z_center,det_edge,inpaint_on,user_name,file_id_cp,dp_version, algorithm_name, submit)
%% Get 3D S1 correction (need about 1,000,000 events)
%bin the detector in 30mm for z

%////////////////Set Up Variables

s2_cut=inrange(s2_phe_both,[200 40000]);%cut on appropriate S2 size

s1=s1_phe_both(s2_cut); %No Corrections. Use Bottom, Top or Both PMT.
s1_bottom=s1_phe_bottom(s2_cut);
xc=s2x(s2_cut);
yc=s2y(s2_cut);
dT=drift_time(s2_cut);

Kr_events=length(s1);


%////////////////

%s1=s1_xyz;

bin=5;
x=40:bin:500; %Set up S1 energy bins
x=x';

cut_fit=x>50 & x<490; %remove bin edges

xy_step=3;% 3 cm bins
r_max=25; % max radius

z_step=ceil( (0.95*det_edge-10)/15 );  %xy_step size about = 30 mm, to create 16 z bins
s1zbins=floor(10):z_step:floor(10)+z_step*15; %20 us


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


Count_S1_3D=zeros(16,16,16);
Count_S1_3D_bottom=zeros(16,16,16);
Mean_S1_3D_both=zeros(16,16,16);
Mean_S1_3D_bottom=zeros(16,16,16);
Sigma_S1_3D_both=zeros(16,16,16);
Sigma_S1_3D_bottom=zeros(16,16,16);

for k = s1zbins; % out to 485 mm % start at about 20 mm bin center. xy_steps of 30 mm
    
    tic;
    
    for i = (-r_max+xy_step):xy_step:r_max %map to rows going down (y_cm)
        for j = (-r_max+xy_step):xy_step:r_max %map columns across (x_cm)
                       
            l=int8(i/xy_step+(r_max)/xy_step); %real y from 1:25
            m=int8(j/xy_step+(r_max)/xy_step); % real x from 1:25
            n=(k-floor(10))/z_step + 1; %z from 1:16
            
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
            
            if (Count_S1_3D(l,m,n) >= 50) % at least 30 counts before fitting. 
                
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
   
            if (Count_S1_3D_bottom(l,m,n) >= 50) % at least 10 counts before fitting. 
                
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


%get the correction for the top PMT array by taking Both-Bottom
Mean_S1_3D_top=Mean_S1_3D_both-Mean_S1_3D_bottom;
Sigma_S1_3D_top=sqrt(Sigma_S1_3D_both.^2+Sigma_S1_3D_bottom.^2);


%% S1 3D correction. Both PMT (normalized to dT=160)
 if inpaint_on==1;  Mean_S1_3D_both(Mean_S1_3D_both==0)=nan; Mean_S1_3D_both=inpaint_nans3(Mean_S1_3D_both,0); end;
center_phe_both_xyz=interp3(s1xbins,s1ybins,s1zbins,Mean_S1_3D_both,x_center,y_center,z_center,'spline');%Normalize to the center (x=y=0. dT=160us)
norm_s1_both_xyz=center_phe_both_xyz./Mean_S1_3D_both; %Normalize to center (x=y=0. dT=160)
norm_s1_both_xyz(isinf(norm_s1_both_xyz))=1;%remove infinity. no correction outside 25cm
norm_s1_both_xyz(isnan(norm_s1_both_xyz))=1;%remove nan. no correction outside 25cm

%Get 1 sigma of the corrections matrix.
 if inpaint_on==1;  Sigma_S1_3D_both(Sigma_S1_3D_both==0)=nan; Sigma_S1_3D_both=inpaint_nans3(Sigma_S1_3D_both,0); end;
sigma_center_phe_both_xyz=interp3(s1xbins,s1ybins,s1zbins,Sigma_S1_3D_both,x_center,y_center,z_center,'spline');%Normalize to the center (x=y=0. dT=160us)
sigma_norm_s1_both_xyz=sqrt((sigma_center_phe_both_xyz./Mean_S1_3D_both).^2+(Sigma_S1_3D_both.*center_phe_both_xyz./Mean_S1_3D_both.^2).^2);
sigma_norm_s1_both_xyz(isinf(sigma_norm_s1_both_xyz))=1;%no infinity. Set Sigma=1, which is 100% uncertainty
sigma_norm_s1_both_xyz(isnan(sigma_norm_s1_both_xyz))=1;%no nan. Set Sigma=1, which is 100% uncertainty

mean_center_3D_both=interp3(s1xbins,s1ybins,s1zbins,Mean_S1_3D_both,x_center,y_center,z_center,'spline');
sigma_mean_center_3D_both=interp3(s1xbins,s1ybins,s1zbins,Sigma_S1_3D_both,x_center,y_center,z_center,'spline');
%to Apply the correction use the following.
%s1_phe_both_xyz=s1_phe_both.*interp3(xx,yy,zz,norm_s1_both_xyz,s2x,s2y,drift_time,'spline');


%% S1 3D correction. Bottom PMT (normalized to dT=160)
 if inpaint_on==1;  Mean_S1_3D_bottom(Mean_S1_3D_bottom==0)=nan; Mean_S1_3D_bottom=inpaint_nans3(Mean_S1_3D_bottom,0); end;
center_phe_bottom_xyz=interp3(s1xbins,s1ybins,s1zbins,Mean_S1_3D_bottom,x_center,y_center,z_center,'spline');%Normalize to the center (x=y=0. dT=160us)
norm_s1_bottom_xyz=center_phe_bottom_xyz./Mean_S1_3D_bottom; %Normalize to center (x=y=0. dT=160)
norm_s1_bottom_xyz(isinf(norm_s1_bottom_xyz))=1;%remove infinity. no correction outside 25cm
norm_s1_bottom_xyz(isnan(norm_s1_bottom_xyz))=1;%remove nan. no correction outside 25cm

%Get 1 sigma of the corrections matrix.
 if inpaint_on==1;  Sigma_S1_3D_bottom(Sigma_S1_3D_bottom==0)=nan; Sigma_S1_3D_bottom=inpaint_nans3(Sigma_S1_3D_bottom,0); end;
sigma_center_phe_bottom_xyz=interp3(s1xbins,s1ybins,s1zbins,Sigma_S1_3D_bottom,x_center,y_center,z_center,'spline');%Normalize to the center (x=y=0. dT=160us)
sigma_norm_s1_bottom_xyz=sqrt((sigma_center_phe_bottom_xyz./Mean_S1_3D_bottom).^2+(Sigma_S1_3D_bottom.*center_phe_bottom_xyz./Mean_S1_3D_bottom.^2).^2);
sigma_norm_s1_bottom_xyz(isinf(sigma_norm_s1_bottom_xyz))=1;%no infinity. Set Sigma=1, which is 100% uncertainty
sigma_norm_s1_bottom_xyz(isnan(sigma_norm_s1_bottom_xyz))=1;%no nan. Set Sigma=1, which is 100% uncertainty

mean_center_3D_bottom=interp3(s1xbins,s1ybins,s1zbins,Mean_S1_3D_bottom,x_center,y_center,z_center,'spline');
sigma_mean_center_3D_bottom=interp3(s1xbins,s1ybins,s1zbins,Sigma_S1_3D_bottom,x_center,y_center,z_center,'spline');
%to Apply the correction use the following.
%s1_phe_bottom_xyz=s1_phe_bottom.*interp3(xx,yy,zz,norm_s1_bottom_xyz,s2x,s2y,drift_time,'spline');

%% S1 3D correction. Top PMT (normalized to dT=160)
 if inpaint_on==1;  Mean_S1_3D_top(Mean_S1_3D_top==0)=nan; Mean_S1_3D_top=inpaint_nans3(Mean_S1_3D_top,0); end;
center_phe_top_xyz=interp3(s1xbins,s1ybins,s1zbins,Mean_S1_3D_top,x_center,y_center,z_center,'spline');%Normalize to the center (x=y=0. dT=160us)
norm_s1_top_xyz=center_phe_top_xyz./Mean_S1_3D_top; %Normalize to center (x=y=0. dT=160)
norm_s1_top_xyz(isinf(norm_s1_top_xyz))=1;%remove infinity. no correction outside 25cm
norm_s1_top_xyz(isnan(norm_s1_top_xyz))=1;%remove nan. no correction outside 25cm

%Get 1 sigma of the corrections matrix.
 if inpaint_on==1;  Sigma_S1_3D_top(Sigma_S1_3D_top==0)=nan; Sigma_S1_3D_top=inpaint_nans3(Sigma_S1_3D_top,0); end;
sigma_center_phe_top_xyz=interp3(s1xbins,s1ybins,s1zbins,Sigma_S1_3D_top,x_center,y_center,z_center,'spline');%Normalize to the center (x=y=0. dT=160us)
sigma_norm_s1_top_xyz=sqrt((sigma_center_phe_top_xyz./Mean_S1_3D_top).^2+(Sigma_S1_3D_top.*center_phe_top_xyz./Mean_S1_3D_top.^2).^2);
sigma_norm_s1_top_xyz(isinf(sigma_norm_s1_top_xyz))=1;%no infinity. Set Sigma=1, which is 100% uncertainty
sigma_norm_s1_top_xyz(isnan(sigma_norm_s1_top_xyz))=1;%no nan. Set Sigma=1, which is 100% uncertainty


mean_center_3D_top=interp3(s1xbins,s1ybins,s1zbins,Mean_S1_3D_top,x_center,y_center,z_center,'spline');
sigma_mean_center_3D_top=interp3(s1xbins,s1ybins,s1zbins,Sigma_S1_3D_top,x_center,y_center,z_center,'spline');

%to Apply the correction use the following.
%s1_phe_top_xyz=s1_phe_top.*interp3(xx,yy,zz,norm_s1_top_xyz,s2x,s2y,drift_time,'spline');


%% Save The data
save(strcat('LUX_corrections/',file_id_cp,'/',file_id_cp,'_3D_IQs'),'norm_s1_both_xyz','norm_s1_bottom_xyz','norm_s1_top_xyz', 'sigma_norm_s1_both_xyz'...
,'sigma_norm_s1_bottom_xyz','sigma_norm_s1_top_xyz','s1xbins','s1ybins','s1zbins'...
,'Count_S1_3D', 'Count_S1_3D_bottom', 'Mean_S1_3D_both', 'Mean_S1_3D_bottom', 'Mean_S1_3D_top', 'Sigma_S1_3D_both', 'Sigma_S1_3D_bottom', 'Sigma_S1_3D_top'...
,'mean_center_3D_both','sigma_mean_center_3D_both','mean_center_3D_bottom','sigma_mean_center_3D_bottom','mean_center_3D_top','sigma_mean_center_3D_top','det_edge');



%% Submit IQ to LUG

if submit == 1
     
    user=user_name;
    source='Kr-83';
    
    file_id=file_id_cp(1:19);
    cp=file_id_cp(23:27);
    
 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %   Check if this is a LUX_SIM or Alternate Data Processing Chain.
        %   Need to modify algoythim version and event selction
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

       

        % Querying the CP uning the cp numbers from the IQ entries to get
        % the gs numbers and so see if the alternate Chain is being used
             query_str = ['select * from lug_complete_process_record where cp="' cp '" ;' ];
             data2 = MySQLQuery_UMD(query_str,'read'); %Modified to only look on master and not sanfordlab mirror -AD
                 if str2num(cp) == data2.cp 
                    gs = data2.gs;
                 end
                    %end get gs number


      %{
                  comment_str=data2.comments{1}; 

             if ( strfind(comment_str,'Alternative Chain')>0 )
                 if strcmp(file_id(1:5),'luxsm') % If SIM, label as SIM
                    algorithm_name='LUXkrypCal_SIM_AC';
                 else
                    algorithm_name='LUXkrypCal_AC'; % label as Alternative Chain
                 end 
                 %get the version number of the Alternative Chain
                 str_index=strfind(comment_str,'Chain v');
                 dp_version=sprintf('%4.3g',sscanf(comment_str(str_index:end),'Chain v%g')) ; %char format
             end  
   
          %}
                    
                    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    version=dp_version;
    clear lug_val Image_paths xmlinput
    
    %S1 XYZ Correction%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    lug_val(1) = {file_id};
    lug_val(2) = {cp}; %Proccessing version number
    lug_val(3) = {gs};
    lug_val(4) = {user};
    lug_val(5) = {algorithm_name};
    lug_val(6) = {version};
    lug_val(7) = {Kr_events};              %total count
    lug_val(8) = {s1xbins};                
    lug_val(9) = {s1ybins};                
    lug_val(10) = {s1zbins};
    lug_val(11) = {norm_s1_both_xyz};        %Normalization Matrix
    lug_val(12) = {sigma_norm_s1_both_xyz}; %1 sigma normalization matric
    lug_val(13) = {norm_s1_bottom_xyz};     %Normalization Matrix
    lug_val(14) = {sigma_norm_s1_bottom_xyz};%1 sigma of the normalization matrix
    lug_val(15) = {norm_s1_top_xyz};        %Normalization Matrix
    lug_val(16) = {sigma_norm_s1_top_xyz};  %1 sigma of the normalization matrix
    lug_val(17) = {mean_center_3D_both};        %Phe Mean matrix
    lug_val(18) = {sigma_mean_center_3D_both};
    lug_val(19) = {mean_center_3D_bottom};      
    lug_val(20) = {sigma_mean_center_3D_bottom};
    lug_val(21) = {mean_center_3D_top};      
    lug_val(22) = {sigma_mean_center_3D_top}; 
    lug_val(23) = {Count_S1_3D};    %Count in each bin
    lug_val(24) = {det_edge};
    lug_val(25) = {x_center};
    lug_val(26) = {y_center};
    lug_val(27) = {z_center};
    
    
     %Image_paths = strcat(file_id_cp,'/S1_3D_both_norm_',file_id_cp,'.jpg');
        %,strcat(file_id_cp,'/S1_3D_Count_',file_id_cp,'.jpg') };
    %strcat(file_id_cp,'/S1_3D_both_',file_id_cp,'.jpg')
     
            xmlinput.iq.global.filename_prefix = lug_val(1);
            xmlinput.iq.global.cp_number = lug_val(2);
            xmlinput.iq.global.gs_number = lug_val(3);
            xmlinput.iq.global.computed_by = lug_val(4);
            xmlinput.iq.global.computed_date = datestr(now,'yyyymmddTHHMM');
            xmlinput.iq.global.source = source;
            xmlinput.iq.global.source_id = 'unknown';
            xmlinput.iq.global.notes = ' ';
            xmlinput.iq.global.algorithm_name = lug_val(5);
            xmlinput.iq.global.algorithm_version = lug_val(6);
            xmlinput.iq.correction.fit.number_of_kr_events= lug_val(7);
            xmlinput.iq.correction.fit.x_bin_center = lug_val(8); 
            xmlinput.iq.correction.fit.y_bin_center = lug_val(9);
            xmlinput.iq.correction.fit.z_bin_center = lug_val(10);
            xmlinput.iq.correction.fit.norm_s1_both_xyz = lug_val(11);
            xmlinput.iq.correction.fit.sigma_norm_s1_both_xyz = lug_val(12);
            xmlinput.iq.correction.fit.norm_s1_bottom_xyz = lug_val(13);
            xmlinput.iq.correction.fit.sigma_norm_s1_bottom_xyz = lug_val(14);
            xmlinput.iq.correction.fit.norm_s1_top_xyz = lug_val(15);
            xmlinput.iq.correction.fit.sigma_norm_s1_top_xyz = lug_val(16);
            xmlinput.iq.correction.fit.mean_center_both_xyz = lug_val(17);
            xmlinput.iq.correction.fit.sigma_mean_center_both_xyz = lug_val(18);
            xmlinput.iq.correction.fit.mean_center_bottom_xyz = lug_val(19);
            xmlinput.iq.correction.fit.sigma_mean_center_bottom_xyz = lug_val(20);
            xmlinput.iq.correction.fit.mean_center_top_xyz = lug_val(21);
            xmlinput.iq.correction.fit.sigma_mean_center_top_xyz = lug_val(22);
            xmlinput.iq.correction.fit.count_s1_both_xyz = lug_val(23);
            xmlinput.iq.correction.fit.detector_edge = lug_val(24);
            xmlinput.iq.correction.fit.x_center = lug_val(25);
            xmlinput.iq.correction.fit.y_center = lug_val(26);
            xmlinput.iq.correction.fit.z_center = lug_val(27);
            
            LUXSubmitIQ(user,'s1_xyz_correction',algorithm_name,version,file_id,cp,xmlinput,'Kr83 Calibration','','');
   
   

    clear lug_val Image_paths xmlinput
       
    

end   %%end submit



fprintf('KrypCal XYZ Finished \n');