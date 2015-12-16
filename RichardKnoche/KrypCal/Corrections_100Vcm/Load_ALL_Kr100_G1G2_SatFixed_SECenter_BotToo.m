dir_path='/scratch4/LUXData/G1G2_Data/Kr_100Vcm/';

potential_folders_to_load={' '};
folders_to_load={' '};

% Check that we are only loading folders with lux10_
% and that they contain at least 10 MB of data
list=dir(dir_path);

    k=1;
for ii=1:length(list)
    
    if length(list(ii).name) > 6
        
        if  strcmp(list(ii).name(1:6),'lux10_') 
            
            folder_contents=dir(strcat(dir_path,'/',list(ii).name));
            dir_size=sum([folder_contents.bytes]);
            
            if dir_size > 10^7 % 10 MB
             
                potential_folders_to_load{k}=list(ii).name;
                k=k+1;
            
            end
                
                
        end
    end
    
end


% Once we determined good folders to load, check their multiplicity
% Keep higher CP unless duplicated folder with a higher CP number contains no data
folder_names_char=char(potential_folders_to_load);
folder_name_no_cp=cellstr([folder_names_char(:,1:19)]); %convert the char back to a cell array of chars

k=1;
good_index=1;
while good_index<=length(potential_folders_to_load)
        
        
       folder_multiplicity=sum(strcmp(potential_folders_to_load{good_index}(1:19),folder_name_no_cp ));
       good_index=good_index+folder_multiplicity-1;
            
       folders_to_load{k}=potential_folders_to_load{good_index};
       
       k=k+1;     
       good_index=good_index+1;
    
end


rqs_to_load = {'pulse_area_phe','event_timestamp_samples'...
   ,'pulse_classification' ...
   ,'z_drift_samples' , 's1s2_pairing','golden'...
   ,'z_corrected_pulse_area_all_phe'...
   ,'x_corrected','y_corrected','selected_s1_s2','top_bottom_ratio',...
   'correction_electron_lifetime','correction_s2_xy_dependence','correction_s1_xyz_dependence',...
   'correction_s2_xy_dependence_bot','correction_s1_xyz_dependence_bot','aft_t0_samples',...
   'correction_s2_xy_dependence','spike_count','xyz_corrected_pulse_area_all_phe','xyz_corrected_pulse_area_bot_phe','s2_rec','x_cm','y_cm'};
    
ii=1;

    path=strcat(dir_path,'/',folders_to_load{ii});
    
    d = LUXLoadMultipleRQMs_framework(path,rqs_to_load);
    d = se_wrt_first_s1_func(d);

%% VUV and Spike Corrections
    % These correction factors are put in here by hand; they are the ratio of
% the two gain versions, per channel.
gainCorr = [ 1.07852885,  1.02147302,  1.14704873,  1.06599796,  1.,1.06565116,  1.17260696,  1.23915269,  1.05237586,  1.04581427,1.19078187,  1.01920447,  1.10806574,  1.04368312,  1.09359576,        1.20713643,  1.01269677,  0.97496447,  1.15328135,  1.04632389,1.09657976,  1.0540953 ,  1.06966811,  1.18564767,  1.03347816,1.12664979,  1.21384238,  1.04435383,  1.096459  ,  1.07002231,1.13197335,  1.        ,  0.99399306,  1.01098228,  1.2029547 ,1.07295668,  1.0696273 ,  1.00686858,  1.08148498,  1.16649075,1.04850371,  1.05249835,  1.14049429,  1.00785073,  1.0943977 ,1.02962349,  1.08931321,  1.18235792,  1.03934669,  1.00709524,1.14178315,  1.04650374,  1.07221255,  1.03081791,  1.1338022 ,1.23879521,  0.9849912 ,  1.02989339,  1.20126104,  1.0419401 ,1.08283583,  1.05111627,  1.16620644,  1.19375442,  1.06260808,1.01510865,  1.13129887,  1.09737434,  1.05748364,  1.08007699,1.07439593,  1.2135306 ,  1.05225111,  1.06547247,  1.13095721,1.07715678,  1.0901655 ,  1.06751444,  1.12705449,  1.17057193,1.09581789,  1.09417162,  1.23750879,  1.0706705 ,  1.18163888,1.10125419,  1.14468798,  1.19880273,  1.12981092,  1.06592446,1.19612131,  1.14547149,  1.        ,  1.06916959,  1.1169941 ,1.15057597,  1.08348327,  1.08020233,  1.13099728,  1.08311348,1.11291194,  1.04754566,  1.10311997,  1.18206575,  1.0658937 ,1.03623549,  1.17802513,  1.07419023,  1.04725227,  1.08861405,1.10805692,  1.19684977,  1.07895852,  1.05339989,  1.15828588,1.08408482,  1.02837954,  1.05048981,  1.09229397,  1.14161773,1.01265111,  1.06664276];

    
n_sig=2;

rcut_min = 0;
rcut_max = 25;%cm
        s1area_bound_min = 0;
        s1area_bound_max = 10^5;

        s2area_bound_min = 50;
        s2area_bound_max = 2*10^6;
   
   
d.z_drift_samples(isnan(d.z_drift_samples)) = 0.0; % get rid of NaN
    d.phe_bottom=d.pulse_area_phe./(1+d.top_bottom_ratio); %bottom PMT pulse area

    events=sort(size(d.pulse_area_phe)); %The first element is the number of pulses. sometimes it's 5, sometimes it's 10

    s1_area_cut= inrange(d.pulse_area_phe,[s1area_bound_min,s1area_bound_max]);
    s2_area_cut= inrange(d.pulse_area_phe,[s2area_bound_min,s2area_bound_max]);

    s1_class=(d.pulse_classification==1 )& s1_area_cut ; %area cut for Kr events
    s2_class=(d.pulse_classification==2) & s2_area_cut ;

    s1_single_cut =logical( (s1_class & d.golden & d.selected_s1_s2).*repmat(sum(s2_class & d.golden & d.selected_s1_s2)==1,events(1),1) ); % was s1_before_s2_cut before using golden
    s2_single_cut =logical( (s2_class & d.golden & d.selected_s1_s2).*repmat(sum(s1_class & d.golden & d.selected_s1_s2)==1,events(1),1) );   

    s1=d.pulse_area_phe(s1_single_cut);
    s2=d.pulse_area_phe(s2_single_cut);
    
    s2_bot=d.phe_bottom(s2_single_cut);
    s2_both=d.pulse_area_phe(s2_single_cut);
    s1_bot=d.phe_bottom(s1_single_cut);
    s1_both=d.pulse_area_phe(s1_single_cut);
    
    s1_phe_both_xyz_VUV=mean(gainCorr).*s1_both;
    s2_phe_both_xyz_VUV=mean(gainCorr).*s2_both;
    s2_phe_bot_xyz_VUV=mean(gainCorr).*s2_bot;
    s1_phe_bot_xyz_VUV=mean(gainCorr).*s1_bot;

    s2x=d.x_cm(s2_single_cut);
    s2y=d.y_cm(s2_single_cut);
    drift_time = d.z_drift_samples(s2_single_cut)/100; %units of us
             
    %% No tracking of XeAct means, sigma, errors, and number of events since we don't do the cuts in this code

      save('/scratch4/LUXData/G1G2_Data/Kr_100Vcm/Run03Kr_100Vcm_BotAndBoth','s1_phe_both_xyz_VUV',...
      's2_phe_both_xyz_VUV','s2x','s2y','drift_time','s1_phe_bot_xyz_VUV','s2_phe_bot_xyz_VUV');        
