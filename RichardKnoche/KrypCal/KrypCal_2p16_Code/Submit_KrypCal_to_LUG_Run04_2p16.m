function [lug_val] = Submit_KrypCal_to_LUG(user_name, dp_version, algorithm_name, dir_path, use_ccv_flag)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Once a directory is populated with Kr83 data sets, loop through all and 
% submit KrypCal IQs to the LUG.
%
% Enter a string for user_name, dp_version, algorithm_name, dir_path
% (directory path), use_ccv_flag (bool)
% Example path = '/scratch2/adobi/DataTransfer/SR_V13/'
% Example dp_version = '1.3';
% Example algorithm_name ='LUXkrypCal'
% Enter 0 or 1 to use the 'ccv_flag' the file processing complete flag
%
% -AD 11/3/13
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
if nargin==4 
    use_ccv_flag=1; % by defult use the ccv_flag
end

if nargin<4 
    error('Invalid number of input, enter three strings \n');
end

if ~ischar(user_name) && ~ischar(dp_version) && ~ischar(dir_path)
    error('Invalid input, enter three strings and ccv_flag (0 or 1) \n');
end



%% Find IQs that have already been submitted without a strike and ignore those files

strg_start='select values_xml from control.lug_iqs where strikeme=0 AND algorithm_version = ';
strg_end=' AND algorithm = "';
group_end=' group by filename_prefix' ; %add this to prevent duplication.
qry_strg=sprintf('%s%s%s%s"%s;',strg_start,dp_version,strg_end,algorithm_name, group_end);

out = MySQLQuery_UMD(qry_strg);
LUG_file_name={' '};


    for j=1:size(out.values_xml,1); % 

    value1 = XMLParser(out.values_xml{j});

    LUG_file_name{j}=char(value1.iq.global.filename_prefix); % create an array of Char
            %Could just selct filename_prefix from the LUG as well
    end


%Get all Kr data filename prefixes
out_Kr = MySQLQuery_UMD('select * from control.lug_acquisitions where source = "Kr-83" ;');
Kr_file_name={' '};

    for j=1:length(out_Kr.filename_prefix) % 

    Kr_file_name{j}=char(out_Kr.filename_prefix(j));

    end



%check that the folder contains data and that the IQ has not already been
%uploaded for that data set (given the IQ version number)
%Also check that it is an actual Kr data set
potential_folders_to_load={' '};

list=dir(dir_path);

    k=1;
for ii=1:length(list)
    
    if length(list(ii).name) > 18
        
        if ( strcmp(list(ii).name(1:6),'lux10_') && sum( strcmp(list(ii).name(1:19),LUG_file_name))==0 ...
            && sum( strcmp(list(ii).name(1:19),Kr_file_name))==1) % A LUX folder & No IQ on LUG & a Kr set
            
            folder_contents=dir(strcat(dir_path,'/',list(ii).name));
            dir_size=sum([folder_contents.bytes]);
            
            if dir_size > 10^9 % 1 GB
                
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



no_folders_to_load=0;

    if ( strcmp(folders_to_load{1},' ') && length(folders_to_load) < 2 )
           no_folders_to_load=1;
    end
    
    
%% Get all 100 V/cm filename prefixes

out_field = MySQLQuery_UMD('select * from control.lug_acquisitions where collection like "%100 V/cm%" ;');
    
    Field_100_names={' '};

    for j=1:length(out_field.filename_prefix) % 

        Field_100_names{j}=char(out_field.filename_prefix(j));

    end
    
% Is this a 100 V/cm data set? (will need to add more fields for Run04).

for ii=1:length(folders_to_load)    
        
        if ( sum( strcmp(folders_to_load{ii}(1:19),Field_100_names))==0 )
                Field_100_flag(ii)=0;
        else
                Field_100_flag(ii)=1;
        end
            
    
end
    
%%    Run KrypCal on each folder one by one. Will need to edit the path


   
rqs_to_load = {'pulse_area_phe','event_timestamp_samples'...
   ,'pulse_classification' ...
   ,'z_drift_samples' , 's1s2_pairing'...
   ,'top_bottom_ratio','x_cm','y_cm'...
   ,'full_evt_area_phe',...
   'event_number','chi2','prompt_fraction','aft_t1_samples','pulse_start_samples',...
   'pulse_end_samples','top_bottom_asymmetry','aft_t0_samples','aft_t2_samples',...
   'full_evt_area_phe','admin','Kr83fit_s1a_area_phe','Kr83fit_s1b_area_phe',...
   'hft_t50r_samples','hft_t50l_samples','Kr83fit_dt_samples'};
    
% removed s1_before_s2 and s1s2_pairing


if length(folders_to_load)>7
    for   ii=1:7 % after 13 submitions Matlab will crash due to java memory leak
            path=strcat(dir_path,'/',folders_to_load{ii});         
                
            files_available=dir(path); %only load if file complete flag is in the directory    
            if ( ismember({'ccv_rqs_done'},{files_available.name}) || use_ccv_flag==0 )
            
                d = LUXLoadMultipleRQMs_framework(path,rqs_to_load);  
                
                temp_name=folders_to_load{ii};
                if datenum(str2num(temp_name(7:10)),str2num(temp_name(11:12)),str2num(temp_name(13:14)),str2num(temp_name(16:17)),str2num(temp_name(18:19)),00) >= datenum(2014,09,03,19,18,00);
                LUXkrypCal_Run04_2p16(d,user_name,dp_version,algorithm_name,1); %rq,use pulse_class,submit IQ to LUG,delay time min     
                else
                LUXkrypCal_Run03(d,user_name,dp_version,algorithm_name,1); %rq,use pulse_class,submit IQ to LUG,delay time min     
                end

            end
            clear files_available
    end
        
elseif no_folders_to_load == 0
    for  ii=1:length(folders_to_load) % after 13 submitions Matlab will crash due to java memory leak
           path=strcat(dir_path,'/',folders_to_load{ii});         
        
            files_available=dir(path); %only load if file complete flag is in the directory    
            if ( ismember({'ccv_rqs_done'},{files_available.name}) || use_ccv_flag==0 )
            
                d = LUXLoadMultipleRQMs_framework(path,rqs_to_load);  

                temp_name=folders_to_load{ii};
                if datenum(str2num(temp_name(7:10)),str2num(temp_name(11:12)),str2num(temp_name(13:14)),str2num(temp_name(16:17)),str2num(temp_name(18:19)),00) >= datenum(2014,09,03,19,18,00);
                LUXkrypCal_Run04_2p16(d,user_name,dp_version,algorithm_name,1); %rq,use pulse_class,submit IQ to LUG,delay time min     
                else
                LUXkrypCal_Run03(d,user_name,dp_version,algorithm_name,1); %rq,use pulse_class,submit IQ to LUG,delay time min     
                end

            end
            clear files_available     

    end
    
end

lug_val=1;    

fprintf('Finished \n');

exit; % quit matlab

