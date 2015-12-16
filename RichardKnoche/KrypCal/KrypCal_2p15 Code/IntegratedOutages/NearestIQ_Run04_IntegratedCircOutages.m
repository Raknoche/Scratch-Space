function [index iq_time]  = NearestIQ_Run04(current_data_set_time, dataset_times, noflow_times)

   %Index iq times
  [a b] = size(dataset_times);
  for i=1:b
    id(i) = i;
  end
  
  %Index circ outage times
  a  = length(noflow_times);
  for i=1:a
    outage_id(i) = i;
  end
 
  %define first and last IQ
  [first_iq_time first_iq_index] = min(dataset_times);
  [last_iq_time last_iq_index] = max(dataset_times);
  
  %define first and last Outage
  [first_outage_time first_outage_index] = min(noflow_times);
  [last_outage_time last_outage_index] = max(noflow_times);  

   no_left_outage=0;
   no_right_outage=0;


   %Find the IQ and outage timestamps that are closest to the current timestamp

  if current_data_set_time > first_iq_time && current_data_set_time < last_iq_time
        [lmin lind left_iq_index]=FindNearestIndex(dataset_times,current_data_set_time,(dataset_times<=current_data_set_time),id);
        [rmin rind right_iq_index]=FindNearestIndex(dataset_times,current_data_set_time,(dataset_times>=current_data_set_time),id);
  
    %Check that there is an outage to the left and right
    if first_outage_time <= current_data_set_time && last_outage_time >= current_data_set_time
        [lmin lind left_outage_index]=FindNearestIndex(noflow_times,current_data_set_time,(noflow_times<=current_data_set_time),outage_id);
        [rmin rind right_outage_index]=FindNearestIndex(noflow_times,current_data_set_time,(noflow_times>=current_data_set_time),outage_id);

    %if there is only an outage to the right
    elseif first_outage_time > current_data_set_time && last_outage_time >= current_data_set_time
        [rmin rind right_outage_index]=FindNearestIndex(noflow_times,current_data_set_time,(noflow_times>=current_data_set_time),outage_id);
        left_outage_index = 1; %Just set it to 1 so the code doesn't break
        no_left_outage=1;
    %if there is only an outage to the left
    elseif  first_outage_time < current_data_set_time && last_outage_time <= current_data_set_time
        [lmin lind left_outage_index]=FindNearestIndex(noflow_times,current_data_set_time,(noflow_times<=current_data_set_time),outage_id);
        right_outage_index = 1; %Just set it to 1 so the code doesn't break
        no_right_outage=1;
    %If the is no outage to the left or the right
    else
       left_outage_index=1;
       right_outage_index=1; %Just set it to 1 so the code doesn't break
       no_left_outage=1;
       no_right_outage=1;
    end
   

     %If there are no outages between the current timestamp and the two nearest IQs, choose the nearest IQ
   if ( no_left_outage==0 && no_right_outage==0 && noflow_times(left_outage_index) <= dataset_times(left_iq_index) && noflow_times(right_outage_index) >= dataset_times(right_iq_index))...
       || ( noflow_times(left_outage_index) <= dataset_times(left_iq_index) && no_left_outage==0 && no_right_outage==1) ...
       || ( no_left_outage==1 && no_right_outage==0 && noflow_times(right_outage_index) >= dataset_times(right_iq_index)) ...
       || ( no_left_outage==1 && no_right_outage==1)

        index = interp1(dataset_times,id,current_data_set_time,'nearest');
        
    %If there is an outage between the current timestamp and the IQ to the right, use the IQ to the left
    elseif ( no_left_outage==0 && no_right_outage==0 && noflow_times(left_outage_index) < dataset_times(left_iq_index) && noflow_times(right_outage_index) < dataset_times(right_iq_index))...
            || (no_left_outage==1 && no_right_outage==0  && noflow_times(right_outage_index) < dataset_times(right_iq_index))
        index = left_iq_index;

                
    %If there is an outage between the current timestamp and the IQ to the left, use the IQ to the right
    elseif ( no_left_outage==0 && no_right_outage==0 && noflow_times(left_outage_index) > dataset_times(left_iq_index) && noflow_times(right_outage_index) > dataset_times(right_iq_index))...
            || (no_left_outage==0 && no_right_outage==1 && noflow_times(left_outage_index) > dataset_times(left_iq_index))

        index = right_iq_index;
        
    %If there is an outage between the current timestamp and every IQ, return an error
    else
        fprintf('Found no acceptable IQs.....problem');
    end
  elseif current_data_set_time <= first_iq_time

      index = first_iq_index;
      
  elseif current_data_set_time >= last_iq_time

      index= last_iq_index;
      
  else
    fprintf('Found no acceptable IQs.....problem');
  end
  
    iq_time = dataset_times(index);

end