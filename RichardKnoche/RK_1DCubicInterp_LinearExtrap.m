function [ first_interp ] = Untitled( X, Y, X0)
%Does cubic interpolation within the bounds, and linear when extrapolation fails
%example: y=s2_both_norm_z_means; Y=y; X0=[0:1:350];

if length(X0)>1
    first_interp= interp1(X,Y,X0,'v5cubic'); %normalize to right below the gate

    %find the working value of X0 closest to the NaN value of X0
    badX0=X0(isnan(first_interp));
    goodX0=X0(~isnan(first_interp));
    goodY0=first_interp(~isnan(first_interp));

    %sort to find the first and last 10% cubic values
    [sorted_goodX0, sorted_index]=sort(goodX0);
    sorted_goodY0=goodY0(sorted_index);

    first20_fit=polyfit(sorted_goodX0(1:floor(length(sorted_goodX0)*0.1)),sorted_goodY0(1:floor(length(sorted_goodX0)*0.1)),1);
    last20_fit=polyfit(sorted_goodX0(ceil(length(sorted_goodX0)*0.9):end),sorted_goodY0(ceil(length(sorted_goodX0)*0.9):end),1);

    %Doing linear extrapolation for events on the edges
    first_interp(isnan(first_interp) & X0<=sorted_goodX0(1))=polyval(first20_fit,X0(isnan(first_interp) & X0<=sorted_goodX0(1)));
    first_interp(isnan(first_interp) & X0>=sorted_goodX0(end))=polyval(last20_fit,X0(isnan(first_interp) & X0>=sorted_goodX0(end)));

    %clean up for anything that is NaN in the middle of the interpolation
    %find the working value of X0 closest to the NaN value of X0
    badX0=X0(isnan(first_interp));
    if length(badX0)>=1;
     goodX0=X0(~isnan(first_interp));
     goodY0=first_interp(~isnan(first_interp));

        for i=1:length(badX0);
            temp=abs(badX0(i)-goodX0);
            [minx idx]=min(temp);
            FixedX0(i)=goodX0(idx);
            FixedX0Index(i)=idx;
        end

     first_interp(isnan(first_interp))=goodY0(FixedX0Index);
    end

else
  
    first_interp= interp1(X,Y,X0,'v5cubic'); %normalize to right below the gate

    if isnan(first_interp)
        
      first_interp= interp1(X,Y,X,'v5cubic'); %normalize to right below the gate

      %find the working value of X0 closest to the NaN value of X0
      badX=X(isnan(first_interp));
      goodX=X(~isnan(first_interp));
      goodY=first_interp(~isnan(first_interp));

      %sort to find the first and last 20% cubic values
      [sorted_goodX, sorted_index]=sort(goodX);
      sorted_goodY=goodY(sorted_index);

      first20_fit=polyfit(sorted_goodX(1:floor(length(sorted_goodX)*0.1)),sorted_goodY(1:floor(length(sorted_goodX)*0.1)),1);
      last20_fit=polyfit(sorted_goodX(ceil(length(sorted_goodX)*0.9):end),sorted_goodY(ceil(length(sorted_goodX)*0.9):end),1);

      if X0<=sorted_goodX(1)
       first_interp=polyval(first20_fit,X0);
      elseif X0>=sorted_goodX(end)
       first_interp=polyval(last20_fit,X0);          
      else     
       first_interp=interp1(X,Y,X0,'nearest');
      end
    end
end

