function [ first_interp ] = Untitled( X, Y, X0)
%Does cubic interpolation within the bounds, and nearest neighbor when extrapolation fails
%example: y=s2_both_norm_z_means; Y=y; X0=[0:1:350];
first_interp= interp1(X,Y,X0,'v5cubic'); %normalize to right below the gate

%find the working value of X0 closest to the NaN value of X0
badX0=X0(isnan(first_interp));
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

