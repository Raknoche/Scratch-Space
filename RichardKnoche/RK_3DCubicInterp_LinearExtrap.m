function [ first_interp ] = Untitled( X, Y, Z, V, X0, Y0, Z0)
%Does cubic interpolation within the bounds, and linear when extrapolation fails
%example: y=s2_both_norm_z_means; Y=y; X0=[0:1:350];

    %Grow the V matrix to accomidate X0,Y0,Z0 bounds
    X0_min=floor(min(X0)); X0_max=ceil(max(X0));
    Y0_min=floor(min(Y0)); Y0_max=ceil(max(Y0));
    Z0_min=floor(min(Z0)); Z0_max=ceil(max(Z0));
    
    %In Z slices, grow the X values
    new_V=zeros(size(V,1)+2,size(V,2)+2,size(V,3)+2);
    new_V(2:end-1,2:end-1,2:end-1)=V;

    %fill in the x boundary values
    for k=2:size(new_V,3)-1;
        for j=2:size(new_V,2)-1;
            lefttempfit=polyfit(Y(1:3).',new_V(2:4,j,k),1);
            righttempfit=polyfit(Y(end-3:end).',new_V(end-4:end-1,j,k),1);
            new_V(1,j,k)=polyval(lefttempfit,Y0_min);
            new_V(end,j,k)=polyval(righttempfit,Y0_max);
        end
    end
    
    %fill in the y boundary values
    for k=2:size(new_V,3)-1;
        for i=2:size(new_V,1)-1;
            lefttempfit=polyfit(X(1:3),new_V(i,2:4,k),1);
            righttempfit=polyfit(X(end-3:end),new_V(i,end-4:end-1,k),1);
            new_V(i,1,k)=polyval(lefttempfit,X0_min);
            new_V(i,end,k)=polyval(righttempfit,X0_max);
        end
    end    
    
        %fill in the y boundary values
    for i=2:size(new_V,1)-1;
        for j=2:size(new_V,2)-1;
            lefttempfit=polyfit(Z(1:3).',squeeze(new_V(i,j,2:4)),1);
            righttempfit=polyfit(Z(end-3:end).',squeeze(new_V(i,j,end-4:end-1)),1);
            new_V(i,j,1)=polyval(lefttempfit,Z0_min);
            new_V(i,j,end)=polyval(righttempfit,Z0_max);
        end
    end    
    
    new_X=zeros(length(X)+2,1);
    new_Y=zeros(length(Y)+2,1);
    new_Z=zeros(length(Z)+2,1);
    
    new_X(2:end-1)=X; new_Y(2:end-1)=Y; new_Z(2:end-1)=Z;
    new_X(1)=X0_min; new_Y(1)=Y0_min; new_Z(1)=Z0_min;
    new_X(end)=X0_max; new_Y(end)=Y0_max; new_Z(end)=Z0_max;
    
    first_interp= interp3(new_X,new_Y,new_Z,new_V,X0,Y0,Z0,'cubic'); %normalize to right below the gate
end

