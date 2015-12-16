function [ minx, idx, miny] = FindNearestIndex( X,Xo, cut, Y )
%Written for apply corrections code
%Takes vector X and value Xo, a cut, and a second vector Y
%Finds value of X closest to Xo
%Returns the value of X, and the corresponding index, and the value of Y for that index (after the cut)
    x= X(cut);
    y= Y(cut);
    temp=abs(Xo-x);
    [minx idx]=min(temp);
    miny=y(idx);
    
end

