function y = myTry(e,ne,liste,listne)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

for i=1:5
    if ismember(i, liste)
        idx = find(liste==i,1);
        if i==1
            a = e(idx);
        elseif i==2
            b = e(idx);
        elseif i==3
            c = e(idx);
        elseif i==4
            d = e(idx);
        else
            f = e(idx);
        end
    else
        idx = find(listne==i,1);
        if i==1
            a = ne(idx);
        elseif i==2
            b = ne(idx);
        elseif i==3
            c = ne(idx);
        elseif i==4
            d = ne(idx);
        else
            f = ne(idx);
        end       
    end
end
y1 = c+d;
y2 = .5*a^2;
y3 = .5*b^2;
y4 = d*a;
y5 = f*a^-1;
y = [y1;y2;y3;y4;y5]
end

