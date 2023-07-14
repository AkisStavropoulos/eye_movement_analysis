function [x,y,b,c] = find_regress(x,y)

if size(x,1) < size(x,2)
    x = x';
end
if size(y,1) < size(y,2)
    y = y';
end

X = [x ones(length(x),1)];
[b]=regress(y,X);
c = b(2);
b = b(1);
