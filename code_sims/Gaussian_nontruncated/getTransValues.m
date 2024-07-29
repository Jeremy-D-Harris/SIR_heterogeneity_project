%get transmissibility values

function [d] = getTransValues(n,ifVar,end_val)
if ifVar
    d = linspace(0,end_val,n);
else 
    d = ones(1,n);
end
