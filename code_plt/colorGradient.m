function [colors] = colorGradient(col1,col2,N)
%interpolate between two colors col1 and col2 and output N colors. All in rgb
%with values between 0  and 1.

%check N valid
mustBeInteger(N) %N must be integer valued
mustBeGreaterThanOrEqual(N,2) %N must produce output of at least two entries  - at minimum col1 and col2 with no interpolation
%check col1 and col2 are valid
validatecolor(col1)
validatecolor(col2)


%calculate output colormap
gradient = (col2 - col1)/(N-1);
colors = zeros([N 3]);

%scale by N
for aa = 1:N
    colors(aa,:) = col1 + (aa-1)*gradient;
end

end