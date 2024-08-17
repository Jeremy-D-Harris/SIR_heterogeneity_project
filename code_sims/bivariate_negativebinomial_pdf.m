function [pdf] = bivariate_negativebinomial_pdf(X,Y,theta1,theta2,m1,m2,lambda)
% bivariate negative binomial pdf

c1 = ((1-theta1)/(1-theta1*exp(-1)))^(1/m1);
c2 = ((1-theta2)/(1-theta2*exp(-1)))^(1/m2);

dim = size(X);
m = dim(1);
n = dim(2);

pdf = zeros(dim);

for xx = 1:m
    for yy = 1:n
        

        term1 = gamma(m1^(-1)+X(xx,yy))/(gamma(X(xx,yy)+1)*gamma(m1^(-1)));
        term2 = gamma(m2^(-1)+Y(xx,yy))/(gamma(Y(xx,yy)+1)*gamma(m2^(-1)));
        
        pdf(xx,yy) = term1*(theta1^X(xx,yy))*(1-theta1)^(m1^(-1))*term2*(theta2^Y(xx,yy))*(1-theta2)^(m2^(-1))*(1 + lambda*(exp(-X(xx,yy)) - c1)*(exp(-Y(xx,yy))-c2));
    end 
end

% pdf = pdf';
end
