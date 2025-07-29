function ta2 = acf_lag1(y)
% ACF_K - Autocorrelation at Lag 1
% acf(y,k)
%
% Inputs:
% y - series to compute acf for

k = 1;
 
ybar = mean(y);
N = max(size(y)) ;

cross_sum = zeros(N-k,1) ;

% Numerator, unscaled covariance
for i = (k+1):N
    cross_sum(i) = (y(i)-ybar)*(y(i-k)-ybar) ;
end

% Denominator, unscaled variance
yvar = (y-ybar)'*(y-ybar) ;

ta2 = sum(cross_sum) / yvar ;