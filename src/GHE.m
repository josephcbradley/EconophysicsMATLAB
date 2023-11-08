function [Hq,SE_Hq,A,B,SE_A,SE_B] = GHE(price,q,tauMIN,tauMAX,plot_choice)
%GHE computes the multiscaling proxies for a financial time series with their standard errors
%scaling region chosen by the user
%   Input variables:
%                   prices:     price time series;
%                   q:           row vector of powers;
%                   tauMIN:      minumum value of aggregation;
%                   tauMAX:      maximum value of aggreation;
%                   plot_choice: set to 1 if you want the scaling shown long with the fit;
%   Output variables:
%                   Hq = Hurst exponent vector H(q) for each moment obtained from linear regression in log-log scale;
%                   SE_Hq = Hurst exponent vector H(q) standard errors;
%                   A: it is 0.5 for uniscaling time series;
%                   B: it is less than zero for multiscaling time series;

price = price(~isnan(price));

if isempty(price) == 1
    fprintf('\nAll values are NaN or the time series is empty\n')
    return;
end

model=fittype('B*x^2+A*x');
%model=fittype('B*x+A');

[Hq,SE_Hq]=genhurst_fast(price,q,tauMIN,tauMAX,plot_choice); 

% SE obtained from the q vs Hq fit
Var = (q.^2)'.*(SE_Hq.^2);
%Var = (SE_Hq.^2);
Weights = sqrt(1./Var); %weights as the inverse standard error
[f]=fit(q',q'.*(Hq),model,'lower',[-inf, -inf ],'upper',[+inf, +inf],'StartPoint',[0 0],'Weight',Weights);
%[f]=fit(q',Hq,model,'lower',[-inf, -inf ],'upper',[+inf, +inf],'StartPoint',[0 0],'Weight',Weights);
B=f.B;
A=f.A;
CI = confint(f);
NN = length(q')-2;
p = 0.975;
xq = tinv(p,NN);
SE_B = (CI(2,1)-CI(1,1))./(2*xq);
SE_A = (CI(2,2)-CI(1,2))./(2*xq);

end
