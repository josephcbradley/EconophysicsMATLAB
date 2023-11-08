function [Hq,SE] = genhurst_fast(price,q,tauMIN,tauMAX,plot_choice)
%GENHURST_FAST This function computes Hurst exponent(s) of a single time-series in a certain
%scaling region chosen by the user
%   Input variables:
%                   price:      price time series;
%                   q:           row vector of powers;
%                   tauMIN:      minumum value of aggregation
%                   tauMAX:      maximum value of aggreation
%                   plot_choice: set to 1 if you want the scaling shown along
%                                with the fit.
%   Output variables:
%                   Hq:          column vector of Hurst exponents of same
%                                length as q.

X=log(price);
lX=length(X);
lq=length(q);
q_mom=zeros(tauMAX-tauMIN+1,lq);

flag=1;
for j=tauMIN:tauMAX
    pricetemp=X(j:j:j*(floor(lX/j))); %time series windows of Length = floor(lX/tau) and Time step =  tau
    %     pricetemp=X(1:j:j*(floor(lX/j))); GB

    pricetemp=detrend(pricetemp); %Y = detrend(X) removes the best straight-line fit linear trend from the data in vector X and returns the residual in vector Y, it makes data stationary
    r=diff(pricetemp);
    q_mom(flag,:)=mean(bsxfun(@power,abs(r),q)); % E[|r|^q] element wise
    k(flag,:)=kurtosis(r);
    flag=flag+1;
end

x=log(tauMIN:1:tauMAX)';
% x=log(1./(tauMIN:1:tauMAX))'; %GB
%  x=((tauMIN:1:tauMAX)./(tauMIN:1:tauMAX))'; %GB
y_emp=log(q_mom);
coefs=zeros(lq,2);
%tic
for i=1:lq
    coefs(i,:)=polyfit(x,y_emp(:,i),1);
    coeff_stat{i}=regstats(y_emp(:,i),x); %addded by GB %By default, regstats uses a linear additive model with a constant term
    %coeff2(i,:)=regress(y_emp(:,i),x); %added by GB
    SE(i,:)=coeff_stat{i}.tstat.se; %standard error on the obtained H(q); 
    pvalues(i,:)=coeff_stat{i}.tstat.pval; %added by GB %p-value for the t-statistic of the two-sided hypothesis test
end
%toc
Hq=coefs(:,1)./q';
SE=SE(:,2);
if(plot_choice==1)
    figure
    hold on
    plot(x,y_emp,'ro')
    for i=1:lq
        plot(x,polyval(coefs(i,:),x),'b-')
        if (i == 1 || i == floor(lq/2) || i == lq)
            text(x(end),y_emp(end,i),["  \leftarrow  q = "+num2str(q(i))])
        end
    end
    xlabel("log(\tau)")
    ylabel("log(K_q(\tau))")
    set(gcf, 'Color', 'w');
    set(findall(gcf,'type','axes'),'fontSize',11)
end

