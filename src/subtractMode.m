function newlogdiff = subtractMode(logdiff) 
    num_stocks = size(logdiff, 2);
    market_mode = mean(logdiff,2);
    newlogdiff= zeros(size(logdiff));
     
    for stock=1:num_stocks
        p = polyfit(market_mode,logdiff(:,stock),1);
        newlogdiff(:,stock)=logdiff(:,stock)-polyval(p,market_mode);
    end
end
 
