function [ pdf ] = lp3pdf(X, gamma, sigma, mu)
%convert skewness, standard deviation, and mean to parameters of PIII distribution
alpha=4./gamma.^2;
beta=(1/2).*sigma.*abs(gamma);
tau=mu-2.*sigma./gamma;
%compute density in log-space
if gamma(1) > 0 
    log_pdf = (alpha - 1) .* log((X - tau )) - (X -tau)./ beta - alpha .* log(beta) - gammaln(alpha);
    for i = 1:length(log_pdf) 
        if isreal(log_pdf(i))
            pdf(i) = exp(log_pdf(i)); 
        else %if imaginary number returned, X(i) outside bounds of distribution, set denisty to 0
            pdf(i) = 0;
        end
    end
elseif gamma(1) < 0 
    log_pdf = (alpha - 1) .* log((tau - X)) - (tau - X)./ beta - alpha .* log(beta) - gammaln(alpha);
    for i = 1:length(log_pdf) 
        if isreal(log_pdf(i)) 
            pdf(i) = exp(log_pdf(i)); 
        else %if imaginary number returned, X(i) outside bounds of distribution, set denisty to 0
            pdf(i) = 0;
        end
    end
end

end

