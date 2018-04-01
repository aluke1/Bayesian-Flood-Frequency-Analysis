function [ post,log_like ] = post_pdf(params,data,cf,Mj,prior)

X=log10(data(:,2).*cf);           %Data record in log10(m^3/s)
t = data(:,1) - min(data(:,1));   %t in years from first observation
vector_ones = ones(length(t),1);  %vector of ones same length as record

%create vectors of parameters
if Mj == 1; 
mu=    vector_ones .* params(3);
sigma= vector_ones .* params(2);
gamma= vector_ones .* params(1);
end
if Mj == 2; 
mu=    vector_ones .* params(3) + (t.*params(4)) ;
sigma= vector_ones .* params(2);
gamma= vector_ones .* params(1);
end

log_like = sum(log(lp3pdf(X,gamma,sigma,mu))); %likelihood of parameters
log_prior = log(prior(params));                %prior density of parameters
log_post = log_like + log_prior;               %unnormalized log-post denisty
post = exp(log_post);                          %unnormalized posterior density

end
