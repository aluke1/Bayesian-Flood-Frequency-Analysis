function [ density ] = prior_pdf(params,d,mp,mpr)
marg_dens = nan(1,d); 
%compute the prior density of a parameter combination
for i = 1:d 
    marg_dens(1,i) = eval([mp{1,i} 'pdf(' num2str(params(1,i)) ',' num2str(mpr(1,i)) ',' num2str(mpr(2,i)) ')']);
end
density = prod(marg_dens);
end

