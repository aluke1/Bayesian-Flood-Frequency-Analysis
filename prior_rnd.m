function [draw] = prior_rnd(r,d,mp,mpr)
draw = nan(r,d); 
%initialize by randomly drawing from the prior
 for i = 1:d
  draw(:,i) = eval([mp{1,i} 'rnd(' num2str(mpr(1,i)) ',' num2str(mpr(2,i)) ',' num2str(r) ',1)']);
 end    
end

