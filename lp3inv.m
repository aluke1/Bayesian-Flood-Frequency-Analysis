function [ Xp ] = lp3inv(p,gamma,sigma,mu)
z=norminv(p,0,1); 
%Compute frequency factor based on Wilson-Hilferty transformation
K = 2./gamma .* (1 + gamma .* z/6 - gamma.^2/36).^3 - 2./gamma;
Xp = mu + sigma.*K;
Xp = 10.^(Xp); 
end

