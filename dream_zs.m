function [x,p_x,Z] = dream_zs(prior,pdf,N,T,d,marg_prior,marg_par)
% DiffeRential Evolution Adaptive Metropolis with sampling from past archive
% and snooker update: DREAM_(ZS) algorithm
% Please study DREAM manual published in Environmental Modeling and Software, 2015:
% 
% Code Written by Jasper A. Vrugt, Oct. 2015
% Full code of DREAM_ZS with all its options available from second author: jasper@uci.edu

[delta,c,c_star,n_CR,p_g] = deal(1,0.05,1e-12,3,0.2);              % Default of algorithmic parameters 
k = 10; p_s = 0.1; m0 = max(N,20 * d); n_d = 3*delta;              % DREAM_ZS algorithmic variables
x = nan(T,d,N); p_x = nan(T,N);                                    % Preallocate chains and density
[J,n_id] = deal(zeros(1,n_CR));                                    % Variables selection prob. crossover
CR = [1:n_CR]/n_CR; p_CR = ones(1,n_CR)/n_CR;                      % Crossover values and select. prob.

Z = prior(m0,d);                                                   % Create initial archive
for i = m0-N+1:m0,  Z(i,d+1) = pdf(Z(i,1:d)); end                  % Compute density initial population
X = Z(m0-N+1:m0,1:d); p_X = Z(m0-N+1:m0,d+1); m = m0;              % Define X population and prepare Z
x(1,1:d,1:N) = reshape(X',1,d,N); p_x(1,1:N) = p_X';               % Store initial states and density
a = 1;
for t = 2:T,    % Dynamic part: Evolution of N chains
    [~,draw] = sort(rand(N-1,N));                                     % Permute [1,...,N-1] N times
    dX = zeros(N,d);                                                  % Set N jump vectors to zero
    lambda = unifrnd(-c,c,N,1);                                       % Draw N lambda values
    std_X = std(X);                                                   % Compute std each dimension
    R = randsample(1:m,N*n_d,'false'); R = reshape(R,n_d,N);          % Sample N*n_d integer values from [1..m]
    method = randsample({'parallel' 'snooker'},1,'true',[1-p_s p_s]); % Mix of parallel and snooker updates
    for i = 1:N,                                                      % Create proposals and accept/reject
        D = randsample([1:delta],1,'true');                           % Select delta (equal select. prob.)
        a = R(1:D,i); b = R(D+1:2*D,i); c_sn = R(2*D+1:3*D,i);        % Define a and b (parallel) + c (snooker)
        if strcmp(method,'parallel'),                                 % PARALLEL DIRECTION update
            id(i) = randsample(1:n_CR,1,'true',p_CR);                 % Select index of crossover value
            z = rand(1,d);                                            % Draw d values from U[0,1]
            A = find(z < CR(id(i)));                                  % Derive subset A selected dimensions
            d_star = numel(A);                                        % How many dimensions sampled?
            if d_star == 0, [~,A] = min(z); d_star = 1; end           % A must contain at least one value
            gamma_d = 2.38/sqrt(2*D*d_star);                          % Calculate jump rate
            g = randsample([gamma_d 1],1,'true',[1-p_g p_g]);         % Select gamma: 80/20 mix [default 1]
            dX(i,A) = c_star*randn(1,d_star) + ...
                (1+lambda(i))*g*sum(Z(a,A)-Z(b,A),1);                 % Compute ith jump diff. evol.
         elseif strcmp(method,'snooker'),                             % SNOOKER update
            id(i) = n_CR;                                             % Full crossover
            F = X(i,1:d)-Z(a,1:d); D_e = max(F*F',1e-300);            % Define projection X(i,1:d) - Z(a,1:d)
            zP = F*( sum((Z(b,1:d)-Z(c_sn,1:d)).*F ) / D_e );         % Orthogonally project zR1 and zR2 onto F
            g = 1.2 + rand;                                           % Determine jump rate
            dX(i,1:d) = c_star*randn(1,d) + (1+lambda(i))*g*zP;       % Calculate ith jump snooker update
        end
        Xp(i,1:d) = X(i,1:d) + dX(i,1:d);                             % Compute ith proposal
        
        for kk = 1:d  
          if strfind('unif',marg_prior{1,kk})                         % Bound handling: If proposal outside
            if Xp(i,kk) < marg_par(1,kk)                              % uniform prior bound, add/subtract distance
             Xp(i,kk) = Xp(i,kk) - marg_par(1,kk) + marg_par(2,kk);   % to opposite boundary (i.e. fold)
            elseif Xp(i,kk) > marg_par(2,kk)
             Xp(i,kk) =  Xp(i,kk) - marg_par(2,kk) + marg_par(1,kk);
            end
          end
        end
    end

    if strcmp(method,'snooker'),                                      % Snooker correction: non-symmetry
        alfa_sn = (sum((Xp - Z(R(1,1:N),1:d)).^2,2)./sum((X - ...     % proposal distribution
                Z(R(1,1:N),1:d)).^2,2)).^((d-1)/2);
    elseif strcmp(method,'parallel'),                                 % otherwise no correction needed
        alfa_sn = ones(N,1);
    end
    for i = 1:N,                                                      % Accept/reject proposals (can use "parfor")
        p_Xp(i,1) = pdf(Xp(i,1:d));                                   % Calculate density ith proposal
        p_acc(i) = min(1,alfa_sn(i)*[p_Xp(i,1)./p_X(i,1)]);           % Compute acceptance probability
    end
    for i = 1:N,                                                      % Accept/reject proposals
        if p_acc(i) > rand,                                           % p_acc(i) larger than U[0,1]?
            X(i,1:d) = Xp(i,1:d); p_X(i,1) = p_Xp(i,1);               % True: Accept proposal
        else
            dX(i,1:d) = 0;                                            % Set jump back to zero for p_CR
        end
        J(id(i)) = J(id(i)) + sum((dX(i,1:d)./std_X).^2);             % Update jump distance crossover idx
        n_id(id(i)) = n_id(id(i)) + 1;                                % How many times idx crossover used
    end
    x(t,1:d,1:N) = reshape(X',1,d,N); p_x(t,1:N) = p_X';              % Append current X and density
    if (mod(t,k) == 0),                                               % Check whether to append X to archive Z
        Z(m+1:m+N,1:d+1) = [X p_X]; m = m + N;                        % Append current values to Z after k generations        
        if t<T/10, p_CR = J./n_id; p_CR = p_CR/sum(p_CR); end         % Update selection prob. crossover
    end
    % --------------------------------------------------------
    % If you want you can compute here convergence diagnostics 
    % --------------------------------------------------------  
    % End dynamic part
    a = a+1;
    if a == 5 || t == T;
    clc
    disp(['DREAM_(ZS) ' num2str(round(t/T*100)) '% complete'])
    a = 1; 
    end
    
end

