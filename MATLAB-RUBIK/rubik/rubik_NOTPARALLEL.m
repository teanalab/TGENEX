function [T, C] = rubik(X, K, guide)
% Rubik computes non-negative sparse tensor factorization with guidance information
%
% INPUT
% X: observed 3-way tensor
% K: rank of CP factorization
% guide: guidance matrix of size: numDiagnosis-by-numGuidance
%
% OUTPUT
% T is an interaction ktensor
% C is a bias tensor
%
% Examples: see RunThisFile.m
%
% Citation: Rubik: Knowledge Guided Tensor Factorization and Completion for Health Data Analytics 
% Yichen Wang, Robert Chen, Joydeep Ghosh, Joshua Denny, Abel Kho, You Chen, Bradley Malin, and Jimeng Sun, KDD 2015.
%
% Dependency: Matlab tensor toolbox v 2.6 by T. Kolda, B. Bader, and others.
% http://www.sandia.gov/~tgkolda/TensorToolbox.


%% set algorithm parameters
printitn = 1;

maxiter= 200;

fitchangetol = 1e-4;

lambda_a = 5; % weight for guide

lambda_q = 1; % weight for pairwise

gamma = [0.01 , 0.01,  0.01]; % sparse threshold

mu = 1e-6;
eta = 1e-6;
mumax = 1e10;
etamax = 1e10;
rho = 1.15;

%% compute statistics for X and guide
dim = size(X);

normX = norm(X);

numGuide = size(guide, 2);

guide = [guide, zeros(dim(2), K - numGuide)];

W = [eye(numGuide), zeros(numGuide, K - numGuide); zeros(K - numGuide,K)]; % K-by-K symetric matrix

%% initialization

% interaction tensor 
A = cell(1,3);
for i = 1:3
    A{i} = rand(dim(i),K);
    A{i}(A{i} <= gamma(i)) = 0; 
end
B = A;

% bias tensor
u = cell(1,3);

for i = 1:3
    u{i} = rand(dim(i),1)/10;
end
v = u;

% Lagrange multipliers
Y = cell(1,3);
for i = 1:3
    Y{i} = zeros(dim(i),K);
end

p = cell(1,3);
for i = 1:3
    p{i} = zeros(dim(i),1);
end

%% main loop

fit = 0;

for iter = 1: maxiter
    
    fitold = fit;
    
    R = X - tensor(ktensor(u)); 
    E = X - tensor(ktensor(A));
    
    id1 = find(R<0);
    id2 = find(E<0);
    if ~isempty(id1)
        R(id1) = 0;
    end
    if ~isempty(id2)
        E(id2) = 0;
    end
    
    for n = 1: 3
        
        % update A{n}        
        pitpi = ones(K,K);
        for i = [1:n-1,n+1:3]
            pitpi = pitpi .* (A{i}' * A{i}); % compute \Pi^t\Pi in Eq.5
        end
        
        if n == 2
            a = lambda_q * B{n} * B{n}';
            b = 2 * pitpi + lambda_a * W + mu * eye(K);
            c = 2 * mttkrp(R,A,n) + lambda_a * guide * W  + lambda_q * B{n} + mu * B{n} + Y{n};
            
            A{n} = sylvester(a,b,c);
  
        else
            term1 = 2 * mttkrp(R,A,n) + mu * B{n} + Y{n};
            
            term2 = 2 * (pitpi) + mu * eye(K);
            
            A{n} = term1/term2;
        end
        
        % update B{n}
        B{n} = A{n} + Y{n}/mu;
        
        % update Y{n}
        Y{n} = Y{n} + mu * (B{n} - A{n});
        
        % update u{n}      
        LambdatLambda = 1;
        for i = [1:n-1,n+1:3]
            LambdatLambda = LambdatLambda .* (u{i}' * u{i});
        end
        
        term1 = 2 * mttkrp(E,u,n) + eta * v{n} + p{n};
        
        term2 = 2 * LambdatLambda + eta;
        
        u{n} = term1/term2;
        
        % update v{n}
        v{n} = u{n} + p{n}/eta;
        
        % update p{n}
        p{n} = p{n} + eta * (v{n} - u{n});
        
    end
    
    % update mu, eta
    mu = min(rho * mu, mumax);
    
    eta = min(rho * eta, etamax);
    
    % compute the fit
    T = ktensor(A);
    
    C = ktensor(u);
    
    normresidual = sqrt( normX^2 + norm(T+C)^2 - 2 * innerprod(X,T+C) );
    
    fit = 1 - (normresidual / normX);
    
    fitchange = abs(fitold - fit);
    
    if mod(iter,printitn)==0
        fprintf(' Iter %2d: fit = %e fitdelta = %7.1e\n', iter, fit, fitchange);
    end
    
    % Check for convergence
    if (iter > 1) && (fitchange < fitchangetol)
        break;
    end
    
end

%% clean up final results
T = arrange(T); % columns are normalized
C = arrange(C);

for i = 1: 3
    T{i}(T{i} <= gamma(i)) = 0; % sparse factor
    
end

if printitn>0
    
    normresidual = sqrt( normX^2 + norm(T+C)^2 - 2 * innerprod(X,T+C) );
    
    fit = 1 - (normresidual / normX);
    
    fprintf(' Final fit = %e \n', fit);
end

end