function [X, guide] = geneData_mod(dim, K)
m = dim(1);
n = dim(2);
k = dim(3);

A1 = rand(m,K);
A3 = rand(k,K);

A2 = zeros(n,K); % A2 has K orthogonal columns
A2(1:3,1) = ones(3,1);
A2(4:6,2) = ones(3,1);
A2(7:8,3) = ones(2,1);

A2(9,4) = 1;
A2(10,5) = 1;

numGuide = 0;
guide = A2(:,1:numGuide); % guidance matrix

u1 = rand(m,1)/10;
u2 = rand(n,1)/10;
u3 = rand(k,1)/10;

A = {A1,A2,A3};
u = {u1,u2,u3};
X = ktensor(A) + ktensor(u);
