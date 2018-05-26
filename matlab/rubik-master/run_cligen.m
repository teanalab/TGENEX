clear all
clc
% 
% addpath(genpath('./tensor_toolbox'));
% 
% %%% commented out yichen's code
% % dim = [10, 10, 10];
% % K = 5;
% % 
% % [X_kt, guide] = geneData_mod_RC(dim, K);
% % X = tensor(X_kt);
% %%%%%%%%%%%
% 
% define the rank
K = 10;


%read binary matrices
read_binaC
read_binaM

 
TenCxM = zeros(2714, 456, 38);
for i = 1:38
    for j = 1:456
        patient_i_j =  binaC(j,i);
        if patient_i_j == 1
            TenCxM(:,j,i) = binaM(:,j);
        end
    end
end

sum(sum(sum(TenCxM)))
%94419

% define the tensor

sprintf('create tensor -- large (32gb)')
tensor_X = tensor(TenCxM);

sptensor_X = sptensor(tensor_X);

sprintf('define the guidance -- all 0')
guide = []; %no guidance
 
sprintf('rank: ')
K % print rank

sprintf('run rubik')
tic  % to time it
%X = sptensor_X;
[T, C]= rubik(sptensor_X, K, guide);
toc % to time it
 



% check the guidance matrix
guide
T{2}

% check orthogonality of diagnosis mode
T{2}'*T{2}


lambda = T.lambda;
u1 = T.u{1};
u2 = T.u{2};
u3 = T.u{3};
 
time_now = datestr(datetime('now'));
save_folder = strcat('R_', num2str(K), '_', time_now, '/');
mkdir(save_folder);
csvwrite(strcat(save_folder, 'lambda.csv'),lambda)
csvwrite(strcat(save_folder, 'u1.csv'),u1)
csvwrite(strcat(save_folder, 'u2.csv'),u2)
csvwrite(strcat(save_folder, 'u3.csv'),u3)


% bias tensor
bias_lambda = C.lambda;
bias_u1 = C.u{1};
bias_u2 = C.u{2};
bias_u3 = C.u{3};
csvwrite(strcat(save_folder, 'bias_lambda.csv'),bias_lambda)
csvwrite(strcat(save_folder, 'bias_u1.csv'),bias_u1)
csvwrite(strcat(save_folder, 'bias_u2.csv'),bias_u2)
csvwrite(strcat(save_folder, 'bias_u3.csv'),bias_u3)
