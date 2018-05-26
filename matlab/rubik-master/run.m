clear all
clc

addpath(genpath('./tensor_toolbox'));

%%% commented out yichen's code
% dim = [10, 10, 10];
% K = 5;
% 
% [X_kt, guide] = geneData_mod_RC(dim, K);
% X = tensor(X_kt);
%%%%%%%%%%%

% define the rank
K = 30;




% define the tensor
sprintf('load the csv indexes')
sptensor_indexes = csvread('../sptensor_sutter_CHFCASES_PATIENT_DATE_Tminus3_Tminus1_forPipeline_TENSORCOORDS.csv');
sptensor_indexes = sptensor_indexes+1; % add one because matlab is 1-indexed
binary_values = ones(numel(sptensor_indexes(:,1)),1);

sprintf('create sptensor')
sptensor_X = sptensor(sptensor_indexes, binary_values);

sprintf('create tensor -- large (32gb)')
tensor_X = tensor(sptensor_X);

sprintf('define the guidance -- all 0')
guide = []; %no guidance

sprintf('rank: ')
K % print rank

sprintf('run rubik')
tic  % to time it
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
