%% NBS demo of recovering and simulating data
% 

%library_path = '/Users/diamac/Desktop/nbs_release_v0.2';
%basedata_path = '/Users/diamac/Desktop/nbs_release_v0.2/data';

library_path = 'C:\\Users\\Diana\\Google Drive\\nbs_release_v0.2_workingC';
basedata_path = 'C:\\Users\\Diana\\Google Drive\\nbs_release_v0.2_workingC\\data';

addpath(genpath(library_path))

%% Load an NBS propagation network

NBSnetwork = load([ basedata_path '/networks/ST90Q_adj_mat.mat' ]);   

%% Load an NBS graph influence measure 

glap = load([ basedata_path '/networks/glap_subnetwork_ST90.mat' ]);
outDeg = 11;
[knnGlap] = sim_to_knn_glap(glap.glap,outDeg);
clear glap;

%% Load somatic mutation cancer data

baseSMDataOld = load([ basedata_path '/TCGA_somatic_mutations/somatic_data_UCEC.mat' ]);

%        gene_id_all: [17968×1 double]
%     gene_indiv_mat: [248×17968 double]
%          sample_id: {248×1 cell}
%            batchIs: [1×1 struct]
%     gene_id_symbol: {17968×1 cell}

%BRCA
load('C:\Users\Diana\Google Drive\nbs_release_v0.2_workingC\BCRAdata\somatic_data_BRCA.mat');
baseSMData = baseSMDataBRCA;
clearvars baseSMDataBRCA;
min_mutation = 1; %mutations already filtered in R
baseSMData.sample_id(sum(baseSMData.gene_indiv_mat,2) < min_mutation) =  [];
baseSMData.gene_indiv_mat((sum(baseSMData.gene_indiv_mat,2) < min_mutation),:) =  [];

for inx = 1:456
		baseSMData.sample_id{inx} = char(baseSMData.sample_id{inx});
end

%regexprep(baseSMData.sample_id,'(TCGA-..-....).*','$1') %already came
%shorten

%% Load phenotype

load([basedata_path '/TCGA_somatic_mutations/UCEC_clinical_phenotype.mat' ]);
oldPheno = UCECppheno;
%These are the variables that the example has:
% UCECppheno = 
%   struct with fields:
%         sample_id: {451×1 cell}
%     days_to_death: [451×1 double]
%      days_to_last: [451×1 double]
%     days_survival: [451×1 double]
%          diag_age: [451×1 double]
%              race: {451×1 cell}
%         ethnicity: {451×1 cell}
%            gender: {451×1 cell}
%            living: {451×1 cell}
%         histICDO3: {451×1 cell}
%       hist_simple: {451×1 cell}
%         ICDO3site: {451×1 cell}
%             ICD10: {451×1 cell}
%      ICD10_simple: [451×1 double]
%             stage: {451×1 cell}
%      stage_simple: {451×1 cell}
%      stage_ismeta: [451×1 double]
%             grade: {451×1 cell}
%      grade_simple: {451×1 cell}
%          residual: {451×1 cell}

UCECppheno = load('C:\Users\Diana\Google Drive\nbs_release_v0.2_workingC\BCRAdata\BRCA_clinical_phenotype.mat' );


%% 

[~,ia,ib] = intersect(baseSMData.sample_id,UCECppheno.sample_id);

%jointPheno = cellfun(@(x,y)([x '-' y]),UCECppheno.hist_simple,UCECppheno.grade_simple,'uniformoutput',0);
jointPheno = cellfun(@(x,y)([x '-' y]),UCECppheno.patienthistological_type,UCECppheno.AJCCStage,'uniformoutput',0);
cphenoJoint.sample_id = UCECppheno.sample_id(ib);
cphenoJoint.phenoMat = [grp2idx(jointPheno(ib)) grp2idx(UCECppheno.patienthistological_type(ib)) grp2idx(UCECppheno.AJCCStage(ib)) ];


%% Standard NMF

cczoptions.iter = 1000;
cczoptions.tof = 1e-3;
cczoptions.dis = true;
cczoptions.dist = 'nnls';
cnum = 15; %

[A,Y]=nmfnnls(baseSMData.gene_indiv_mat',cnum,cczoptions);

indClust_ctrl=NMFCluster(Y);
[table, chi2, p] = crosstab (indClust_ctrl(ia),cphenoJoint.phenoMat(:,1));
[AR, RI, MI, HI] = RandIndex(indClust_ctrl(ia),cphenoJoint.phenoMat(:,1));
fprintf('Sim_matrix: P:%e AR:%f\n', p, AR);
% output
% NNLS based NMF successes! 
%  # of iterations is 220. 
%  The final residual is 1.0601e+02.
% Sim_matrix: P:2.549643e-01 AR:0.003821


csvwrite('BCRAdata/output_standardNMF_k15.csv',Y);

%% Run NBS w/o consensous clustering

% Type of NMF to use:
propnmf_options.nmf_type = 'netnmf';
% Number of clusters 
propnmf_options.K = [ 15 ];
% Minimum number of mutations persample to be includeded 
propnmf_options.min_mutations = 1; %already filtered
% Propagation value
propnmf_options.proV = NBSnetwork.network.propVal;
% Verbose output
propnmf_options.dis = true;
% Show progress bar (requries X11)
propnmf_options.pb = true;
% Quantile normzlize
propnmf_options.normalize_rows = 1;
% Use netwpork
propnmf_options.isNetwork = 1;


[indClust] = ...    
    nbs_single_main(baseSMData.gene_indiv_mat,baseSMData.gene_id_all,NBSnetwork.network,propnmf_options,knnGlap);

% Using quick propagate
% Fast prop
% 1.270600
% Converged 23
% Normalize
% Net NMF


csvwrite('BCRAdata/output_clustersNBS_k15.csv',indClust);

%%
[table, chi2, p] = crosstab(indClust(ia),cphenoJoint.phenoMat(:,1));
[AR,RI,MI,HI]=RandIndex(indClust(ia),cphenoJoint.phenoMat(:,1));
fprintf('NBS: P:%e AR:%f\n',p,AR);

%NBS: P:1.383783e-02 AR:0.010385

%for 10
%NBS: P:2.724142e-03 AR:0.003950

%for 15
%NBS: P:3.219363e-04 AR:0.005643


%% Run consensous clustering (CC) NBS

% Type of NMF to use:
propnmf_options.nmf_type = 'netnmf';
% Number of times to perform NMF as part of CC
% (100 is probably the bare minimum)
propnmf_options.nsample = 1000;
% Proportion of rows (patient samples) to include in CC
propnmf_options.smp_ind = 0.8;
% Proportion of columns (patient samples) to include in CC
propnmf_options.smp_feat = 0.8;
% Number of clusters 
propnmf_options.K = [ 2 3 4 5 6 7 8 9 10 11 12 13 14 15];
% Minimum number of mutations persample to be includeded 
propnmf_options.min_mutations = 1; %already taken into account
% Propagation value
propnmf_options.proV = NBSnetwork.network.propVal;
% Verbose output
propnmf_options.dis = true;
% Show progress bar (requries X11)
propnmf_options.pb = true;
% Quantile normzlize
propnmf_options.normalize_rows = 1;
% Use netwpork
propnmf_options.isNetwork = 1;

[NBS_cc_cluster] = ...    
    nbs_cc_main(baseSMData.gene_indiv_mat,baseSMData.gene_id_all,NBSnetwork.network,propnmf_options,knnGlap);
%%
[NBS_cc_label] = calc_label(NBS_cc_cluster.clust_hard.cc_matrix{cnum},cnum,'av',1);
[table, chi2, p] = crosstab(NBS_cc_label(ia),cphenoJoint.phenoMat(:,1));
[AR,RI,MI,HI]=RandIndex(NBS_cc_label(ia),cphenoJoint.phenoMat(:,1));
fprintf('CC NBS: P:%e AR:%f\n',p,AR);


% Normalize
% Net NMF
% Net NMF
% Round 36 time - 88.429341
% Using quick propagate
% Fast prop
% 0.268793
% Converged 17
% Normalize
% Net NMF
% Net NMF
% Round 37 time - 85.172626
% Using quick propagate
% Fast prop
% 0.268840
% Converged 15
% Normalize
% Net NMF
% Net NMF
% Round 38 time - 78.254978
% Using quick propagate
% Fast prop
% 0.272440
% Converged 17
% Normalize
% Net NMF
% Net NMF
% Round 39 time - 90.609277
% Using quick propagate
% Fast prop
% 0.265763
% Converged 17
% Normalize
% Net NMF
% Net NMF
% Round 40 time - 89.063606
% Using quick propagate
% Fast prop
% 0.247059
% Converged 17
% Normalize
% Net NMF
% Net NMF
% Round 41 time - 89.358525
% Using quick propagate
% Fast prop
% 0.263747
% Converged 17
% Normalize
% Net NMF
% Net NMF
% Round 42 time - 78.545026
% Using quick propagate
% Fast prop
% 0.276953
% Converged 17
% Normalize
% Net NMF
% Net NMF
% Round 43 time - 87.901241
% Using quick propagate
% Fast prop
% 0.267930
% Converged 17
% Normalize
% Net NMF
% Net NMF
% Round 44 time - 83.025532
% Using quick propagate
% Fast prop
% 0.269265
% Converged 17
% Normalize
% Net NMF
% Net NMF
% Round 45 time - 90.224702
% Using quick propagate
% Fast prop
% 0.275821
% Converged 17
% Normalize
% Net NMF
% Net NMF
% Round 46 time - 92.648659
% Using quick propagate
% Fast prop
% 0.270385
% Converged 17
% Normalize
% Net NMF
% Net NMF
% Round 47 time - 90.069617
% Using quick propagate
% Fast prop
% 0.267965
% Converged 17
% Normalize
% Net NMF
% Net NMF
% Round 48 time - 93.856065
% Using quick propagate
% Fast prop
% 0.274726
% Converged 17
% Normalize
% Net NMF
% Net NMF
% Round 49 time - 87.474521
% Using quick propagate
% Fast prop
% 0.270868
% Converged 17
% Normalize
% Net NMF
% Net NMF
% Round 50 time - 94.418060
% Using quick propagate
% Fast prop
% 0.264664
% Converged 17
% Normalize
% Net NMF
% Net NMF
% Round 51 time - 88.819274
% Using quick propagate
% Fast prop
% 0.250738
% Converged 17
% Normalize
% Net NMF
% Net NMF
% Round 52 time - 85.930492
% Using quick propagate
% Fast prop
% 0.265768
% Converged 17
% Normalize
% Net NMF
% Net NMF
% Round 53 time - 88.458115
% Using quick propagate
% Fast prop
% 0.269076
% Converged 17
% Normalize
% Net NMF
% Net NMF
% Round 54 time - 92.949137
% Using quick propagate
% Fast prop
% 0.274031
% Converged 16
% Normalize
% Net NMF
% Net NMF
% Round 55 time - 87.025897
% Using quick propagate
% Fast prop
% 0.269091
% Converged 17
% Normalize
% Net NMF
% Net NMF
% Round 56 time - 88.372516
% Using quick propagate
% Fast prop
% 0.268119
% Converged 16
% Normalize
% Net NMF
% Net NMF
% Round 57 time - 89.926677
% Using quick propagate
% Fast prop
% 0.268697
% Converged 17
% Normalize
% Net NMF
% Net NMF
% Round 58 time - 89.397841
% Using quick propagate
% Fast prop
% 0.266946
% Converged 17
% Normalize
% Net NMF
% Net NMF
% Round 59 time - 88.886274
% Using quick propagate
% Fast prop
% 0.274392
% Converged 17
% Normalize
% Net NMF
% Net NMF
% Round 60 time - 90.313408
% Using quick propagate
% Fast prop
% 0.261103
% Converged 17
% Normalize
% Net NMF
% Net NMF
% Round 61 time - 90.486625
% Using quick propagate
% Fast prop
% 0.249941
% Converged 16
% Normalize
% Net NMF
% Net NMF
% Round 62 time - 82.126236
% Using quick propagate
% Fast prop
% 0.276649
% Converged 16
% Normalize
% Net NMF
% Net NMF
% Round 63 time - 92.994990
% Using quick propagate
% Fast prop
% 0.270371
% Converged 17
% Normalize
% Net NMF
% Net NMF
% Round 64 time - 95.542806
% Using quick propagate
% Fast prop
% 0.267490
% Converged 17
% Normalize
% Net NMF
% Net NMF
% Round 65 time - 91.536484
% Using quick propagate
% Fast prop
% 0.272190
% Converged 17
% Normalize
% Net NMF
% Net NMF
% Round 66 time - 92.500646
% Using quick propagate
% Fast prop
% 0.264021
% Converged 17
% Normalize
% Net NMF
% Net NMF
% Round 67 time - 86.655124
% Using quick propagate
% Fast prop
% 0.267773
% Converged 17
% Normalize
% Net NMF
% Net NMF
% Round 68 time - 85.002445
% Using quick propagate
% Fast prop
% 0.265825
% Converged 17
% Normalize
% Net NMF
% Net NMF
% Round 69 time - 92.914176
% Using quick propagate
% Fast prop
% 0.269835
% Converged 17
% Normalize
% Net NMF
% Net NMF
% Round 70 time - 92.248765
% Using quick propagate
% Fast prop
% 0.269583
% Converged 17
% Normalize
% Net NMF
% Net NMF
% Round 71 time - 91.735192
% Using quick propagate
% Fast prop
% 0.262149
% Converged 17
% Normalize
% Net NMF
% Net NMF
% Round 72 time - 90.431345
% Using quick propagate
% Fast prop
% 0.271893
% Converged 17
% Normalize
% Net NMF
% Net NMF
% Round 73 time - 93.204114
% Using quick propagate
% Fast prop
% 0.257554
% Converged 17
% Normalize
% Net NMF
% Net NMF
% Round 74 time - 88.060999
% Using quick propagate
% Fast prop
% 0.250917
% Converged 17
% Normalize
% Net NMF
% Net NMF
% Round 75 time - 91.277038
% Using quick propagate
% Fast prop
% 0.250185
% Converged 16
% Normalize
% Net NMF
% Net NMF
% Round 76 time - 87.106299
% Using quick propagate
% Fast prop
% 0.243447
% Converged 17
% Normalize
% Net NMF
% Net NMF
% Round 77 time - 91.335137
% Using quick propagate
% Fast prop
% 0.270694
% Converged 17
% Normalize
% Net NMF
% Net NMF
% Round 78 time - 95.025571
% Using quick propagate
% Fast prop
% 0.276222
% Converged 16
% Normalize
% Net NMF
% Net NMF
% Round 79 time - 88.454535
% Using quick propagate
% Fast prop
% 0.267072
% Converged 17
% Normalize
% Net NMF
% Net NMF
% Round 80 time - 86.052881
% Using quick propagate
% Fast prop
% 0.270255
% Converged 17
% Normalize
% Net NMF
% Net NMF
% Round 81 time - 88.334851
% Using quick propagate
% Fast prop
% 0.275695
% Converged 17
% Normalize
% Net NMF
% Net NMF
% Round 82 time - 87.811294
% Using quick propagate
% Fast prop
% 0.275822
% Converged 17
% Normalize
% Net NMF
% Net NMF
% Round 83 time - 92.147016
% Using quick propagate
% Fast prop
% 0.259308
% Converged 17
% Normalize
% Net NMF
% Net NMF
% Round 84 time - 93.073754
% Using quick propagate
% Fast prop
% 0.265695
% Converged 16
% Normalize
% Net NMF
% Net NMF
% Round 85 time - 84.551936
% Using quick propagate
% Fast prop
% 0.235955
% Converged 16
% Normalize
% Net NMF
% Net NMF
% Round 86 time - 87.703584
% Using quick propagate
% Fast prop
% 0.258163
% Converged 17
% Normalize
% Net NMF
% Net NMF
% Round 87 time - 86.421909
% Using quick propagate
% Fast prop
% 0.266257
% Converged 17
% Normalize
% Net NMF
% Net NMF
% Round 88 time - 91.881203
% Using quick propagate
% Fast prop
% 0.272639
% Converged 17
% Normalize
% Net NMF
% Net NMF
% Round 89 time - 92.959388
% Using quick propagate
% Fast prop
% 0.273778
% Converged 17
% Normalize
% Net NMF
% Net NMF
% Round 90 time - 87.285914
% Using quick propagate
% Fast prop
% 0.267394
% Converged 17
% Normalize
% Net NMF
% Net NMF
% Round 91 time - 89.228854
% Using quick propagate
% Fast prop
% 0.278440
% Converged 17
% Normalize
% Net NMF
% Net NMF
% Round 92 time - 92.405676
% Using quick propagate
% Fast prop
% 0.277200
% Converged 17
% Normalize
% Net NMF
% Net NMF
% Round 93 time - 92.083420
% Using quick propagate
% Fast prop
% 0.256370
% Converged 17
% Normalize
% Net NMF
% Net NMF
% Round 94 time - 90.334567
% Using quick propagate
% Fast prop
% 0.254270
% Converged 17
% Normalize
% Net NMF
% Net NMF
% Round 95 time - 90.585575
% Using quick propagate
% Fast prop
% 0.273945
% Converged 17
% Normalize
% Net NMF
% Net NMF
% Round 96 time - 84.426440
% Using quick propagate
% Fast prop
% 0.270031
% Converged 17
% Normalize
% Net NMF
% Net NMF
% Round 97 time - 81.715675
% Using quick propagate
% Fast prop
% 0.262260
% Converged 17
% Normalize
% Net NMF
% Net NMF
% Round 98 time - 86.123044
% Using quick propagate
% Fast prop
% 0.275417
% Converged 17
% Normalize
% Net NMF
% Net NMF
% Round 99 time - 92.679643
% Using quick propagate
% Fast prop
% 0.272601
% Converged 16
% Normalize
% Net NMF
% Net NMF
% Round 100 time - 86.956908
% Done clustering compiling results
% Notify: Transforming into distance metric
% Subscript indices must either be real positive integers or logicals.
% 
% Error in Contingency (line 20)
%    Cont(Mem1(i),Mem2(i))=Cont(Mem1(i),Mem2(i))+1;
% 
% Error in RandIndex (line 18)
% C=Contingency(c1,c2);	%form contingency matrix



%% Run standard consensous clustering NMF

propnmf_options.nmf_type = 'nmf';
propnmf_options.isNetwork = 0;

[ctrl_cc_cluster] = ...    
    nbs_cc_main(baseSMData.gene_indiv_mat,baseSMData.gene_id_all,[],propnmf_options);

[ctrl_cc_label] = calc_label(ctrl_cc_cluster.clust_hard.cc_matrix{cnum},cnum,'av',1);
[table, chi2, p] = crosstab(ctrl_cc_label(ia),cphenoJoint.phenoMat(:,1));
[AR,RI,MI,HI]=RandIndex(ctrl_cc_label(ia),cphenoJoint.phenoMat(:,1));
fprintf('Standard CC NMF: P:%e AR:%f\n',p,AR);

