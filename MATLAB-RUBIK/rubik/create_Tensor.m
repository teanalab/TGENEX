%read binary files as matrices
readClinical
readMutation

binaC = MatlabbinaC';
binaM = MatlabbinaM;

dC = size(binaC)
dM = size(binaM)

sum(sum(binaM)) %checking
sum(sum(binaC)) %checking
 
TenCxM = zeros(dM(1), dC(1), dC(2) ); %genes x patients x clinical
for i = 1:dC(2) %num of clinical
    for j = 1:dC(1) %num of patients
        patient_i_j =  binaC(j,i);  
        if patient_i_j == 1
            TenCxM(:,j,i) = binaM(:,j);
        end
    end
end

save("TenCxM.mat",'TenCxM')