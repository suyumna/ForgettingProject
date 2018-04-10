% Script to analyze RS in healthy subjects stimulated with cTBS en Left M1.
% Each subject goes to the scanner for a RS measurement during 10 min before the stimulation
% After that the subject is stimulated during 40s with cTBS over the left Motor cortex
% Immediately after the stimulation, they go to the scanner again for a RS measurement during 10 min 




cd '/data/projects/NeFF_B5/scripts_connectivity/'
disp('Loading AAL data...')
Y=spm_vol('roi_mnifsl_222.img'); aal=spm_read_vols(Y);
disp('AAL data loaded...')

% Before stimulation
disp('Processing Before left M1 stimulation...');


% DR, DF, Neu, Neg

for i = 1
        
    i
    subj_path = strcat('/data/projects/NeFF_B5/Real_Study/Forgetting/Basic/S', num2str(i), '/analysis/DF_GT/SPM.mat');
    dir_path= strcat('/data/projects/NeFF_B5/Real_Study/Forgetting/Basic/S', num2str(i), '/analysis/DF_GT/');
    load(subj_path)     % load spm.mat file to get the conditions
    nbeta=numel(SPM.Vbeta);  % here I get the different parameters of each beta estimate
    
    c1=1; c2=1; c3=1; c4=1; % need a counter to add the time dimension to the 3D beta data
    
    for i_beta=1:nbeta
        
        data_path{i_beta} = strcat(dir_path, 'beta_', num2str(i_beta, '%04.f'),'.nii'); % the instruction '%04.f' indicates that there will be 4 digits before the decimal point padded with zero depending on the number i_beta
        
        Y = spm_vol(data_path{i_beta}); data=spm_read_vols(Y);
        
        name_beta=SPM.Vbeta(i_beta).descrip; process_set=name_beta(29:30); % Here I get the name of the condition of each beta estimate:NR, NF, DR, DF
        switch process_set            %SPM.Vbeta(1) from here get the conditions
        case 'NR'
        data_1(:,:,:,c1)=spm_read_vols(Y); c1=c1+1; % creates the 4D data and increment counter when condition found
        case 'NF'
        data_2(:,:,:,c2)=spm_read_vols(Y); c2=c2+1;
        case 'DR'
        data_3(:,:,:,c3)=spm_read_vols(Y); c3=c3+1;
        case 'DF'
        data_4(:,:,:,c4)=spm_read_vols(Y); c4=c4+1;
        end
    end
    
    c1
    c2
    c3
    c4
    % ======    
    for i_data=1:4
        if i_data==1   
         data=data_1;
         save_path = strcat('/data/projects/NeFF_B5/Real_Study/Forgetting/Basic/S', num2str(i), '/analysis/graph_theory/NR_results_AAL.mat') ;  
         elseif i_data==2
         data=data_2;
         save_path = strcat('/data/projects/NeFF_B5/Real_Study/Forgetting/Basic/S', num2str(i), '/analysis/graph_theory/NF_results_AAL.mat') ;
         elseif i_data==3
         data=data_3;
         save_path = strcat('/data/projects/NeFF_B5/Real_Study/Forgetting/Basic/S', num2str(i), '/analysis/graph_theory/DR_results_AAL.mat') ;
         elseif i_data==4
         data=data_4;
         save_path = strcat('/data/projects/NeFF_B5/Real_Study/Forgetting/Basic/S', num2str(i), '/analysis/graph_theory/DF_results_AAL.mat') ; 
         end
         
    % now data is [x,y,z,time] matrix and we shouldn't need to modify any of these scripts
     [adj_matrix, c,l,e,Q,Ci, c_rand, l_rand, e_rand, Q_rands, Ci_rand, degree, bet, thr, degree_rand, bet_rand] = network_measures_aal_newman_Ben(data, aal);
     clearvars -except adj_matrix i c l e Q Ci c_rand l_rand e_rand Q_rands Ci_rand degree bet thr degree_rand bet_rand aal ids num id
     
     
    save(save_path) 
    end
       
end

        
        
  