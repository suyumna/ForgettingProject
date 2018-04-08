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
   %if not  (i ==7 ||i ==18)   % para evitar el 7 and 18
    subj_path = strcat('/data/projects/NeFF_B5/Real_Study/Forgetting/Basic/S', num2str(i), '/analysis/DF_GT/SPM.mat');
    load(subj_path)
    % load functional volumes  -- defino el path e inserto el numero de
    % sujeto pasandolo a string con num2str y concateno los strings con strcat
    % different paths for different trial types
    
    % ====== I added this, according to what I understood from you.....
    nbeta=numel(SPM.Vbeta);
    % so nbeta is your 'time' dimension, we gotta make your data 4d before we can do anything
    c=1; % need a counter for this style
    for i_beta=1:nbeta
        % see n2sp(number,padding) is simpler imo e.g. n2sp(1,4) produces the string "0001"
       if i_beta < 10
          data_path{i_beta} = strcat('/data/projects/NeFF_B5/Real_Study/Forgetting/Basic/S', num2str(i), '/analysis/DF_GT/beta_000', num2str(i_beta),'.nii');
       else
          data_path{i_beta} = strcat('/data/projects/NeFF_B5/Real_Study/Forgetting/Basic/S', num2str(i), '/analysis/DF_GT/beta_00', num2str(i_beta),'.nii'); 
       end
       % data_path should be a cell structure imo, depends how you wanna implement and what knowledge you wanna save
       
        name_beta=SPM.Vbeta(i_beta).descrip;
        process_set=name_beta(29:30); % explicit references are safer imo e.g. some type of str matching
    % ======    
    if strcmp(process_set,condition) % condition needs defined in 2nd loop
    
       Y = spm_vol(data_path{i_beta}); data=spm_read_vols(Y);
        data(:,:,:,c)=spm_read_vols(Y); c=c+1; % increment counter when condition found
        % This loop is essential to making the 4d file! its gotta stay in
        % for jj=1:length(data_path)
         %    Y = spm_vol(data_path{jj}); 
         %     data(:,:,:,i_beta)=spm_read_vols(Y);
        %end
        
        % if you're looping through all nbeta you're gonna get constants and motion parameters too, which I don't think you want
        
    end % of nbeta loop! we need a 4d file before we can do any network stuff

    % now data is [x,y,z,time] matrix and we shouldn't need to modify any of these scripts
        [adj_matrix, c,l,e,Q,Ci, c_rand, l_rand, e_rand, Q_rands, Ci_rand, degree, bet, thr, degree_rand, bet_rand] = network_measures_aal_newman_Ben(data, aal);
        clearvars -except adj_matrix i c l e Q Ci c_rand l_rand e_rand Q_rands Ci_rand degree bet thr degree_rand bet_rand aal ids num id
        
        % ======
        % ====== I added this, according to what I understood from you.....
        % this confuses me? i think this loop should be 2nd e.g. subject loop, condition loop, vbeta loop
        
        switch process_set
            %SPM.Vbeta(1) from here get the conditions
        case 'NR'
        save_path = strcat('/data/projects/NeFF_B5/Real_Study/Forgetting/Basic/S', num2str(i), '/analysis/graph_theory/NR_results_AAL_', num2str(i_beta),'.mat') ;  
        case 'NF'
        save_path = strcat('/data/projects/NeFF_B5/Real_Study/Forgetting/Basic/S', num2str(i), '/analysis/graph_theory/NF_results_AAL_', num2str(i_beta),'.mat') ;
        case 'DR'
        save_path = strcat('/data/projects/NeFF_B5/Real_Study/Forgetting/Basic/S', num2str(i), '/analysis/graph_theory/DR_results_AAL_', num2str(i_beta),'.mat') ;
        case 'DF'
        save_path = strcat('/data/projects/NeFF_B5/Real_Study/Forgetting/Basic/S', num2str(i), '/analysis/graph_theory/DF_results_AAL_', num2str(i_beta),'.mat') ;
             
        end
        % ======
        
        save(save_path)
        
    % end should be sooner, must be sooner
        
  
 
    
end
  
% cd '/data/projects/NeFF_B5/scripts_connectivity/'
% disp('Loading AAL data...')
% Y=spm_vol('roi_mnifsl_444.img'); aal=spm_read_vols(Y);
% disp('AAL data loaded...')
% 
% % Before stimulation
% disp('Processing Before left M1 stimulation...');
% 
% for i = 9
%         
%     i
% 
%  %if not(i==17) % para evitar el 4
% 
%     % load functional volumes  -- defino el path e inserto el numero de
%     % sujeto pasandolo a string con num2str y concateno los strings con
%     % strcat
%     data_path = strcat('/data/projects/NeFF_B5/Real_Study/MRI_Session_1/H/S', num2str(i), '/bic_epi_rest_state_physio/fM/preprocessed/preprocessed_444_regression_mov_peaks_erased.img')
%     Y = spm_vol(data_path); data=spm_read_vols(Y);
%     [adj_matrix, c,l,e,Q,Ci, c_rand, l_rand, e_rand, Q_rands, Ci_rand, degree, bet, thr, degree_rand, bet_rand] = network_measures_aal2(data, aal);
%     clearvars -except i c l e Q Ci c_rand l_rand e_rand Q_rands Ci_rand degree bet thr degree_rand bet_rand aal ids num id
%     % defino el directorio que voy a crear (como defini data_path) y lo
%     % creo con mkdir
% %    dir_path = strcat('/data/projects/NeFF_B5/Real_StudyMRI_Session_1/P/P', num2str(i), '/analysis/graph_theory');
% %     mkdir(dir_path);
%     save_path = strcat('/data/projects/NeFF_B5/Real_Study/MRI_Session_1/H/S', num2str(i), '/analysis/graph_theory/results_AAL_RS_reg.mat')
%     save(save_path)
%  %end
%  
%     
% end
