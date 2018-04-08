function [adj_matrix , c,l,e,Q,Ci, c_rand, l_rand, e_rand, Q_rands, Ci_rands, degree, bet, thr, degree_rand, bet_rand] = network_measures_aal_newman_Ben(data, template)
c=[];l=[];e=[];Q=[];Ci=[];
c_rand=[];l_rand=[];e_rand=[];Q_rands=[];Ci_rands=[];
degree=[];
degree_rand=[];
bet=[];
bet_rand=[];
thr=[];
thr_rand=[];

% computes network measures for networks constructed in a density range
% between 0.2:0.02:0.50. ***Requires Brain Connectivity Toolbox***

%%% inputs %%%
% data file (4D)
% template = AAL template (or other) to define ROIs (same resolution as data) %

% outputs %
% c = average clustering coefficient %
% l = average shortest path length %
% e = average global efficiency %
% Q = network modularity (Louvain algorithm)
% Ci = module membership
% X_rand = same with degree preserving randomization
randsteps = 100; % steps for randomization
trials_num=100; % trials for the Louvain algorithm

disp('Extracting time series...')
dim=size(data);
time_series = extract_timeseries_aal_Ben(data, template);
%time_series = time_series(:,1:round(size(time_series,2)/2));%%% Borrar....time series first half
%time_series = time_series(:,round(size(time_series,2)/2)+1:end);%%% Borrar....time series second half
dim_time_series=size(time_series);
disp('Computing correlation matrices...')
adj_matrix  = adjacency_matrix_lincorr(time_series, 0.01);

adjdim=size(adj_matrix);
disp('Computing network measures...')
for n = 1:adjdim(3)
    adjdim(3)-n
    adj=reshape(adj_matrix(:,:,n), adjdim(1), adjdim(2) );

    
    %%%randomice networks
    disp('Randomizing...')
    adj_rand =  randomizer_bin_und(adj,randsteps);
    %%%
     
    % centrality measures
    disp('Centrality measures...');
      for node=1:1:dim_time_series(1);
        degree(node,n) =sum(adj(node,:));
        degree_rand(node,n) = sum(adj_rand(node,:));
      end
    %
    b=betweenness_bin( adj); bet(:,n) = b;
    b=betweenness_bin( adj_rand); bet_rand(:,n) = b;
    
    
    
    disp('C,L,E...')
    C=clustering_coef_bu(adj); % clustering coefficient
    D=distance_bin(adj); % distance matrix
    [L,E] = charpath(D); % average path length & efficiency
    c(n) = mean(C);
    l(n) = L;
    e(n) = E;
    
    disp('C,L,E (randomized)...')
    C_rand=clustering_coef_bu(adj_rand); % clustering coefficient
    D_rand=distance_bin(adj_rand); % distance matrix
    [L_rand,E_rand] = charpath(D_rand); % average path length & efficiency
    c_rand(n) = mean(C_rand);
    l_rand(n) = L_rand;
    e_rand(n) = E_rand;
    
    
    
    disp('Community analysis...')                   
    
    mods_1=0;mods_rand=0;
    Ci_final=[];  Ci_final_rand=[];
    % performs hierarchical decomposition into communities,
    % Louvain Algorithm
        for trials=1:trials_num
                        
        [Ci_1, Q_1] = modularity_und( adj ); 
        [Ci_rand,Q_rand] =  modularity_und( adj_rand );
               if Q_1 > mods_1
              mods_1 = Q_1;
              Ci_final_1=Ci_1;
        end
        if Q_rand > mods_rand
              mods_rand = Q_rand;
              Ci_final_rand=Ci_rand;
        end
        
        end
   
     Q(n) = mods_1; 
     
  
     Ci(:,n) = Ci_final_1; 
     
     mods_rand;
     Q_rands(n) = mods_rand;
     
    
     Ci_rands(:,n) = Ci_final_rand';
                    
         
end
    

