set(0,'DefaultTextInterpreter','none');

%% Load Data %%
        load('./Data/Subject_1/fMRI_BOLD_FC.mat');
        data=FC_cc_68;
        roi_names= csvimport('./Data/Subject_1/freesurfer_regions_68_sort_full.txt');
        N=length(data);   % Number of ROIs
        issymmetric(data) % Check if FC is symmetric
%% Clean data %%
        fc=atanh(data);                         % Z-transform original data to avoid false-positives
        fc=weight_conversion(fc,'autofix');     % Set main diagonal to zero, remove any NaN (or Inf values) , correct round-off errors  
        fc=weight_conversion(fc,'normalize');   % rescale the values between 0 to 1 (Got read of negative weights!)
%% Visualize adjacency Matrix %%
        imagesc(fc); colorbar; %Notice negative weights

%% Community Detection %%
         [Ci0 Q0]= community_louvain(fc,[],[],'negative_sym');
         
% Optional-- for better community partition %
%                           for i=1:100
%                                  [Ci0 Q0]= community_louvain(fc,[],[],'negative_sym');
%                                  Ci_iter(i,:)=Ci0;
%                           end
%                           Ci=transpose(Ci_iter);
%                           D= agreement(Ci);
%                           ci0=consensus_und(D,80,100);
 %% Visualize Modular Structure of Adjacency Matrix %%
        figure;
        [X,Y,INDSORT] = grid_communities(Ci0); % call function
        imagesc(fc(INDSORT,INDSORT));           % plot ordered adjacency matrix
        hold on;                                 % hold on to overlay community visualization
        plot(X,Y,'y','linewidth',3);
        title('Community Structure');
        hold off;
 %% Participation Coefficient        
        [Ppos,Pneg]=participation_coef_sign(fc,Ci0);
        [Ppos,I]=sort(Ppos,1,'descend'); 
	       ROIs_sorted=roi_names(I);
           figure;
           bar(Ppos);
          title('Positive Clustering Coefficient');
          set(0,'defaulttextinterpreter','none')
          set(gca, 'XTickLabel',ROIs_sorted, 'XTick',1:numel(ROIs_sorted));
          rotateXLabels(gca(), 90 );
          
          [Pneg,I]=sort(Pneg,1,'descend'); 
	       ROIs_sorted=roi_names(I);
           figure;
           bar(Pneg);
          title('Negative Clustering Coefficient');
          set(0,'defaulttextinterpreter','none')
          set(gca, 'XTickLabel',ROIs_sorted, 'XTick',1:numel(ROIs_sorted));
          rotateXLabels(gca(), 90 );
    
 %% Effect of thresholding & binarizing
    figure;
      for i=0.01:0.01:0.1
        th_fc=threshold_proportional(fc,i);
        bin_fc=weight_conversion(th_fc,'binarize');
        d=sum(bin_fc);
        histogram(d,10); title([i]);
        pause(0.5);
      end  
      
      th_fc=threshold_proportional(fc,0.01);
      bin_fc=weight_conversion(th_fc,'binarize');

 %% Centrality Measures %%
      % Betweenness Centrality %
           b=betweenness_bin(bin_fc);
           norm_b=b./((N-1)*(N-2));
	       [B,I]=sort(norm_b,2,'descend'); 
	       ROIs_sorted=roi_names(I);
           Bplot= bar(B);
          title('Betweenness Centrality');
          set(0,'defaulttextinterpreter','none')
          set(gca, 'XTickLabel',ROIs_sorted, 'XTick',1:numel(ROIs_sorted));
          rotateXLabels(gca(), 90 );
         % setbarcolor(Bplot,B,y);           


      % Eigen Centrality %
          eig=eigenvector_centrality_und(th_fc);
          [E,I]=sort(eig,1,'descend'); 
	      ROIs_sorted=roi_names(I);
          figure;
          Eplot= bar(E);
          title('Eigenvector Centrality');
          set(0,'defaulttextinterpreter','none')
          set(gca, 'XTickLabel',ROIs_sorted, 'XTick',1:numel(ROIs_sorted));
          rotateXLabels( gca(), 90 );

%%%% See tutorial for visualization of network %%%%
  %% Length Matrix - inverse of adj weights
           L=weight_conversion(th_fc,'lengths');
        
  %% Distance Matrix - Shortest path for each node pair
           D= distance_bin(bin_fc);
        
  %% Efficiency 
          Eglob=efficiency_bin(bin_fc);
          Eloc= efficiency_bin(bin_fc,1);
          [SEloc,I]=sort(Eloc,1,'descend'); 
	       ROIs_sorted=roi_names(I);
           Lplot= bar(SEloc);
          title('Local Efficiency');
          set(0,'defaulttextinterpreter','none')
          set(gca, 'XTickLabel',ROIs_sorted, 'XTick',1:numel(ROIs_sorted));
          rotateXLabels(gca(), 90 );
 
  %% Small World Index
         C=   clustering_coef_bu(bin_fc); %Clustering Coefficient
         cpl= charpath(D,0,0); % Characteristic Path Length
         small_worldness= C./cpl;
