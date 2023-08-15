clear;
cd /Users/bhanumamillapalli/Desktop/Models/
addpath('/Users/bhanumamillapalli/Desktop/Models/')
addpath(genpath('/Users/bhanumamillapalli/Desktop/CellNetAnalyzer/'))
startcna(1)
% Cplex Path
addpath('/Applications/CPLEX_Studio1210/cplex/matlab/x86-64_osx')

load MCS_iAF_iAF

%% Build Rxn Gene Matrix
idx.bm = find(ismember(Com.c,1));
% Change this to prevent issues with genes
isSameOrganism = true;

Com.rxnGeneMat = zeros(length(Com.rxns),length(Com.genes));

org1_rxns = find(contains(Com.rxns,'_org1'));
org2_rxns = find(contains(Com.rxns,'_org2'));

if isSameOrganism
    dup = find(ismember(Com.genes,Com.genes(1)));
    org1_genes = 1:dup(2)-1;
    org2_genes = length(org1_genes)+1:size(Com.rxnGeneMat,2);
else
% If the organisms aren't the same, their genes will have different names
% so we don't need to differentiate between the 2 while making the gene
% matrix
    org1_genes = 1:size(Com.genes);
    org2_genes = org1_genes;
end
% Organism 1
for i = org1_rxns
    for j = org1_genes 
        Com.rxnGeneMat(i,j) = contains(Com.grRules(i),Com.genes(j));
    end
end

% Organism 2
for i = org2_rxns
    for j = org2_genes  
        Com.rxnGeneMat(i,j) = contains(Com.grRules(i),Com.genes(j));
    end
end

% Share Reactions
share_rxns = find(contains(Com.rxns,'_s_org'));
share_genes = find(contains(Com.genes,Com.grRules(share_rxns)));

for i = share_rxns
    for j = 1:length(share_genes)
        gene = share_genes(j);
        Com.rxnGeneMat(i,gene) = contains(Com.grRules(i),Com.genes(gene));
    end
 end 

%% Basic Variables
idx.bm = find(ismember(Com.c,1));
eliminated = [];
cnap_original = cnap;

%% Things to test
% Genes
% Whether there is actually codependence via transfer reactions

%% Test if the genes work
for i = 1:size(mcs,2)
    cnap = cnap_original;
    
    %  Find what rxns are KO'd
    mcs_test =  mcs(:,i);
    rxnKOs = find(ismember(mcs_test,-1));
    
    % Find which genes must be KO'd for this
    [~,geneKOs] = find(Com.rxnGeneMat(rxnKOs,:));
    geneKOs = unique(geneKOs);
    
    % Apply those gene KO's and find what reactions get KO'd
    KOrxnGeneMat = Com.rxnGeneMat;
    KOrxnGeneMat(:,geneKOs) = 0;
    
    % Find which reactions already had no genes
    non_targetable = find(cellfun(@isempty,Com.grRules));
    
    % Find which reactions now have no genes
    rxn_genes = sum(full(KOrxnGeneMat),2);
    actual_KOs = find(ismember(rxn_genes,0));
    
    % Determine which ones will end up KO'd
    actual_KOs = setdiff(actual_KOs,non_targetable);
  
    % Apply the KO's
    cnap.reacMin(actual_KOs) = 0;
    cnap.reacMax(actual_KOs) = 0;
    
    % Check for growth
    FBA = CNAoptimizeFlux(cnap,[],[],2,0,0,[],[]); 
    
    % Get rid of options that don't work
    if sum(FBA(idx.bm)) == 0
        eliminated(end+1) = i;
    else
        % GET THE ACTUAL KO's
        mcs(:,i) = zeros(1,size(mcs,1));
        mcs(actual_KOs,i) = -1;
    end
    
    if rem(i,1000) == 0
        i
    end
end

mcs(:,eliminated) = [];
cnap = cnap_original;

save ('GeneFiltered_iAF_iAF.mat','mcs', 'Com', 'cnap');


clear ans dup i j non_targetable org1_genes org1_rxns org2_rxns


    
