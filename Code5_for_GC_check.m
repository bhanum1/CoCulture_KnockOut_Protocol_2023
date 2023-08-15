clear;
cd /Users/bhanumamillapalli/Desktop/Models/
addpath('/Users/bhanumamillapalli/Desktop/Models/')
addpath(genpath('/Users/bhanumamillapalli/Desktop/CellNetAnalyzer/'))
startcna(1)
% Cplex Path
addpath('/Applications/CPLEX_Studio1210/cplex/matlab/x86-64_osx')

load GeneFiltered_iAF_BS_total
%% Form 2 separate models
eliminated2 = [];
idx.bm = find(ismember(Com.c,1));

cnap_original = cnap;
org1_rxns = find(contains(Com.rxns,'_org1'));
org2_rxns = find(contains(Com.rxns,'_org2'));
share_rxns = find(contains(Com.rxns,'_s_org'));

%Apply KO's

for i = 1:size(mcs,2)
    cnap = cnap_original;
    mcs_test = mcs(:,i);
    cnap.reacMin(mcs_test == -1) = 0;
    cnap.reacMax(mcs_test == -1) = 0;

    Org1 = cnap;
    Org2 = cnap;

    Org1.reacMin(org2_rxns) = 0;
    Org1.reacMax(org2_rxns) = 0;

    Org2.reacMin(org1_rxns) = 0;
    Org2.reacMax(org1_rxns) = 0;
    
    % Find what metabolites are coming in with Com FBA
    Com_FBA = CNAoptimizeFlux(cnap,[],[],2,0,0,[],[]);
    Org1_inputs = Com_FBA(intersect(org1_rxns,share_rxns)) < 0;
    Org2_inputs = Com_FBA(intersect(org2_rxns,share_rxns)) < 0;
    
    %Find Ex rxns for these mets
    [mets,~] = find(Com.S(:,share_rxns) == 1);
    mets = Com.mets(mets);
    Ex_rxns = find(contains(Com.rxns,strcat('EX_',mets)));
    
    Org1_Ex_In = Ex_rxns(Org1_inputs);
    Org2_Ex_In = Ex_rxns(Org2_inputs);
    
    % Open all output EX rxns
    Org1.reacMax(Ex_rxns) = 10;
    Org2.reacMax(Ex_rxns) = 10;
    
    % Check if they can survive without any mets
    FBA1 = CNAoptimizeFlux(Org1,[],[],2,0,0,[],[]);
    FBA2 = CNAoptimizeFlux(Org2,[],[],2,0,0,[],[]);
    
    % If they can't survive without additional mets
    if ~sum(FBA1(idx.bm)) && ~sum(FBA2(idx.bm))
        % Check which mets are needed
        
        % Open input EX for the mets coming in
        Org1.reacMin(Org1_Ex_In) = -1;
        Org2.reacMin(Org2_Ex_In) = -1;
        
        % Close each met one at a time and check
        Org1_needed = [];
        Org2_needed = [];
        
        for j = 1:length(Org1_Ex_In)
            test = Org1;
            test.reacMin(Org1_Ex_In(j)) = 0;
            FBA_test = CNAoptimizeFlux(test,[],[],2,0,0,[],[]);
            if ~sum(FBA_test(idx.bm))
                Org1_needed(end+1) = Org1_Ex_In(j);
            end
        end
        
        for k = 1:length(Org2_Ex_In)
            test = Org2;
            test.reacMin(Org2_Ex_In(k)) = 0;
            FBA_test = CNAoptimizeFlux(test,[],[],2,0,0,[],[]);
            if ~sum(FBA_test(idx.bm))
                Org2_needed(end+1) = Org2_Ex_In(k);
            end
        end
        
        % Run individual FBA
        Org1.reacMin(Org1_needed) = -1;
        Org2.reacMin(Org2_needed) = -1;

        FBA1 = CNAoptimizeFlux(Org1,[],[],2,0,0,[],[]);
        FBA2 = CNAoptimizeFlux(Org2,[],[],2,0,0,[],[]);

        % Check GC mets and compare with needed mets of other organism
        GC1 = FBA1(Ex_rxns) > 0;
        GC2 = FBA2(Ex_rxns) > 0;

        % This checks if there are any mets that are required by one,
        % and growth coupled in the other.
        check1 = intersect(Org1_needed,Ex_rxns(GC2));
        check2 = intersect(Org2_needed,Ex_rxns(GC1));

        if isempty(check1) && isempty(check2)
            eliminated2(end+1) = i;
        end
        
    else
        eliminated2(end+1) = i;
    end
end

mcs(:,eliminated2) = [];
cnap = cnap_original;
save ('Final_iAF_BS.mat','mcs', 'Com', 'cnap');
