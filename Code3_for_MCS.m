clear;
cd /Users/bhanumamillapalli/Desktop/Models/
addpath('/Users/bhanumamillapalli/Desktop/Models/')
addpath(genpath('/Users/bhanumamillapalli/Desktop/CellNetAnalyzer/'))
startcna(1)
% Cplex Path
addpath('/Applications/CPLEX_Studio1210/cplex/matlab/x86-64_osx')

load MCS_Com_iAF_PP

%% Make a list of all Ex_reactions that remains open
idx.bm = find(ismember(Com.c,1));
idx.o2 = find(contains(Com.rxns,{'EX_o2_e'}));
% Setting the Cplex problem
    Cplex_Com.Aeq=Com.S;
    Cplex_Com.lb=Com.lb;
    Cplex_Com.ub=Com.ub;
    Cplex_Com.f=-Com.c;
    Cplex_Com.beq=Com.b;
% Run FBA with both organisms growth as objective
    Cplex_Com.f(idx.bm(2)) = 0;
    FBA_org1=cplexlp(Cplex_Com.f,[],[],Cplex_Com.Aeq,Cplex_Com.beq,Cplex_Com.lb,Cplex_Com.ub);
    Cplex_Com.f(idx.bm(2)) = -1;
    Cplex_Com.f(idx.bm(1)) = 0;
    FBA_org2=cplexlp(Cplex_Com.f,[],[],Cplex_Com.Aeq,Cplex_Com.beq,Cplex_Com.lb,Cplex_Com.ub);
% Set and run pFBA
    Nogenes=cellfun(@isempty,Com.grRules);
    Nogenes(find(ismember(Com.grRules,strcat(strcat('x(',num2str(find(contains(Com.genes,'shared')))),')')))) = 1; %Do not include the artificial genes
    Genes=Nogenes==0;
    Cplex_Com.lb(idx.bm) = FBA_org1(idx.bm);
    Cplex_Com.f(:,1) = 0;
    Cplex_Com.f(Genes,1)=1;
    pFBA_org1=cplexlp(Cplex_Com.f,[],[],Cplex_Com.Aeq,Cplex_Com.beq,Cplex_Com.lb,Cplex_Com.ub);
    Cplex_Com.lb(idx.bm) = FBA_org2(idx.bm);
    pFBA_org2=cplexlp(Cplex_Com.f,[],[],Cplex_Com.Aeq,Cplex_Com.beq,Cplex_Com.lb,Cplex_Com.ub);
    
    % Find all rxns that are always open in both cases
    EX_rxns = find(startsWith(Com.rxns,'EX_'));
    EX_Open_aero = [Com.rxns(EX_rxns(find(pFBA_org1(EX_rxns)~=0)));...
        Com.rxns(EX_rxns(find(pFBA_org2(EX_rxns)~=0)))];
    EX_Open_FBA = [Com.rxns(EX_rxns(find(FBA_org1(EX_rxns)~=0)));...
        Com.rxns(EX_rxns(find(FBA_org2(EX_rxns)~=0)))];
    Keep = [EX_Open_aero;...
        Com.rxns(find(ismember(Com.rxns,strrep(strrep(Com.rxns(find(startsWith(Com.rxns,'S_'))),'S_','EX_'),'_s','_e_org1'))));...
        Com.rxns(find(ismember(Com.rxns,strrep(strrep(Com.rxns(find(startsWith(Com.rxns,'S_'))),'S_','EX_'),'_s','_e_org2'))))];
  
%% Make the model

%if ~exist('cnan','var')
%    startcna(1);
%    wait(parfevalOnAll(@startcna,0,1));
%end
max = 100;
Com.lb(EX_rxns) = 0;
Com.ub(EX_rxns) = 0;
Com.lb(find(ismember(Com.rxns,Keep))) = Cplex_Com.lb(find(ismember(Com.rxns,Keep)));
Com.ub(find(ismember(Com.rxns,Keep))) = Cplex_Com.ub(find(ismember(Com.rxns,Keep)));
% Com.ub(contains(Com.rxns,'EX_') & contains(Com.rxns,'_s') & ~contains(Com.rxns,'_e_org')) = 100;
Com.lb(find(contains(Com.rxns,'tr_s_'))) = -100;
cnap = CNAcobra2cna(Com);

% Add gprRules from Cobra to CNA
cnap = CNAsetGenericReactionData_with_array(cnap,'geneProductAssociation',Com.grRules);

cnap.reacMin(cnap.reacMin==-max) = -inf;
cnap.reacMax(cnap.reacMax== max) =  inf;

idx.glc  = find(ismember(cnap.reacID,{'EX_glc__D_e_com'}));
idx.atpm_org1 = find(ismember(cnap.reacID,{'ATPM_org1'}));
idx.atpm_org2 = find(ismember(cnap.reacID,{'ATPM_org2'}));

cnap.reacMax(idx.atpm_org1) = inf;
cnap.reacMax(idx.atpm_org2) = inf;

% Uncomment to test feasibility 
% FBA = CNAoptimizeFlux(cnap,[],[],2,0,0,[],[]); 
% FBA(idx.bm)
% CNAplotPhasePlane(cnap,idx.bm)

%% Set up cMCS input
modules{1}.sense = 'target';
modules{1}.type = 'lin_constraints';
modules{1}.V = zeros(2,cnap.numr);
modules{1}.V(1,idx.bm(1)) = 1;
modules{1}.V(2,idx.bm(2)) = -1;
modules{1}.v = zeros(2,1);
modules{1}.v(2,1) = -0.05;
% modules{1}.v(1,1) = 0.01;

modules{2}.sense = 'target';
modules{2}.type = 'lin_constraints';
modules{2}.V = zeros(2,cnap.numr);
modules{2}.V(1,idx.bm(2)) = 1;
modules{2}.V(2,idx.bm(1)) = -1;
modules{2}.v = zeros(2,1);
modules{2}.v(2,1) = -0.05;
% modules{1}.v(1,1) = 0.01;

modules{3}.sense = 'desired';
modules{3}.type = 'lin_constraints';
modules{3}.V = zeros(1,cnap.numr);
modules{3}.V(1,idx.bm(1)) = -1;
modules{3}.v = -0.05;

modules{4}.sense = 'desired';
modules{4}.type = 'lin_constraints';
modules{4}.V = zeros(1,cnap.numr);
modules{4}.V(1,idx.bm(2)) = -1;
modules{4}.v = -0.05;

% modules{5}.sense = 'target';
% modules{5}.type = 'lin_constraints';
% modules{5}.V = zeros(2,cnap.numr);
% modules{5}.V(1,find(ismember(Com.rxns,'AKGtr_s_org1'))) = 1;
% modules{5}.V(2,find(ismember(Com.rxns,'GLU__Ltr_s_org2'))) = -1;
% modules{5}.v = zeros(2,1);
% modules{5}.v(1,:) = 0.01;
% modules{5}.v(2,:) = -0.01;

%% Rxns should not feasible ko
rkoCost = ones(cnap.numr,1);
rkoCost(find(ismember(Genes,0))) = nan;
rkoCost(find(contains(Com.rxns,'H2Ot_'))) = nan;
rkoCost(find(contains(Com.rxns,'tr_s_'))) = nan;
rkiCost = nan(cnap.numr,1);
% rkiCost(find(contains(Com.rxns,'_s_org'))) = 1;
% specifying gene knockout costs
% Generate list of genes as a template to define k/o genes
% [~,~,genes,gpr_rules] = CNAgenerateGPRrules(cnap);
% gkoCost = ones(numel(genes),1);
% gkiCost = nan(numel(genes),1);
% gkiCost(find(ismember(genes,strcat('x',num2str(find(contains(Com.genes,'shared'))))))) = 0;

maxSolutions = 10;
maxCost = 100;
options.milp_solver = 'cplex';
options.milp_time_limit  = 43200;
% options.milp_bigM = 10000; %Try true, options.seed (change could help)
options.mcs_search_mode = 1;
verbose = 1;

%% Run cMCS and save outcome 

[mcs, status, cmp_mcs, cmp_cnap, mcs_idx_cmp_full, obj] = ...
    CNAMCSEnumerator3(cnap,modules,rkoCost,rkiCost,maxSolutions,maxCost,options,verbose);

% [mcs, status, cmp_mcs, cmp_cnap, mcs_idx_cmp_full, obj] = ...
%     CNAgeneMCSEnççççumerator3(cnap,modules,{},{},maxSolutions,maxCost,gkoCost,gkiCost,{},options,verbose);
Com.keep = Keep;
save ('ComStrategyiAF_PP_test.mat','Com','cnap','mcs');
%save ('MCS_ComCore_test.mat','Com','cnap','idx')  
mcs_test = mcs(:,1);
cnap_original = cnap;
cnap.reacMin(mcs_test == -1) = 0;
cnap.reacMax(mcs_test == -1) = 0;
CNAplotPhasePlane(cnap,idx.bm)
cnap.objFunc(idx.bm(1)) = 0;
FBA = CNAoptimizeFlux(cnap,[],[],2,0,0,[],[]); 
FBA(idx.bm)