%% This code makes a community model and rename the reactions and metabolites in alphabetical order to ease search
% The output is a combined matrix where nothing can be exchanged between
% the two organisms
%clear;
cd /Users/bhanumamillapalli/Desktop/Models/
addpath('/Users/bhanumamillapalli/Desktop/Models/')
addpath(genpath('/Users/bhanumamillapalli/Desktop/CellNetAnalyzer/'))
% Cplex Path
addpath('/Applications/CPLEX_Studio1210/cplex/matlab/x86-64_osx')
%% Making the two models 

% Enter the model of your first organism here and rename
    load 'iAF1260'
    model1=iAF1260;
    Ratio1=0.5;
    
% Enter the model of your second organism here and rename
    load 'iJN1463';
    model2=iJN1463;
    Ratio2=1-Ratio1;
    
%% Rename models to fit standard notation
%find if there is [e]
  [n,~]=size(model1.mets(contains(model1.mets,'[e]'))); 
%if there is then change it
  if n>0
    model1.mets=strrep(model1.mets,'_','__');
    model1.mets=strrep(model1.mets,'-','__');
    model1.rxns=strrep(model1.rxns,'-','__');
    model1.mets=strrep(model1.mets,'[e]','_e');
    model1.mets=strrep(model1.mets,'[c]','_c');
    model1.rxns=strrep(model1.rxns,'(e)','_e');
    model1.rxns=strrep(model1.rxns,'EX_glc_e','EX_glc__D_e');
  end
%model 2
  [n,~]=size(model2.mets(contains(model2.mets,'[e]'))); 
%if there is then change it
  if n>0
    model2.mets=strrep(model2.mets,'_','__');
    model2.mets=strrep(model2.mets,'-','__');
    model2.rxns=strrep(model2.rxns,'-','__');
    model2.mets=strrep(model2.mets,'[e]','_e');
    model2.mets=strrep(model2.mets,'[c]','_c');
    model2.rxns=strrep(model2.rxns,'(e)','_e');
    model2.rxns=strrep(model2.rxns,'EX_glc_e','EX_glc__D_e');
  end
% Sometimes, there is an s at MetCharge, sometimes not
if isfield(model1,'metCharges') == 1
    model1.metCharge=model1.metCharges;
end
if isfield(model2,'metCharges') == 1
    model2.metCharge=model2.metCharges;
end
% Sometimes, grRules is just rules
if isfield(model1,'rules') == 1
    model1.grRules=model1.rules;
end
if isfield(model2,'rules') == 1
    model2.grRules=model2.rules;
end
% Sometimes it does not exist but do not worry
if isfield(model1,'metCharge') == 0
    model1.metCharge=zeros(size(model1.mets,1),1);
end
if isfield(model2,'metCharge') == 0
    model2.metCharge=zeros(size(model2.mets,1),1);
end
%% Order models in alphabetical order

[model1.mets,sortMets]=sort(model1.mets);
model1.metCharge=model1.metCharge(sortMets);
model1.metNames=model1.metNames(sortMets);
model1.metFormulas=model1.metFormulas(sortMets);
model1.S=model1.S(sortMets,1:end);
[model1.rxns,sortRxns]=sort(model1.rxns);
model1.S=model1.S(1:end,sortRxns);
if isfield(model1,'rxnGeneMat') == 1
    if isfield(model2,'rxnGeneMat') == 1
        model1.rxnGeneMat=model1.rxnGeneMat(sortRxns,1:end);
    end
end
model1.grRules=model1.grRules(sortRxns);
model1.rxnNames=model1.rxnNames(sortRxns);
model1.lb=model1.lb(sortRxns);
model1.ub=model1.ub(sortRxns);
model1.c=model1.c(sortRxns);
if isfield(model1,'rev') == 1
    if isfield(model2,'rev') == 1
        model1.rev=model1.rev(sortRxns);
    end
end
[model2.mets,sortMets]=sort(model2.mets);
model2.metCharge=model2.metCharge(sortMets);
model2.metNames=model2.metNames(sortMets);
model2.metFormulas=model2.metFormulas(sortMets);
model2.S=model2.S(sortMets,1:end);
[model2.rxns,sortRxns]=sort(model2.rxns);
model2.S=model2.S(1:end,sortRxns);
if isfield(model1,'rxnGeneMat') == 1
    if isfield(model2,'rxnGeneMat') == 1
        model2.rxnGeneMat=model2.rxnGeneMat(sortRxns,1:end);
    end
end
model2.grRules=model2.grRules(sortRxns);
model2.rxnNames=model2.rxnNames(sortRxns);
model2.lb=model2.lb(sortRxns);
model2.ub=model2.ub(sortRxns);
model2.c=model2.c(sortRxns);
if isfield(model1,'rev') == 1
    if isfield(model2,'rev') == 1
        model2.rev=model2.rev(sortRxns);
    end
end

%% limit ub and lb to have the same
    Max = 100;
    model1.lb(model1.lb<-Max)=-Max;
    model1.ub(model1.ub>Max)=Max;
    model2.lb(model2.lb<-Max)=-Max;
    model2.ub(model2.ub>Max)=Max;

%% identify EX reaction and metabolites
%Find where are the exchanged metabolites and if their name is alright
for n1 = 1:size(model1.S,2)
    nbz(n1) = sum(model1.S(:,n1)~=0); 
    if nbz(n1) == 1 & isempty(find(contains(model1.rxnNames(n1),"exchange"))) ~= 1;
        EX(n1) = 1; % there will be a 1 if this rxns is an exchange
        emets(n1) = find(model1.S(:,n1)~=0); %location of ex_metabolites
    else
        EX(n1) = 0;
        emets(n1) = 0;
    end
end
emets_org1 = emets(find(~ismember(emets,0)));
EX_org1 = find(ismember(EX,1));
if isempty(find(~contains(model1.mets(emets_org1),'_e'))) ~=1
    eloc_org1 = find(~contains(model1.mets(emets_org1),'_e'));  
    extra_org1 = find(contains(model1.metNames(emets_org1(eloc_org1)),'extracellular'));
    model1.mets(emets_org1(eloc_org1(extra_org1))) = strcat(model1.mets(emets_org1(eloc_org1(extra_org1))),'_e');
end
model1.rxns(EX_org1(find(~contains(model1.rxns(EX_org1),'EX'))))=...
    strcat(strcat('EX_',model1.rxns(EX_org1(find(~contains(model1.rxns(EX_org1),'EX'))))),'_e');

% Identify and change structure
    pattern_mets ="_e";
    Exit_mets = contains(model1.mets,pattern_mets);
    A=find(ismember(Exit_mets,1)); 
    pattern_rxns ="EX_";
    Exit_rxns=contains(model1.rxns,pattern_rxns);
    B=find(ismember(Exit_rxns,1));
    [EX_mets_m,~]=size(A);
    [EX_rxns_m,~]=size(B);
%% Manipulate the matrix to have external metabolites first
    [m,n]=size(model1.S);
    C=find(ismember(Exit_mets,0)); % get all the internal and transport
    D=find(ismember(Exit_rxns,0)); % get all the internal and transport
    S1=sparse(m,n);
    S1(1:EX_mets_m,1:n)=model1.S(A,1:n);
    S1(EX_mets_m+1:m,1:n)=model1.S(C,1:n);
    S2=sparse(m,n);
    S2(1:m,1:EX_rxns_m)=S1(1:m,B);
    S2(1:m,EX_rxns_m+1:n)=S1(1:m,D);
    
%% Manipulate the rest of the variables in the model to do the same
% Now change order
Org1.mets(1:EX_mets_m,1)=model1.mets(A);
Org1.mets(EX_mets_m+1:m,1)=model1.mets(C);
Org1.metNames(1:EX_mets_m,1)=model1.metNames(A);
Org1.metNames(EX_mets_m+1:m,1)=model1.metNames(C);
Org1.metFormulas(1:EX_mets_m,1)=model1.metFormulas(A);
Org1.metFormulas(EX_mets_m+1:m,1)=model1.metFormulas(C);
Org1.metCharge(1:EX_mets_m,1)=model1.metCharge(A);
Org1.metCharge(EX_mets_m+1:m,1)=model1.metCharge(C);
Org1.genes=model1.genes;
if isfield(model1,'rxnGeneMat') == 1
    if isfield(model2,'rxnGeneMat') == 1
        Org1.rxnGeneMat=model1.rxnGeneMat;
        Org1.rxnGeneMat(1:EX_rxns_m,1:end)=model1.rxnGeneMat(B,1:end);
        Org1.rxnGeneMat(EX_rxns_m+1:n,1:end)=model1.rxnGeneMat(D,1:end);
    end
end
Org1.grRules(1:EX_rxns_m,1)=model1.grRules(B);
Org1.grRules(EX_rxns_m+1:n,1)=model1.grRules(D);
Org1.rxns(1:EX_rxns_m,1)=model1.rxns(B);
Org1.rxns(EX_rxns_m+1:n,1)=model1.rxns(D);
Org1.rxnNames(1:EX_rxns_m,1)=model1.rxnNames(B);
Org1.rxnNames(EX_rxns_m+1:n,1)=model1.rxnNames(D);
Org1.S=S2;
clear S1
clear S2
Org1.lb(1:EX_rxns_m,1)=model1.lb(B);
Org1.lb(EX_rxns_m+1:n,1)=model1.lb(D);
Org1.ub(1:EX_rxns_m,1)=model1.ub(B);
Org1.ub(EX_rxns_m+1:n,1)=model1.ub(D);
Org1.c(1:EX_rxns_m,1)=model1.c(B);
Org1.c(EX_rxns_m+1:n,1)=model1.c(D);
Org1.b=zeros(m,1);
if isfield(model1,'rev') == 1
    if isfield(model2,'rev') == 1
        Org1.rev(1:EX_rxns_m,1)=model1.rev(B);
        Org1.rev(EX_rxns_m+1:n,1)=model1.rev(D);
    end
end
Org1.description=model1.description;

%% Second organism
clearvars EX emets n1 nbz
for n1 = 1:size(model2.S,2)
    nbz(n1) = sum(model2.S(:,n1)~=0); 
    if nbz(n1) == 1 & isempty(find(contains(model2.rxnNames(n1),"exchange"))) ~= 1;
        EX(n1) = 1; % there will be a 1 if this rxns is an exchange
        emets(n1) = find(model2.S(:,n1)~=0); %location of ex_metabolites
    else
        EX(n1) = 0;
        emets(n1) = 0;
    end
end
emets_org2 = emets(find(~ismember(emets,0)));
EX_org2 = find(ismember(EX,1));
% Add _e to metabolites
if isempty(find(~contains(model2.mets(emets_org2),'_e'))) ~=1
    eloc_org2 = find(~contains(model2.mets(emets_org2),'_e'));  
    extra_org2 = find(contains(model2.metNames(emets_org2(eloc_org2)),'extracellular'));
    model2.mets(emets_org2(eloc_org2(extra_org2))) = strcat(model2.mets(emets_org2(eloc_org2(extra_org2))),'_e');
end
% Add EX_ to rxns
model2.rxns(EX_org2(find(~contains(model2.rxns(EX_org2),'EX'))))=...
    strcat(strcat('EX_',model2.rxns(EX_org2(find(~contains(model2.rxns(EX_org2),'EX'))))),'_e');
%% identify EX reaction and metabolites

pattern_mets ="_e";
Exit_mets = contains(model2.mets,pattern_mets);
A=find(ismember(Exit_mets,1)); 
pattern_rxns ="EX_";
Exit_rxns=contains(model2.rxns,pattern_rxns);
B=find(ismember(Exit_rxns,1));
[EX2_mets_m,~]=size(A);
[EX2_rxns_m,~]=size(B);

%% Manipulate the matrix to have external metabolites first
[m2,n2]=size(model2.S);
C=find(ismember(Exit_mets,0)); % get all the internal and transport
D=find(ismember(Exit_rxns,0)); % get all the internal and transport
S1=sparse(m2,n2);
S1(1:EX2_mets_m,1:n2)=model2.S(A,1:n2);
S1(EX2_mets_m+1:m2,1:n2)=model2.S(C,1:n2);
S2=sparse(m2,n2);
S2(1:m2,1:EX2_rxns_m)=S1(1:m2,B);
S2(1:m2,EX2_rxns_m+1:n2)=S1(1:m2,D);

%% Manipulate the rest of the variables in the model2 to do the same

Org2.mets(1:EX2_mets_m,1)=model2.mets(A);
Org2.mets(EX2_mets_m+1:m2,1)=model2.mets(C);
Org2.metNames(1:EX2_mets_m,1)=model2.metNames(A);
Org2.metNames(EX2_mets_m+1:m2,1)=model2.metNames(C);
Org2.metFormulas(1:EX2_mets_m,1)=model2.metFormulas(A);
Org2.metFormulas(EX2_mets_m+1:m2,1)=model2.metFormulas(C);
Org2.metCharge(1:EX2_mets_m,1)=model2.metCharge(A);
Org2.metCharge(EX2_mets_m+1:m2,1)=model2.metCharge(C);
Org2.genes=model2.genes;
if isfield(model1,'rxnGeneMat') == 1
    if isfield(model2,'rxnGeneMat') == 1
        Org2.rxnGeneMat=model2.rxnGeneMat;
        Org2.rxnGeneMat(1:EX2_rxns_m,1:end)=model2.rxnGeneMat(B,1:end);
        Org2.rxnGeneMat(EX2_rxns_m+1:n2,1:end)=model2.rxnGeneMat(D,1:end);
    end
end
Org2.grRules(1:EX2_rxns_m,1)=model2.grRules(B);
Org2.grRules(EX2_rxns_m+1:n2,1)=model2.grRules(D);
Org2.rxns(1:EX2_rxns_m,1)=model2.rxns(B);
Org2.rxns(EX2_rxns_m+1:n2,1)=model2.rxns(D);
Org2.rxnNames(1:EX2_rxns_m,1)=model2.rxnNames(B);
Org2.rxnNames(EX2_rxns_m+1:n2,1)=model2.rxnNames(D);
Org2.S=S2;
clear S1
clear S2
Org2.lb(1:EX2_rxns_m,1)=model2.lb(B);
Org2.lb(EX2_rxns_m+1:n2,1)=model2.lb(D);
Org2.ub(1:EX2_rxns_m,1)=model2.ub(B);
Org2.ub(EX2_rxns_m+1:n2,1)=model2.ub(D);
Org2.b=zeros(m2,1);
Org2.c(1:EX2_rxns_m,1)=model2.c(B);
Org2.c(EX2_rxns_m+1:n2,1)=model2.c(D);
if isfield(model1,'rev') == 1
    if isfield(model2,'rev') == 1
        Org2.rev(1:EX2_rxns_m,1)=model2.rev(B);
        Org2.rev(EX2_rxns_m+1:n2,1)=model2.rev(D);
    end
end
Org2.description=model2.description;

clearvars -except Org1 Org2

%% Make the matrix
S1 = zeros(size(Org1.S,1)+size(Org2.S,1),size(Org1.S,2)+size(Org2.S,2));
S1(1:size(Org1.S,1),1:size(Org1.S,2)) = Org1.S;
S1(size(Org1.S,1)+1:end,size(Org1.S,2)+1:end) = Org2.S;
%% Make the model
Com.mets(1:size(Org1.S,1),1) = strcat(Org1.mets,'_org1');
Com.mets(size(Org1.S,1)+1:size(S1,1)) = strcat(Org2.mets,'_org2');
Com.metNames(1:size(Org1.S,1),1) = strcat(Org1.metNames,'_org1');
Com.metNames(size(Org1.S,1)+1:size(S1,1)) = strcat(Org2.metNames,'_org2');
Com.metFormulas(1:size(Org1.S,1),1) = Org1.metFormulas;
Com.metFormulas(size(Org1.S,1)+1:size(S1,1)) = Org2.metFormulas;
Com.metCharge(1:size(Org1.S,1),1) = Org1.metCharge;
Com.metCharge(size(Org1.S,1)+1:size(S1,1)) = Org2.metCharge;
Com.genes(1:size(Org1.genes,1),1) = Org1.genes;
Com.genes(size(Org1.genes,1)+1:size(Org1.genes,1)+size(Org2.genes,1),1) =...
    Org2.genes;

Com.grRules(1:size(Org1.S,2),1) = Org1.grRules;
Com.grRules(size(Org1.S,2)+1:size(S1,2)) = Org2.grRules;
if isempty(find(contains(Com.grRules,'|'))) == 0
    Com.grRules = strrep(Com.grRules,'|','or');
    Com.grRules = strrep(Com.grRules,'&','and');
end
Com.rxns(1:size(Org1.S,2),1) = strcat(Org1.rxns,'_org1');
Com.rxns(size(Org1.S,2)+1:size(S1,2)) = strcat(Org2.rxns,'_org2');
Com.rxnNames(1:size(Org1.S,2),1) = strcat(Org1.rxnNames,'_org1');
Com.rxnNames(size(Org1.S,2)+1:size(S1,2)) = strcat(Org2.rxnNames,'_org2');
Com.S = sparse(S1);
Com.lb(1:size(Org1.S,2),1) = Org1.lb;
Com.lb(size(Org1.S,2)+1:size(S1,2)) = Org2.lb;
Com.ub(1:size(Org1.S,2),1) = Org1.ub;
Com.ub(size(Org1.S,2)+1:size(S1,2)) = Org2.ub;
Com.c(1:size(Org1.S,2),1) = Org1.c;
Com.c(size(Org1.S,2)+1:size(S1,2)) = Org2.c;
Com.b(1:size(Org1.S,1),1) = Org1.b;
Com.b(size(Org1.S,1)+1:size(S1,1)) = Org2.b;
Com.rev(1:size(Org1.S,2),1) = Org1.rev;
Com.rev(size(Org1.S,2)+1:size(S1,2)) = Org2.rev;
Com.description = strcat(strcat(Org1.description,'_&_'),Org2.description);

save ('4MCSCom.mat','Com')

