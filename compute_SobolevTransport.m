%
% compute Sobolev transport distance matrix
%
% Choose:
% (1) typeGG = 'RandLLE' (G_Log) or typeGG = 'RandSLE' (G_Sqrt)
% (2) pp=1 or pp=2 (the p-order parameter of Sobolev transport)
%

clear all
clc

typeGG = 'RandLLE'; % log-linear #edges (G_Log)
% typeGG = 'RandSLE'; % sqrt-linear #edges (G_Sqrt)

dsName = 'twitter';
maxKC = 100;
nSS = 20; % #tree (average for Sobolev)

pp = 1;
% pp = 2;

% DD_SS1, 5, 10, 20
load([dsName '_' num2str(maxKC) '_' typeGG '_Graph.mat']);

randSArray = randperm(nGG);
wwGG = GG.Edges.Weight;
DD_SS = cell(nSS, 1);

runTime_Prep = zeros(nSS, 1);
runTime_Dist = zeros(nSS, 1);

for idSS = 1:nSS

    % ------- FOR EACH S0 (randomly choose) ---------
    s0 = randSArray(idSS);
    
    tic
    disp(['...[' num2str(idSS) '] compute the tree path']);
    % tree path!!!
    [trPP, trDD, trEP] = shortestpathtree(GG, s0, 'OutputForm', 'cell');
    
    disp(['...[' num2str(idSS) '] vector representation for each vertex']);
    
    % ---------------
    % ===For GRAPH===
    % vector representation for each vertex 1 --> nGG
    
    disp('......vector representation for each vertex');
    % length(wwGG): #edges in graph GG (can be reduced into #edges in tree)
    vecGG_VV = zeros(nGG, length(wwGG));
    for ii = 1:nGG % each vertex in graph
        vecGG_VV(ii, trEP{ii}) = 1;
    end
     
    % V2: extract ---> TREE
    sumEdgeVal = sum(vecGG_VV, 1);
    idNZ = find(sumEdgeVal>0);
    vecGG_VV_TR = vecGG_VV(:, idNZ); % spare version of vecGG_VV
    wwGG_TR = wwGG(idNZ);
    
    disp('......vector representation for each distribution');
    % ===For Data===
    % N: #samples (input data)
    % Input: WW, 

    % V2: --> spare version
    XX_SI = zeros(N, length(idNZ));
    
    for ii = 1:N % each distribution
        tmpWW = WW{ii}/sum(WW{ii}); % normalization for weight!!!
        tmpXX = XX_ID{ii};
        
        tmpXX_GG_TR = vecGG_VV_TR(tmpXX, :);
        tmpWW_GG_TR = repmat(tmpWW, 1, length(idNZ));

        tmpWWXX = tmpXX_GG_TR .* tmpWW_GG_TR;
        XX_SI(ii, :) = sum(tmpWWXX, 1);
    end
    runTime_Prep(idSS) = toc;
   
    tic
    % compute the Lp distance matrix
    DD_SS_II = zeros(N, N);
    for ii = 1:(N-1)
        % ii --> (ii+1):N        
        if mod(ii, 20) == 0
            disp(['...' num2str(ii)]);
        end
    
        tmpII_vec = XX_SI(ii, :);
        
        tmpJJ_mat = XX_SI((ii+1):N, :);
        tmpII_mat = repmat(tmpII_vec, N-ii, 1);
        
        tmpAbsDD_mat = abs(tmpII_mat - tmpJJ_mat);
        
        if pp > 1
            tmpPP_AbsDD_mat = tmpAbsDD_mat.^pp;
        else
            tmpPP_AbsDD_mat = tmpAbsDD_mat;
        end
        
        wwGG_TR_mat = repmat(wwGG_TR', N-ii, 1); 
        % --
        tmpWWPP_AbsDD_mat = wwGG_TR_mat .* tmpPP_AbsDD_mat;
        
        tmpPP_DD_vec = sum(tmpWWPP_AbsDD_mat, 2); % sum over rows --> column
        
        if pp > 1
            tmpDD_vec = tmpPP_DD_vec.^(1/pp);
        else
            tmpDD_vec = tmpPP_DD_vec;
        end
        
        DD_SS_II(ii, (ii+1):N) = tmpDD_vec';
        DD_SS_II((ii+1):N, ii) = tmpDD_vec;
    end
    runTime_Dist(idSS) = toc;
    
    % save distance matrix
    DD_SS{idSS} = DD_SS_II;
end

runTime_Prep_Avg = sum(runTime_Prep) / nSS;
runTime_Dist_Avg = sum(runTime_Dist) / nSS;

runTime_Dist_ALL = runTime_Prep + runTime_Dist;
runTime_Dist_ALL_Avg = sum(runTime_Dist_ALL) / nSS;

% Average
tmpNN = [1, 5, 10, 20];
tmpDDSS_Cell = cell(length(tmpNN), 1);

for iiRR = 1:length(tmpNN)
    
    tmpDDSS = zeros(N, N);
    for ii = 1:tmpNN(iiRR)
        tmpDDSS = tmpDDSS + DD_SS{ii};
    end
    tmpDDSS = tmpDDSS / tmpNN(iiRR);
    
    tmpDDSS_Cell{iiRR} = tmpDDSS;
end

DD_SS1 = tmpDDSS_Cell{1};
DD_SS5 = tmpDDSS_Cell{2};
DD_SS10 = tmpDDSS_Cell{3};
DD_SS20 = tmpDDSS_Cell{4};

outName = [dsName '_Sobolev_V2_' num2str(maxKC) '_' typeGG '_S' num2str(nSS) 'P' num2str(pp) '.mat'];
    
save(outName, 'DD_SS1', 'DD_SS5', 'DD_SS10', 'DD_SS20', ...
     'runTime_Dist', 'runTime_Prep', 'runTime_Dist_ALL', ...
     'runTime_Dist_Avg', 'runTime_Prep_Avg', 'runTime_Dist_ALL_Avg', ...
     'randSArray', 'pp', 'nSS', ...
     'YY');

disp('FINISH !!!');



