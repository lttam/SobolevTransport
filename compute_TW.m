%
% compute tree-Wasserstein distance matrix (with random tree metric from graph structure)
%
% Choose:
% (1) typeGG = 'RandLLE' (G_Log) or typeGG = 'RandSLE' (G_Sqrt)
%

clear all
clc

typeGG = 'RandLLE'; % log-linear #edges
% typeGG = 'RandSLE'; % sqrt-linear #edges

dsName = 'twitter';
maxKC = 100;
nSS = 20; % #tree (average for TW)

load([dsName '_' num2str(maxKC) '_' typeGG '_Graph.mat']);

DD_SS = cell(nSS, 1);

runTime_Prep = zeros(nSS, 1);
runTime_Dist = zeros(nSS, 1);

for idSS = 1:nSS

    disp('... randomly sampling tree');
    [minTR_EdgeID, minTR_EdgeWW] = RandomlySamplingTree(nGG, GG.Edges.EndNodes, GG.Edges.Weight);

    minTRGG = graph(minTR_EdgeID(:, 1), minTR_EdgeID(:, 2), minTR_EdgeWW);

    s0=1;

    tic
    disp(['... compute the tree path']);
    % tree path!!!
    [trPP, trDD, trEP] = shortestpathtree(minTRGG, s0, 'OutputForm', 'cell');

    disp(['...vector representation for each vertex']);

    % ---------------
    % ===For TREE===
    % vector representation for each vertex 1 --> nGG

    disp('......vector representation for each vertex');
    % length(wwGG): #edges in graph GG (can be reduced into #edges in tree)
    vecGG_VV = zeros(nGG, length(minTR_EdgeWW));
    for ii = 1:nGG % each vertex in graph/tree
        vecGG_VV(ii, trEP{ii}) = 1;
    end

    disp('......vector representation for each distribution');
    % ===For Data===
    % N: #samples (input data)
    XX_SI = zeros(N, length(minTR_EdgeWW));
    for ii = 1:N % each distribution
        tmpWW = WW{ii}/sum(WW{ii}); % normalization for weight!!!
        tmpXX = XX_ID{ii};

        tmpXX_GG = vecGG_VV(tmpXX, :);
        tmpWW_GG = repmat(tmpWW, 1, length(minTR_EdgeWW));

        tmpWWXX = tmpXX_GG .* tmpWW_GG;
        XX_SI(ii, :) = sum(tmpWWXX, 1);
    end
    runTime_Prep_II = toc;


    tic
    % compute the Lp distance matrix
    DD_TW = zeros(N, N);
    for ii = 1:(N-1)
        % ii --> (ii+1):N
        tmpII_vec = XX_SI(ii, :);

        tmpJJ_mat = XX_SI((ii+1):N, :);
        tmpII_mat = repmat(tmpII_vec, N-ii, 1);

        tmpAbsDD_mat = abs(tmpII_mat - tmpJJ_mat);
        wwGG_mat = repmat(minTR_EdgeWW', N-ii, 1); 
        tmpWW_AbsDD_mat = wwGG_mat .* tmpAbsDD_mat;

        tmpDD_vec = sum(tmpWW_AbsDD_mat, 2); % sum over rows --> column

        DD_TW(ii, (ii+1):N) = tmpDD_vec';
        DD_TW((ii+1):N, ii) = tmpDD_vec;
    end
    runTime_Dist_II = toc;

    % saving !!!
    runTime_Dist(idSS) = runTime_Dist_II;
    runTime_Prep(idSS) = runTime_Dist_II;
    DD_SS{idSS} = DD_TW;
    
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

DD_TW1 = tmpDDSS_Cell{1};
DD_TW5 = tmpDDSS_Cell{2};
DD_TW10 = tmpDDSS_Cell{3};
DD_TW20 = tmpDDSS_Cell{4};

DD_TW = DD_SS;

outName = [dsName '_TW_Random_' num2str(maxKC) '_' typeGG '_S' num2str(nSS) '.mat'];
    
save(outName, 'DD_TW1', 'DD_TW5', 'DD_TW10', 'DD_TW20', ...
     'runTime_Dist', 'runTime_Prep', 'runTime_Dist_ALL', ...
     'runTime_Dist_Avg', 'runTime_Prep_Avg', 'runTime_Dist_ALL_Avg', ...
     'nSS', ...
     'YY');

disp('FINISH !!!');


