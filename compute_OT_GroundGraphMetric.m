%
% compute optimal transport distance matrix (with ground graph metric)
%
% Choose:
% (1) typeGG = 'RandLLE' (G_Log) or typeGG = 'RandSLE' (G_Sqrt)
%

clear all
clc

typeGG = 'RandLLE'; % log-linear #edges (G_Log)
% typeGG = 'RandSLE'; % sqrt-linear #edges (G_Sqrt)

dsName = 'twitter';
maxKC = 100;

load([dsName '_' num2str(maxKC) '_' typeGG '_Graph.mat']);

GM = zeros(nGG, nGG);
disp('compute the ground graph metric');
tic
for ii = 1:(nGG-1)
    [~, TRD_II, ~] = shortestpathtree(GG, ii, [(ii+1):nGG], 'OutputForm', 'cell');  
    GM(ii, (ii+1):nGG) = TRD_II; 
    GM((ii+1):nGG, ii) = TRD_II';
end
runTime_GroundGM = toc;
 
% histogram
XX_ID_vec = zeros(N, nGG);
tic
for ii = 1:N
    % WW{ii}
    tmpWW = WW{ii}/sum(WW{ii}); % normalization
    tmpXX_ID = XX_ID{ii};
    
    XX_ID_vec(ii, tmpXX_ID) = tmpWW'; 
end
runTime_Hist = toc;

tic
% compute the OT
DD_OT = zeros(N, N);
for ii = 1:(N-1)
    
    if mod(ii, 20) == 0
        disp(['...' num2str(ii)]);
    end

    for jj = (ii+1):N
        % preprocessing
        tmpALL = XX_ID_vec(ii, :) + XX_ID_vec(jj, :); 
        idNZ = find(tmpALL > 0);
        
        tmpII = XX_ID_vec(ii, idNZ);
        tmpJJ = XX_ID_vec(jj, idNZ);
        GMIJ = GM(idNZ, idNZ);
        
        DD_OT(ii, jj) = mexEMD(tmpII', tmpJJ', GMIJ);
        DD_OT(jj, ii) = DD_OT(ii, jj);
    end
end
runTime_Dist = toc;

runTime_Dist_ALL = runTime_Dist + runTime_GroundGM + runTime_Hist;

outName = [dsName '_OT_' num2str(maxKC) '_' typeGG '.mat'];

save(outName, 'DD_OT', ...
     'runTime_Dist', 'runTime_GroundGM', 'runTime_Hist', ...
     'runTime_Dist_ALL', ...
     'YY');
    

disp('FINISH !!!');

