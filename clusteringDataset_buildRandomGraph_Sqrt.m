%
% build G_Sqrt graph from supports
%

clear all
clc

dsName = 'twitter';
maxKC = 100;
constEPS = 1e-10;

load([dsName '.mat']);
% N: scalar
% WW: Nx1 cell
% XX: Nx1 cell

% start & end ID for each XX{ii}
XX_sID = zeros(N, 1);
XX_eID = zeros(N, 1);

sID = 0;
for ii = 1:N
    nnII = size(XX{ii}, 1);
    
    XX_sID(ii) = sID + 1;
    XX_eID(ii) = sID + nnII;

    sID = XX_eID(ii);
end
sumNN = sID;
disp(['Sample size: ' num2str(sumNN)]);

XX_ALL = zeros(sumNN, size(XX{1}, 2));

for ii = 1:N
    XX_ALL(XX_sID(ii):XX_eID(ii), :) = XX{ii};
end

[NN, dim] = size(XX_ALL);

disp('Clustering by the farthest-point clustering !!!');
tic
[KC, ~, XX_ALL_ID, XX_CC, ~, ~] = ...
    figtreeKCenterClustering(dim, NN, XX_ALL', maxKC);

% note: XX_ALL_ID: 0 --> (KC-1)
XX_ALL_ID = XX_ALL_ID + 1;
runTime_Clustering = toc

% XX_IdCC --> ID
% XX_CC: value

XX_ALL_ID = XX_ALL_ID'; % NN x 1
XX_CC = XX_CC'; % KC x dim

XX_ID = cell(N, 1);
for ii = 1:N
    tmp = XX_ALL_ID(XX_sID(ii):XX_eID(ii));
    XX_ID{ii} = tmp;
end

% XX_CC: data (cluster centers) --> KC x dim
% XX_ID: cell Nx1 --> each cell: vector of ID (of clusters)

% -----------
% Build graph
nGG = KC;

% --- PRUNNING -- (kk*N)
kGG = ceil(sqrt(nGG));

disp('Compute distance matrix');
tic
sqDD = sqdistance(XX_CC');
runTime_DM = toc

% -- correct machine precision
sqDD(sqDD<constEPS) = 0;
DD = sqrt(sqDD);

% #edges before prunning
nEE = nGG*(nGG-1)/2;
% only kept #edges (after prunning)
mEE = nGG*kGG;

EdgeID = zeros(nEE, 2);
EdgeWW = zeros(nEE, 1);

id = 1;
for ii = 1:(nGG-1)
    for jj = (ii+1):nGG
        EdgeID(id, :) = [ii, jj];
        EdgeWW(id) = DD(ii, jj);
        id = id+1;
    end
end

% save original
EdgeID_org = EdgeID;
EdgeWW_org = EdgeWW;
Key_org = EdgeID_org(:, 1)*(nGG+1) + EdgeID_org(:, 2); % Key for searching !!!

% --- random perturbation----
tmpRand = randperm(length(EdgeWW_org));
EdgeWW = EdgeWW_org(tmpRand);
EdgeID = EdgeID_org(tmpRand, :);

% only keep mEE = kGG*nGG
EdgeWW = EdgeWW(1:mEE);
EdgeID = EdgeID(1:mEE, :);

GG = graph(EdgeID(:, 1), EdgeID(:, 2), EdgeWW);

% Test connected components
disp('check connected components');
tic
[bGG, bsGG] = conncomp(GG, 'OutputForm','cell');
runTime_Components = toc

numCC = length(bGG);
disp(num2str(numCC));

runTime_Connecting = 0;
% ############################
if numCC > 1
    
    disp('connecting the graph');
    
    tic
    addEdgeID = zeros(numCC-1, 2);
    addEdgeWW = zeros(numCC-1, 1);
    for ii = 1:(numCC-1)
        
        % random II
        allEdgeII = bGG{ii};
        rID1 = randperm(length(allEdgeII));
        idII = allEdgeII(rID1(1));
        
        % random JJ
        allEdgeJJ = bGG{ii+1};
        rID2 = randperm(length(allEdgeJJ));
        idJJ = allEdgeJJ(rID2(1));
        
        % edge & weight
        addEdgeIJ = sort([idII, idJJ]);
        key_addIJ = addEdgeIJ(:, 1)*(nGG+1) + addEdgeIJ(:, 2);
        addEdgeWWIJ = EdgeWW_org(Key_org==key_addIJ);
        
        % add
        addEdgeID(ii, :) = addEdgeIJ;
        addEdgeWW(ii) = addEdgeWWIJ;
    end
    
    % update graph as connected one
    EdgeID = [EdgeID; addEdgeID];
    EdgeWW = [EdgeWW; addEdgeWW];

    GG = graph(EdgeID(:, 1), EdgeID(:, 2), EdgeWW);
    runTime_Connecting = toc
end

disp('...saving data');

outName = [dsName '_' num2str(maxKC) '_RandSLE_Graph.mat'];
save(outName, 'GG', 'EdgeID', 'EdgeWW', 'XX_ID', 'XX_CC', ...
    'runTime_Connecting', 'runTime_Components', ...
    'maxKC', 'EdgeID_org', 'EdgeWW_org', 'Key_org', ...
    'nEE', 'mEE', 'DD', 'nGG', 'kGG', 'XX_sID', 'XX_eID', ...
    'N', 'XX', 'YY', 'WW');


disp('FINISH !!!');



