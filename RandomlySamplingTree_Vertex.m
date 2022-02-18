function [sEdgeID, sEdgeWW] = RandomlySamplingTree_Vertex(N, EdgeID, EdgeWW)

% input:
% N: number of vertices
% EdgeID: mx2 (m edges -- IDs)
% EdgeWW: mx1 (weight)

% output
% sEdgeID: sampling edges
% sEdgeWW: corresponding weights of sampled edges

% (Assumption EdgeID(1) < EdgeID(2)

setID1 = [1];
setID0 = [2:N];

sEdgeID = zeros(N-1, 2);
sEdgeWW = zeros(N-1, 1);

KeyID = EdgeID(:, 1)*(N+1) + EdgeID(:, 2);

for ii = 1:(N-1)

    disp(['...' num2str(ii)]);
    randID1 = randperm(length(setID1));
    randID2 = randperm(length(setID0));
    
    rsetID1 = setID1(randID1);
    rsetID0 = setID0(randID2);
    
    flagEE = 0;
    for id1 = 1:length(rsetID1)
        val1 = rsetID1(id1);
        for id2 = 1:length(rsetID0)
            val2 = rsetID0(id2);
            
            % check [val1 val2]
            edge12 = sort([val1 val2]);
            key12 = val1*(N+1) + val2;
            
            fID = find(KeyID == key12);
            if length(fID) > 0 % found
                % save edge
                sEdgeID(ii, :) = edge12;
                sEdgeWW(ii) = EdgeWW(fID);
                % update
                setID1 = [setID1 val2]; % add val2 --> setID1
                setID0(setID0 == val2) = []; % delete val2 from setID0
                % update flag (found!)
                flagEE = 1;
                break;
            end
        end
        if flagEE > 0
            break;
        end
    end
end

end


