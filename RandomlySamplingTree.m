function [sEdgeID, sEdgeWW] = RandomlySamplingTree(N, EdgeID, EdgeWW)

% input:
% N: number of vertices
% EdgeID: mx2 (m edges -- IDs)
% EdgeWW: mx1 (weight)

% output
% sEdgeID: sampling edges
% sEdgeWW: corresponding weights of sampled edges

% (Assumption EdgeID(1) < EdgeID(2)

sEdgeID = zeros(N-1, 2);
sEdgeWW = zeros(N-1, 1);
flagEdge = zeros(N, 1);

rid = randperm(length(EdgeWW));

sEdgeID(1, :) = EdgeID(rid(1), :);
sEdgeWW(1) = EdgeWW(rid(1));
flagEdge(sEdgeID(1, :)) = 1;

sID = 1;
while sID < (N-1)
    for ii = 2:length(EdgeWW)
        id = rid(ii);
        val1 = EdgeID(id, 1);
        val2 = EdgeID(id, 2);

        % "add" edge
        if (flagEdge(val1) + flagEdge(val2)) == 1
            sID = sID + 1;
            sEdgeID(sID, :) = EdgeID(id, :);
            sEdgeWW(sID) = EdgeWW(id);
            % flag !!!
            flagEdge([val1 val2]) = 1;
            break;
        end
    end
end

end


