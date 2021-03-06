% Matrix of the squares of the mutual distances of 3D points
function dist2 = SUB_mat_dist2_3D(XY)
    N = length(XY(:,1));
    x = XY(:,1);
    y = XY(:,2);
    z = XY(:,3);
    
    dist2 = zeros(N,N); % dist2(i,i) = 0
    for i=1:N
        for j=1:N
            if i~=j
                dist2(i,j) = (x(i)-x(j))^2 + (y(i)-y(j))^2 + (z(i)-z(j))^2;
            end
        end
    end    
end
