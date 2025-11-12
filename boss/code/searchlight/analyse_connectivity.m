function [mappedPoints, mappingMatrix] = analyse_connectivity(r, minDist)
    S = zeros(size(r, 1),size(r, 1));
    for a = 1: size(r, 1)
        for b = 1: size(r, 1)
            num = 0;
            den = 0;
            mu = mean((r(a, :)-r(b, :))/2);
            for j = 1: size(r, 2)
                muj = (r(a, j)-r(b, j))/2;
                num = num + (r(a, j)-muj).^2-(r(b, j)-muj).^2;
                den = den + (r(a, j)-mu).^2-(r(b, j)-mu).^2;
            end
            S(a, b) = abs(1 - num/den);
        end
    end
    
    W = zeros(size(r, 1),size(r, 1));
    for a = 1:size(r, 1)
        for b = 1:size(r, 1)
            W(a,b) = norm(S(a,~isnan(S(a,:)))-S(b,~isnan(S(b,:))))^2;
        end
    end
    
    W(W<minDist) = 0;
    W(W~=0) = 1;
    L = eye(size(W, 1), size(W, 1)).*sum(W,2) - W;
    [V, ~] = eig(L);
    
    mappingMatrix = V(:, 1:2);
    mappedPoints = r' * mappingMatrix;
end