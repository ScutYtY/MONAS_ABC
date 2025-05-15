function PopObj = Normalization_Obj(PopObj)    
%% Normalization
    % Detect the extreme points
    [N,M]  = size(PopObj);
    Zmin = min(PopObj,[],1);
    Zmax = max(PopObj,[],1);
    % Normalization
    for i = 1:N
        for j = 1:M
            PopObj(i,j) = (PopObj(i,j) - Zmin(j))/(Zmax(j) - Zmin(j));
        end
    end
end