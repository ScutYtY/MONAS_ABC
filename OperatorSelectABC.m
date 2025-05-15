function Offspring = OperatorSelectABC(Population,N,M)
%OperatorGAhalf - Crossover and mutation operators of genetic algorithm.
%
%   This function is the same to OperatorGA, while only the first half of
%   the offsprings are evaluated and returned.
%
%   See also OperatorGA

%------------------------------- Copyright --------------------------------
% Copyright (c) 2021 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
%j = randi(D);
% %j = randperm(D,randi(D));
% x = Parent.dec(j);
% k = Prandom.dec(j);
%b = best(j);
% r = rand*2 - 1;
% v = Parent.dec;
% %v(j) = b + r*(b-k);
% v(j) = normrnd((x+b)/2,abs(x-b));
% Offspring = SOLUTION(v);
%  

[FrontNo,MaxFNo] = NDSort(Population.objs,Population.cons,N);
Next = FrontNo < MaxFNo; 

% for i=1:2*N
%     allfunction(i) = sum(Population(i).obj,2);
% end
b = [Population.obj];
c = [Population.best.obj];
% b = reshape(b,3,[]);
% e= sum(b);
c = reshape(c,[],M);
nobj = Normalization_Obj(c);
e= sum(nobj,2);
[q,index] = sort(e);
AI = index(1:0.5*length(q));
Next(AI) = true;

if sum(Next) < N      %保证Next为种群大小 91
%     Temp = find(FrontNo==MaxFNo & KneePoints==0);
%     [~,Rank] = sort(Distance(Temp),'descend');
%     Next(Temp(Rank(1:(K-sum(Next))))) = true;
%     for i=1:182
%         allfunction(i) = sum(Population(i).obj,2);
%     end
    b = [Population.obj];
    b = reshape(b,[],M);
    nobj = Normalization_Obj(b);
    e= sum(nobj,2);
    %Nextsum = sum(b) > median(allfunction,2);
    %[c,d]      = sort(sum(b));
    %c = find(Next==1);
    Temp = find(FrontNo==MaxFNo); %c(Temp) 除最后一层以外的解的索引
    [q,index] = sort(e(Temp));
%     w = q(1:(N-sum(Next)));
    Next(Temp(index(1:(N-sum(Next))))) = true;
   

elseif sum(Next) > N
%     Temp = find(FrontNo==MaxFNo & KneePoints==1);
%     [~,Rank] = sort(Distance(Temp));
%     Next(Temp(Rank(1:(sum(Next)-K)))) = false;
%     for i=1:182
%         allfunction(i) = sum(Population(i).obj,2);
%     end
    b = [Population.obj];
    b = reshape(b,M,[]);
    %Nextsum = sum(b) > median(allfunction,2);
    
    Temp = find(FrontNo~=1); %c(Temp) 除最后一层以外的解的索引
    e= sum(b);
    [q,index] = sort(e(Temp),'descend');
    index(1:(sum(Next)-N));
    
    Next(Temp(index(1:(sum(Next)-N)))) = false;

end
% for i=1:N
%     Population(Next).dec = v;
% 
% end


Offspring = Population(Next);
% for i=1:182
%     a = sum(Population(i).obj,2);
% end
% b = sort(a);




    