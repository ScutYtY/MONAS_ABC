function Offspring = OperatorEmployedABC(Parent,D,N,M,Problem)
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
[FrontNo,MaxFNo] = NDSort(Parent.objs,Parent.cons,N);
Next = FrontNo < 2;
%Problem = Problem.Current();
%多项式变异参数
proM=1;
disM=20;

Offspring = Parent;


% for i=1:91
%     allfunction(i) = sum(Parent(i).obj,2);
% end
% [x1,i1] = min(allfunction);
% [x2,i2] = max(allfunction);
b = [Parent.obj];
c = [Parent.best.obj];
%b = reshape(b,3,[]);
b = reshape(b,[],M);
c = reshape(c,[],M);
nobj = Normalization_Obj(b);
bnobj = Normalization_Obj(c);
e= sum(nobj,2);
f = sum(bnobj,2);
[x1,i1] = min(e);
[x2,i2] = max(e);
[q,popindex] = sort(e);
[o,bindex] = sort(f);
% AI = index(1:0.1*sum(Next));
% AI2 = index(11:sum(Next));
% %AI3 = index(sum(Next)+1:91);
% RE = randi(10);


% if(0.2*sum(Next)<=1)
%     NON1 = bindex;
%     NON2 = bindex;
%     RE = randi(length(bindex));
%     RE1 = randi(length(NON2));
% else
%     NON1 = bindex(1:0.2*sum(Next));
%     NON2 = bindex(0.2*sum(Next):sum(Next));
%     RE = randi(length(NON1));
%     RE1 = randi(length(NON2));
%     RE2 = randi(length(bindex));
% 
% end

    NON1 = bindex;
    NON2 = bindex;
    RE = randi(length(bindex));
    RE1 = randi(length(NON2));


%for i=1:0.1*sum(Next)
for i=1:N 
    K = [1:i-1 i+1:N];
    k = K(randi([1 numel(K)]));
    r1 = K(randi([1 numel(K)]));
    r2 = K(randi([1 numel(K)]));
    r3 = K(randi([1 numel(K)]));
    while 1
        r1 = K(randi([1 numel(K)]));
        r2 = K(randi([1 numel(K)]));
        r3 = K(randi([1 numel(K)]));
        if(r1~=r2 && r2~=r3 && r1~=r3 && i~= r1 && i~= r2 && i~= r3)
        %if(r1~=r2 && i~= r1 && i~= r2)
            break
        end
    end
    %j = randi(D);
    %R = rand*2 - 1;
    j = randi(D);
    j1 = randperm(D,randi(D));
    R = 2*rand(1,length(j))-1;
    j2 = R;
    m=1;
    for m = 1:length(R)
        j2(m) = randi(D);
    end
    
    %j = randperm(D,3);
    R1 = rand*2 - 1;
    x = Parent(i).dec(j);
    v = Parent(i).dec;
%     if rand<0.5
%         v(j1) = normrnd((Parent(NON1(RE)).dec(j1)+Parent(NON2(RE1)).dec(j1))/2,abs(Parent(NON2(RE1)).dec(j1)-Parent(NON1(RE)).dec(j1)));
%         
%     else
%         v(j1) = Parent(NON2(RE1)).dec(j1);
%     end
    %v(j) = Parent(AI(RE)).dec(j) + R*(Parent(AI(RE)).dec(j)-Parent(k).dec(j));
    %v(j) = Parent(AI(RE)).dec(j) + R*(Parent(r1).dec(j)-Parent(r2).dec(j));
    %v(j) = normrnd((Parent(i).dec(j)+Parent(AI(RE)).dec(j))/2,abs(Parent(i).dec(j)-Parent(AI(RE)).dec(j)));
    %v(j) = Parent(i1).dec(j) + Parent(i2).dec(j)-Parent(i).dec(j);
    %v(j)= Parent(i).dec(j) + R1*(Parent(i).dec(j2)-Parent(r2).dec(j));
    %v(j)= Parent(i).dec(j) + R*(Parent(i).dec(j)-Parent(r2).dec(j));
    v(j1) = Parent(NON1(RE)).dec(j1) + R.*(Parent(NON1(RE)).dec(j1)-Parent(r2).dec(j1));
    %v(j) = Parent(NON1(RE)).dec(j) + R*(Parent(k).dec(j)-Parent(r2).dec(j));
    %v(j) = Parent(i).dec(j) + R*(Parent(i).dec(j)-Parent(k).dec(j));
    %进行多项式编译
    Lower = repmat(Problem.lower,1,1);
    Upper = repmat(Problem.upper,1,1);
    v       = min(max(v,Lower),Upper);
%     Site  = rand(1,length(j1)) < proM/D;
%     mu    = rand(1,length(j1));
% 
%     temp  = Site & mu<=0.5;
%     v(temp) = v(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
%                       (1-(v(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
%     temp = Site & mu>0.5; 
%     v(temp) = v(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
%                       (1-(v(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));
%     if rand<0.1
%         v(j1) = Parent(r3).dec(j2) + R.*(Parent(r1).dec(j2)-Parent(r2).dec(j2));
%         %v(j1) = Parent(i1).dec(j2) + Parent(i2).dec(j2)-Parent(NON1(RE)).dec(j2);
%     end

% MaOACB-TA
     if rand<1-(i/N)
         v(j1) = Parent(r3).dec(j2) + R.*(Parent(r1).dec(j2)-Parent(r2).dec(j2));
         v(j1) = Parent(i1).dec(j2) + Parent(i2).dec(j2)-Parent(NON1(RE)).dec(j2);
     end

    Offspring(i) = Problem.Evaluation(v);
    
  
    
end


% for i = 1:length(AI2)
%         K = [1:i-1 i+1:91];
%         k = K(randi([1 numel(K)]));
%         r1 = K(randi([1 numel(K)]));
%         r2 = K(randi([1 numel(K)]));
%         b = best(randi(numel(best))).dec;
%         j = randi(D);
%         R = rand*2 - 1;
%         R1 = rand*1.5;
%         %j = randperm(D,randi(D));
%         %j = randperm(D,3);
%         a = randi(length(AI));
%         v = Parent(AI2(i)).dec;
%         v(j) = Parent(AI2(i)).dec(j) + R*(Parent(AI2(i)).dec(j)-b(j)) + R1*(Parent(AI(a)).dec(j)-Parent(AI2(i)).dec(j));
%         Offspring(AI2(i)) = SOLUTION(v);
% end
% 
% for i = 1:length(AI3)
%         K = [1:i-1 i+1:91];
%         k = K(randi([1 numel(K)]));
%         r1 = K(randi([1 numel(K)]));
%         r2 = K(randi([1 numel(K)]));
%         
%         b = best(randi(numel(best))).dec;
%         %j = randi(D);
%         %R = rand*2 - 1;
%         j = randperm(D,randi(D));
%         R = 2*rand(1,length(j))-1;
% 
%         %j = randperm(D,3);
%   
%         x = Parent(AI3(i)).dec(j);
%         v = Parent(AI3(i)).dec;
%         a = randi(length(AI));
%         v(j) = (Parent(AI(a)).dec(j)+ b(j))/2 + R.*((Parent(AI(a)).dec(j)+ b(j))/2 -Parent(AI3(i)).dec(j));           
%         Offspring(AI3(i)) = SOLUTION(v);
% 
% end
    