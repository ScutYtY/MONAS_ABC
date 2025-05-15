function Offspring = OperatorOnlookerABC(Parent,D,N,M,Problem)
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

%allfunction = 0;

[FrontNo,MaxFNo] = NDSort(Parent.objs,Parent.cons,N);
Next = FrontNo < 2;

%Problem = PROBLEM.Current();
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
%e= sum(b);
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
% AI = index(1:sum(Next));
% RE = randi(sum(Next));
% NON1 = bindex(1:0.2*sum(Next));
% %disp(length(NON1));
% NON2 = bindex(0.8*sum(Next):sum(Next));
% 
% if(length(NON1)==0) 
%     
%     
%     RE = randi(length(bindex));
%     RE1 = randi(length(NON2));
% else
%     NON1 = bindex(1:0.2*sum(Next));
%     disp(length(NON1));
%     NON2 = bindex(0.8*sum(Next):sum(Next));
%     RE = randi(length(NON1));
%     RE1 = randi(length(NON2));
% end

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
    
t=1;
i=1;
while t<=N
    P = (max(e)-q(i))/(max(e)-min(e));
    %Efunction = sum(Parent(i).obj,2);
    
    if P > rand
        t = t + 1;
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
       
        j = randi(D); %单维
        j1 = randperm(D,randi(D));%多维
        %j = randperm(D,3);
        %R = rand*2 - 1;
       
        R = 2*rand(1,length(j))-1;
        j2 = R;
        m=1;
        for m = 1:length(R)
            j2(m) = randi(D);
        end       %异维
        
        x = Parent(i).dec(j);
        v = Parent(i).dec;
        %v(j) = normrnd((x+b)/2,abs(x-b));
        %v(j) = Parent(i1).dec(j) + r*(Parent(i1).dec(j)-Parent(k).dec(j));
%         if rand<0.5
%             v(j1) = normrnd((Parent(NON1(RE)).dec(j1)+Parent(NON2(RE1)).dec(j1))/2,abs(Parent(NON2(RE1)).dec(j1)-Parent(NON1(RE)).dec(j1)));
%         else
%             v(j1) = Parent(NON1(RE)).dec(j1) ;
%         end
        %v(j) = Parent(AI(RE)).dec(j) + R*(Parent(AI(RE)).dec(j)-Parent(k).dec(j));     
        %v(j) = Parent(AI(RE)).dec(j) + R*(Parent(r1).dec(j)-Parent(r2).dec(j));
        %v(j) = Parent(i).dec(j) + R*(Parent(i).dec(j)-Parent(k).dec(j));
        v(j1) = Parent(NON1(RE)).dec(j1) + R.*(Parent(NON1(RE)).dec(j1)-Parent(r2).dec(j1));  
        %v(j) = Parent(NON1(RE)).dec(j) + R*(Parent(k).dec(j)-Parent(r2).dec(j));  

%       MaOACB-TA
        
         if rand<1-(t/N)
             v(j1) = Parent(NON2(RE1)).dec(j2) + R.*(Parent(r1).dec(j2)-Parent(r2).dec(j2));
             v(j1) = Parent(i1).dec(j2) + Parent(i2).dec(j2)-Parent(NON1(RE)).dec(j2);
         end
%         if rand<0.05
%             v(j1) = Parent(NON1(RE)).dec(j2) + R.*(Parent(r1).dec(j2)-Parent(r2).dec(j2));
%             %v(j1) = Parent(i1).dec(j2) + Parent(i2).dec(j2)-Parent(NON1(RE)).dec(j2);
%         end
         Lower = repmat(Problem.lower,1,1);
         Upper = repmat(Problem.upper,1,1);
         v       = min(max(v,Lower),Upper);
%         Site  = rand(1,length(j1)) < proM/D;
%         mu    = rand(1,length(j1));
% 
%         temp  = Site & mu<=0.5;
%         v(temp) = v(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
%                       (1-(v(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
%         temp = Site & mu>0.5; 
%         v(temp) = v(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
%                       (1-(v(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));
%         if rand<0.1
%             v(j1) = Parent(r3).dec(j2) + R.*(Parent(r1).dec(j2)-Parent(r2).dec(j2));
%             %v(j1) = Parent(i1).dec(j2) + Parent(i2).dec(j2)-Parent(NON1(RE)).dec(j2);
%         end
%          if rand<0.1
%              v(j1) = Parent(r3).dec(j2) + R.*(Parent(r1).dec(j2)-Parent(r2).dec(j2));
%              %v(j1) = Parent(i1).dec(j2) + Parent(i2).dec(j2)-Parent(NON1(RE)).dec(j2);
%          end
        %Offspring(NON1(i)) = SOLUTION(v);  
        Offspring(i) = Problem.Evaluation(v);
        i = i +1;
    end
    i =  mod(i,length(NON1)) + 1;
    %i =  mod(i,sum(Next)) + 1;
% for i=1:50
    
    %end
    
end
 
            