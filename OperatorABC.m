function Offspring = OperatorABC(Parent1,Parent2,Parent3,Parameter)


    %% Parameter setting
    if isa(Parent1(1),'SOLUTION')
        calObj  = true;
        Parent1 = Parent1.decs;
        Parent2 = Parent2.decs;
        Parent3 = Parent3.decs;
    else
        calObj = false;
    end
    [N,D]   = size(Parent1);
    Problem = PROBLEM.Current();

    %% Differental evolution
    Site = rand(N,D) < 0.45;
    F=rand*2-1;
    Offspring       = Parent1;
    Offspring(Site) = Offspring(Site) + F*(Parent2(Site)-Parent3(Site));

    %% Polynomial mutation

    if calObj
        Offspring = SOLUTION(Offspring);
    end
end