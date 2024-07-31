function [obj] = RunStochasticGillepsieEnsemble(obj,filename,Trajectories,displayProgress)
    %check inputs
    arguments
        obj
        filename string
        Trajectories double
        displayProgress logical
    end

    %% Load/set parameters
    if isempty(filename)
        %if not file, run an SIR model with no variation in susceptibility
        %or transmissibility. Keeping a structure here such that one could
        %extend the number of susceptibility and transmissibility
        %variables.
        susceptibility = 1;
        transmissibility = 1;
        RES_s = length(susceptibility); %number of susceptibility values
        RES_t = length(transmissibility); %number of transmissibility values
        S0 = zeros(RES_t,RES_s);
        I0 = zeros(RES_t,RES_s);
        R0 = zeros(RES_t,RES_s);
        for aa = 1:obj.Population_Size
            x = randi([1 RES_t],1,1);
            y = randi([1 RES_s],1,1);
            S0(x,y) = S0(x,y)+1;
        end

        beta = obj.beta;
        gamma = obj.gamma;
        obj.init_mean_transmissibility = 1;
        obj.init_mean_susceptibility = 1;
        obj.init_correlation = NaN;

    else

        [beta,gamma,susceptibility,transmissibility,initjoint] = stoch_loadgrab(filename);
        %seeding initial population
        RES_s = length(susceptibility); %number of susceptibility values
        RES_t = length(transmissibility); %number of transmissibility values
        S0 = zeros(RES_t,RES_s);
        I0 = zeros(RES_t,RES_s);
        R0 = zeros(RES_t,RES_s);

        relativized = initjoint./sum(sum(initjoint));
        cumulative = cumsum(relativized(:));
        rands = rand(1,obj.Population_Size);
        Susvec = 0.*rands;
        Travec = Susvec;
        for aa = 1:length(rands)
            indexs = find(cumulative>=rands(aa),1);
            [row,col] = ind2sub([RES_t RES_s],indexs);
            S0(row,col) = S0(row,col)+1;
            Susvec(aa) = susceptibility(col);
            Travec(aa) = transmissibility(row);
        end

        %% check  average sus + tra %%
        tra = 0;
        sus = 0;
        for xx = 1:size(S0,1)
            for yy = 1:size(S0,2)
                sus = sus + S0(xx,yy)*susceptibility(yy);
                tra = tra + S0(xx,yy)*transmissibility(xx);
            end
        end

        rho = corrcoef(Susvec,Travec); %correlation coefficient
        obj.init_mean_transmissibility = tra/obj.Population_Size;
        obj.init_mean_susceptibility = sus/obj.Population_Size;
        obj.init_correlation = rho(1,2);
        obj.beta = beta;
        obj.gamma=gamma;
        clear relativized cumulative rands Travec Susvec rho;
    end

    %% Initialize model architecture
    %architecture to run/save stochastic ensemble
    saveOutbreakSize = zeros(1,Trajectories);
    saveOutbreakTime = zeros(1,Trajectories);
    %index cases for each trajectory:
    randIperson = randi([1 obj.Population_Size],1,Trajectories);
    r = zeros(1,RES_t*RES_s+RES_t*RES_s);  %should be length of infection routes and recovery routes
    act1 = r; %first index of action
    act2 = r; %second index of action

    %set action indices (according to rules in stochasticReactions.m)
    act1 = act1.*0;
    act2 = act2.*0;
    rind=1;
    %S->I
    for s_a = 1:RES_s
        for s_b = 1:RES_t
            act1(rind) = s_a;
            act2(rind) = s_b;
            rind = rind+1;
        end
    end
    %I->R
    for aa = 1:RES_s
        for bb = 1:RES_t
            act1(rind) = aa;
            act2(rind) = bb;
            rind=rind+1;
        end
    end

    %% simulate each trajectory and store relevant information
    for aa = 1:Trajectories
        S = S0;
        I = I0;
        R = R0;
        %seed the potential outbreak
        INF = randIperson(aa);
        JJ = cumsum(S(:));
        indexi = find(JJ>=INF,1);
        [row,col] = ind2sub([RES_t RES_s],indexi);
        S(row,col) = S(row,col)-1;
        I(row,col) = 1;
        Itot=1;
        Rtot = 0;

        %start the reaction loop
        time = 0;
        %Inext = 1;
    
        while Itot > 0  %while more than one person infected
            % total reaction rate
            r = stochasticReactions(S,I,susceptibility,transmissibility,RES_s,RES_t,beta,gamma,obj.Population_Size,r);
            if sum(r<0)>0
                disp("reaction rate less than 0. exiting.")
                break
            end
            
            r_prob = cumsum(r)/sum(r);
            %update time step (not computing as not needed, so commented out)
                tau = 1/sum(r) *log(1/rand());
                %time = [time time(end)+tau];
                time = time+tau;
            %choose reaction
            choose = rand();
            index = find(r_prob>choose,1); %choose reaction
            if (index > RES_t*RES_s)  % recovery
                I(act1(index),act2(index)) = I(act1(index),act2(index))-1;
                R(act1(index),act2(index)) = R(act1(index),act2(index))+1;
                Itot = Itot-1;
                Rtot = Rtot+1;
            else  %infection
                S(act1(index),act2(index)) = S(act1(index),act2(index))-1;
                I(act1(index),act2(index)) = I(act1(index),act2(index))+1;
                Itot = Itot+1;
            end
            %disp(Itot)
            %not recording so commenting out
            %Inext = [Inext Itot];
        end
    
        %record output
        saveOutbreakSize(aa) = Rtot;
        saveOutbreakTime(aa) = time;
        if displayProgress==true
            disp(aa)
        end
    end

    obj.Final_Size = saveOutbreakSize;
    obj.Final_Time = saveOutbreakTime;

end