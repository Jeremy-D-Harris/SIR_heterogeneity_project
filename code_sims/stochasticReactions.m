function [r] = stochasticReactions(S,I,susceptibility,transmissibility,RES_s,RES_t,beta,gamma,PopulationSize,r)
    %calculate total force of infection
    totinfforce = 0;
    for i_a = 1:RES_t
        for i_b = 1:RES_s
            totinfforce = totinfforce + transmissibility(i_a)*I(i_a,i_b)/PopulationSize;
        end
    end
    %S->I
    rind=1;
    for s_a = 1:RES_t
        for s_b = 1:RES_s
            r(rind) = beta*susceptibility(s_b)*S(s_a,s_b)*totinfforce;
            rind = rind+1;
        end
    end
    %I->R
    for aa = 1:RES_t
        for bb = 1:RES_s
            r(rind) =  gamma*I(aa,bb);
            rind=rind+1;
        end
    end
end