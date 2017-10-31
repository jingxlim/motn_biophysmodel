function U = make_U(t,inj_cmprt,N,radius,length,...
                    tp_AMPA,tp_GABA,Gs_AMPA,Gs_GABA,D_AMPA,D_GABA,...
                    AMPA_inputt,GABA_inputt)

    Gp_AMPA = @(r,dl) Gs_AMPA * D_AMPA * 2*pi*r*dl;
    Gp_GABA = @(r,dl) Gs_GABA * D_GABA * 2*pi*r*dl;

    g_AMPA = @(t,inputt,r,dl) ((t-inputt)/tp_AMPA).*exp(1-((t-inputt)/tp_AMPA))...
                               .*Gp_AMPA(r,dl);
    g_GABA = @(t,inputt,r,dl) ((t-inputt)/tp_GABA).*exp(1-((t-inputt)/tp_GABA))...
                                .*Gp_GABA(r,dl);

    U = zeros(2*N,1);  % initialize U matrix
    
    AMPA_cond = 0;  % conductnace is 0 unless there's an input
    for i=1:numel(AMPA_inputt)
        inputt = AMPA_inputt(i);
        % time shifted conductance
        AMPA_g = @(t) subplus(g_AMPA(t,inputt,radius(inj_cmprt),length(inj_cmprt)));
        AMPA_cond = AMPA_cond + AMPA_g(t);
    end
    
    U(2*inj_cmprt-1,1) = AMPA_cond;  % excitatory
    
    GABA_cond = 0;  % conductnace is 0 unless there's an input
    
    for i=1:numel(GABA_inputt)
        inputt = GABA_inputt(i);
        % time shifted conductance
        GABA_g = @(t) subplus(g_GABA(t,inputt,radius(inj_cmprt),length(inj_cmprt)));
        GABA_cond = GABA_cond + GABA_g(t);
    end
    
    U(2*inj_cmprt,1) = GABA_cond;  % inhibitory
end