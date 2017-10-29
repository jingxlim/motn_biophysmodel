function G_ = make_G(t,inj_cmprt,N,AMPA_cond,sim_time)
    G_ = zeros(N);
    G_(inj_cmprt,inj_cmprt) = -AMPA_cond(find(sim_time==t));
end