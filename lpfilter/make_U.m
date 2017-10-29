function U = make_U(t,inj_cmprt,N,AMPA_cond,sim_time)
    U = zeros(2*N,1);
    time = t
    index = find(sim_time==t)
    U(2*inj_cmprt-1,1) = AMPA_cond(index);  % excitatory
end