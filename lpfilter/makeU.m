function U = makeU(t,inj_cmprt,N,mu,sigma,seed)

U = zeros(N,1);
U(inj_cmprt,1) = mu + sigma.*rand(1);

end