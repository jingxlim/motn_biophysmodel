figure(3); clf; hold on;
set(gcf,'units','points','position',[100,100,1000,400])
subplot(1,2,1); hold on;
subplot(1,2,2); hold on;

injects = round(linspace(142,173,5));

for j=1:numel(injects)
    inject = injects(j);

    % construct U vector
    U_mat = @(t) makeU(t,inject,N,mu,sigma,seed);
    
    % check changes in conductance in U matrix
    figure(); clf; hold on;
    set(gcf,'units','points','position',[100,100,1000,400])
    
    subplot(1,2,1); hold on;
    plot(sim_time, Iapp_rand);
    ylabel('I_{app} [mA]'); xlabel('Time [us]');
    
%     for i=1:numel(sim_time)
%         t = sim_time(i)
%         mat_U = U_mat(t);
%         u = mat_U(inject,1);
%         plot(t, u,'.');
%     end
    
    % Solve for voltage over time
    % system of differential equations
    dvdt = @(t,v) A*v + B_mat*U_mat(t);

    sim_T = sim_time/(Rm*Cm);

    [t,v] = ode23(@(t,v) dvdt(t,v), sim_time, zeros(N,1), options);
    
    v_store{end+1} = v;  % save data    

    subplot(1,2,2); hold on;
    plot(sim_T,v(:,inject),'DisplayName','V_{injection}');
    plot(sim_T,v(:,1),'DisplayName','V_{soma}');
    title('Time evolution of V(X,T) in response to I_{app}');
    xlabel('T'); ylabel('V [mV]'); legend('show');    
    saveas(gcf, strcat('IO_Iapp',outdate,'.png'));

    % spectral analysis for I_app input
    fs = size(sim_time,2);
    psa = [inject,1];
    
    figure(); clf; hold on;
    set(gcf,'units','points','position',[100,100,1000,400])
    
    subplot(1,2,1); hold on;
    x = Iapp_rand';
    n = size(x,1);
    xdft = fft(x);
    xdft = xdft(1:n/2+1);
    psdx = (1/(fs*n)) * abs(xdft).^2;
    psdx(2:end-1) = 2*psdx(2:end-1);
    freq = 0:fs/n:fs/2;
    plot(freq,10*log10(psdx),'DisplayName', 'Input current')
    xlabel('Frequency (Hz)'); ylabel('Power/Frequency (dB/Hz)')    

    subplot(1,2,2); hold on;
    powers = {};
    for i=1:numel(psa)
        ccmprt = psa(i);
        x = v(:,ccmprt);
        n = size(x,1);
        xdft = fft(x);
        xdft = xdft(1:n/2+1);
        psdx = (1/(fs*n)) * abs(xdft).^2;
        psdx(2:end-1) = 2*psdx(2:end-1);
        freq = 0:fs/n:fs/2;
        power = 10*log10(psdx);
        powers{end+1} = power;
        plot(freq, power, 'DisplayName', 'cmprt '+string(ccmprt))
    end
    xlabel('Frequency (Hz)'); ylabel('Power/Frequency (dB/Hz)');
    grid on; legend('show')
    saveas(gcf, strcat('spectra_',num2str(inject),'_Iapp',outdate,'.png'));

    figure(3); subplot(1,2,1);
    plot(freq, powers{2}-powers{1});

    figure(3); subplot(1,2,2);
    plot(freq, powers{2}-powers{1});
    
end

figure(3);
subplot(1,2,1);
xlabel('Frequency (Hz)')
ylabel('Power/Frequency difference (dB/Hz)')
subplot(1,2,2);
xlabel('Frequency (Hz)')
ylabel('P/F diff (dB/Hz)')
xlim([0 50]);
suptitle('Periodogram Using FFT')
saveas(gcf, strcat('spectradiff_Iapp',outdate,'.png'));