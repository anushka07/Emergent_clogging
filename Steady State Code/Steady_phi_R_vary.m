% plots steady state for imposed phim(z)

clear all;

z = linspace(0,1,1000);

mu1 = 0.3;
p_arr = 0.2;        % Inlet particle fraction
phim0 = 0.65;       % Max packing fraction

Rc = 0.5;           % Set the (1-R_min)

arr = -2;
Da_arr = 10.^arr;   % Da

Da = Da_arr; 
p0 = p_arr;
pm = phim0;

d=10;
b = Rc/2;
a = 1-b;
RR = a + b*(erf(-(z-0.5).*d));

muw1 = linspace((1+1e-1)*mu1,50*mu1,1e6);
muw2 = linspace((1+1e-15)*mu1,(1+1e-1)*mu1,1e3);
muw = cat(2,muw2,muw1);
m = muw;
mw=m;
mu=mu1;
rc = mu1./muw;
f1 = (2.*rc)./m - 2.*((2.*rc)./(3.*m.^(1./2)) - (2.*(m.*rc + 1))./m.^(3./2)).*(1 - rc).^(1./2) + (4.*(1 - rc).^(1./2))./(3.*m.^(1./2)) - 2./m + rc.^2 - (pi.*(m.*rc + 1).*2i)./m.^2 + (2.*log(-(m - m.*rc - 2.*m.^(1./2).*(1 - rc).^(1./2) + 1)./(m.*rc - m + 1)).*(m.*rc + 1))./m.^2 - (2.*log(1 - (m.*rc + 1)./m).*(m.*rc + 1))./m.^2 + (2.*log(rc - (m.*rc + 1)./m).*(m.*rc + 1))./m.^2;
f3 = (rc.^2.*(rc - 1./2))./2 - rc./6 + (rc.^2.*(rc - 1).^2)./4 - (5.*rc.^4)./24 + 1./8;
f2 = (rc.^2.*(rc - 1).^2)./4 - (rc - 1).^2./(2.*m) - (((2.*rc.*(rc - 1).^2)./(3.*m.^(1./2)) - (2.*(m.*rc + 1).*(rc - 1).^2)./m.^(3./2)).*(1 - rc).^(1./2))./2 + (m.*(rc./m.^3 + 1./m.^4))./2 + (m - m.*rc).^3./(6.*m.^4) - (m - m.*rc).^(7./2)./(7.*m.^4) + (1 - rc).^(5./2)./(3.*m.^(1./2)) - ((m - m.*rc).^(1./2).*((2.*rc)./m.^3 + 2./m.^4))./2 + ((m - m.*rc).^2.*(rc./(2.*m.^3) + 1./(2.*m.^4)))./2 - ((m - m.*rc).^(3./2).*((2.*rc)./(3.*m.^3) + 2./(3.*m.^4)))./2 - ((m - m.*rc).^(5./2).*((2.*rc)./(5.*m.^3) + 2./(5.*m.^4)))./2 + (log((m - m.*rc).^(1./2) + 1).*(2.*m.*rc + 2))./(2.*m.^4) - (m.*rc.*(rc./m.^3 + 1./m.^4))./2 + (rc.*(rc - 1).^2)./(2.*m) + (log(rc - (m.*rc + 1)./m).*(m.*rc + 1).*(rc - 1).^2)./(2.*m.^2) - (pi.*(m.*rc + 1).*(rc - 1).^2.*1i)./(2.*m.^2) + (log(-(m - m.*rc - 2.*m.^(1./2).*(1 - rc).^(1./2) + 1)./(m.*rc - m + 1)).*(m.*rc + 1).*(rc - 1).^2)./(2.*m.^2) - (log(1 - (m.*rc + 1)./m).*(m.*rc + 1).*(rc - 1).^2)./(2.*m.^2);
f4 = @(phim)-(16.*mu.^2.*(mw - mu).^(1./2) - 24.*mw.^2.*(mw - mu).^(1./2) + 60.*phim.^3.*(mw - mu).^(1./2) + 15.*mu.*mw.^2 + 30.*mu.*phim.^3 + 45.*mw.^2.*phim - 30.*mw.*phim.^3 - 5.*mu.^3 - 15.*mw.^2 - 10.*mw.^3 - 30.*phim.^3.*log((mw - mu + 2.*(mw - mu).^(1./2) + 1)) + 15.*mu.^2.*phim.^3 - 45.*mw.^2.*phim.^2 + 8.*mu.*mw.*(mw - mu).^(1./2) - 30.*mu.*phim.^3.*log((mw - mu + 2.*(mw - mu).^(1./2) + 1)) - 24.*mu.^2.*phim.*(mw - mu).^(1./2) + 40.*mu.*phim.^3.*(mw - mu).^(1./2) + 36.*mw.^2.*phim.*(mw - mu).^(1./2) + 20.*mw.*phim.^3.*(mw - mu).^(1./2) - 12.*mu.*mw.*phim.*(mw - mu).^(1./2))./(15.*mw.^2.*phim.^2);

minR = min(RR);
bot = (f3 + Da.*f4(pm)./minR.^2);
top = pm.*f2;
flux1 = top./bot;

[maxflux I] = max(flux1);

%find mu0
[~,k0] = min(abs(phim0*f1-p0)); % Find index of muw for phi_in
%find initial flux
bot = (f3 + Da*f4(pm)./RR(1).^2);
top = pm*f2;
flux = top./bot;                % Inlet flux function for all muw
flux0 = flux(k0);               % Inlet flux corresponding to inlet muw m(k0) 
Sin(l) = flux0;


%NO CLOG
if(maxflux>flux0)
    for i=1:length(z)
        Rval = RR(i);
        bot = (f3 + Da*f4(pm)./Rval.^2);
        top = pm*f2;
        flux = top./bot;
        [~,j] = max(flux);                  % Get the lower phi branch/ higher mu branch closer to inital condition
        if(k0<j)
        [~,k] = min(abs(flux(1:j)-flux0));  % gives index relative to j
        pbar(i) = pm*f1(k);
        mw_s(i) = muw(k);
        In(i) = k;
        u = k;
        else
        [~,k] = min(abs(flux(j:end)-flux0));% gives index relative to j
        pbar(i) = pm*f1(k+j-1);             % to get accurate index 
        mw_s(i) = muw(k+j-1);
        In(i) = k+j-1;
        u = k+j-1;
        end
        G(i) = 1./bot(u)/Rval.^4;
    end

else %CLOG
    for i=1:length(z)
       Rval = RR(i);
       bot = (f3 + Da*f4(pm)./Rval.^2);
       top = pm*f2;
       flux = top./bot;
       flux(1)=0;
      [~,j] = max(flux);
           [~,k] = min(abs(flux(1:j)-maxflux));% higher phi vs flux branch
           pbar(i) = pm*f1(k);
           mw_s(i) = muw(k);
           In(i) = k;
           G(i) = 1./bot(k)/Rval.^4;
    end
end


%output total pressure drop relative to unchanged phim and resistances:
F4 = f4(pm(1));
DP(l,h) = trapz(z,G);
bot0 = ((f3(k0) + Da*F4(k0)).*RR.^4);
DP0 = 1./bot0;
Pmin(l,h) = min(pm - pbar);


%clf
figure;
yyaxis right
plot(z,RR,'r-','linewidth',1);
ylim([-0.2 1]);
hold on;
yyaxis left
plot(z,pbar,'b-','linewidth',1.5)
xlabel('z');ylabel('\phi')


% Save for every new Da
r_min = 1-Rc;
save('Steady_Rvary_rc_0.5_phim0_0.65_mu1_0.3.mat','z','Da_arr','DP','p_arr','r_min','Pmin');




    