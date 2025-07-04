% plots steady state for imposed phim(z)
clear all

z = linspace(0,1,1000);

mu1 = 0.3;
p_arr = 0.3;
phim0 = 0.8;
plen = length(p_arr);

c = 0.2;
Da_arr = 10^-(1); 
Da = Da_arr(h); 

p0 = p_arr;
d=10;
b = c/2; 
a = phim0-b;
pm = a + b*(erf(-(z-0.5).*d));


RR = 1;
muw = linspace((1+1e-10)*mu1,50*mu1,1e5);
m = muw;
mw=m;
mu=mu1;
rc = mu1./muw;
f1 = (2.*rc)./m - 2.*((2.*rc)./(3.*m.^(1./2)) - (2.*(m.*rc + 1))./m.^(3./2)).*(1 - rc).^(1./2) + (4.*(1 - rc).^(1./2))./(3.*m.^(1./2)) - 2./m + rc.^2 - (pi.*(m.*rc + 1).*2i)./m.^2 + (2.*log(-(m - m.*rc - 2.*m.^(1./2).*(1 - rc).^(1./2) + 1)./(m.*rc - m + 1)).*(m.*rc + 1))./m.^2 - (2.*log(1 - (m.*rc + 1)./m).*(m.*rc + 1))./m.^2 + (2.*log(rc - (m.*rc + 1)./m).*(m.*rc + 1))./m.^2;
f3 = (rc.^2.*(rc - 1./2))./2 - rc./6 + (rc.^2.*(rc - 1).^2)./4 - (5.*rc.^4)./24 + 1./8;
f2 = (rc.^2.*(rc - 1).^2)./4 - (rc - 1).^2./(2.*m) - (((2.*rc.*(rc - 1).^2)./(3.*m.^(1./2)) - (2.*(m.*rc + 1).*(rc - 1).^2)./m.^(3./2)).*(1 - rc).^(1./2))./2 + (m.*(rc./m.^3 + 1./m.^4))./2 + (m - m.*rc).^3./(6.*m.^4) - (m - m.*rc).^(7./2)./(7.*m.^4) + (1 - rc).^(5./2)./(3.*m.^(1./2)) - ((m - m.*rc).^(1./2).*((2.*rc)./m.^3 + 2./m.^4))./2 + ((m - m.*rc).^2.*(rc./(2.*m.^3) + 1./(2.*m.^4)))./2 - ((m - m.*rc).^(3./2).*((2.*rc)./(3.*m.^3) + 2./(3.*m.^4)))./2 - ((m - m.*rc).^(5./2).*((2.*rc)./(5.*m.^3) + 2./(5.*m.^4)))./2 + (log((m - m.*rc).^(1./2) + 1).*(2.*m.*rc + 2))./(2.*m.^4) - (m.*rc.*(rc./m.^3 + 1./m.^4))./2 + (rc.*(rc - 1).^2)./(2.*m) + (log(rc - (m.*rc + 1)./m).*(m.*rc + 1).*(rc - 1).^2)./(2.*m.^2) - (pi.*(m.*rc + 1).*(rc - 1).^2.*1i)./(2.*m.^2) + (log(-(m - m.*rc - 2.*m.^(1./2).*(1 - rc).^(1./2) + 1)./(m.*rc - m + 1)).*(m.*rc + 1).*(rc - 1).^2)./(2.*m.^2) - (log(1 - (m.*rc + 1)./m).*(m.*rc + 1).*(rc - 1).^2)./(2.*m.^2);
f4 = @(phim)-(16.*mu.^2.*(mw - mu).^(1./2) - 24.*mw.^2.*(mw - mu).^(1./2) + 60.*phim.^3.*(mw - mu).^(1./2) + 15.*mu.*mw.^2 + 30.*mu.*phim.^3 + 45.*mw.^2.*phim - 30.*mw.*phim.^3 - 5.*mu.^3 - 15.*mw.^2 - 10.*mw.^3 - 30.*phim.^3.*log((mw - mu + 2.*(mw - mu).^(1./2) + 1)) + 15.*mu.^2.*phim.^3 - 45.*mw.^2.*phim.^2 + 8.*mu.*mw.*(mw - mu).^(1./2) - 30.*mu.*phim.^3.*log((mw - mu + 2.*(mw - mu).^(1./2) + 1)) - 24.*mu.^2.*phim.*(mw - mu).^(1./2) + 40.*mu.*phim.^3.*(mw - mu).^(1./2) + 36.*mw.^2.*phim.*(mw - mu).^(1./2) + 20.*mw.*phim.^3.*(mw - mu).^(1./2) - 12.*mu.*mw.*phim.*(mw - mu).^(1./2))./(15.*mw.^2.*phim.^2);
 
minpm = phim0 - c;
bot = (f3 + Da*f4(minpm));
top = minpm*f2;
flux1 = top./bot;

maxflux = max(flux1);

%find mu0
[~,k0] = min(abs(phim0*f1-p0));
%find initial flux
bot = (f3 + Da*f4(pm(1)));
top = pm(1)*f2;
flux = top./bot;
flux0 = flux(k0);

%NO CLOG
if(maxflux>flux0)
    for i=1:length(z)
        pmm = pm(i);
        bot = (f3 + Da*f4(pmm));
        top = pmm*f2;
        flux = top./bot;
        [~,j] = max(flux);                  % Get the lower phi branch/ higher mu branch closer to inital condition
        if(k0<j)
        [~,k] = min(abs(flux(1:j)-flux0));  % gives index relative to j
        pbar(i) = pmm*f1(k);
        u=k;
        else
        [~,k] = min(abs(flux(j:end)-flux0));% gives index relative to j
        pbar(i) = pmm*f1(k+j-1);            % to get accurate index 
        u=k+j-1;
        end
        G(i) = 1./bot(u)/RR.^4;
    end

else %CLOG
    for i=1:length(z)
        pmm = pm(i);
        bot = (f3 + Da*f4(pmm));
        top = pmm*f2;
        flux = top./bot;
        flux(1)=0;
        [~,j] = max(flux);
        [~,k] = min(abs(flux(1:j)-maxflux));% higher phi vs flux branch
        pbar(i) = pmm*f1(k);
        G(i) = 1./bot(k)/RR.^4;
    end
end


figure;
plot(z,pm,'r--','linewidth',1);
hold on;
plot(z,pbar,'b-','linewidth',1.5)
xlabel('z');ylabel('\phi')

%output total pressure drop relative to unchanged phim:
DP(l,h) = trapz(z,G);
F4 = f4(pm(1));
bot0 = (f3(k0) + Da*F4(k0)).*RR.^4;
DP0 = 1./bot0;
Pmin(l,h) = min(pm - pbar);


phim_min = 0.9-c_arr;
save('Steady_phimvary_c_0.2_phim0_0.8_mu1_0.3.mat','z','Da_arr','DP','c_arr','phim_min','Pmin');



  