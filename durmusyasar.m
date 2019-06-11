
%*************************************************
%  SISTEM SIMULASYONU - bm_sist_2d_01.m
%  Copyright : Ismail H. ALTAS, 2002
%*************************************************
clear;  % Hafizayi ve ekrani temizle
A=[ -1.1 -1;1 0 ];
B=[1 0; 1 0]; C=[ 1 0]; D=0;
% initials for simulation
x10=0; x20=0; u1=1; u2=1; dt=0.01;
tend=5; t0=0;
k=1; % k is dimension counter

% initials for simulation
x100=0; x20=1; u1=10.0; u2=10; dt=0.01;
tend=5; t0=0; k=1; % k is dimension counter.

%Baslangiç degerleri
U0=[u1;u2];    BOY=size(A);  LS=BOY(1); LK=BOY(2);
for n=1:LS
    x0(n)=0;
end
% -----------------Denetim için gerekli veriler
r0=input('Referans girisi degeri ==> ');
tend=input('Simülasyon bitis zamani ==> ');
%************************************************
% Bulanik denetleyici girisleri
% e, de ve du kesin uzaylarinin sinirlari
	EMAX  = r0;        EMIN  = -EMAX;
	DEMAX = EMAX/10;   DEMIN = -DEMAX;
	DUMAX = 1;         DUMIN = -1;
% Üyelik fonksiyonlari için veriler
   NLe=EMIN;   NTe=NLe;    NRe=0; 
   SLe=NTe;    STe=0;      SRe=EMAX;
   PLe=STe;    PTe=EMAX;   PRe=PTe;

   NLde=DEMIN;     NTde=NLde;     NRde=0;
   SLde=NTde;      STde=0;        SRde=DEMAX;
   PLde=STde;      PTde=DEMAX;    PRde=PTde;

   NLdu=DUMIN;       NTdu=DUMIN;    NRdu=0;
   SLdu=NTdu;        STdu=0;        SRdu=DUMAX;
   PLdu=STdu;        PTdu=DUMAX;    PRdu=PTdu;
%*************************************************
%  Hata ve denetleyici için baslangiç degerleri
   ee=EMAX;   dee=0;  E=EMAX;  dE=0;    
   e0=EMAX;   e(1)=0; e(2)=0;  de=e(2)-e(1); 
   C(1)=0;
%*************************************************
% membership matrix 
     DU=[ NTdu NTdu STdu
          NTdu STdu PTdu
          STdu PTdu PTdu ];
%*************************************************
% -------------- iteratif cözüm
while t0<tend-dt
      E=limiter(EMIN,EMAX,ee);%---------- limit E 
      % --------------------------- Bulaniklastirma
      EN=[NLe NTe NRe];
      ES=[SLe STe SRe];
      EP=[PLe PTe PRe];
      FSE=uyelik3(EN,ES,EP,E);% E için üyelik Degeri 
      DE=limiter(DEMIN,DEMAX,dee); %----- limit DE 
      DEN=[NLde NTde NRde];
      DES=[SLde STde SRde];
      DEP=[PLde PTde PRde];
      FSDE=uyelik3(DEN,DES,DEP,DE); % DE için ÜD 
     % ------------------------------Durulastirma
     N1=length(FSE);
     N2=length(FSDE);
     NM=N1*N2;
     nn=1;
     for mm=1:N1
        for qq=1:N2
            FSDU(nn)=min( [FSE(mm) FSDE(qq)] );
            nn=nn+1;
        end
     end
     nn=1;
     for mm=1:N1
         for qq=1:N2
             DDU(nn)=FSDU(nn)*DU(mm,qq);
             nn=nn+1;
         end
      end
      DUTOP1 = sum(DDU) ;            
      DUTOP2 = sum(FSDU);
      DV = (DUTOP1/DUTOP2);

     % ------------------------------- PI etki
     C(k+1) = C(k) + DV;
     CC=limiter(0,DUMAX,C(k+1));                    
     UU0 = CC*U0;
     [x]=runge(A,B,U0,x0,dt);%Runge ile denklem çöz.
     UU(k)=UU0(1);
     t(k)=t0+dt;    t0=t(k);
     r(k)=r0;       y(k)=x0(1);
     x1(k)=x(1); x2(k)=x(2);
     e(k)=r(k)-y(k);    de(k)=e(k)-e0;
     ee=e(k);  dee=de(k);   e0=e(k);   duty(k)=CC;
      for n=1:LS
           x0(n)=x(n);  XX(k,n)=x(n); 
      end
      k=k+1;
end;
% -------------------------- Grafikler
subplot(311)
plot(t,x1); xlabel('Zaman (sn)'); ylabel('x1');grid

subplot(312)
plot(t,x2);  xlabel('Zaman (sn)');  ylabel('x2'); grid

subplot(313)
plot(t,y,t,r,t,e); xlabel('Zaman (sn)'); ylabel('y');grid
