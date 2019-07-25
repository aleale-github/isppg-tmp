function [Stat,ProbTIR,ProbAuto,RendPerc,MeanYearTIR,BackDetTIR,ProbAutoBack1,ProbAutoBack2,RatioST1,OT,Moment] = PT_Standard_Long_Autocallable_Barrier6(Indice,Info_Cedole,Info,Tipo,Rf,Underline,NumScen)

% Inizializzazione Campionamento
    reset(RandStream.getGlobalStream) 
    
    Data = Indice(:,2);
    Today = x2mdate(Info_Cedole(1,1));
    Today5y = datenum(datetime(year(Today)-5,month(Today),day(Today)));
    DateHist = x2mdate(Indice(:,1));
    Check5yTime = sum(DateHist < Today5y)+1;
    [~,nIndex] = size(Data);
    PxRatio = Data(2:end,:)./Data(1:end-1,:);
    LogRend = log(PxRatio(Check5yTime:end,:));
    LogRend = LogRend-repmat(mean(LogRend),size(LogRend,1),1);
    NumSteps = wrkdydif(Today,x2mdate(Info_Cedole(end,1)))-1;
    OT = yearfrac(Today,x2mdate(Info_Cedole(Info_Cedole(:,2) == 1,1)),3);
    if Tipo == 'B'
        LogRendScen = 0;
        for t = 1:NumSteps
            LogRendScen_t = Bootstrap(LogRend,NumScen);            
            LogRendScen = LogRendScen+LogRendScen_t-mean(LogRendScen_t);
        end
    elseif Tipo == 'M'
        Rend = log(Data(2:end,:)./Data(1:end-1,:));
        S0 = Data(end,:);
        Mu = repmat(mean(Rend,1)',1,NumSteps);
        Q = zeros(nIndex,NumSteps);
        Sigma = repmat((std(Rend)*sqrt(252))',1,NumSteps);
        Rho = corrcoef(Rend);
        DT = repmat(1/252,nIndex,NumSteps);
        ST1 = Monte_Carlo1(S0,Mu,Q,Sigma,Rho,DT,NumScen,'Y',NumSteps);
        if (size(ST1,2) == 1) && (ndims(ST1) == 3)
            ST = permute(ST1,[1 3 2]);
        else
            ST = ST1;
        end        
        LogRendScen = log(ST(:,2)./ST(:,1));
    end
    PxRatioSim = exp(LogRendScen);
    Prc = 1;
%    meanRendSim = (mean(PxRatioSim)).^(1/OT);
    RendPerc = prctile(PxRatioSim,[99 95 75 50 25 5 1])';
    RendPerc = [RendPerc;mean(PxRatioSim)];    
    RendPercSim = prctile(PxRatioSim,[75 25 (100-Prc) Prc])';
    RendPercSim = [RendPercSim;mean(PxRatioSim)];
    %Annual = (RendPercSim).^(1/OT)+Rf-1;  
    R_tT = (RendPercSim).^(1/(OT(end)-OT(1)))+Rf-1;    
    Drift_R = [R_tT(1);Rf;R_tT(2);R_tT(3);R_tT(4)];        
    Stat = zeros(size(Drift_R,1),8);
 
 % Calcolo TIR deterministico
    Date_Bond = Info_Cedole(Info_Cedole(:,2) == 1,1);
    [DetTIR] = DetTirProd(Date_Bond,Info,Drift_R,OT);
    DetTIR = round(DetTIR,4);
    ProbTIR = zeros(size(Drift_R,1),8);      
    ProbAuto = zeros(1,length(Date_Bond)-Info(end)-1);
    Moment = zeros(6,size(Drift_R,1));
    DriftOT = prctile(LogRendScen,[75 50 25 (100-Prc) Prc])';    
    DriftOT(2) = Rf*(OT(end)-OT(1));
    
    for i = 1:size(DriftOT,1)
        [TIR,~,~,~,Px_Scenari,Autocall] = TIR_Standard_Long_Autocallable_Barrier(Indice(Check5yTime-1:end,:),Info_Cedole,Info,DriftOT(i),'B','G',NumScen);
        for j = 1:size(TIR,1)
            TIRxY = round(TIR{j,1}/0.74,4);
            ProbAuto(i,j) = size(TIRxY,1)./NumScen;
        end
        TIREndY = round(TIR{end,1}/0.74,4);       
        Stat(i,1) = prctile(TIREndY,50);
        Stat(i,2) = prctile(Px_Scenari(:,end).^(1/(OT(end)-OT(1)))-1,50); 
        PxScenari = Px_Scenari(Autocall(:,end) > 0,:);
        Stat(i,3) = mean(PxScenari(:,end).^(1/(OT(end)-OT(1)))-1,1);
        Stat(i,4) = mean(TIREndY(TIREndY < DetTIR(i),1),1);        
        Stat(i,5) = mean(TIREndY(TIREndY == DetTIR(i),1),1); 
        Stat(i,6) = mean(TIREndY(TIREndY > DetTIR(i),1),1);
        Stat(i,7) = mean(TIREndY(TIREndY < 0,1),1);        
        Stat(i,8) = mean(TIREndY(TIREndY >= 0,1),1);        
        ProbTIR(i,end-4) = sum(TIREndY < DetTIR(i))/size(TIREndY,1);
        ProbTIR(i,end-3) = sum(TIREndY == DetTIR(i))/size(TIREndY,1);
        ProbTIR(i,end-2) = sum(TIREndY > DetTIR(i))/size(TIREndY,1);
        ProbTIR(i,end-1) = sum(TIREndY < 0)/size(TIREndY,1);
        ProbTIR(i,end) = sum(TIREndY >= 0)/size(TIREndY,1);
        Moment(1,i) = mean(TIREndY,1);
        Moment(2,i) = std(TIREndY,1);
        Moment(3,i) = (Moment(1,i)-Rf)/Moment(2,i);
        Moment(4,i) = mean(Px_Scenari(:,end)).^(1/OT(end))-1;
        Moment(5,i) = std(Px_Scenari(:,end).^(1/OT(end))-1,1); 
        Moment(6,i) = (Moment(4,i)-Rf)/Moment(5,i);        
    end
    Stat(isnan(Stat)) = 0;
    ProbTIR(isnan(ProbTIR)) = 0;
    
% Calcolo del BackTesting
    [MeanYearTIR,BackDetTIR,ProbAutoBack1,ProbAutoBack2,RatioST1] = BackTestTirProd(Indice,Info,Info_Cedole,DetTIR);
      
    figure(1)
    subplot(2,1,1)
    %step1 = 1.1*round(min(RatioST1(:,2)*100,[],1):(abs(min(RatioST1(:,2)*100,[],1))+max(RatioST1(:,2)*100,[],1))/5:max(RatioST1(:,2)*100,[],1));
    step1 = round(linspace(min(RatioST1(:,2)*100,[],1)*0.9,max(RatioST1(:,2)*100,[],1)*1.1,6));
    %step2 = round(min(Indice(1:size(RatioST1,1),2),[],1)*0.9:(min(Indice(1:size(RatioST1,1),2),[],1)*0.9+max(Indice(1:size(RatioST1,1),2),[],1)*1.1)/15:max(Indice(1:size(RatioST1,1),2),[],1)*1.1);
    step2 = round(linspace(min(Indice(1:size(RatioST1,1),2),[],1)*0.9,max(Indice(1:size(RatioST1,1),2),[],1)*1.1,6));    [x] = plotyy(RatioST1(:,1),RatioST1(:,2:3)*100,RatioST1(:,1),Indice(1:size(RatioST1,1),2));
    a = x(1).XTick;
    a(1) = RatioST1(1,1); a(end) = RatioST1(end,1);
    x(1).XTick = a;
    x(2).XTick = a;
    xlim(x(1),[RatioST1(1,1) RatioST1(end,1)]);
    xlim(x(2),[RatioST1(1,1) RatioST1(end,1)]);
    ylim(x(1),[min(RatioST1(:,2)*100,[],1)*1.1 max(RatioST1(:,2)*100,[],1)*1.1]);
    ylim(x(2),[min(Indice(1:size(RatioST1,1),2),[],1)*0.9 max(Indice(1:size(RatioST1,1),2),[],1)*1.1]);    
    x(1).YTick = step1;
    x(2).YTick = step2;    
    datetick('x','mmmyy','keepticks') 
    xlabel('Date'); ylabel(x(1),'IRR (%)'),ylabel(x(2),'Und Level')
    grid(x(1),'on')
    title('Back Testing of Certificate')
    legend(['IRR Underline' 'IRR Certificate' Underline],'Location','Best')
    hold on
    plot(RatioST1(:,1),zeros(size(RatioST1,1),1),'black')    
    hold off  
    subplot(2,1,2)
    bar(RatioST1(:,1),RatioST1(:,4:end));
    ax = gca;
    ax.XLim = [RatioST1(1,1) RatioST1(end,1)];    
    datetick('x','mmmyy','keepticks') 
    xlabel('Date'); ylabel('Months')
    grid on
    title('Autocallabilty Time')
    
function [MeanYearTIR,BackDetTIR,ProbAutoBack1,ProbAutoBack2,RatioST1] = BackTestTirProd(Indice,Info,Info_Cedole,DetTIR)    
    
% Calcolo di Backtesting    
    HistDate = x2mdate(Indice(:,1));
    Info_Cedole(:,1) = x2mdate(Info_Cedole(:,1));
    Date = Info_Cedole(Info_Cedole(:,2) == 1,1);
    NumDays = cumsum(daysdif(Date(1:end-1),Date(2:end)));  
% Variabili ad uso interno
    C1 = Info(1,1);
    C2 = Info(2,1);
    C3 = Info(3,1);
    B1 = Info(4,1);
    B2 = Info(5,1);
    B3 = Info(7,1);
    NumCedoleFixed = Info(8,1); 
% Costruzione del payoff
    Und = zeros(size(NumDays,1)+1,1);
    Und(1) = -1;
    %if NumCedoleFixed > 0
    %    Und(2:1+NumCedoleFixed,1) = C1;
    %end
    BackTestTIR = zeros(1,2);
    ST = zeros(1,1);
    NumAutoBack = size(NumDays,1)-NumCedoleFixed;         
    ProbAutoBack = zeros(1,NumAutoBack);
    ProbAutoBack1 = zeros(1,NumAutoBack);
    ProbAutoBack2 = zeros(1,NumAutoBack);
    PointerAutoBack = 1:1:NumAutoBack;
    RatioST1 = zeros(1,3+NumAutoBack);    
    for t = 1:size(HistDate,1)
        Cedole = zeros(1,1);
        Cedole(1,1) = -1;     
        %Cedole(2,1) = C1;
        HistStartDate = HistDate(t,1);
        ST0 = Indice(t,2);
        HistDateBond = HistStartDate+NumDays;
        if HistDateBond(end) > HistDate(end) 
            break
        end        
        for i = 1:size(HistDateBond,1)
            check = find(HistDate >= HistDateBond(i),1);
            ST(i,1) = Indice(check(1),2);
        end
        RatioST = ST./ST0;
        
        if NumCedoleFixed > 0
            for k = 1:NumCedoleFixed 
                Cedole(1+k,1)=C1*(RatioST(k)>=B2);
            end    
           Und(2:1+NumCedoleFixed,1) = C1;
        end
        
        
        for A = NumCedoleFixed+1:size(NumDays,1)
            if (RatioST(A) >= B1) && A < size(RatioST,1)
                Cedole(1+A,1) = C3+1;
                %Date = [HistStartDate;HistDateBond(1:A,1)];
                BackTestTIR(t,1) = HistStartDate;
                BackTestTIR(t,2) = xirr(Cedole,Date(1:1+A,1),0.03,1000,3);
                
                %Und(A+1+NumCedoleFixed) = RatioST(A);
                Und(A+1) = RatioST(A);
                Und1 = [Date,Und];
                Und1(Und1(:,2) == 0,:) = [];                
                RatioST1(t,1) = HistStartDate;                
                RatioST1(t,2) = xirr(Und1(:,2),Und1(:,1),0.03,1000,3);
                RatioST1(t,3) = BackTestTIR(t,2);          
                
                for j = 1:NumAutoBack
                    if A == j+NumCedoleFixed
                        ProbAutoBack(t,j) = 1;
                    end
                end
                RatioST1(t,4:end) = ProbAutoBack(t,:).*PointerAutoBack;
                Und(2:end,:) = 0;                
                break
            elseif (RatioST(A) < B1) && (RatioST(A) >= B2) && A < size(RatioST,1)
                if (C2 == 0) || (C1 > 0)
                    Cedole(1+A,1) = C1;
                else
                    Cedole(1+A,1) = C2*(RatioST(A) >= B2);
                end
            elseif (RatioST(A) < B2) && A < size(RatioST,1)
                Cedole(1+A,1) = 0;
            end
            if (A == size(RatioST,1))
                if RatioST(A) < B3
                    ProbAutoBack(t,end) = 1;                    
                    Cedole(1+A,1) = RatioST(A,1);
                    %Date = [HistStartDate;HistDateBond(1:A,1)];
                    BackTestTIR(t,1) = HistStartDate;
                    BackTestTIR(t,2) = xirr(Cedole,Date(1:1+A,1),0.03,1000,3); 

                    %Und(A+1+NumCedoleFixed) = RatioST(A);
                    Und(A+1) = RatioST(A);
                    Und1 = [Date,Und];
                    Und1(Und1(:,2) == 0,:) = [];                
                    RatioST1(t,1) = HistStartDate;                
                    RatioST1(t,2) = xirr(Und1(:,2),Und1(:,1),0.03,1000,3);
                    RatioST1(t,3) = BackTestTIR(t,2);
                    RatioST1(t,4:end) = ProbAutoBack(t,:).*PointerAutoBack;
                    Und(2:end,:) = 0;                    
                elseif (RatioST(A) >= B2)
                    ProbAutoBack(t,end) = 1+C3;                     
                    Cedole(1+A,1) = 1;
                    BackTestTIR(t,1) = HistStartDate;
                    BackTestTIR(t,2) = xirr(Cedole,Date(1:1+A,1),0.03,1000,3);    

                    %Und(A+1+NumCedoleFixed) = RatioST(A);
                    Und(A+1) = RatioST(A);
                    Und1 = [Date,Und];
                    Und1(Und1(:,2) == 0,:) = [];                
                    RatioST1(t,1) = HistStartDate;                
                    RatioST1(t,2) = xirr(Und1(:,2),Und1(:,1),0.03,1000,3);
                    RatioST1(t,3) = BackTestTIR(t,2);
                    RatioST1(t,4:end) = ProbAutoBack(t,:).*PointerAutoBack;
                    Und(2:end,:) = 0;
                end
            end
        end
    end
    BackTestTIR(:,2) = round(BackTestTIR(:,2),4);
    YearStart = year(BackTestTIR(1,1));
    YearEnd = year(HistStartDate);
    NumYear = YearEnd-YearStart;
    MeanYearTIR = zeros(1,7);
    a = 1;
    for Y = 0:NumYear
        YearTIR = round(BackTestTIR(year(BackTestTIR(:,1)) == YearStart+Y,2),4);
        MeanYearTIR(a,1) = round(mean(YearTIR,1),4);
        MeanYearTIR(a,2) = round(min(YearTIR,[],1),4);
        MeanYearTIR(a,3) = round(max(YearTIR,[],1),4);
        MeanYearTIR(a,4) = sum(YearTIR < MeanYearTIR(a,1),1)/size(YearTIR,1);        
        MeanYearTIR(a,5) = sum(YearTIR == MeanYearTIR(a,1),1)/size(YearTIR,1);
        MeanYearTIR(a,6) = sum(YearTIR > MeanYearTIR(a,1),1)/size(YearTIR,1);
        MeanYearTIR(a,7) = sum(YearTIR < 0,1)/size(YearTIR,1);
        MeanYearTIR(a,8) = sum(YearTIR >= 0,1)/size(YearTIR,1);
        MeanYearTIR(a,9) = size(YearTIR,1);
        TotAutoBack = sum(sum(ProbAutoBack(year(BackTestTIR(:,1)) == YearStart+Y,:),1),2);
        ProbAutoBack1(a,:) = sum(ProbAutoBack(year(BackTestTIR(:,1)) == YearStart+Y,:),1)./TotAutoBack;
        a = a+1;
    end
    BackDetTIR = zeros(size(DetTIR,1),2);
    BackDetTIR(:,1) = DetTIR;
    TotAutoBack1 = sum(sum(ProbAutoBack,1),2);
    for D = 1:size(DetTIR,1)
        BackDetTIR(D,2) = sum(BackTestTIR(:,2) < DetTIR(D,1),1)/size(BackTestTIR,1);          
        BackDetTIR(D,3) = sum(BackTestTIR(:,2) == DetTIR(D,1),1)/size(BackTestTIR,1);  
        BackDetTIR(D,4) = sum(BackTestTIR(:,2) > DetTIR(D,1),1)/size(BackTestTIR,1);
        BackDetTIR(D,5) = sum(BackTestTIR(:,2) < 0,1)/size(BackTestTIR,1);
        BackDetTIR(D,6) = sum(BackTestTIR(:,2) >= 0,1)/size(BackTestTIR,1);        
        ProbAutoBack2(D,:) = sum(ProbAutoBack,1)./TotAutoBack1;        
    end
        
function [TIR] = DetTirProd(Date,Para,ST,OT)
% Date Matlab
    Date(:,1) = x2mdate(Date(:,1),0);
% Variabili ad uso interno
    C1 = Para(1,1);
    C2 = Para(2,1);
    C3 = Para(3,1);
    B1 = Para(4,1);
    B2 = Para(5,1);
    B3 = Para(7,1);
    NumCouponFix = Para(8,1);
    T = OT(2:end,1)-OT(1,1);
    ST = (repmat(ST',size(T,1),1)+1).^repmat(T,1,size(ST',2));
% Costruzione del payoff
    Cedole = zeros(size(Date,1),1);
    Cedole(1) = -1;
    %if (NumCouponFix > 0)
    %    Cedole(2:NumCouponFix+1) = C1;
    %end
% Calcolo del TIR
    TIR = zeros(size(ST,2),1);    
% Calcolo del TIR        
    for i = 1:size(ST,2)
        if (C2 == 0) && (NumCouponFix == 0)
            Cedole(2:end,:) = C1;
        else 
            Cedole(2:end,:) = (ST(:,i) >= B2).*C1;
        end
        if (ST(end,i) >= B3)
            Cedole(end) = 1+Cedole(end);
        elseif ST(end,i) < B3
            Cedole(end) = ST(end,i);
        end
        TIR(i) = xirr(Cedole,Date,0.03,1000,3);
    end
    
function ST = Monte_Carlo1(S0,Mu,Q,Sigma,Rho,DT,Scenari,Opzioni,Steps)

% Funzione che genara percorsi possibili per i prezzi di un titolo.

% INPUT
%        S0 = (N,1), valore in t0 degli Asset
%        Mu = (N,P), N drift degli Asset per P tenor (Struttura a termine)
%         Q = (N,P), N Dividendi degli Asset per P tenor (Struttura a
%             termine)
%     Sigma = (N,P), matrice delle volatilità dove N corrisponde al numero
%             di Asset e P al numero di Tenor
%       Rho = (N,N), matrice delle correlazioni tra gli Asset
%       dtT = (N,P), tempo tra P passi (in frazione di anno)
%   Scenari = (1,1), numero di Scenari
%   Opzioni = (1,1), se 'Y' aggiunge ad ST il passo t0, se 'N' no

% OUTPUT
%        ST = (S,T,N), matrice degli (S) Scenari, per i (T) passi, per gli
%        N Titoli

% Inizializzazione Campionamento
    reset(RandStream.getGlobalStream)    

% Dimensionamento delle variabili
    [Num_Titoli] = size(S0,1);
    [Num_Passi] = size(DT,2);
    [Num_Scenari] = Scenari;
% Crezione dello Spazio di memoria
    ST = zeros(Num_Scenari,size(Steps,1),Num_Titoli);
    S0 = repmat(S0',Num_Scenari,1);
    if Opzioni(1,1) == 'Y'    
        S01 = S0;
    end
    Ch = chol(Rho)';
    %dtT_tot = repmat(DT',Num_Titoli,1);
%    dtT_tot = DT';
% Simulazione Monte Carlo
    t = 1;
    for P = 1:Num_Passi
    % Drift    
        mu = Mu(:,P)';
        mu = repmat(mu,Num_Scenari,1);
    % Dividendi    
        q = Q(:,P)';
        q = repmat(q,Num_Scenari,1);
    % Passi    
        dtT = repmat(DT(:,P)',Num_Scenari,1);
    % Sigma
        sigma = Sigma(:,P)';
        sigma = repmat(sigma,Num_Scenari,1);
    % Campionamento su asset con Correlazione Rho;
        %W = randn(Num_Scenari,Num_Titoli);
        W = norminv(LHSSample(Num_Scenari,Num_Titoli));
        W = Ch*W';            
    % Brownian Motion
 %       St = S0.*exp((mu-q-sigma.^2.*0.5).*dtT+sigma.*sqrt(dtT).*W');
        St = S0.*exp(+sigma.*sqrt(dtT).*W');
        S0 = St;
        if P == Steps(t)
            ST(:,t,:) = St;
            t = t+1;
        end
    end 
    
% Aggiunge lo scenario di partenza    
    if Opzioni(1,1) == 'Y'
        ST = permute(ST,[1 3 2]);
        ST = cat(3,S01,ST);
        ST = permute(ST,[1 2 3]);
    end

function [TIR,TIR_Perc,Stat,VFR_Barr,Px_Scenari,Autocall] = TIR_Standard_Long_Autocallable_Barrier(Indice,Info_Cedole,Info,Drift,Modo,Frequenza,Scenari)

% Viene calcolato contemporaneamente il TIR del contratto ed 
% effettuato il confronto con un BTP di pari scadenza.
%
% Cedola
% if
%   1Y    se  B2*Index_0 <= Index_t1 < B1*Index_0     => C2
%   2Y    se B2*Index_0 <= Index_t1 < B1*Index_0      => C2
% else
%   1Y    C1
%   2Y    C1
%   3Y    se Index_t3 <= B2*Index_0                   => Index_t3
%         se Index_t3 >= B2*Index_0                   => 1 + C3

% Rimborso anticipato
% 1Y    se Index_t1 >= B1*Index_0 => 1
% 2Y    se Index_t2 >= B1*Index_0  => 1


% Index1_t = valore dell'Index1 alla data cedolare
% Index1_0 = valore dell'Index1 alla data di partenza del contratto

% INPUT:
%           Indice = (N,M), storico degli M titoli dell'Indice
%      Date_Cedole = (P,2), riporta le date corrispondenti alle cedole del
%                    contratto. La seconda colonna indica quali sono date di pagamento (1) e
%                    quali date di fixing (2).
%            Drift = (1,1), se 0 si usa il drift storico, altrimenti il drift è il
%                    valore che si inserisce
%             modo = (1,1), se 'A' campionamento semplice, se 'H' Latin
%                    Hypercube
%             Info = (5,1), dove (1,1) tasso fisso C1, (2,1) tasso fisso C2,
%                    (3,1) Partecipazione P1, (4,1) Partecipazione P2, (5,1) TIR del BTP

% OUTPUT:
%     Payoff = (1,1), Struttura di Payoff del contratto
%        TIR = (S,1), TIR del contratto simulati per ogni scenario
%   TIR_Perc = (7,1), percentili della PDF dei TIR del contratto all'1%, 5%, 25%, 50%, 75%, 95%, 99%
% Statistica = (4,1), dove (1,1) TIR max, (2,1) TIR medio, (3,1) TIR minimo, (4,1) Probabilità di battere il BTP



tic
% Inizializzazione Campionamento
    reset(RandStream.getGlobalStream)

% Conversione Date da Excel -> Matlab
    Indice(:,1) = x2mdate(Indice(:,1),0);
    
% Variabili ad uso interno
    C1 = Info(1,1);
    C2 = Info(2,1);
    C3 = Info(3,1);
    B1 = Info(4,1);
    B2 = Info(5,1);
    B3 = Info(7,1);
    Num_Cedole_Fix = Info(8,1);
    Px_t0 = Indice(end,2:end);
    
    Date_Bond = x2mdate(Info_Cedole(Info_Cedole(:,2) == 1,1));
    %Date_Bond1 = Info_Cedole(Info_Cedole(:,2) == 1,1);
    Date_Opz = x2mdate(Info_Cedole(Info_Cedole(:,3) == 2,1));

    TIR_BTP = zeros(size(Date_Bond,1)-1-Num_Cedole_Fix,1);    
% Calcolo delle Date_Cedole
    if Frequenza == 'G'
        Passi = wrkdydif(Date_Opz(1:end-1,1),Date_Opz(2:end,1))-1;
        Passi(Passi(:,1) < 1) = 1;
        Passi = cumsum(Passi);        
     elseif Frequenza == 'M'
    %    Passi = months(Date_Opz(1:end-1,1),Date_Opz(2:end,1));
        Passi = (12:12:12*5)';
    end
    
 % Costruzione Date per drift Cambi
    if size(Drift,1) >= 1
        mesi = (12:12:12*(size(Date_Bond,1)-1))';
        Date_Drift = [Date_Bond(1);datemnth(Date_Bond(1),mesi)];
        Passi_Drift = wrkdydif(Date_Drift(1:end-1,1),Date_Drift(2:end,1))-1;
        Passi_Drift(Passi_Drift <= 0) = 1;
        Passi_Drift_Cum = cumsum(Passi_Drift);
        if Passi_Drift_Cum(end) < Passi(end) 
            Passi_Drift_Cum(end) = Passi(end);
        end
    end   
    Passi_Tot = Passi(end,1);
    [Num_Cedole] = size(Date_Bond,1);
    [Num_Cedole_Opz] = size(Date_Opz,1)-1;
    [Num_Titoli] = size(Indice,2)-1;
    
    
% Creazione dello spazio di memoria    
    Px = zeros(Scenari,Num_Titoli,Num_Cedole_Opz);
    Cedole = zeros(Scenari,Num_Cedole);
%    check = zeros(Scenari,Num_Cedole-1);
%    check_VRI = zeros(Scenari,1);
%    check_VRF = zeros(Scenari,1);
% Dimensionamento delle variabili
    [r_Indice] = size(Indice,1);
% Vettore del rapporto tra i prezzi (giornaliero)
    LogRend = log(Indice(2:end,2:end)./Indice(1:end-1,2:end));
    LogRend = LogRend-repmat(mean(LogRend),size(LogRend,1),1); 
    
%    I_nd = cumprod(Px_Ratio);
    if Drift(1,1) ~= 0
        if size(Drift,1) == 1
            if Frequenza == 'G'
                %Px_Ratio = Px_Ratio./repmat(mean(Px_Ratio),r_Indice-1,1).*(1+Drift).^(1/250);
                LogRend = LogRend-repmat(mean(LogRend),r_Indice-1,1)+(1+Drift).^(1/250);
            elseif Frequenza == 'M'
                LogRend = LogRend./repmat(mean(LogRend),r_Indice-1,1).*(1+Drift).^(1/12);
            end
        end
    end
    c = 1; 
    LogRendScen = 0;
    DriftS = Drift/Passi_Drift_Cum(end);   
    if Modo == 'B'      % bootstrap semplice
        for i = 1:Passi_Tot
            LogRendScen_t = Bootstrap(LogRend,Scenari);
            LogRendScen = LogRendScen+LogRendScen_t-mean(LogRendScen_t);            
            LogRendScen = LogRendScen-repmat(mean(LogRendScen),Scenari,1)+DriftS*i; 
            if i == Passi(c)
                Px(:,:,c) = Px_t0.*exp(LogRendScen);                
                c = c+1;     
            end
            Px_Scenari = exp(LogRendScen);
        end
    elseif Modo == 'H'       % ipercubi latini
        Hyper = zeros(Scenari,Passi_Tot);
        Hyper(:,:) = LHSSample(Scenari,Passi_Tot);
        [r_sample] = size(LogRend,1);
        Hyper = 1+floor(Hyper*r_sample);
        for i = 1:Passi_Tot
            LogRendScen_t = LogRend(Hyper(:,i),:);
            LogRendScen = LogRendScen+LogRendScen_t-mean(LogRendScen_t);            
            LogRendScen = LogRendScen-repmat(mean(LogRendScen),Scenari,1)+DriftS*i;          
            if i == Passi(c)
                Px(:,:,c) = Px_t0.*exp(LogRendScen);
                c = c+1;     
            end
            Px_Scenari = exp(LogRendScen);            
        end
    end

% Conversione a 2 dimensioni
    if Num_Titoli == 1
        Sim(:,:) = Px(:,1,:);
    else
        Sim(:,:,:) = Px;
    end
    
    VRI = Sim(:,1);
    VRF_T = Sim(:,end-size(Passi,1)+2:end);
   
% Cedole del Bond e calcolo del TIR
    Cedole(:,1) = -1;
    if Num_Cedole_Fix > 0 
        for p = 1:Num_Cedole_Fix 
            Cedole(:,1+p)= C1 * (VRF_T(:,p) >= B2*VRI); 
        end  
    end 
    
    TIR = cell(Num_Cedole-Num_Cedole_Fix-1,1);
    TIR_Perc = zeros(7,Num_Cedole-Num_Cedole_Fix-1);    
    TIR_Max = zeros(1,Num_Cedole-Num_Cedole_Fix-1);        
    TIR_Media = zeros(1,Num_Cedole-Num_Cedole_Fix-1);            
    TIR_Min = zeros(1,Num_Cedole-Num_Cedole_Fix-1);
    Prob_BTP = zeros(1,Num_Cedole-Num_Cedole_Fix-1);
    Prob_Y = zeros(1,Num_Cedole-Num_Cedole_Fix-1);
    TIR_Mediana = zeros(1,Num_Cedole-Num_Cedole_Fix-1);    
    check =zeros(Scenari,Num_Cedole-Num_Cedole_Fix-1);      
    check1 =zeros(Scenari,Num_Cedole-Num_Cedole_Fix-1);  
    
    for j = 1:Num_Cedole-2-Num_Cedole_Fix
        check(VRF_T(:,j+Num_Cedole_Fix) < B2*VRI,j) = 1;
        check1(VRF_T(:,j+Num_Cedole_Fix) < B1*VRI,j) = 1;
        Cedole(:,j+1+Num_Cedole_Fix) = C3.*(VRF_T(:,j+Num_Cedole_Fix) >= B2*VRI);
        Cedole(:,j+1+Num_Cedole_Fix) = Cedole(:,j+1+Num_Cedole_Fix)+ (1-check1(:,j));
        if j > 1
            Cedole(sum(Cedole(:,Num_Cedole_Fix+1:j+1+Num_Cedole_Fix),2) > 1,j+1+Num_Cedole_Fix) = 0;
        end
        
    % Calcolo del TIR
        TIR_C = [Cedole(:,1),Cedole(:,2:j+1+Num_Cedole_Fix)];
        %TIR_C(TIR_C(:,end) ~= 1+C3,:) = [];
        TIR_C(TIR_C(:,1+j+Num_Cedole_Fix) < 1+C1,:) = [];
        if isempty(TIR_C) == 1
            TIR_Perc(:,j) = zeros(7,1);
            TIR_Max(1,j) = 0;
            TIR_Media(1,j) = 0;
            TIR_Min(1,j) = 0;
            Prob_BTP(1,j) = 0;
            Prob_Y(1,j) = 0;
            TIR_Mediana(1,j) = 0;
            TIR{j,1} = [];
        elseif isempty(TIR_C) == 0
            %TIR_C(:,end) = TIR_C(:,end)+1;
            TIR1 = zeros(size(TIR_C,1),1);         
            for i = 1:size(TIR_C,1)
                TIR1(i,1) = xirr(TIR_C(i,:)',[Date_Bond(1,:);Date_Bond(2:j+1+Num_Cedole_Fix,:)],0.03,1000,3)*0.74;
            end
            TIR{j,1} = TIR1;
        % Statistica
            TIR_Perc(:,j) = prctile(TIR1,[1 5 25 50 75 95 99],1)';
            TIR_Max(1,j) = TIR_Perc(end,j);
            TIR_Media(1,j) = mean(TIR1,1);
            TIR_Min(1,j) = min(TIR1,[],1);
            Prob_BTP(1,j) = sum(TIR1 > TIR_BTP(j,1),1)./size(TIR1,1);
            Prob_Y(1,j) = size(TIR1,1)./Scenari; 
            TIR_Mediana(1,j) = prctile(TIR1,50,1);
        end
    end
    
    if j == Num_Cedole-2-Num_Cedole_Fix
        check_VRI = zeros(Scenari,1);
        check_VRF = zeros(Scenari,1);
        check_End = sum(check1,2);
        check_VRI(check_End == Num_Cedole-2-Num_Cedole_Fix,1) = VRF_T(check_End == Num_Cedole-2-Num_Cedole_Fix,end) >= B2.*VRI(check_End == Num_Cedole-2-Num_Cedole_Fix);
        check_VRF(check_End == Num_Cedole-2-Num_Cedole_Fix,1) = VRF_T(check_End == Num_Cedole-2-Num_Cedole_Fix,end) < B2.*VRI(check_End == Num_Cedole-2-Num_Cedole_Fix);
        Fix_C3 = repmat(1+C3,Scenari,1);
        Cedole(check_VRI == 1,j+2+Num_Cedole_Fix) = Fix_C3(check_VRI == 1);
        Cedole(check_VRF == 1,j+2+Num_Cedole_Fix) = VRF_T(check_VRF == 1,end)./VRI(check_VRF == 1);
        Autocall = Cedole(:,Num_Cedole_Fix+2:end);
        
    % Calcolo del TIR
        TIR_C = Cedole(:,1:end);
        TIR_C(TIR_C(:,end) == 0,:) = [];
        VFR_Barr = sum(TIR_C(:,end) < B2)/Scenari;
        if isempty(TIR_C) == 1
            TIR_Perc(:,j+1) = zeros(7,1);
            TIR_Max(1,j+1) = 0;
            TIR_Media(1,j+1) = 0;
            TIR_Min(1,j+1) = 0;
            Prob_BTP(1,j+1) = 0;
            Prob_Y(1,j+1) = 0;
            TIR_Mediana(1,j+1) = 0;
            TIR{j+1,1} = [];
        elseif isempty(TIR_C) == 0
            TIR1 = zeros(size(TIR_C,1),1);         
            for i = 1:size(TIR_C,1)
                TIR1(i,1) = xirr(TIR_C(i,:)',[Date_Bond(1,:);Date_Bond(2:end,:)],-0.24,1000,3)*0.74;
            end
            TIR{j+1,1} = TIR1;
        % Statistica
            TIR_Perc(:,j+1) = prctile(TIR1,[1 5 25 50 75 95 99],1)';
            TIR_Max(1,j+1) = TIR_Perc(end,j+1);
            TIR_Media(1,j+1) = mean(TIR1,1);
            TIR_Min(1,j+1) = TIR_Perc(1,j+1);
            Prob_Y(1,j+1) = size(TIR1,1)./Scenari;            
            Prob_BTP(1,j+1) = sum(TIR1 > TIR_BTP(j+1,1),1)./size(TIR1,1);
        % Statistica su tutta la distribuzione a scadenza
            TIR2 = [TIR{j+1,1}];
            TIR_Mediana(1,j+1) = prctile(TIR2,50,1);            
            Prob_0(1,j+1) = sum(TIR2 > TIR_BTP(j+1,1),1)./size(TIR2,1);
            Prob_0(1,j+1) = sum(TIR2 < -1.53222779501094e-10,1)./size(TIR2,1);                
        end        
    end
    Stat = [TIR_Max; TIR_Media; TIR_Min; Prob_BTP; Prob_Y; Prob_0; TIR_Mediana];    
%    hist({TIR(j+1,1)})

toc 