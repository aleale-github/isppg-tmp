function [TIR,TIR_Perc,Stat,VFR_Barr] = TIR_Certificate_Digital_Autocallable_Barrier_v1(Indice,Info_Cedole,Info,Drift,Modo,Frequenza,Scenari)

% Viene calcolato contemporaneamente il TIR del contratto ed 
% effettuato il confronto con un BTP di pari scadenza.
%
% Cedola
% 1Y    C1
% 2Y    se Index_t2 < B1*Index_0 => C2
% 3Y    se Index_t3 < B1*Index_0 => C2

% Rimborso anticipato
% 1Y    se Index_t1 >= B1*Index_0 => 1+C3 
% 2Y    se Index_t2 >= B1*Index_0  => 1+C3*1
% 3Y    se Index_t3 >= B2*Index_0  => 1+C3*2
% 3Y    se Index_t3 < B2*Index_0 => Index_t3/Index_0

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
    Info_Cedole(:,1) = x2mdate(Info_Cedole(:,1),0);
    
% Variabili ad uso interno
    C1 = Info(1,1);
    C2 = Info(2,1);
    C3 = Info(3,1);
    B1 = Info(4,1);
    B2 = Info(5,1);
    B3 = Info(6,1);
    Num_Cedole_Fix = Info(7,1);
    TIR_BTP = Info(8:end,1);
    Px_t0 = Indice(end,2:end);
    
    Date_Bond = Info_Cedole(Info_Cedole(:,2) == 1,1);
    Date_Opz = Info_Cedole(Info_Cedole(:,3) == 2,1);
    
% Calcolo delle Date_Cedole
    if Frequenza == 'G'
        Passi = wrkdydif(Date_Opz(1:end-1,1),Date_Opz(2:end,1))-1;
        Passi(Passi(:,1) < 1) = 1;
        Passi = cumsum(Passi);        
     elseif Frequenza == 'M'
    %    Passi = months(Date_Opz(1:end-1,1),Date_Opz(2:end,1));
        Passi = (12:12:12*5)';
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
    Px_Ratio = Indice(2:end,2:end)./Indice(1:end-1,2:end);
    I_nd = cumprod(Px_Ratio);
    if Drift(1,1) ~= 0
        Drift = repmat(Drift,r_Indice-1,1);
%        Px_Ratio = Px_Ratio./repmat(mean(Px_Ratio),r_Indice-1,1).*(Drift/12+1);
        if Frequenza == 'G'
            Px_Ratio = Px_Ratio./repmat(mean(Px_Ratio),r_Indice-1,1).*(1+Drift).^(1/250);
            Px_Ratio = Px_Ratio-repmat(mean(Px_Ratio-1),r_Indice-1,1)+Drift/250;
            I_d = cumprod(Px_Ratio);
            plot([[1,1];[I_nd,I_d]]);
            grid on;
            grid off;
            %Px_Ratio = Px_Ratio./repmat(mean(Px_Ratio),r_Indice-1,1).*(1+Drift/250);
            
        elseif Frequenza == 'M'
            Px_Ratio = Px_Ratio./repmat(mean(Px_Ratio),r_Indice-1,1).*(1+Drift).^(1/12);
        end       
    end
    c = 1;
    Px_Scenari = repmat(Px_t0,Scenari,1);
    if Modo == 'B'      % bootstrap semplice
        for i = 1:Passi_Tot
            Px_Scenari_t = Bootstrap(Px_Ratio,Scenari);
            Px_Scenari = Px_Scenari.*Px_Scenari_t;
            if i == Passi(c)
                Px(:,:,c) = Px_Scenari;                
                c = c+1;     
            end    
        end
    elseif Modo == 'H'       % ipercubi latini
        Hyper = zeros(Scenari,Passi_Tot);
        Hyper(:,:) = LHSSample(Scenari,Passi_Tot);
        [r_sample] = size(Px_Ratio,1);
        Hyper = 1+floor(Hyper*r_sample);
        for i = 1:Passi_Tot
            Px_Scenari_t = Px_Ratio(Hyper(:,i),:);
            Px_Scenari = Px_Scenari.*Px_Scenari_t;
            if i == Passi(c)
                Px(:,:,c) = Px_Scenari;
                c = c+1;     
            end
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
            Cedole(:,1+p) = C1*(VRF_T(:,p) >= B3*VRI); 
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
    verifica = VRF_T./VRI;
    for j = 1:Num_Cedole-2-Num_Cedole_Fix
        check(VRF_T(:,j+Num_Cedole_Fix) < B3*VRI,j) = 1;
        check1(VRF_T(:,j+Num_Cedole_Fix) < B1*VRI,j) = 1;        
        Cedole(:,j+1+Num_Cedole_Fix) = C1.*(VRF_T(:,j+Num_Cedole_Fix) >= B3*VRI);
        Cedole(:,j+1+Num_Cedole_Fix) = Cedole(:,j+1+Num_Cedole_Fix)+ (1-check1(:,j));
        if j > 1
            Cedole(sum(Cedole(:,Num_Cedole_Fix+1:j+Num_Cedole_Fix),2) > 1,j+1+Num_Cedole_Fix) = 0;
        end
        
    % Calcolo del TIR
        TIR_C = [Cedole(:,1),Cedole(:,2:j+1+Num_Cedole_Fix)];
        TIR_C(TIR_C(:,1+j+Num_Cedole_Fix) < 1+C1,:) = [];
        %TIR_C(sum(check(:,1:j),2) ~= j-1,:) = [];
        %TIR_C(TIR_C(:,end) ~= 1+C1,:) = [];
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
            k = max(ceil((size(TIR_C,2)-1)/6)-1,1);
            Prob_BTP(1,j) = sum(TIR1 > TIR_BTP(k,1),1)./size(TIR1,1);
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
        Fix_C1 = repmat(1+C1,Scenari,1);
        Cedole(check_VRI == 1,j+2+Num_Cedole_Fix) = Fix_C1(check_VRI == 1);
        Cedole(check_VRF == 1,j+2+Num_Cedole_Fix) = VRF_T(check_VRF == 1,end)./VRI(check_VRF == 1);
                
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
                TIR1(i,1) = xirr(TIR_C(i,:)',[Date_Bond(1,:);Date_Bond(2:end,:)],-0.24,1000,3).*0.74;
            end
            TIR{j+1,1} = TIR1;
        % Statistica
            TIR_Perc(:,j+1) = prctile(TIR1,[1 5 25 50 75 95 99],1)';
            TIR_Max(1,j+1) = TIR_Perc(end,j+1);
            TIR_Media(1,j+1) = mean(TIR1,1);
            TIR_Min(1,j+1) = TIR_Perc(1,j+1);
            k = max(ceil((size(TIR_C,2)-1)/6)-1,1);
            Prob_Y(1,j+1) = size(TIR1,1)./Scenari;            
            Prob_BTP(1,j+1) = sum(TIR1 > TIR_BTP(k,1),1)./size(TIR1,1);
        % Statistica su tutta la distribuzione a scadenza
            TIR2 = [TIR{j+1,1}];
            TIR_Mediana(1,j+1) = prctile(TIR2,50,1);            
            Prob_0(1,j+1) = sum(TIR2 > TIR_BTP(k,1),1)./size(TIR2,1);
            Prob_0(1,j+1) = sum(TIR2 < -1.53222779501094e-10,1)./size(TIR2,1);                
        end        
    end
    Stat = [TIR_Max; TIR_Media; TIR_Min; Prob_BTP; Prob_Y; Prob_0; TIR_Mediana];    
%    hist({TIR(j+1,1)})

toc