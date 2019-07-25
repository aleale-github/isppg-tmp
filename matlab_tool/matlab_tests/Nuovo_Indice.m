function [Indice_New,Info] = Nuovo_Indice(Indice_True,Indice_False,Position,Option)

% Funzione che per qualsiasi Basket di dati, riallinea i dati del basket
% ed eventualmente allunga le serie con nuove serie (Indice_False). Position
% indica quali serie di Indice_True vanno allungate.
% Eventuali buchi finali vengono interpolati linearmente

% NB1: Nuovo_Indice(Indice_True): il codice gira senza allungare le serie ma
% allineando le originali... default: Option = 1
% Nuovo_Indice(Indice_True,[],[],[]): come sopra
% Nuovo_Indice(Indice_True,Indice_False,[],[]): allunga le/la serie
% originale (Indice_False può anche avere meno indici di quelli della serie
% originale).... default: Position = [1;2;3;...], Option = 1


% INPUT:
%  Indice_True = (M,N), dove M sono i dati storici degli Indici "Veri". La 1°, 3°, 5°, ecc... colonna le date
%                della serie e la 2°, 4°, 6°, ecc... i valori di mercato.
% Indice_False = (P,Q), dove M sono i dati storici delle "Proxy". La 1°, 3°, 5°, ecc... colonna le date
%                della serie e la 2°, 4°, 6°, ecc... i valori di mercato.
%     Position = (P,1), vettore che indicata quale Indice_True associare
%                all'Indice_False (es. [2;1] significa che l'Indice_True 1
%                è allungato con l'Indice_False 2 e l'Indice_True 2 con
%                Indice_False 1
%       Option = (1,1) se 2 allinea sulla serie originale con più date 
%                (possono capitare degli NaN) se 1 allinea sulla serie
%                originale con meno date. Se 3 allinea sulla serie con +
%                date ed estrapola linearmente per le serie più corte.
% OUTPUT:
%   Indice_New = (D,N/2), matrice che riporta nella prima colonna le date e
%                nelle successive i valori degli indici allineati
%         Info = (N/2,2), (N/2,1) indica quanti zeri sono stati
%               interpolati, (N/2,2) indica quanti NaN (da estrapolazione sono
%               stati cancellati

% Separazione Tassi di Cambio
    if (nargin > 1) && (isempty(Indice_False) == 0)
        Num_Cambi = size(find(Position(:,end) > 0),1);
        if Num_Cambi > 0
            Cambi = Indice_False(:,end-Num_Cambi-1:end);
            c = 0;
            C = cell(Num_Cambi,1);
            for j = 1:Num_Cambi
                Cambio = Cambi(:,1+c:2+c);
                Cambio(Cambio(:,1) == 0,:) = [];
                C{j,1} = Cambio;
                c = 2+c;
            end
            Indice_False(:,end-Num_Cambi-1:end) = [];        
        end
    end
  


% Variabili ad uso interno    
    Num_IT = size(Indice_True,2)/2;
    if (nargin == 1) || (isempty(Indice_False) == 1)
        Num_IF = 1;
    elseif nargin > 1
        Num_IF = size(Indice_False,2)/2;
    end
    A = cell(Num_IT,2);
    Num_Dati = zeros(Num_IT,Num_IF);
    Date_Start = zeros(Num_IT,Num_IF);
    
% Creazione Cell Array per Indici True ed eliminazione degli zeri   
    c = 0;
    for i = 1:Num_IT
        Indice = Indice_True(:,1+c:2+c);
        Indice(Indice(:,1) == 0,:) = [];
        A{i,1} = Indice;
        Num_Dati(i,1) = size(Indice,1);
        Date_Start(i,1) = Indice(1,1);
        c = 2+c;
    end
    B = A(:,1);    
    if (nargin > 1) && (isempty(Indice_False) == 0)
    % Creazione Cell Array per Indici False ed eliminazione degli zeri   
        c = 1;
        for i = 0:Num_IF-1
            Indice = Indice_False(:,i+c:i+1+c);
            Indice(Indice(:,1) == 0,:) = [];
            A{i+1,2} = Indice;
            Date_Start(i+1,2) = Indice(1,1);
            c = c+1;
        end    

    % Allungamento delle serie storiche dove necessario
        if (nargin == 2) || (isempty(Position) == 1)
            Position = (1:1:Num_IF)';
        end
        p = 1;
        for i = 1:Num_IF
            Indice_T = A{Position(i,1),1};
            Indice_F = A{Position(i,2),2};
            Data_T = Date_Start(Position(i,1),1);
            Indice_F(Indice_F(:,1) > Data_T,:) = [];

            Num_Date_F = size(Indice_F,1);
            Indice_N = zeros(Num_Date_F-1,2);
            Indice_N(1:Num_Date_F-1) = Indice_F(1:end-1,1);
        % Distingue se serie di tasso o Indice
            if Position(i,4) > 0
                Valuta = C{p,1};
                Valuta_A = interp1(Valuta(:,1),Valuta(:,2),Indice_F(:,1));
                Indice_F(:,2) = Indice_F(:,2)./Valuta_A;
                p = p+1;
            end
            if Position(i,3) == 0
                Price_Ratio = Indice_F(2:end,2)./Indice_F(1:end-1,2);
            elseif Position(i,3) == 1
                Price_Ratio = 1./(1+(Indice_F(2:end,2)-Indice_F(1:end-1,2))./100);
            end
            for d = 0:Num_Date_F-2
                if d == 0
                    Indice_N(end-d,2) = Indice_T(1,2)./Price_Ratio(end-d);
                else
                    Indice_N(end-d,2) = Indice_N(end-d+1,2)./Price_Ratio(end-d);
                end
            end
            Indice_N = [Indice_N; Indice_T];
            Num_Dati(Position(i,1),1) = size(Indice_N,1);        
            B{Position(i,1),1} = Indice_N;
        end
    end
    
% Allineamento Indici New
    if (nargin <= 3) || (isempty(Option) == 1)
        Option(1,1) = 1;
    end
    if Option(1,1) == 1
        [r,c] = find(Num_Dati == min(Num_Dati(:,1)));
        if size(r,1) > 1
            r = r(1,1);
            c = c(1,1);
        end
    elseif Option(1,1) >= 2
        [r,c] = find(Num_Dati == max(Num_Dati(:,1)));
        if size(r,1) > 1
            r = r(1,1);
            c = c(1,1);
        end
    end
    Date_Min = Num_Dati(r,c);
    Indice = B{r,c};
    Allineo_New = zeros(Date_Min,Num_IT+1);    
    Allineo_New(:,1) = Indice(:,1);
    for i = 1:Num_IT
        d = 1;
        while d <= Num_Dati(i,1)
            Indice = B{i,1};
            Allineo_New(Allineo_New(:,1) == Indice(d,1),i+1) = Indice(d,2);
            d = d+1;
        end
    end

% Interpolazione per eventuali Buchi nell'allineamento
    Indice_New = Allineo_New(:,1);
    Info = zeros(Num_IT,2);    
    for j = 1:Num_IT
        Allineo = [Allineo_New(:,1), Allineo_New(:,j+1)];
        Date_Zero = Allineo(Allineo(:,2) == 0,1);
        if isempty(Date_Zero) == 0
            Allineo(Allineo(:,2) == 0,:) = [];
            Indice_Zero = interp1(Allineo(:,1),Allineo(:,2),Date_Zero);
            Allineo = [Allineo;[Date_Zero,Indice_Zero]];
            Allineo = sortrows(Allineo,1);
        end
            Indice_New = [Indice_New, Allineo(:,2)];
    % Info
        Info(j,1) = size(Date_Zero,1); 
    end
    
% Cancellazione eventuali NaN ottenuti causa interpolazione senza
% estrapolazione (output di interp1)
    for j = 1:Num_IT
        if Option == 1
            Num_NaN = sum(isnan(Indice_New(:,j+1)));        
            Indice_New(isnan(Indice_New(:,j+1)) == 1,:) = [];
    % Info
        Info(j,2) = Num_NaN;            
        elseif Option == 2
            Num_NaN = sum(isnan(Indice_New(:,j+1)));        
            Indice_New(isnan(Indice_New(:,j+1)) == 1,:) = [];
        elseif Option == 3
            Num_NaN = sum(isnan(Indice_New(:,j+1)));            
            r = find(isnan(Indice_New(:,j+1)));
            if isempty(r) == 0
                Indice_New(r,j+1) = Indice_New(r(end,1)+1,j+1);
            end
            Info(j,2) = Num_NaN;            
        elseif Option == 4
            Num_NaN = sum(isnan(Indice_New(:,j+1)));            
            r = find(isnan(Indice_New(:,j+1)));
            if isempty(r) == 0
                Indice_New(r,:) = [];
            end
            Info(j,2) = Num_NaN;            
        end
    end
    
    
    
    
            
        
    
    
    
    