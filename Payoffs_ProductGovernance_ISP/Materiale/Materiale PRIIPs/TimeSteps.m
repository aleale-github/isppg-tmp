function [Steps,NumSteps,T] = TimeSteps(StartDate,EndDate,TimeFreq,ProductType)



% Time Horizon
    if strcmp(ProductType,'TC') || strcmp(ProductType,'TCA') || ... 
        strcmp(ProductType,'TCS') || strcmp(ProductType,'TCSA') || ...
        strcmp(ProductType,'TPPF') || strcmp(ProductType,'TPPFA') || ...
        strcmp(ProductType,'TPPU') || strcmp(ProductType,'TPPUA') || ... 
        strcmp(ProductType,'TPPVI') || ...
        strcmp(ProductType,'CPPVI') || strcmp(ProductType,'CPPVE') || ...
        strcmp(ProductType,'CCTI') || strcmp(ProductType,'CCTE') || ...
        strcmp(ProductType,'CPFKII') || strcmp(ProductType,'CPFKIE') || ...
        strcmp(ProductType,'CPCI') || strcmp(ProductType,'CPCE') || ...
        strcmp(ProductType,'CPCA') || strcmp(ProductType,'CPCV')
        T = round(yearfrac(StartDate,EndDate),1); 
    else
        T = round(yearfrac(StartDate,EndDate),12); 
    end

% Simulation Steps
    DateSteps = [StartDate;EndDate];
    a = 1; t = 1;
    Steps = zeros(size(DateSteps,1)-1,1);
    if strcmp(TimeFreq{1},'D')
        while StartDate <= DateSteps(end)
            StartDate = busdate(StartDate,1);
            if StartDate == DateSteps(t+1)
                Steps(t) = a;
                t = t+1;
                if t == size(DateSteps,1)
                    break
                end
            end
            a = a+1;
        end
    elseif strcmp(TimeFreq{1},'W')

    elseif strcmp(TimeFreq{1},'M')
        %calcolo delle date a cui fare riferimento per il numero di mesi:
        %lascia ferme le date a fine mese e porta le date intermedie alla
        %fine del mese precedente.
        
        DateSteps1 = datenum(dateshift(datetime(DateSteps,'ConvertFrom','datenum')+1,'end','month','previous')); 
        Steps = months(DateSteps1(1:end-1),DateSteps1(2:end));
        Steps = cumsum(Steps,1);        
    end
    NumSteps = size(Steps,1);
 
  