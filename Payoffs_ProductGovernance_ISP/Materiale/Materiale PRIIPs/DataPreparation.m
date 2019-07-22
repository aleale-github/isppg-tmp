function [Underlying,Div,DataType,Freq,UnderlyingName,IRCnew,LogRend,LogRendAsset] = DataPreparation(Asset,Dividends,DataType,Freq,UnderlyingName,IRC,PricingDate)

% Internal Variable
NumAsset = size(Asset,2)/2;
Div=zeros(NumAsset,1);
% Adjust Matricies
a = 2:2:NumAsset*2;
col = 1:2:(NumAsset*2)-1;

 for i = 1:NumAsset
        if strcmp(DataType{1,col(i)},'Rate')
            Asset(:,a(i)) = Asset(:,a(i))./100;
        end
 end
    
% From Excel to Matlab Dates
for i=1:NumAsset
    Asset(Asset(:,a(i)-1)>0,a(i)-1)=x2mdate(Asset(Asset(:,a(i)-1)>0,a(i)-1));
end

% Only dates coherent with Pricing date and series lenght

LogRendAsset=zeros(size(Asset,1),NumAsset*2);

MaxSeriesLenght=0;

for i=1:NumAsset
    LastUndDate=Asset(find(Asset(:,a(i)-1),1,'last'),a(i)-1);
    MinDataSim = datenum(year(LastUndDate)-5,month(LastUndDate),day(LastUndDate));
    Asset1(:,:)=Asset(Asset(:,a(i)-1)>=MinDataSim & Asset(:,a(i)-1)<=LastUndDate(1),a(i)-1:a(i));
    LogRendAsset(1:size(Asset1,1),(i*2)-1:i*2)=Asset1;
    MaxSeriesLenght=max(MaxSeriesLenght,size(Asset1,1));
    Asset1=[];
end

LogRendAsset(MaxSeriesLenght+1:end,:)=[];

% Last Dividends value

for NA = 1:NumAsset
    if size(Dividends(:,a(NA)-1),1) > 1
        CheckDataDiv = Dividends(Dividends(:,a(NA)-1) > 0,a(NA)-1);
        CheckDataDiv2=Dividends(Dividends(:,a(NA)-1) > 0,a(NA));
        Div(NA)= CheckDataDiv2(end,1)./100;
    else
        Div(NA)=0;
    end
end

% Common dates between series

if NumAsset>1
    ND = size(LogRendAsset(LogRendAsset(:,1) > 0,1),1);
    h=0;
    for i = 1:ND
        ch=1;
        for j=2:NumAsset
            if isempty(find(LogRendAsset(:,a(j)-1)==LogRendAsset(i,1))) == 0
                ch=ch+1;
            end
        end
        if ch == NumAsset
            h=h+1;
            CommDates(h,:)= LogRendAsset(i,1);
        end
    end
    CommDates(CommDates(:,:)==0,:)=[];
    
    % Different dates between series
    
    DiffDates=zeros(1,1);
    h=0;
    for i=1:NumAsset
        ND = size(LogRendAsset(LogRendAsset(:,a(i)-1) > 0,1),1);
        for j = 1:ND
            if isempty(find(CommDates(:,:)==LogRendAsset(j,a(i)-1))) && isempty(find(DiffDates(:,:)==LogRendAsset(j,a(i)-1)))
                h=h+1;
                DiffDates(h,:)= LogRendAsset(j,a(i)-1);
            end
        end
    end
    
    % Create DataType New
    DataType = DataType(:,col);

    % LogRend calculation
    DataType = DataType';   
    for i = 1:size(DataType,1)
        if strcmp(DataType{i,1},'Rate')
            if strcmp(DataType{i,end},'MinRate')    
                LogRendAsset(LogRendAsset(:,a(i)-1)>0,a(i)) = LogRendAsset(LogRendAsset(:,a(i)-1)>0,a(i))+DataType{i,end-1};
            end
        end
    end
 
    DataType = DataType';  
    
    for i=1:NumAsset
       
        h=0;
        ND = size(LogRendAsset(LogRendAsset(:,a(i)-1) > 0,1),1);
        j=1;
        while j<=ND
            if isempty(find(CommDates(:,1) == LogRendAsset(j,a(i)-1))) == 0
                if  j==1
                    j=j+1;
                elseif isempty(find(CommDates(:,1)==LogRendAsset(j,a(i)-1)))==0 && isempty(find(CommDates(:,1)==LogRendAsset(j-1,a(i)-1)))==0
                    if DiffDates(LogRendAsset(j,a(i)-1)>DiffDates & DiffDates>LogRendAsset(j-1,a(i)-1),:)
                        j=j+1;
                    else
                        h=h+1;
                        LogRend(h,i)=log(LogRendAsset(j,a(i))/LogRendAsset(j-1,a(i)));
                        j=j+1;
                    end
                else
                j=j+1;
                end
            else
                j=j+2;
            end
        end
    end
else

LogRendAsset=LogRendAsset(LogRendAsset(:,1)>0,:);
LogRend = log(LogRendAsset(2:end,2)./LogRendAsset(1:end-1,2));
end

% Per le Currency????

%     Currency1 = zeros(size(Asset,1),size(Asset,2));
%         if size(Currency,a(NA)-1) > 1
%             CheckDataCurr = Currency(Currency(:,a(NA)-1) > 0,a(NA)-1);
%             CheckDataCurr(1:end-1) = [];
%             CheckDataAsset = Asset(Asset(:,a(NA)-1) > 0,a(NA)-1);
%             NumDataAsset = size(CheckDataAsset,1);
%             CheckDataAsset(1:end-1) = [];
%             if CheckDataCurr < CheckDataAsset
%                 Currency1(:,a(NA)-1) = Asset(:,a(NA)-1);
%                 for ND = 1:NumDataAsset
%                     if sum(Currency(ND,a(NA)-1) == Asset(:,a(NA)-1),1) == 1
%                         Currency1(ND,a(NA)) = Currency(ND,a(NA));
%                     else
%                         Currency1(ND,a(NA)) = Currency1(ND-1,a(NA));
%                     end
%                 end
%             end
%         end


LogRend=[LogRend ones(size(LogRend,1),size(LogRend,2))];

% Create Underlying matrix
Underlying = Asset;
Today = Underlying(end,1);


% Create Freq New
Freq = Freq(:,col);
% Create UnderlyingName New
UnderlyingName = UnderlyingName(:,col);
% Create IRC variable
NumIRC = size(IRC,1);
IRCnew = zeros(size(IRC{1},1),NumIRC+1);
IRC1 = IRC{1};
IRCnew(:,1:2) = IRC1;
if NumIRC > 1
    for r = 1:NumIRC-1
        IRC1 = IRC{1+r};
        IRCnew(:,2+r) = IRC1(:,2);
    end
    
end
end