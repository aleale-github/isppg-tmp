function [varargout] = DateMatWrkDays(varargin)

    nVarargs = length(varargin);
    varargout = cell(1,nargout);
% Data Conversion from Excel to Matlab and Next business days
    for k = 1:nVarargs
        d = varargin{k};
        if sum(d) ~= 0         
            Date = x2mdate(d(:,1),0);
            check1 = busdate(busdate(Date,1),-1);
            check1(check1 < Date,1) = busdate(check1(check1 < Date,1),+1);
            d(:,1) = check1; 
        end
        varargout{k} = d;
    end
