function [DataSorted,varargout] = DataSorting(varargin)

    nVarargs = length(varargin);
    varargout = cell(1,nargout-1);
    d1 = [];
% Data Conversion from Excel to Matlab and Next business days
    for k = 1:nVarargs
        d = varargin{:,k};
        d2 = sort([d1;d]);
    % Check for Double
        check0 = diff(d2);
        check1 = [1;check0];
        d2(check1 == 0,:) = [];
        DataSorted = d2;
        d1 = d2;
    end

% Check for position
    for k = 1:nVarargs
        if isempty(varargin{:,k}) == 1
            varargout{:,k} = [];
        else
            Date = varargin{:,k};
            T = zeros(size(Date,1),1);
            for t = 1:size(Date,1)
                T(t) = find(Date(t) == d2);
            end
            varargout{:,k} = T;
        end
    end