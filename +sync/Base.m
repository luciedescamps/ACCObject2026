classdef (Abstract) Base
    % BASE abstract superclass for all sync classes. Just provides
    % versioning functionality.
    
    properties (SetAccess = protected, Hidden)
        version = sync.version();
        datenum = now();
        upgraded = false;
    end 
    
    properties (Dependent)
        date
    end
    
    methods
        function val = get.date(self)
            if self.version >= 1
                val = datestr(self.datenum);
            else
                val = "<unknown>";
            end
        end
    end
end