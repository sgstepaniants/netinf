classdef progress <handle
    properties
        N;
        c;
    end
    methods 
        function obj=progress(N)
            obj.N = N;
            obj.c = 0;
        end
        function count(obj)
           obj.c = obj.c + 1;
           fprintf('%f%%\n',100 * obj.c / obj.N);
        end
    end
end
