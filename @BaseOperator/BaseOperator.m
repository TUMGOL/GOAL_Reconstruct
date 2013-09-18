classdef BaseOperator 
    properties
        adjoint = 0;
    end
    methods        
        function res = eq(obj,b)
            res = false;
            if strcmp(class(obj),class(b)) && (obj.adjoint == b.adjoint)
                res = true;
            end
        end
        
        function obj = ctranspose(obj)
            obj.adjoint = ~obj.adjoint;
        end
        
        function res = times(obj,b)
            res = mtimes(obj,b);
        end

    end
    
    methods (Abstract)
        res = mtimes(obj,b)
    end
    
    
end