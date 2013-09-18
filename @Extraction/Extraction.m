classdef Extraction < BaseOperator
    properties
        M
        sz
    end
    methods
        function obj = Extraction(M,sz)
            obj.M = unique(M);
            obj.sz = sz;
        end
        function res = mtimes(obj,b)
            if obj.adjoint
                res = zeros(obj.sz);
                res(obj.M) = b;       
            else
                res = b(obj.M);
            end
        end
    end
end


