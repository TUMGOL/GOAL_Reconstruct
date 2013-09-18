% CSLasso Stands for Compressed Sensing based on the Noiselt transformation
% and the Wavelet transform as the sparsifying transformation
classdef SparseSolver < Optimizationbase
    properties
        p_sm = 2;   % p_sm <= 1 for l0 approximation by a smoothed funciton
        grad_g
    end
    methods
        %% Constructor
        function obj   = SparseSolver(y, x0, M_in, AOP_in, l1Smooth, p, range)
            
            if exist('AOP_in','var') && ~isempty(AOP_in)
                obj.AOP      = AOP_in;
            end
            
            if exist('M_in','var') && ~isempty(M_in)
                obj.M      = M_in;
            end
            
            % Creating Measurements
            obj.y      = y;
            
            % Inital Guess of the signal to be reconstructed
            obj.x      = x0;
            
            % Setting the Lagrange alpha_smltiplyer
            obj.lambda      = 0;
            
            % Setting the smoothing value for the derivative
            obj.l1Smooth    =  l1Smooth;
            
            % The value for the Positivity Constraint
            obj.p = p;
            
            % The range in which the positivity constraint lies
            obj.range = range;
        end
        
        %% Init function for CG method or even GD method
        function obj = init(obj)
            % First we set t = 0 to evalute this function at the current
            % poistion
            obj.gx_old  = obj.gx_old*0;
            t_s         = obj.t;
            obj.t       = 0;
            
            % These variable are computed only once, and consecutively
            % updated. In that way we can save some matrix multiplications
            obj.c_y    = obj.M*obj.x - obj.y;
            obj.sp_x   = (obj.AOP * obj.x);
            
            obj.Norm_lambda = numel(obj.y)/numel(obj.sp_x);
            
            obj.grad_g = 1 + obj.l1Smooth*obj.sp_x.^2;
            
            % Computing the gradient
            obj.c_dy    = zeros(size(obj.c_y));
            obj.sp_dx   = zeros(size(obj.sp_x));
            obj.dx      = zeros(size(obj.x));
            Evalute(obj);
            obj.g       = gradient(obj);
            %Setting the descent direction
            obj.dx      = -obj.g;
            
            preobjective(obj);
            Evalute(obj);
            % Set it back to the inital values
            obj.t   = t_s;%min(1,1/obj.beta^2*t_s);
            obj.k   = 0;
            
        end
        
        %% Gradient of the objective function
        function grad_f = df_x(obj)    % Objective   Function
            % We store this, because this values is always needed for
            % evaluating the objective function.  In this way it is
            % computed only once
            if obj.p_sm == 1
                grad_f          =  obj.M'*(sign(obj.c_y));
            else
                grad_f          =  obj.M'*(obj.c_y.^(obj.p_sm-1));
            end
        end
        
        %% Computing the entire gradient
        function grad = gradient(obj)
            grad_f = df_x(obj) + obj.d_Inrange(obj.x);
            obj.grad_g = obj.l1Smooth*(obj.AOP'*((obj.sp_x)./(obj.grad_g)));
            grad = (grad_f + obj.Norm_lambda*obj.lambda*obj.grad_g);
        end
        
        %% Preobjective function required for accelerating the algorithmn
        function obj = preobjective(obj)
            obj.c_dy  = obj.M*(obj.dx);
            obj.sp_dx = obj.AOP*obj.dx;
        end
        
        %% Objective Function
        function val = f_x(obj)
            val =  (obj.c_y + obj.t*obj.c_dy);
            if obj.p_sm == 1
                val =  sum(abs(val(:))) + max(1,obj.Norm_lambda*obj.lambda)*((obj.Inrange(obj.x+obj.t*obj.dx)));
            else
                val =  1/obj.p_sm*sum(abs(val(:)).^obj.p_sm) + max(1,obj.Norm_lambda*obj.lambda)*((obj.Inrange(obj.x+obj.t*obj.dx)));
            end
        end
        
        %% Regularizer Function
        function val = g_x(obj)
            obj.grad_g = (obj.sp_x + obj.sp_dx*obj.t);
            obj.grad_g = 1+obj.l1Smooth*obj.grad_g.^2;
            val        = 0.5*sum((log(obj.grad_g(:))));
        end
    end
end