classdef Optimizationbase < handle & hgsetget
    %% Accessible from outside
    properties
        t                   = 1;    % Stepsize
        alpha               = 1e-1; % First linesearch parameter
        beta                = .9;   % Second linesearch parameter
        Ck                  = 0;
        Qk                  = 0;
        nuk                 = .0;
        maxLineSearchSteps  = 150;
        p_norm              = 1;    % ThatÄs currently just for testing
        delta_conv          = 1e-9; % Parameter for stopping iterations
        max_iterations      = 500;  % Maximum number of iterations
        verbose             = 1;    % If == 0 no output, otherwise every veboseth frame 1 output
        side_eval_function  = {};   % If this is a function handle, this function will be evaluate when verbose is true
        lambda              = 1;    % Lagrange Multiplier
        lambda_min                  % Minimum Lagrange Multiplier if continuation
        lambda_st           = 1;
        l1Smooth                    % l1-norm smoothing term
        p                   = 0;    % For positivity constraint, if == 0 no pos contraint
        sigma
        Norm_lambda
    end
    
    properties %(SetAccess = protected)
        %% Variable Required for the optimization
        k = 0;          % Current iteration
        f0              % Last value ob objective function
        
        %% Vectors used throughout the optimization
        x           % Variable to be optimized
        dx          % Stepdirection
        y           % Measured Data
        g           % Previous Gradient
        c_y         % Store A*x-y = c_y
        c_dy        % Store A*dx  = c_dy
        sp_x = 1;       % Store W*x   = sp_x
        sp_dx = 1;      % Store W*dx  = sp_dx
        range = [0,255]      % The range where the positivity constraint lies in
        
        %% Matrices, perhaps this must move into the methods section but I
        %  don't know it yet. But these are rather operators then
        %  propoerties
        AOP     = 1;  % Sparsifying transformation
        M           % This extracts the required points
        
        %% Things to control the stopping
        gx_old = zeros(11,1);
    end
    
    methods
        %% This function evalutes the total cost function to be minimized.
        %  As we are dealing with unconstrainted lagrange functions, this
        %  always has the form objective + lambda* regularizer
        function obj = Evalute(obj)
            fx = f_x(obj);
          %  fprintf('FIDELITY SHIT %e\n',fx)
            obj.f0 = obj.f0 - fx;
            if obj.f0 < 0
                obj.f0 = fx;
                return;
            end
            gx = g_x(obj);

            %fprintf('Value of Inf Norm: %f\n',obj.lambda_st);
            %obj.lambda = min(1,min(obj.lambda,obj.lambda_st));
            % Required for the stopping criterion
            obj.gx_old = [obj.gx_old(2:end);gx];
            obj.f0 = fx + obj.Norm_lambda*obj.lambda*gx;
        end
        
        %% linesearch updating the direction
        function [obj,worked] = linesearch(obj)
            if obj.k  == 1
                obj.Ck = 0;
                obj.Qk = 0;
            end
            
            obj.Ck = (obj.nuk*obj.Qk*obj.Ck+obj.f0)/(obj.nuk*obj.Qk+1);
            obj.Qk =  obj.nuk*obj.Qk+1;
            f_c = obj.f0;
            ls_iter     = 0;
            worked      = 1;
            val = obj.alpha * ((obj.g(:)'*obj.dx(:)));
            
            %Performing the linesearch
            while (ls_iter == 0) || ... % Quasi Do while loop
                    (obj.f0 > obj.Ck + obj.t * val) && ... % Check Wolfe Condition
                    (ls_iter < obj.maxLineSearchSteps) % Check Maximum number of Iterations
                
                
                % Update the step size
                obj.t  = obj.t*obj.beta;
                % HACK
                obj.f0 = f_c + obj.t * val;% + eps;
                
                % Evalute the Objective Function at the taken step
                Evalute(obj);
                
                % Increase the number of linesearch iterates as we only
                % allow a certain amount of steps
                ls_iter = ls_iter + 1;
            end
            
            if ls_iter >= obj.maxLineSearchSteps
                worked = 0;
                return;
            end
            if obj.verbose
                fprintf('ls_iter %d, t: %f\n',ls_iter, obj.t);
            end
            update_variables(obj);

            obj.t  =  obj.t/obj.beta^2;  
        end

        
        %% Update everything for the next iteration
        function [obj, stop] = next_iteration(obj, cg_beta, g_new)
            stop = 0;
            
            if obj.k >= numel(obj.gx_old)
                mean_gx = mean(obj.gx_old(1:end-1));
                delta_gx = abs(obj.gx_old(end)-mean_gx)/mean_gx;
                if delta_gx <= obj.delta_conv
                    stop = 1; return;
                end
            end
            
            % Update the iteration counter
            obj.k  = obj.k  + 1;
            
            if obj.k >= obj.max_iterations
                stop = 1; return;
            end
            
            if obj.verbose && ~mod(obj.k,obj.verbose)
                for func_id = 1:numel(obj.side_eval_function)
                    if isa(obj.side_eval_function{func_id},'function_handle')
                        obj.side_eval_function{func_id}(obj.x);
                    end
                end
            end
            
            obj.dx = -g_new + cg_beta*obj.dx;
            % Store the Gradient
            obj.g  = g_new;
            preobjective(obj);
        end
        
        %% This function must be overloaded if other variable are added
        %  that should be updated
        function obj = update_variables(obj)
            obj.x    = obj.x     + obj.t * obj.dx;
            obj.c_y  = obj.c_y   + obj.t * obj.c_dy;
            obj.sp_x = obj.sp_x  + obj.t * obj.sp_dx;
        end
        
        % Calculates the positivity penalty value for xp
        function val = Inrange(obj, xp)
            if obj.p == 0
                val = 0; return;
            end

            xp = (xp);
          
            sel1     = xp > obj.range(2);
            sel2     = xp < obj.range(1);
            
            xp(sel1) = (xp(sel1) - obj.range(2));
            xp(sel2) = (xp(sel2) - obj.range(1));
            
            xp = xp*1;
            
            val     = (sum(abs(xp((sel1 | sel2))).^obj.p));
        end
        
        % This is the gradient of the inrange function
        function grad = d_Inrange(obj, xp)
            if obj.p == 0
                grad = 0; return;
            end
            %grad = xp;
            grad = zeros(size(xp));
            
            %grad(grad>=obj.range(1) & grad<=obj.range(2))=0;
            
            sel1     = xp >= obj.range(2);
            sel2     = xp <= obj.range(1);
            
            v1 = xp(sel1) - obj.range(2);
            v2 = xp(sel2) - obj.range(1);
            
            mul = obj.p*(1)^obj.p;
            
            grad(sel1) = mul*sign(v1).*abs(v1).^(obj.p-1);
            grad(sel2) = mul*sign(v2).*abs(v2).^(obj.p-1);
        end
    end
    
    
    %% Abstract methods that must be reimplemented by inheriting classes
    methods (Abstract)
        val = f_x(obj)    % Objective   Function
        val = g_x(obj)    % Regularizer Function
        obj = preobjective(obj)
    end
    
end
