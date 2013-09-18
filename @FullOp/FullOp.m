% classdef FullOp < BaseOperator
%     properties
%         kernel
%         dim
%         %         large_scale_mode = 0;
%         h_step = 1;
%         v_step = 1;
%         sz
%         sz2
%         half_kernel
%     end
%     methods
%         function obj = FullOp(kernel,patch_size,h,v,sz)
%             obj.dim = size(kernel,1);
%             obj.kernel=cell(obj.dim,1);
%             for i=1:obj.dim
%                 obj.kernel{i} = reshape(kernel(i,:),patch_size,patch_size);
%             end
%             obj.sz2 = sz;
%             obj.sz  = sz-size(obj.kernel{1})+1;
%             obj.sz  = sz;
%             
%             obj.half_kernel = floor(patch_size/2);
%             
%             obj.v_step = unique([1:v:obj.sz(1),obj.sz(1)]);
%             obj.h_step = unique([1:h:obj.sz(2),obj.sz(2)]);
%             %             obj.large_scale_mode = ls_mode;
%             %obj.sz = ceil(size(obj.kernel{1})./2)-1;
%         end
%         function res = mtimes(obj,b)
%             if obj.adjoint
%                 %                 if obj.large_scale_mode
%                 %                     tmp = b.Data;
%                 %                     b.Data = [];
%                 %                     parfor i=1:obj.dim
%                 %                         tmp{i} = conv2(tmp{i},obj.kernel{i},'full');
%                 %                     end
%                 %                     res = zeros(size(tmp{1}));
%                 %                     parfor i = 1:obj.dim
%                 %                         res = res +  tmp{i};
%                 %                     end
%                 %                     b.Data = tmp;
%                 %                     tmp = [];
%                 %                 else
%                 b = mat2cell(b,ones(obj.dim,1)*size(b,1)/obj.dim);
%                 res = 0;
%                 for i = 1:obj.dim
%                     x    = zeros(obj.sz);
%                     x(obj.v_step,obj.h_step) = b{i};
%                     %b{i} = conv2(x,obj.kernel{i},'full');
%                     b{i} = conv2(x,obj.kernel{i},'same');
%                     %b{i} = imfilter(x,obj.kernel{i},'conv','same','replicate');
%                 end
%                 parfor i = 1:obj.dim
%                     res = res +  b{i};
%                 end
%                 %                 for outer = 1:16:obj.dim
%                 %                     crnt_sz = min(16,obj.dim-outer);
%                 %                     q = cell(crnt_sz,1);
%                 %                     for i=1:crnt_sz
%                 %                         x    = zeros(obj.sz);
%                 %                         x(obj.v_step,obj.h_step) = b{i};
%                 %                         q{i} = conv2(x,obj.kernel{i+outer-1},'full');
%                 %                     end
%                 %                     b(1:crnt_sz) = [];
%                 %                     for i = 1:crnt_sz
%                 %                         res = res +  q{i};
%                 %                     end
%                 %
%                 %                 end
%                 %                end
%             else
%                 %                 if obj.large_scale_mode
%                 %                     res = LargeScale(numel(obj.kernel));
%                 %                     tmp = res.Data;
%                 %                     res.Data = [];
%                 %                     parfor i=1:obj.dim
%                 %                         tmp{i} = filter2(obj.kernel{i},b,'valid');
%                 %                     end
%                 %                     res.Data = tmp;
%                 %                 else
%                 res = cell(numel(obj.kernel),1);
%                 b=b([ones(1,obj.half_kernel),1:end,end*ones(1,obj.half_kernel)],:);
%                 b=b(:,[ones(1,obj.half_kernel),1:end,end*ones(1,obj.half_kernel)]);
%                 parfor i=1:obj.dim
%                     x = filter2(obj.kernel{i},b,'valid');
%                     %x = filter2(obj.kernel{i},b,'same');
%                     %x = imfilter(b,obj.kernel{i},'corr','same','replicate');
%                     % x = x(obj.v_step,obj.h_step);
%                     res{i} = x(obj.v_step,obj.h_step);
%                 end
%                 res = cell2mat(res);
%                 %                 end
%             end
%         end
%     end
%     
%     
% end


classdef FullOp < BaseOperator
    properties
        kernel
        dim
        %         large_scale_mode = 0;
        h_step = 1;
        v_step = 1;
        full
        sz
        sz2
        half_kernel
    end
    methods
        function obj = FullOp(kernel,patch_size,h,v,sz)
            obj.dim = size(kernel,1);
            obj.kernel=cell(obj.dim,1);
            for i=1:obj.dim
                obj.kernel{i} = reshape(kernel(i,:),patch_size,patch_size);
            end
            obj.sz2 = sz;
            obj.sz  = sz-size(obj.kernel{1})+1;
            obj.sz  = sz;
            obj.full = 1;
            obj.half_kernel = floor(patch_size/2);
            
            if v==1 && h == 1
                obj.v_step = 1;
                obj.h_step = 1;
            else
                obj.v_step = unique([1:v:obj.sz(1),obj.sz(1)]);
                obj.sz2(1) = numel(obj.v_step);
                obj.h_step = unique([1:h:obj.sz(2),obj.sz(2)]);
                obj.sz2(2) = numel(obj.h_step);
                obj.full = 0;
            end
            
            
        end
        function res = mtimes(obj,b)
            if obj.adjoint
                res = 0;
                for i = 1:obj.dim
                    if obj.full
                        res  = res + conv2(b(:,:,i),obj.kernel{i},'same');
                     else
                         x    = zeros(obj.sz);
                         x(obj.v_step,obj.h_step) = b(:,:,i);
                         res  = res + conv2(x,obj.kernel{i},'same');
                     end
                end
            else
                res = zeros(obj.sz2(1),obj.sz2(2),numel(obj.kernel));
                sub=1;%double(~mod(obj.half_kernel,2));
                b=b([ones(1,obj.half_kernel),1:end,end*ones(1,obj.half_kernel-sub)],:);
                b=b(:,[ones(1,obj.half_kernel),1:end,end*ones(1,obj.half_kernel-sub)]);
                for i=1:obj.dim
                    x = filter2(obj.kernel{i},b,'valid');
                    if obj.full
                        res(:,:,i) = x;
                    else
                        res(:,:,i) = x(obj.v_step,obj.h_step);
                    end
                end
            end
        end
    end
    
    
end