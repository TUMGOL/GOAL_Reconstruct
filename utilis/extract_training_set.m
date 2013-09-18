%% This function extracts n training samples of size sz from Signal S.
% The input S can be either a 1-D signal or 2-D signal (image). For a 1-D
% signal n training sequences having sz samples will be extracted from
% random positions. For a 2-D Signal n training patches of size
% [sz(1),sz(2)] or [sz,sz] will be extracted from n random positions.
% (c) Simon Hawe, Lehrstuhl fuer Datenverarbeitung Technische Universitaet
% Muenchen, 2012. Contact: simon.hawe@tum.de
function [X] = extract_training_set(S, n, sz)

if isscalar(sz) && min(size(S)) > 1
    sz = [sz,sz];
end

patches = im2col(S,sz,'sliding');

patches_ext = bsxfun(@minus,patches,mean(patches));
patches_ext = bsxfun(@times,patches_ext,1./sqrt(sum(patches_ext.^2)));
patches(:,isnan(1./sum(patches_ext)))=[];
patches_ext(:,isnan(1./sum(patches_ext)))=[];

TV_Mat  = create_tvmat(sz(1));
%regions = sum(abs(TV_Mat*patches_ext));
regions = (abs(TV_Mat*patches_ext));
regions = sum(bsxfun(@times,regions,median(regions)));
regions = regions < (1.3*mean(regions));
patches_ext(:,regions) = [];
patches(:,regions) = [];
%[~,p] = unique(patches_ext','rows');
p = 1:length(patches_ext);
sel = randperm(numel(p));
X = patches(:,p(sel(1:min(end,n))));

% sz_h = floor(sz./2);
% % Extraction window. awesomely coded :)
% Extractor = bsxfun(@plus,[-sz_h(2):(-sz_h(2)+sz(2)-1)]',(-sz_h(1):(-sz_h(1)+sz(1)-1))*size(S,1));
% Extractor = Extractor(:);
% 
% if ~isempty(op) && ~isscalar(op)
%     op  = FullOp(op, sz(1), 1,1, size(S));
%     res = sum(log(1+(op*double(S)).^2),3);
%     res(1:sz_h(1),:)=0;
%     res(end-sz_h(1)+1:end,:)=0;
%     res(:,1:sz_h(2))=0;
%     res(:,end-sz_h(2)+1:end)=0;
%     res(IdxSet) = 0;
%     [~,idx]=sort(res(:),'descend');
% else
%     pos = reshape(1:numel(S),size(S,1),size(S,2));
%     pos(IdxSet)= 0;
%     pos(1:sz_h(1),:)=[];
%     pos(end-sz_h(1)+0:end,:)=[];
%     pos(:,1:sz_h(2))=[];
%     pos(:,end-sz_h(2)+0:end)=[];
%     if ~isempty(op) && op == 1
%         idx = pos(1:sz(1):end,1:sz(2):end);
%         idx = idx(:);
%         n = numel(idx);
%     else
%         pos = pos(:);
%         pos(pos==0)=[];
%         rand('state',0);
%         sel = randperm(numel(pos));
%         idx = pos(sel);
%         n = min(numel(idx),n);
%     end
% end
% C_pos = bsxfun(@plus,idx(1:n)',Extractor);
% X = S(C_pos);
% IdxSet = [IdxSet,idx(1:n)'];

end


% if ~isempty(op)
%
%     sz_h = floor(sz./2);
%
%     % Extraction window. awesomely coded :)
%     Extractor = bsxfun(@plus,[-sz_h(2):(-sz_h(2)+sz(2)-1)]',(-sz_h(1):(-sz_h(1)+sz(1)-1))*size(S,1));
%     Extractor = Extractor(:);
%
%     op  = FullOp(op, sz(1), 1,1, size(S));
%     res = sum(abs(op*double(S)).^.2,3);
%     res(1:sz_h(1),:)=0;
%     res(end-sz_h(1)+1:end,:)=0;
%     res(:,1:sz_h(2))=0;
%     res(:,end-sz_h(2)+1:end)=0;
%
%     [~,idx]=sort(res(:),'descend');
%
%     C_pos = bsxfun(@plus,idx(1:n)',Extractor);
%     X = S(C_pos);
%
% else
%     [i1,i2] = reggrid(size(S)-sz+1, n);
%     X = sampgrid(S, sz, i1, i2);
% end


%X = X(:,randperm(n));
%X = bsxfun(@minus,X,mean(X));


% %% This function extracts n training samples of size sz from Signal S.
% % IdxSet contains indizes of patches that should not be extracted.
% % If the Analysis operator op is also provided, this function will extract
% % traning patches which are not cosparse through the given operator.
% function [X, IdxSet] = extract_training_set(S, n, sz, IdxSet, op, overlap)
%     SZ = size(S);
%     sz_h = floor(sz./2);
%     % Extraction window. awesomely coded :)
%     Extractor = bsxfun(@plus,[-sz_h(1):(-sz_h(1)+sz(1)-1)]',(-sz_h(2):(-sz_h(2)+sz(2)-1))*size(S,1));
%     if numel(sz) == 3
%         Extractor = KronSum([0:numel(S(:,:,1:sz(end)))/sz(end):numel(S(:,:,1:sz(end)))-1],Extractor);
%     end
%     Extractor = Extractor(:);
%     
%     
%     if ~isempty(overlap)
%         sz_h = sz;
%     end
%     
%     if ~isempty(op) && ~isscalar(op)
%         op  = FullOp(op, sz(1), 1,1, size(S));
%         res = sum(log(1+(op*double(S)).^2),3);
%         res(1:sz_h(1),:)=0;
%         res(end-sz_h(1)+1:end,:)=0;
%         res(:,1:sz_h(2))=0;
%         res(:,end-sz_h(2)+1:end)=0;
%         res(IdxSet) = 0;
%         [~,idx]=sort(res(:),'descend');
%     else
%         if sz(end) == SZ(end)
%             pos = reshape(1:numel(S(:,:,1)),size(S,1),size(S,2));
%         else
%             pos = reshape(1:numel(S),size(S));
%         end
%         pos(IdxSet)                 = 0; 
%         pos(1:sz_h(1),:,:)          = [];
%         pos(end-sz_h(1)+1:end,:,:)  = [];
%         pos(:,1:sz_h(2),:)          = [];
%         pos(:,end-sz_h(2)+1:end,:)  = [];
%         
%         if numel(sz)>2 && sz(end) ~= SZ(end)
%             mul = 1;
%             if ~isempty(overlap)
%                 mul = 2;
%             end
%             pos(:,:,end-sz(3)*mul+1:end)=[];
%         end
%         
%         if ~isempty(op) && op == 1
%             idx = pos(1:sz(1):end,1:sz(2):end);
%             idx = idx(:);
%             n = numel(idx);
%         else
%             pos = pos(:);
%             pos(pos==0)=[];
%             sel = randperm(numel(pos));
%             idx = pos(sel);
%         end
%     end
%     n = min(n,numel(idx));
%     
%     if ~isempty(overlap)
%         idx_2 = idx;
%         laps = ceil(n/prod(sz));
%         step = numel(Extractor);
%         for i=1:laps
%             ix = randi(numel(idx_2),1);
%             idx(1+(i-1)*step:i*step) = idx_2(ix)+Extractor;
%             idx_2 = setdiff(idx_2,idx(1+(i-1)*step:i*step));
%         end
%         n = laps*prod(sz);
%     end
%     
%     C_pos = bsxfun(@plus,idx(1:n)',Extractor);
%     X = S(C_pos);
%     IdxSet = [IdxSet,idx(1:n)'];
% 
% end
% 
% 
% % if ~isempty(op)
% %     
% %     sz_h = floor(sz./2);
% %     
% %     % Extraction window. awesomely coded :)
% %     Extractor = bsxfun(@plus,[-sz_h(2):(-sz_h(2)+sz(2)-1)]',(-sz_h(1):(-sz_h(1)+sz(1)-1))*size(S,1));
% %     Extractor = Extractor(:);
% %     
% %     op  = FullOp(op, sz(1), 1,1, size(S));
% %     res = sum(abs(op*double(S)).^.2,3);
% %     res(1:sz_h(1),:)=0;
% %     res(end-sz_h(1)+1:end,:)=0;
% %     res(:,1:sz_h(2))=0;
% %     res(:,end-sz_h(2)+1:end)=0;
% %     
% %     [~,idx]=sort(res(:),'descend');
% %     
% %     C_pos = bsxfun(@plus,idx(1:n)',Extractor);
% %     X = S(C_pos);
% %    
% % else
% %     [i1,i2] = reggrid(size(S)-sz+1, n);
% %     X = sampgrid(S, sz, i1, i2);
% % end
% 
% 
% %X = X(:,randperm(n));
% %X = bsxfun(@minus,X,mean(X));
