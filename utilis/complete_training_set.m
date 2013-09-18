%(c) Simon Hawe, Lehrstuhl fuer Datenverarbeitung Technische Universitaet
%Muenchen, 2012. Contact: simon.hawe@tum.de
%Images    => Cell of images used to extract patches or path to folder
%dim       => dimension of single patch
%n_patches => Total Number of patches
function [S] = complete_training_set(Images, n_patches, Patch_size)
S = [];
if  ischar(Images{1})
    Image_list = dir(Images{1});
    Image_names = {};
    k = 0;
    for i=1:numel(Image_list)
        if ~Image_list(i).isdir && ~strcmp(Image_list(i).name,'Thumbs.db') ...
                && ~strcmp(Image_list(i).name, Images{2})
            k = k + 1;
            Image_names{k} = [Images{1},Image_list(i).name];
        end
    end
    n_patches = ceil(n_patches/k);
    for i = 1:k
        Im = imread(Image_names{i});
        if ndims(Im) == 3
            Im = rgb2gray(Im);
        end
        Im = double(Im);
        [Sc] = extract_training_set(Im, n_patches, Patch_size);
        S  = [S,Sc];
        %S  = [S,extract_training_set(Im, n_patches, Patch_size, method)];
    end
else
    n_patches = ceil(n_patches/numel(Images));
    for i = 1:numel(Images)
        S = [S,extract_training_set(Images{i}, n_patches, Patch_size)];
    end
end





% %(c) Simon Hawe, Lehrstuhl fuer Datenverarbeitung Technische Universitaet
% %Muenchen, 2012. Contact: simon.hawe@tum.de
% %Images    => Cell of images used to extract patches or path to folder
% %dim       => dimension of single patch
% %n_patches => Total Number of patches
% function [S,IdxSet] = complete_training_set(Images, n_patches, Patch_size, IdxSet, method)
% if nargin == 3
%     method = 0;
% end
% 
% overlap = 0;
% 
% S = [];
% if  ischar(Images{1})
%     Image_list = dir(Images{1});
%     Image_names = {};
%     k = 0;
%     for i=1:numel(Image_list)
%         if ~Image_list(i).isdir && ~strcmp(Image_list(i).name,'Thumbs.db') ...
%                 && ~strcmp(Image_list(i).name, Images{2})
%             k = k + 1;
%             Image_names{k} = [Images{1},Image_list(i).name];
%         end
%     end
%     n_patches = ceil(n_patches/k);
%     for i = 1:k
%         if strcmp(Image_names{i}(end-2:end),'mat')
%             load(Image_names{i});
%             Im = DI;
%         elseif strcmp(Image_names{i}(end-2:end),'mnc')
%             Im = loadminc(Image_names{i});
%             overlap = 1;
%         else
%             Im = imread(Image_names{i});
%             if ndims(Im) == 3 && numel(Patch_size) ~= 3
%                 Im = rgb2gray(Im);
%             end
%             Im = double(Im);
%         end
%         if length(IdxSet) < i
%             IdxSet{i}=[];
%         end
%         
%         [Sc,IdxSet{i}] = extract_training_set(Im, n_patches, Patch_size, IdxSet{i}, method, overlap);
%         
%         S  = [S,Sc];
%         %S  = [S,extract_training_set(Im, n_patches, Patch_size, method)];
%     end
% else
%     n_patches = ceil(n_patches/numel(Images));
%     for i = 1:numel(Images)
%         S = [S,extract_training_set(Images{i}, n_patches, Patch_size, method)];
%     end
% end
% 
% 
% 
