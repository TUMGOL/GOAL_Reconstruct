%(c) Simon Hawe, Lehrstuhl fuer Datenverarbeitung Technische Universitaet
%Muenchen, 2012. Contact: simon.hawe@tum.de
%Images    => Cell of images used to extract patches or path to folder
%dim       => dimension of single patch
%n_patches => Total Number of patches
function [S] = load_faces(Images, Patch_size)
if nargin == 3
    method = 0;
end

overlap = 0;

S = [];

    Image_list = dir(Images{1});
    Image_names = {};
    k = 0;
    S = zeros(prod(Patch_size),numel(Image_list));
    for i=1:numel(Image_list)
        if ~Image_list(i).isdir && ~strcmp(Image_list(i).name,'Thumbs.db') ...
                && ~strcmp(Image_list(i).name, Images{2})
            k = k + 1;
            Im = imread([Images{1},Image_list(i).name]);
            if ndims(Im) == 3 && numel(Patch_size) ~= 3
                Im = double(rgb2gray(Im));
            end
            %Im = imresize(Im,Patch_size,'bicubic');
            S(:,k)= Im(:);
        end
    end
    S = S(:,1:k);




