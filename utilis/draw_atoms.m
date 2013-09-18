function A_matrix = draw_atoms(M, sz, per_row)

if iscell(M)
    Mn = M{end};
    for i=numel(M)-1:-1:1
        Mn = kron(Mn,M{i});
    end
    M=Mn;
end


if size(M,1) == prod(sz)
    M=M';
end


 NOP = create_tvmat(sz(1));
 NOP(sum(abs(NOP),2)==1,:)=[];
 NOP = bsxfun(@times,NOP,1./sqrt(sum(NOP.^2,2))); 
 [~,bb] = sort(sum(abs(NOP*M')));
 M = M(bb,:);

cols = ceil(size(M,1)/per_row);
RAND = 1;
if numel(sz)==3
    c = 3;
else
    c = 1;
end
A_matrix = zeros(cols*sz(1)+(cols+1)*RAND,per_row*sz(2)+(per_row+1)*RAND,c)+255;
% screen_size = get(0, 'ScreenSize');
% f1=figure(123456);
% set(f1, 'Position', [0 0 screen_size(3)/2 screen_size(4)/2 ] );
i = 1;
for y=1:cols
        if i > size(M,1)
            break;
        end
    for x=1:per_row
        crnt = reshape(M(i,:),sz);
        
        ys = (y-1)*(sz(1)+RAND)+1+RAND;
        xs = (x-1)*(sz(2)+RAND)+1+RAND;
        value = (crnt/(2*max(max(abs(crnt(:)))))+0.5)*255;

        A_matrix(ys:ys+sz(1)-1,xs:xs+sz(2)-1,:) = value;
        i=i+1;
        if i > size(M,1)
            break;
        end
    end
end


figure(123456)
imshow(uint8(A_matrix))
drawnow