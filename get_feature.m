function feature = get_feature(S,patch_size)
       
[M, N] = size(S);
Mp = ceil(M/patch_size);
Np = ceil(N/patch_size);


feature = zeros(Mp*Np, 2);
i = 1;
for m = 1:Mp
    for n = 1:Np
        idx1 = [1,patch_size]+(m-1)*patch_size;
        idx2 = [1,patch_size]+(n-1)*patch_size;
        patch = S(idx1(1):min(idx1(2),M), idx2(1):min(idx2(2),N));
        
         
        
        [val_min,~] = min(patch(:));
        [val_max,~] = max(patch(:));
        

         
        feature(i,1)=val_min;
        feature(i,2)=val_max;

        i=i+1;
        
    end

end

%% normalization
% for j=1:size(feature,1)
%     feature(j,1)=feature(j,1)/sum(feature(:,1));
%     feature(j,2)=feature(j,2)/sum(feature(:,2));
% 
% end

end

