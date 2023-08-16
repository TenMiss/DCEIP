function S = deblur_tv_emp(Im, kernel, lambda, Blambda,mu, opts)

S = Im;
alphamax = 1e5;
[M,N,D] = size(Im);
otfFh = psf2otf([1, -1], [M,N]);
otfFv = psf2otf([1; -1],[M,N]);
otfKER = psf2otf(kernel,[M,N]);

denKER  = abs(otfKER).^2;
denGrad = abs(otfFh).^2 + abs(otfFv ).^2;
Fk_FI = conj(otfKER).*fft2(Im);
alpha = 2.0*mu;

K=3;
kappa=2;

while alpha<alphamax

    for k=1:K
        
          feature = get_feature(S, opts.r);  
            
           %% minimal 
           [~,U_min,obj_min] = fcm(feature(:,1),2,[2,100,1e-5,0]);
           [~,label_min] = max(U_min); 
            
           %% maximal
           [~,U_max,obj_max] = fcm(feature(:,2),2,[2,100,1e-5,0]);
           [~,label_max] = max(U_max); 

    
            joint_max=max(U_min,U_max);
            [~,index] = max(joint_max);   
            
       [Z, Md, thresh_split] = find_optimal_pixels(S,opts.r,index);
   
        x=Z(Md>0);
        z = Z(Md>0 & thresh_split == 1);
        zB= Z(Md>0 & thresh_split == 2);
        
        if opts.s<opts.scales/2
            lambdat = min(max(lambda,mean(abs(z))),0.1);
            lambdatB = max(min(Blambda,mean(abs(zB))),0.9);
            Z(abs(Z)<lambdat & thresh_split == 1) = 0;
            Z(abs(Z)>lambdatB & thresh_split == 2) = 1;
        else
   
            for m = 1:size(S,1)
                for n = 1:size(S,2)
                    if thresh_split(m,n) == 1 & Md(m, n) > 0
                        Z(m, n) = sign(Z(m,n)).*max(Z(m,n)-lambda,0);
                    elseif thresh_split(m,n) == 2 & Md(m, n) > 0
                         Z(m, n) = sign(Z(m,n)).*min(Z(m,n)+Blambda,1);
                    end
                end
            end
%           
        end

  
        S = S.*(1-Md) + Z.*Md;
       

       % g  (Gradient) sub-problem 
        Gh = [diff(S,1,2), S(:,1,:) - S(:,end,:)]; % gh
        Gv = [diff(S,1,1); S(1,:,:) - S(end,:,:)]; % gv
        t = (Gh.^2 + Gv.^2) < mu/alpha;
        Gh(t)=0;    Gv(t)=0;
        
       % I subproblem
        gh = [Gh(:,end,:) - Gh(:, 1,:), -diff(Gh,1,2)];
        gv = [Gv(end,:,:) - Gv(1, :,:); -diff(Gv,1,1)];
        Fs = (Fk_FI + alpha*fft2(gh+gv))./(denKER + alpha*denGrad);
        S = real(ifft2(Fs));
        
    end
        alpha = alpha*kappa;

end
end
