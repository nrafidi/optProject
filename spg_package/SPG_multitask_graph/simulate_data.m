function [Y, X, B]=simulate_data(N,J,K)

         rand('seed', sum(clock));         
                 
         %% generating X and te_X
         X=randn(N,J);             
         X=zscore(X);
         
         %% generate B matrix %10 of 1 for each group, 5% for crossing two
   
         B=zeros(J,K);
         if (K<10)
             B=randn(J,K);
             B(B(:)<median(B(:)))=0;
         else
            g_size=10;
            NG=round(K/g_size);
            nj_1=round(J*0.1);
            for g=1:NG
                j_idx=randsample(J, nj_1);
                B(j_idx, (g-1)*g_size+1: g*g_size)=1;
            end 
             nj_2=round(J*0.05);
             for g=1:NG-1
                 j_idx=randsample(J, nj_2);
                B(j_idx, (g-1)*g_size+1: (g+1)*g_size)=1;
             end         
            nj_3=round(J*0.01);
            for g=1:NG-2
                 j_idx=randsample(J, nj_3);
                 B(j_idx, (g-1)*g_size+1: (g+2)*g_size)=1;
            end  
         end 
        
         %% Generate output Y
         Y=X*B+randn(N,K);   
end 