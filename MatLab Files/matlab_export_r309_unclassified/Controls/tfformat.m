function tfformat(num,den,k,varname)


cnt = 1;
for i = 1:2*N

    if i == 1
        fprintf('%s\n',['      call initialize(',varname,',',...
                num2str(N),num2str(k),'.0D0,'])
    end

    
    if cnt <= N

        if cnt == 1
            fprintf('%s',['     ',num2str(1),'       ']);            
            fprintf('%s','(/');
        end        
        if cnt ~= N
            fprintf('%14.12f%s',num(cnt),'D0, ');
        else
            fprintf('%14.12f%s\n',num(cnt),'D0/)');     
        end        
    else
        if cnt == N+1
            fprintf('%s',['     ',num2str(2),'       ']);            
            fprintf('%s','(/');
        end              
        if cnt ~= 2*N
            fprintf('%14.12f%s',den(cnt-N),'D0, ');
        else
            fprintf('%14.12f%s\n',den(cnt-N),'D0/))');     
        end                
    end
    
    cnt = cnt+1;
    
end






end