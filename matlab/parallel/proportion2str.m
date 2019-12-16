function str = proportion2str(prob)
%prob belongs to [0.00, 0.01, ... ,1]. 
%str is a 4-mark presentation of proportion.

if abs(prob)<1e-3
    str = '0.00';
elseif abs(prob-1) < 1e-3;
    str = '1.00';
else
    prob = round(100*prob);
    if prob<10
        str = ['0.0' num2str(prob)];    
    else
        str = ['0.' num2str(prob)];
    end;        
end;

%-------------------------------------------------------