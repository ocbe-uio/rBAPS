function str = ownNum2Str(number)

absolute = abs(number);

if absolute < 1000
    str = num2str(number);
elseif absolute < 10000000
    first_three = rem(number,1000);
    next_four = (number - first_three) /1000;
    first_three = abs(first_three);
    if first_three<10
        first_three = ['00' num2str(first_three)];
    elseif first_three<100
        first_three = ['0' num2str(first_three)];
    else
        first_three = num2str(first_three);
    end;
    str = [num2str(next_four) first_three]; 
elseif absolute < 100000000
    first_four = rem(number,10000);
    next_four = (number - first_four) /10000;
    first_four = abs(first_four);
    if first_four<10
        first_four = ['000' num2str(first_four)];
    elseif first_four<100
        first_four = ['00' num2str(first_four)];
    elseif first_four<1000
        first_four = ['0' num2str(first_four)];
    else
        first_four = num2str(first_four);
    end;
    str = [num2str(next_four) first_four];     
else
    str = num2str(number);
end;   