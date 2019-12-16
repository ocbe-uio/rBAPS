function mjono = logml2String(logml)
% Palauttaa logml:n string-esityksen.

mjono = '       ';

if isequal(logml,-Inf)
    mjono(7) = '-';
    return
end

if abs(logml)<10000
    %Ei tarvita e-muotoa
    mjono(7) = palautaYks(abs(logml),-1);
    mjono(6) = '.';
    mjono(5) = palautaYks(abs(logml),0);
    mjono(4) = palautaYks(abs(logml),1);
    mjono(3) = palautaYks(abs(logml),2);
    mjono(2) = palautaYks(abs(logml),3);
    pointer = 2;
    while mjono(pointer)=='0' & pointer<7
        mjono(pointer) = ' ';
        pointer=pointer+1;
    end
    if logml<0
        mjono(pointer-1) = '-';
    end 
else
    suurinYks = 4;
    while abs(logml)/(10^(suurinYks+1)) >= 1
        suurinYks = suurinYks+1;
    end
    if suurinYks<10
        mjono(7) = num2str(suurinYks);
        mjono(6) = 'e';
        mjono(5) = palautaYks(abs(logml),suurinYks-1);
        mjono(4) = '.';
        mjono(3) = palautaYks(abs(logml),suurinYks);
        if logml<0
            mjono(2) = '-';
        end
    elseif suurinYks>=10
        mjono(6:7) = num2str(suurinYks);
        mjono(5) = 'e';
        mjono(4) = palautaYks(abs(logml),suurinYks-1);
        mjono(3) = '.';
        mjono(2) = palautaYks(abs(logml),suurinYks);
        if logml<0
            mjono(1) = '-';
        end
    end
end

function digit = palautaYks(num,yks)
% palauttaa luvun num 10^yks termin kertoimen
% string:inä 
% yks täytyy olla kokonaisluku, joka on 
% vähintään -1:n suuruinen. Pienemmillä
% luvuilla tapahtuu jokin pyöristysvirhe.

if yks>=0
    digit = rem(num, 10^(yks+1));
    digit = floor(digit/(10^yks));
else
    digit = num*10;
    digit = floor(rem(digit,10));
end
digit = num2str(digit);