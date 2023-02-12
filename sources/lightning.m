function out = lightning(time,pp)
    if time<=1.4860e-05
        out(1) = ppval(pp, time)*3.5d7;
    else 
        out(1) = 0;
    end
    out(2:12) = 0;
end