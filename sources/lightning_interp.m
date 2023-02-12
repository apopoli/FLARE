function out = lightning_interp(time,Vfulminazione)
    out(1) = interp1(Vfulminazione(:,1),Vfulminazione(:,2),time) * 3.5d7;
    out(2:12) = 0;
end