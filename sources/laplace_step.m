function out = laplace_step(s)

sigma_phase = 3.5d7;
out(1) = 1/s * sigma_phase; % L{gradino unitario}
out(2:12) = 0;

end