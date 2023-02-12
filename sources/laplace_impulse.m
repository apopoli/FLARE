function out = laplace_impulse(s)

sigma_phase = 3.5d7;
out(1) = 1 * sigma_phase; % L{\delta}
out(2:12) = 0;

end