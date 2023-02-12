function out = laplace_sine(s_lapl,omega)

sigma_phase = 3.5d7;
out(1) = omega/(s_lapl^2+omega^2) * sigma_phase; % L{sin(t)}
out(2:12) = 0;

end