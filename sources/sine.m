function out = sine(time,w)

sigma_phase = 3.5d7;
out(1) = sin(w*time) * sigma_phase;
out(2:12) = 0;

end