% guided image filter
% p: the image to be filterd
% I: the guide image
% q_i = a_k * I_i + b_k , for i in a window w_k centered at pixel k
% a_k = {(1/|w|) * Sum(I_i * p_i) - mean(I_w) * mean(p_w)} / (var_I_w + eps)
% b_k = mean(p_w) - a_k * mean(I_w)


function q = guidedFilterMono(I, p, r, eps)
s = 2 * r + 1;
mean_I = imboxfilt(I, s);
mean_p = imboxfilt(p, s);
corr_Ip = imboxfilt(I.*p, s);
cov_Ip = corr_Ip - mean_I.*mean_p;
corr_II = imboxfilt(I.^2, s);
var_I = corr_II - mean_I.^2;
a = cov_Ip ./ (var_I + eps);
b = mean_p - a.*mean_I;
mean_a = imboxfilt(a, s);
mean_b = imboxfilt(b, s);
q = mean_a.*I + mean_b;
end