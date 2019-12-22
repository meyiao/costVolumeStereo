% guided image filter
% p: the image to be filterd
% I: the guide image
% q_i = a_k * I_i + b_k , for i in a window w_k centered at pixel k
% a_k is a 3x1 vector
% I_i is a 3x1 vector
% a_k = (var_I_w + eps * eye(3))^(-1) * {1/|w| * sum(I_i * p_i) - mean(I_w) * mean(p_w)}
% b_k = mean(p_w) - a_k' * mean(I_w)

function q = guidedFilterColor(I, p, r, eps)

s = 2 * r + 1;

Ir = I(:, :, 1);
Ig = I(:, :, 2);
Ib = I(:, :, 3);

mean_Ir = imboxfilt(Ir, s);
mean_Ig = imboxfilt(Ig, s);
mean_Ib = imboxfilt(Ib, s);

mean_p = imboxfilt(p, s);
cov_Irp = imboxfilt(Ir.*p, s) - mean_Ir.*mean_p;
cov_Igp = imboxfilt(Ig.*p, s) - mean_Ig.*mean_p;
cov_Ibp = imboxfilt(Ib.*p, s) - mean_Ib.*mean_p;

var_Irr = imboxfilt(Ir.*Ir, s) - mean_Ir.*mean_Ir + eps;
var_Irg = imboxfilt(Ir.*Ig, s) - mean_Ir.*mean_Ig;
var_Irb = imboxfilt(Ir.*Ib, s) - mean_Ir.*mean_Ib;
var_Igg = imboxfilt(Ig.*Ig, s) - mean_Ig.*mean_Ig + eps;
var_Igb = imboxfilt(Ig.*Ib, s) - mean_Ig.*mean_Ib;
var_Ibb = imboxfilt(Ib.*Ib, s) - mean_Ib.*mean_Ib + eps;

[H, W] = size(p);
a = zeros(H, W, 3);
for i = 1 : H
    for j = 1 : W
        Sigma = [var_Irr(i,j), var_Irg(i,j), var_Irb(i,j);...
                 var_Irg(i,j), var_Igg(i,j), var_Igb(i,j);...
                 var_Irb(i,j), var_Igb(i,j), var_Ibb(i,j)];
        cov_Ip = [cov_Irp(i,j), cov_Igp(i,j), cov_Ibp(i,j)];
        a(i,j,:) = cov_Ip * inv(Sigma);
    end
end

b = mean_p - a(:,:,1).*mean_Ir - a(:,:,2).*mean_Ig - a(:,:,3).*mean_Ib;
mean_a = imboxfilt(a, s);
mean_b = imboxfilt(b, s);
q = mean_a(:,:,1).*Ir + mean_a(:,:,2).*Ig + mean_a(:,:,3).*Ib + mean_b;

return





