

function u = soft_thresholding(v, lambda)

sign_v = sign(v);

u = max(abs(v)-lambda,0).*sign_v;