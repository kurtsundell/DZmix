function [ R2_n ] = weight_cutoffs(xx)

global x
global N
global pdp_sink 
global pdp_source

W1 = xx(1);
W2 = xx(2);

if N >= 3 
W3 = xx(3);
end
if N == 3
weights = [W1,W2,W3];
end

if N >= 4
W4 = xx(4);
end
if N == 4
weights = [W1,W2,W3,W4];
end

if N >= 5
W5 = xx(5);
end
if N == 5
weights = [W1,W2,W3,W4,W5];
end

if N >= 6
W6 = xx(6);
end
if N == 6
weights = [W1,W2,W3,W4,W5,W6];
end

if N >= 7
W7 = xx(7);
end
if N == 7
weights = [W1,W2,W3,W4,W5,W6,W7];
end

if N >= 8
W8 = xx(8);
end
if N == 8
weights = [W1,W2,W3,W4,W5,W6,W7,W8];
end

if N >= 9
W9 = xx(9);
end
if N == 9
weights = [W1,W2,W3,W4,W5,W6,W7,W8,W9];
end

if N >= 10
W10 = xx(10);
end
if N == 10
weights = [W1,W2,W3,W4,W5,W6,W7,W8,W9,W10];
end

w3 = weights.*0.01;

pdp = sum(pdp_source.*(repmat(w3,length(x),1)),2);

R2 = ((sum((pdp_sink - mean(pdp_sink)).*(pdp - mean(pdp))))/(sqrt((sum((pdp_sink - mean(pdp_sink)).*(pdp_sink - mean(pdp_sink))))*(sum((pdp - mean(pdp)).*(pdp - mean(pdp)))))))*...
	((sum((pdp_sink - mean(pdp_sink)).*(pdp - mean(pdp))))/(sqrt((sum((pdp_sink - mean(pdp_sink)).*(pdp_sink - mean(pdp_sink))))*(sum((pdp - mean(pdp)).*(pdp - mean(pdp)))))));

R2_n = 1 - R2;