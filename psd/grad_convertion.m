function G_per_cm = grad_convertion(mT_per_m)
% THis function converts grad units in mT/m to G/cm
G_per_T = 10000;
G_per_mT = G_per_T * 1e-3;
cm_per_m = 100;

G_per_cm = (mT_per_m * G_per_mT)/cm_per_m;
end