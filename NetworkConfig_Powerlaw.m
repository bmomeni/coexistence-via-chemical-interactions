function C = NetworkConfig_Powerlaw(Nc,Nm,q)

C = zeros(Nc,Nm);
for nm = 1:Nm,
    nl = Nc+1;
    while nl>Nc,
        nl = randraw('zeta',q,1);
    end
    prm = randperm(Nc);
    C(prm(1:nl),nm) = 1;
end

return;