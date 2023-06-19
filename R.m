function Rw = R(p,dt,nDOF,f,Sdf,Cor,sigmaw)
    for i = 1:nDOF
        for j = 1:nDOF
            Rw(i,j) = sum(Cor(i,j)*sigmaw(i)*sigmaw(j).*Sdf.*cos(2*pi.*f*p*dt)*dt);
        end
    end
end