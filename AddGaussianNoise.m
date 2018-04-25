function un = AddGaussianNoise(u, strength)

% add gaussian noise on top of u
un = u + normrnd(0, strength, size(u));

end
