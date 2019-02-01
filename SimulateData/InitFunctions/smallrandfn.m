function p = smallrandfn(n, bc)
    p = uniffn(n, bc);
    p = p + (p(2) - p(1)) / 2 * randn([n, 1]);
end
