function BB = load_sets(lp)
    BB = [];

    s = 2^lp;
    BB = dec2bin(0:s-1)-'0';
end