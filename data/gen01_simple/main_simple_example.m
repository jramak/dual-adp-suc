% maxDemand_vec = [1500, 2500, 3500, 4500];
% percent_sig = [0.15, 0.20, 0.25];
maxDemand_vec = [500];
percent_sig = [0.25];

for md = maxDemand_vec
    for psig = percent_sig
        gen_simple_example(md, psig);
    end
end
