maxDemand = 4660
maxDemand_mult = [0.4, 0.6, 0.8]
percent_sig = [0.15, 0.20, 0.25];
% maxDemand_vec = [1500];
% percent_sig = [0.15];

for md = maxDemand_mult
    for psig = percent_sig
        gen_dem_data_15(md, maxDemand, psig);
    end
end
