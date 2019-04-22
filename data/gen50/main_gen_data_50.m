maxDemand = 9870.5
maxDemand_mult = [0.4, 0.6, 0.8]
% maxDemand_vec = [3000];
% percent_sig = [0.15, 0.20, 0.25];
percent_sig = [0.15, 0.2, 0.25];

for md = maxDemand_mult
    for psig = percent_sig
        gen_dem_data_50(md, maxDemand, psig);
    end
end
