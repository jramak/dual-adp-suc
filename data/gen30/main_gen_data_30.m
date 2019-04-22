maxDemand = 6783.8
maxDemand_mult = [0.4, 0.6, 0.8];
% maxDemand_vec = [3000];
percent_sig = [0.15, 0.20, 0.25];
% maxDemand_vec = [1500];
% percent_sig = [0.15];

for md = maxDemand_mult
    for psig = percent_sig
        gen_dem_data_30(md, maxDemand, psig);
    end
end
