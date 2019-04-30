function [gen_state, z_prev, u_prev] = init_state_calc(filename,numGens)

T = 72;

irgdx('init_expdem_mip.gdx', 'u_s', 'z_s');
irgdx(filename, 'Lu', 'Ld')

z_prev = z_s(1:numGens, T);
u_prev = u_s(1:numGens, T);
gen_state = ones(numGens, 1) .* (Lu + Ld); % initial start state

for i = 1:numGens
    for t = 1:T
        if (gen_state(i) < Lu(i))
            gen_state(i) = gen_state(i) + 1;
        elseif ((gen_state(i) == Lu(i)) && (u_s(i, t) == 0))
            gen_state(i) = gen_state(i) + 1;
        elseif ((Lu(i) < gen_state(i)) && (gen_state(i) < (Lu(i) + Ld(i))))
            gen_state(i) = gen_state(i) + 1;
        elseif ((gen_state(i) == (Lu(i) + Ld(i))) && (u_s(i, t) == 1))
            gen_state(i) = 1;
        end
    end
end
