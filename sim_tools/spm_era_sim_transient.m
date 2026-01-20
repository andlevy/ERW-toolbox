function Y = spm_era_sim_transient(trans)

Y = [];
for i = 1:numel(trans)

    pst  = trans(i).pst(1):trans(i).dt:trans(i).pst(2);
    ns   = length(pst);
    y    = zeros(ns,1);

    for e = 1:numel(trans(i).elem)

        ons  = trans(i).elem(e).set.onset;
        t    = pst-ons;
        stim = spm_vec(t>=trans(i).elem(e).set.twin(1)&t<=trans(i).elem(e).set.twin(2));

        switch trans(i).elem(e).type
            case 'sin'

                x  = spm_vec(sin(2*pi*trans(i).elem(e).set.hz*t));
                x  = x.*trans(i).elem(e).set.amp;
                x1 = y(find(stim,1,'first')-1);
                x  = x + x1;

            case 'exp'
                x = spm_vec(exp(-t/trans(i).elem(e).set.tau))-1;
                x = x.*trans(i).elem(e).set.amp;

            case 'sig'
                x  = spm_vec(1./(1+exp(-trans(i).elem(e).set.k*(t-trans(i).elem(e).set.t0))));
                x1 = y(find(stim,1,'first')-1);
                x2 = y(find(stim,1,'last')) + trans(i).elem(e).set.amp;
                x  = x1 + (x2 - x1) * x;
                
        end
        x = x.*stim;
        y = y.*~stim;
        y = y + x;
    end
    Y(:,i) = y;
end
